/*
 * GenomeTrackImportParquet.cpp
 *
 * gtrack.import_parquet: the value columns of Parquet files as dense tracks, one track per
 * column, binned exactly as gtrackimportwig bins a tab-delimited file. With multitasking on,
 * the chromosomes and columns are split between forked kids.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <tuple>

#include "DataFrameUtils.h"
#include "GenomeTrack.h"
#include "GenomeTrackFixedBin.h"
#include "GTrack2DImport.h"
#include "Parquet.h"
#include "UnknownChroms.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

// GenomeTrackImportWig.cpp
void verify_chroms_matched(const GenomeChromKey &chromkey, const char *fname, const UnknownChroms &unknown, uint64_t num_matched, SEXP envir);
void report_skipped_chroms(const char *fname, const UnknownChroms &unknown, uint64_t num_matched, SEXP envir);

namespace {

// A unit of work is a chromosome and a batch of at most this many columns
const int MAX_UNIT_COLS = 16;

// A process imports its units of a batch together, at most this many rows at a time (or one
// chromosome, if it has more): a row group that holds rows of several of the chromosomes is
// decoded once for all of them. Their rows are held in memory together, about 24 bytes a row.
const uint64_t MAX_ROWS_AT_ONCE = 1 << 22;

// The value of one bin, as gtrackimportwig computes it: sum and covered_bases from the
// non-NaN values that overlap the bin, the rest of the bin filled with defval unless it is NaN
inline float bin_value(double sum, uint64_t covered_bases, uint64_t bin_size_actual, double defval)
{
	if (covered_bases < bin_size_actual && !std::isnan(defval)) {
		uint64_t uncovered = bin_size_actual - covered_bases;
		sum += defval * uncovered;
		covered_bases += uncovered;
	}
	return covered_bases ? (sum / covered_bases) : numeric_limits<float>::quiet_NaN();
}

// Rows [begin, end) of a row group, all on one chromosome
struct Run {
	int     file;
	int     row_group;
	int64_t begin;
	int64_t end;
};

// Rows that run past the end of their chromosome: clipped to it, or dropped when they start at or
// past it. Kept in shared memory, one per kid, so that the parent reports them once.
struct ClipReport {
	uint64_t num_clipped;
	uint64_t num_dropped;
	// the first such row, in the order of the files and their rows; file < 0 if there is none
	int      file;
	int      row_group;
	int64_t  row;  // in the row group
	int      chromid;
	int64_t  start;
	int64_t  end;

	void clear() { num_clipped = num_dropped = 0; file = -1; }

	void note(bool dropped, int f, int rg, int64_t r, int c, int64_t s, int64_t e)
	{
		++(dropped ? num_dropped : num_clipped);
		if (file < 0 || make_tuple(f, rg, r) < make_tuple(file, row_group, row)) {
			file = f;
			row_group = rg;
			row = r;
			chromid = c;
			start = s;
			end = e;
		}
	}

	void merge(const ClipReport &other)
	{
		if (other.file >= 0) {
			uint64_t clipped = num_clipped + other.num_clipped;
			uint64_t dropped = num_dropped + other.num_dropped;
			note(false, other.file, other.row_group, other.row, other.chromid, other.start, other.end);
			num_clipped = clipped;
			num_dropped = dropped;
		}
	}
};

// Imports columns of a set of Parquet files. index() finds which rows of which row groups lie
// on each chromosome (reading only the chrom columns); import_chroms() then writes the track
// files of some chromosomes for any subset of the columns.
class ParquetImport {
public:
	enum { CHROM, START, END, NUM_COORDS };

	// names: what messages call the files
	ParquetImport(const GenomeChromKey &chromkey, const vector<string> &files, const vector<string> &names, const vector<string> &columns,
	              const vector<string> &dirs, unsigned binsize, double defval) :
		m_chromkey(chromkey), m_files(files), m_names(names), m_columns(columns), m_dirs(dirs), m_binsize(binsize), m_defval(defval) {}

	// Opens the files and finds the columns in each
	void open();

	// Fills runs, unknown chromosome names and the files they came from
	void index();

	// Writes the track files of the chromosomes for the columns. Each column chunk that holds
	// their rows is decoded once, however many of the chromosomes it holds. Rows that run past
	// the end of their chromosome are clipped (or dropped) and, unless report is NULL, noted in it.
	void import_chroms(const vector<int> &chromids, const vector<int> &cols, ClipReport *report);

	// Row row of row group row_group of the file, counted from 0 in the whole file
	int64_t file_row(int file, int row_group, int64_t row) const;

	// The warning about the rows of report, "" if there were none
	string clip_message(const ClipReport &report) const;

	const UnknownChroms &unknown() const { return m_unknown; }
	uint64_t num_matched_chroms() const;
	uint64_t num_rows(int chromid) const;

	// Where the first row of the chromosome is: file, row group, row; chromosomes without rows
	// come after all the others
	tuple<int, int, int64_t> first_row(int chromid) const;

	// The files that named chromosomes the database does not have, for messages
	string unknown_chroms_files() const;

private:
	const GenomeChromKey             &m_chromkey;
	vector<string>                    m_files;
	vector<string>                    m_names;
	vector<string>                    m_columns;
	vector<string>                    m_dirs;
	unsigned                          m_binsize;
	double                            m_defval;
	vector<unique_ptr<ParquetReader>> m_readers;
	vector<vector<int>>               m_col_idx;  // per file: chrom, start, end, then the value columns
	vector<vector<Run>>               m_runs;     // per chromid
	UnknownChroms                     m_unknown;
	vector<bool>                      m_file_has_unknown;

	// The rows of a chromosome, in the order of its runs
	struct ChromRows {
		vector<int64_t>  starts;
		vector<int64_t>  ends;     // clipped to the chromosome
		vector<float>    vals;     // of the column being written
		vector<uint32_t> order;    // the rows sorted by start, unless aligned; dropped rows left out
		vector<char>     dropped;  // the rows that start at or past the chromosome end; empty if none
		bool             aligned;  // every row (not dropped) is exactly a bin

		bool is_dropped(size_t i) const { return !dropped.empty() && dropped[i]; }
	};

	// A run of chromosome chromids[chrom] of import_chroms()
	struct ChromRun {
		size_t     chrom;
		const Run *run;
	};

	[[noreturn]] void rows_error(const char *what, const Run &run1, const Run &run2) const;
	[[noreturn]] void row_error(const char *what, const Run &run, int64_t row, int chromid, int64_t start, int64_t end) const;

	template <typename T, typename D>
	void gather(const vector<ChromRun> &runs, int col_idx, void (ParquetReader::*read)(int, int, vector<T> &), D dest);
	void prepare_rows(int chromid, ChromRows &rows, ClipReport *report) const;
	void write_bins(int chromid, const ChromRows &rows, int col, vector<float> &bins) const;
};

void ParquetImport::open()
{
	uint64_t genome_size = 0;
	for (int chromid = 0; chromid < (int)m_chromkey.get_num_chroms(); ++chromid)
		genome_size += m_chromkey.get_chrom_size(chromid);

	m_readers.clear();
	m_col_idx.clear();
	int64_t total_rows = 0;
	for (size_t f = 0; f < m_files.size(); ++f) {
		const string &file = m_names[f];
		m_readers.emplace_back(new ParquetReader());
		ParquetReader &reader = *m_readers.back();
		reader.open(m_files[f], file);
		total_rows += reader.metadata().num_rows;

		// Rows that do not overlap cannot outnumber the bases (and overlaps are an error). Checked
		// before decoding, which allocates by the row count: a small corrupt file could ask for
		// tens of GB.
		for (int rg = 0; rg < reader.num_row_groups(); ++rg) {
			if ((uint64_t)reader.num_rows(rg) > genome_size)
				verror("File %s, row group %d: %lld rows, more than the %llu bases of the genome", file.c_str(), rg,
				       (long long)reader.num_rows(rg), (unsigned long long)genome_size);
		}

		vector<string> names = { "chrom", "start", "end" };
		names.insert(names.end(), m_columns.begin(), m_columns.end());
		m_col_idx.emplace_back();
		for (const string &name : names) {
			int idx = reader.column_index(name);
			if (idx < 0)
				verror("File %s has no column \"%s\"", file.c_str(), name.c_str());
			m_col_idx.back().push_back(idx);
		}
	}

	// Only the total counts: an empty shard among non-empty ones is fine
	if (!total_rows) {
		string names;
		for (const string &file : m_names)
			names += (names.empty() ? "" : ", ") + file;
		verror("%s %s no rows", names.c_str(), m_names.size() == 1 ? "has" : "have");
	}
}

void ParquetImport::index()
{
	vector<string> chroms;

	m_runs.assign(m_chromkey.get_num_chroms(), vector<Run>());
	m_unknown.clear();
	m_file_has_unknown.assign(m_files.size(), false);

	for (size_t f = 0; f < m_files.size(); ++f) {
		for (int rg = 0; rg < m_readers[f]->num_row_groups(); ++rg) {
			// (an empty row group may type an all-missing column as null, e.g. pandas' empty chrom column)
			if (!m_readers[f]->num_rows(rg))
				continue;
			m_readers[f]->read_strings(rg, m_col_idx[f][CHROM], chroms);
			for (size_t i = 0, j; i < chroms.size(); i = j) {
				for (j = i + 1; j < chroms.size() && chroms[j] == chroms[i]; ++j)
					;
				int chromid = m_chromkey.find_chromid(chroms[i]);
				if (chromid < 0) {
					m_unknown.record(chroms[i]);
					m_file_has_unknown[f] = true;
				} else
					m_runs[chromid].push_back(Run{ (int)f, rg, (int64_t)i, (int64_t)j });
			}
			check_interrupt();
		}
	}
}

uint64_t ParquetImport::num_rows(int chromid) const
{
	uint64_t num = 0;
	for (const Run &run : m_runs[chromid])
		num += run.end - run.begin;
	return num;
}

tuple<int, int, int64_t> ParquetImport::first_row(int chromid) const
{
	if (m_runs[chromid].empty())
		return make_tuple(numeric_limits<int>::max(), 0, (int64_t)0);
	const Run &run = m_runs[chromid].front();
	return make_tuple(run.file, run.row_group, run.begin);
}

uint64_t ParquetImport::num_matched_chroms() const
{
	uint64_t num = 0;
	for (const vector<Run> &runs : m_runs)
		num += !runs.empty();
	return num;
}

string ParquetImport::unknown_chroms_files() const
{
	string res;
	size_t num = 0;

	for (size_t f = 0; f < m_files.size(); ++f) {
		if (!m_file_has_unknown[f])
			continue;
		if (num < (size_t)UnknownChroms::MAX_REPORTED)
			res += (num ? ", " : "") + m_names[f];
		++num;
	}
	if (num > (size_t)UnknownChroms::MAX_REPORTED)
		res += " (and " + to_string(num - UnknownChroms::MAX_REPORTED) + " more files)";
	return res;
}

void ParquetImport::rows_error(const char *what, const Run &run1, const Run &run2) const
{
	if (run1.file == run2.file && run1.row_group == run2.row_group)
		verror("File %s, row group %d: %s", m_names[run1.file].c_str(), run1.row_group, what);
	verror("File %s, row group %d and file %s, row group %d: %s", m_names[run1.file].c_str(), run1.row_group,
	       m_names[run2.file].c_str(), run2.row_group, what);
}

int64_t ParquetImport::file_row(int file, int row_group, int64_t row) const
{
	for (int rg = 0; rg < row_group; ++rg)
		row += m_readers[file]->num_rows(rg);
	return row;
}

// Rows are numbered from 0 in the whole file, as pyarrow and pandas number them
void ParquetImport::row_error(const char *what, const Run &run, int64_t row, int chromid, int64_t start, int64_t end) const
{
	verror("File %s, row %lld (row group %d): %s (%s, %lld, %lld)", m_names[run.file].c_str(), (long long)file_row(run.file, run.row_group, row),
	       run.row_group, what, m_chromkey.id2chrom(chromid).c_str(), (long long)start, (long long)end);
}

string ParquetImport::clip_message(const ClipReport &report) const
{
	if (report.file < 0)
		return string();

	char msg[10000];
	snprintf(msg, sizeof(msg),
	         "%llu row(s) running past the end of their chromosome were clipped to it, and %llu row(s) starting at or past "
	         "the chromosome end were dropped; the first is row %lld (row group %d) of %s: (%s, %lld, %lld), chromosome size %llu",
	         (unsigned long long)report.num_clipped, (unsigned long long)report.num_dropped,
	         (long long)file_row(report.file, report.row_group, report.row), report.row_group, m_names[report.file].c_str(),
	         m_chromkey.id2chrom(report.chromid).c_str(), (long long)report.start, (long long)report.end,
	         (unsigned long long)m_chromkey.get_chrom_size(report.chromid));
	return msg;
}

// Appends the values of column col_idx (an index into m_col_idx) of each run to dest(run.chrom).
// runs are in the order of the files and row groups, so each column chunk is decoded once.
template <typename T, typename D>
void ParquetImport::gather(const vector<ChromRun> &runs, int col_idx, void (ParquetReader::*read)(int, int, vector<T> &), D dest)
{
	vector<T> buf;

	for (size_t r = 0; r < runs.size(); ++r) {
		const Run &run = *runs[r].run;

		if (!r || run.file != runs[r - 1].run->file || run.row_group != runs[r - 1].run->row_group) {
			((*m_readers[run.file]).*read)(run.row_group, m_col_idx[run.file][col_idx], buf);
			check_interrupt();
		}
		if ((uint64_t)run.end > buf.size())
			rows_error("fewer rows than when it was indexed", run, run);

		auto &out = dest(runs[r].chrom);
		out.insert(out.end(), buf.begin() + run.begin, buf.begin() + run.end);
	}
}

// Checks the coordinates of the rows, clips the rows that run past the end of the chromosome (or
// drops them if they start at or past it), sets aligned and, unless the rows are aligned, sorts
// them
void ParquetImport::prepare_rows(int chromid, ChromRows &rows, ClipReport *report) const
{
	const vector<Run> &runs = m_runs[chromid];
	uint64_t chromsize = m_chromkey.get_chrom_size(chromid);
	uint64_t num_bins = (chromsize + m_binsize - 1) / m_binsize;
	const vector<int64_t> &starts = rows.starts;
	vector<int64_t> &ends = rows.ends;
	vector<size_t> run_of_row;  // the first row of each run, for messages

	rows.aligned = true;
	rows.dropped.clear();
	for (size_t r = 0, i = 0; r < runs.size(); ++r) {
		const Run &run = runs[r];

		run_of_row.push_back(i);
		for (int64_t k = run.begin; k < run.end; ++k, ++i) {
			int64_t start = starts[i];
			int64_t end = ends[i];

			if (start < 0)
				row_error("start coordinate is negative", run, k, chromid, start, end);
			if (start >= end)
				row_error("start coordinate exceeds or equals the end coordinate", run, k, chromid, start, end);
			if ((uint64_t)end > chromsize) {
				bool drop = (uint64_t)start >= chromsize;

				if (report)
					report->note(drop, run.file, run.row_group, k, chromid, start, end);
				if (drop) {
					if (rows.dropped.empty())
						rows.dropped.resize(starts.size(), 0);
					rows.dropped[i] = 1;
					continue;
				}
				ends[i] = end = chromsize;
			}
			rows.aligned = rows.aligned && start % m_binsize == 0 && (uint64_t)end == min((uint64_t)start + m_binsize, chromsize);
		}
	}

	auto run_of = [&](size_t row) -> const Run & {
		return runs[upper_bound(run_of_row.begin(), run_of_row.end(), row) - run_of_row.begin() - 1];
	};

	// Rows that are exactly the bins (the common case) are placed by their start; any other rows
	// are sorted by start, as GenomeArraysCsv sorts them. Overlapping or duplicate rows are errors
	// either way.
	if (rows.aligned) {
		vector<uint32_t> row_of_bin(num_bins, UINT32_MAX);
		for (size_t i = 0; i < starts.size(); ++i) {
			if (rows.is_dropped(i))
				continue;
			uint32_t &row = row_of_bin[starts[i] / m_binsize];
			if (row != UINT32_MAX)
				rows_error("rows overlap", run_of(row), run_of(i));
			row = i;
		}
	} else {
		vector<uint32_t> &order = rows.order;
		order.clear();
		for (size_t i = 0; i < starts.size(); ++i) {
			if (!rows.is_dropped(i))
				order.push_back(i);
		}
		sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) { return starts[a] < starts[b] || (starts[a] == starts[b] && ends[a] < ends[b]); });
		for (size_t k = 1; k < order.size(); ++k) {
			if (ends[order[k - 1]] > starts[order[k]])
				rows_error("rows overlap", run_of(order[k - 1]), run_of(order[k]));
		}
	}
}

// Bins rows.vals and writes the track file of column col
void ParquetImport::write_bins(int chromid, const ChromRows &rows, int col, vector<float> &bins) const
{
	uint64_t chromsize = m_chromkey.get_chrom_size(chromid);
	uint64_t num_bins = (chromsize + m_binsize - 1) / m_binsize;
	const vector<int64_t> &starts = rows.starts;
	const vector<int64_t> &ends = rows.ends;
	const vector<float> &vals = rows.vals;
	const vector<uint32_t> &order = rows.order;
	char filename[FILENAME_MAX];

	bins.resize(num_bins);
	if (rows.aligned) {
		for (uint64_t bin = 0; bin < num_bins; ++bin)
			bins[bin] = bin_value(0, 0, min((bin + 1) * m_binsize, chromsize) - bin * m_binsize, m_defval);
		for (size_t i = 0; i < starts.size(); ++i) {
			if (rows.is_dropped(i))
				continue;
			uint64_t len = ends[i] - starts[i];
			double sum = 0;
			uint64_t covered_bases = 0;
			float v = vals[i];

			if (!std::isnan(v)) {
				sum += (double)v * len;  // in double, as in gtrackimportwig: a float product can overflow
				covered_bases += len;
			}
			bins[starts[i] / m_binsize] = bin_value(sum, covered_bases, len, m_defval);
		}
	} else {
		// gtrackimportwig's dense loop
		size_t k = 0;
		for (uint64_t start_coord = 0, bin = 0; start_coord < chromsize; start_coord += m_binsize, ++bin) {
			double sum = 0;
			uint64_t covered_bases = 0;
			uint64_t end_coord = min(start_coord + m_binsize, chromsize);

			while (k < order.size() && (uint64_t)ends[order[k]] <= start_coord)
				++k;

			for (size_t c = k; c < order.size() && (uint64_t)starts[order[c]] < end_coord; ++c) {
				uint64_t overlap_start = max(start_coord, (uint64_t)starts[order[c]]);
				uint64_t overlap_end = min(end_coord, (uint64_t)ends[order[c]]);

				if (overlap_end > overlap_start) {
					float v = vals[order[c]];
					if (!std::isnan(v)) {
						uint64_t overlap_len = overlap_end - overlap_start;
						sum += (double)v * overlap_len;
						covered_bases += overlap_len;
					}
				}
			}
			bins[bin] = bin_value(sum, covered_bases, end_coord - start_coord, m_defval);
		}
	}

	snprintf(filename, sizeof(filename), "%s/%s", m_dirs[col].c_str(), GenomeTrack::get_1d_filename(m_chromkey, chromid).c_str());
	GenomeTrackFixedBin gtrack;
	gtrack.init_write(filename, m_binsize, chromid);
	gtrack.write_next_bins(bins.data(), bins.size());
	gtrack.flush_writes();
}

void ParquetImport::import_chroms(const vector<int> &chromids, const vector<int> &cols, ClipReport *report)
{
	vector<ChromRows> rows(chromids.size());
	vector<ChromRun> runs;  // the runs of all the chromosomes, in the order of the files and row groups
	vector<float> bins;

	for (size_t c = 0; c < chromids.size(); ++c) {
		for (const Run &run : m_runs[chromids[c]])
			runs.push_back(ChromRun{ c, &run });
	}
	// index() made the runs of each chromosome in this order, so the sort keeps it
	sort(runs.begin(), runs.end(), [](const ChromRun &a, const ChromRun &b) {
		return make_tuple(a.run->file, a.run->row_group, a.run->begin) < make_tuple(b.run->file, b.run->row_group, b.run->begin);
	});

	gather(runs, START, &ParquetReader::read_int64s, [&](size_t c) -> vector<int64_t> & { return rows[c].starts; });
	gather(runs, END, &ParquetReader::read_int64s, [&](size_t c) -> vector<int64_t> & { return rows[c].ends; });
	for (size_t c = 0; c < chromids.size(); ++c)
		prepare_rows(chromids[c], rows[c], report);

	for (int col : cols) {
		for (ChromRows &chrom_rows : rows)
			chrom_rows.vals.clear();
		gather(runs, NUM_COORDS + col, &ParquetReader::read_doubles, [&](size_t c) -> vector<float> & { return rows[c].vals; });
		for (size_t c = 0; c < chromids.size(); ++c) {
			write_bins(chromids[c], rows[c], col, bins);
			check_interrupt();
		}
	}
}

// Unsigned integers are reported by their physical type: they are numbers too
const char *parquet_type_name(const ParquetColumn &col, const parquet::FileMetaData &meta, size_t idx)
{
	if (meta.schema[idx + 1].repetition_type == parquet::FieldRepetitionType::REPEATED || !col.other_type.empty())
		return "other";
	switch (col.physical_type) {
	case parquet::Type::INT32: return "int32";
	case parquet::Type::INT64: return "int64";
	case parquet::Type::FLOAT: return "float";
	case parquet::Type::DOUBLE: return "double";
	case parquet::Type::BYTE_ARRAY: return "string";
	}
	return "other";
}

}

extern "C" {

// The columns of a Parquet file: data.frame(name, type), type one of int32, int64, float,
// double, string, other. name is what messages call the file.
SEXP parquet_columns(SEXP _file, SEXP _name, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_file) || Rf_length(_file) != 1 || !Rf_isString(_name) || Rf_length(_name) != 1)
			verror("File argument is not a string");

		ParquetReader reader;
		reader.open(CHAR(STRING_ELT(_file, 0)), CHAR(STRING_ELT(_name, 0)));

		const vector<ParquetColumn> &cols = reader.columns();
		SEXP answer = DataFrameUtils::create_data_frame(cols.size(), 2);
		SEXP names = rprotect_ptr(RSaneAllocVector(STRSXP, cols.size()));
		SEXP types = rprotect_ptr(RSaneAllocVector(STRSXP, cols.size()));
		SEXP colnames = Rf_getAttrib(answer, R_NamesSymbol);

		for (size_t i = 0; i < cols.size(); ++i) {
			SET_STRING_ELT(names, i, RSaneMkChar(cols[i].name.c_str()));
			SET_STRING_ELT(types, i, RSaneMkChar(parquet_type_name(cols[i], reader.metadata(), i)));
		}
		SET_VECTOR_ELT(answer, 0, names);
		SET_VECTOR_ELT(answer, 1, types);
		SET_STRING_ELT(colnames, 0, RSaneMkChar("name"));
		SET_STRING_ELT(colnames, 1, RSaneMkChar("type"));
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	} catch (const exception &e) {
		rerror("%s", e.what());
	}
	return R_NilValue;
}

// The Arrow field metadata of the columns of a Parquet file (see arrow_field_metadata): a list named
// by column of named character vectors (key = value), only for the columns of _columns (NULL: all)
// that have some. name is what messages call the file.
SEXP parquet_field_metadata(SEXP _file, SEXP _name, SEXP _columns, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_file) || Rf_length(_file) != 1 || !Rf_isString(_name) || Rf_length(_name) != 1)
			verror("File argument is not a string");
		if (!Rf_isNull(_columns) && !Rf_isString(_columns))
			verror("Columns argument is not a character vector");

		ParquetReader reader;
		reader.open(CHAR(STRING_ELT(_file, 0)), CHAR(STRING_ELT(_name, 0)));

		map<string, vector<pair<string, string>>> md = arrow_field_metadata(reader.metadata());
		const char *name = CHAR(STRING_ELT(_name, 0));

		if (!Rf_isNull(_columns)) {
			set<string> wanted;
			for (int i = 0; i < Rf_length(_columns); ++i)
				wanted.insert(CHAR(STRING_ELT(_columns, i)));
			for (auto icol = md.begin(); icol != md.end(); )
				icol = wanted.count(icol->first) ? next(icol) : md.erase(icol);
		}

		// An R string cannot hold a NUL byte: such a column name, key or value is left out, with a
		// warning, rather than failing the import
		auto has_nul = [](const string &s) { return s.find('\0') != string::npos; };
		auto shown = [](const string &s) {
			string res;
			for (char c : s)
				res += c ? string(1, c) : string("\\0");
			return res;
		};
		for (auto icol = md.begin(); icol != md.end(); ) {
			if (has_nul(icol->first)) {
				add_pending_diagnostic(_envir, "warning",
				                       ("The field metadata of column \"" + shown(icol->first) + "\" of " + name +
				                        " was left out: the column name holds a NUL byte").c_str());
				icol = md.erase(icol);
				continue;
			}
			vector<pair<string, string>> kept;
			for (const auto &kv : icol->second) {
				if (has_nul(kv.first) || has_nul(kv.second))
					add_pending_diagnostic(_envir, "warning",
					                       ("Field metadata \"" + shown(kv.first) + "\" of column \"" + icol->first + "\" of " + name +
					                        " was left out: it holds a NUL byte, which an attribute cannot").c_str());
				else
					kept.push_back(kv);
			}
			icol->second.swap(kept);
			++icol;
		}

		SEXP answer = rprotect_ptr(RSaneAllocVector(VECSXP, md.size()));
		SEXP cols = rprotect_ptr(RSaneAllocVector(STRSXP, md.size()));
		size_t i = 0;

		for (const auto &col : md) {
			SEXP values = rprotect_ptr(RSaneAllocVector(STRSXP, col.second.size()));
			SEXP keys = rprotect_ptr(RSaneAllocVector(STRSXP, col.second.size()));
			for (size_t k = 0; k < col.second.size(); ++k) {
				SET_STRING_ELT(keys, k, RSaneMkCharLenCE(col.second[k].first.data(), col.second[k].first.size(), CE_NATIVE));
				SET_STRING_ELT(values, k, RSaneMkCharLenCE(col.second[k].second.data(), col.second[k].second.size(), CE_NATIVE));
			}
			Rf_setAttrib(values, R_NamesSymbol, keys);
			SET_VECTOR_ELT(answer, i, values);
			SET_STRING_ELT(cols, i, RSaneMkCharLenCE(col.first.data(), col.first.size(), CE_NATIVE));
			++i;
		}
		Rf_setAttrib(answer, R_NamesSymbol, cols);
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	} catch (const exception &e) {
		rerror("%s", e.what());
	}
	return R_NilValue;
}

// dirs[i] is an existing empty directory that receives the per-chromosome files of columns[i],
// for every chromosome of the database. names[i] is what messages call files[i]: gtrack.import
// reads a gzipped file from a temporary copy.
SEXP gtrackimportparquet(SEXP _dirs, SEXP _files, SEXP _names, SEXP _columns, SEXP _binsize, SEXP _defval, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_dirs) || !Rf_isString(_columns) || Rf_length(_dirs) != Rf_length(_columns) || !Rf_length(_columns))
			verror("dirs and columns must be character vectors of the same, positive, length");
		if (!Rf_isString(_files) || !Rf_length(_files))
			verror("files must be a non-empty character vector");
		if (!Rf_isString(_names) || Rf_length(_names) != Rf_length(_files))
			verror("names must be a character vector of the length of files");
		if ((!Rf_isReal(_binsize) && !Rf_isInteger(_binsize)) || Rf_length(_binsize) != 1)
			verror("Binsize argument is not a number");
		if ((!Rf_isReal(_defval) && !Rf_isInteger(_defval)) || Rf_length(_defval) != 1)
			verror("Defval argument is not a number");

		double dbinsize = Rf_isReal(_binsize) ? REAL(_binsize)[0] : INTEGER(_binsize)[0];
		double defval = Rf_isReal(_defval) ? REAL(_defval)[0] : INTEGER(_defval)[0];

		// before the cast: converting NaN, Inf or anything >= 2^32 to unsigned is undefined
		if (!(dbinsize > 0) || dbinsize >= 4294967296.0 || dbinsize != floor(dbinsize))
			verror("Invalid value of binsize argument: %g", dbinsize);
		unsigned binsize = (unsigned)dbinsize;

		vector<string> dirs, files, names, columns;
		for (int i = 0; i < Rf_length(_dirs); ++i) {
			dirs.push_back(CHAR(STRING_ELT(_dirs, i)));
			columns.push_back(CHAR(STRING_ELT(_columns, i)));
		}
		for (int i = 0; i < Rf_length(_files); ++i) {
			files.push_back(CHAR(STRING_ELT(_files, i)));
			names.push_back(CHAR(STRING_ELT(_names, i)));
		}

		IntervUtils iu(_envir);
		const GenomeChromKey &chromkey = iu.get_chromkey();
		ParquetImport import(chromkey, files, names, columns, dirs, binsize, defval);

		import.open();
		import.index();

		string srcname = import.unknown_chroms_files();
		uint64_t num_matched = import.num_matched_chroms();
		verify_chroms_matched(chromkey, srcname.c_str(), import.unknown(), num_matched, _envir);

		// Units: each chromosome times batches of columns of nearly equal size, ordered by batch and
		// then by where the chromosome's rows start in the files. Kid i takes the i-th of num_kids
		// contiguous blocks of units of about equal work (rows plus bins) and imports the units of
		// a batch together (see MAX_ROWS_AT_ONCE): a row group that holds many chromosomes is
		// decoded once by each kid that needs it, not once per chromosome. A batch decodes the
		// coordinates again, so in one process all the columns are one batch.
		int num_cols = columns.size();
		int num_batches = (num_cols + MAX_UNIT_COLS - 1) / MAX_UNIT_COLS;
		int num_kids = import_num_kids(iu, chromkey.get_num_chroms() * num_batches);
		if (num_kids == 1)
			num_batches = 1;

		vector<vector<int>> batch_cols(num_batches);
		for (int batch = 0; batch < num_batches; ++batch) {
			batch_cols[batch].resize(num_cols * (batch + 1) / num_batches - num_cols * batch / num_batches);
			iota(batch_cols[batch].begin(), batch_cols[batch].end(), num_cols * batch / num_batches);
		}

		vector<int> chromids(chromkey.get_num_chroms());
		iota(chromids.begin(), chromids.end(), 0);
		stable_sort(chromids.begin(), chromids.end(), [&](int a, int b) { return import.first_row(a) < import.first_row(b); });

		struct Unit {
			int      chromid;
			int      batch;
			uint64_t rows;
			int      kid;
		};
		vector<Unit> units;
		vector<uint64_t> work;
		uint64_t total_work = 0;

		for (int batch = 0; batch < num_batches; ++batch) {
			for (int chromid : chromids) {
				uint64_t rows = import.num_rows(chromid);
				units.push_back(Unit{ chromid, batch, rows, 0 });
				work.push_back(rows + (chromkey.get_chrom_size(chromid) + binsize - 1) / binsize);
				total_work += work.back();
			}
		}
		for (size_t u = 0, before = 0; u < units.size(); before += work[u++])
			units[u].kid = min(num_kids - 1, (int)((double)before / total_work * num_kids));

		// Each kid notes the rows it clips or drops in its own report; only units of the first batch
		// note them, as every batch reads the same rows
		SharedMem shm(num_kids * sizeof(ClipReport));
		ClipReport *reports = (ClipReport *)shm.addr;
		for (int kid = 0; kid < num_kids; ++kid)
			reports[kid].clear();

		// Kids read through the file descriptors they inherit, not by opening the files again: a file
		// replaced since index() would no longer match it
		run_kids(iu, num_kids, [&](int kid) {
			auto by_kid = [](const Unit &unit, int k) { return unit.kid < k; };
			size_t begin = lower_bound(units.begin(), units.end(), kid, by_kid) - units.begin();
			size_t end = lower_bound(units.begin(), units.end(), kid + 1, by_kid) - units.begin();
			Progress_reporter progress;

			progress.init(end - begin, 1);
			for (size_t u = begin, next; u < end; u = next) {
				vector<int> chroms;
				uint64_t rows = 0;

				for (next = u; next < end && units[next].batch == units[u].batch && (next == u || rows + units[next].rows <= MAX_ROWS_AT_ONCE); ++next) {
					chroms.push_back(units[next].chromid);
					rows += units[next].rows;
				}
				check_interrupt();
				import.import_chroms(chroms, batch_cols[units[u].batch], units[u].batch ? NULL : &reports[kid]);
				progress.report(next - u);
			}
			progress.report_last();
			return (char)0;
		});

		// Only now, when nothing can fail, is it safe to leave a pending message for .gcall()
		ClipReport clipped;
		clipped.clear();
		for (int kid = 0; kid < num_kids; ++kid)
			clipped.merge(reports[kid]);
		string clip_msg = import.clip_message(clipped);
		if (!clip_msg.empty())
			add_pending_diagnostic(_envir, "warning", clip_msg.c_str());
		report_skipped_chroms(srcname.c_str(), import.unknown(), num_matched, _envir);
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	} catch (const exception &e) {
		rerror("%s", e.what());
	}
	return R_NilValue;
}

}
