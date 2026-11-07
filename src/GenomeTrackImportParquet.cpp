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
#include <memory>
#include <numeric>
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
	// their rows is decoded once, however many of the chromosomes it holds.
	void import_chroms(const vector<int> &chromids, const vector<int> &cols);

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
		vector<int64_t>  ends;
		vector<float>    vals;     // of the column being written
		vector<uint32_t> order;    // the rows sorted by start, unless aligned
		bool             aligned;  // every row is exactly a bin
	};

	// A run of chromosome chromids[chrom] of import_chroms()
	struct ChromRun {
		size_t     chrom;
		const Run *run;
	};

	[[noreturn]] void rows_error(const char *what, const Run &run1, const Run &run2) const;

	template <typename T, typename D>
	void gather(const vector<ChromRun> &runs, int col_idx, void (ParquetReader::*read)(int, int, vector<T> &), D dest);
	void prepare_rows(int chromid, ChromRows &rows) const;
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

// Checks the coordinates of the rows, sets aligned and, unless they are aligned, sorts them
void ParquetImport::prepare_rows(int chromid, ChromRows &rows) const
{
	const vector<Run> &runs = m_runs[chromid];
	uint64_t chromsize = m_chromkey.get_chrom_size(chromid);
	uint64_t num_bins = (chromsize + m_binsize - 1) / m_binsize;
	const vector<int64_t> &starts = rows.starts;
	const vector<int64_t> &ends = rows.ends;
	vector<size_t> run_of_row;  // the first row of each run, for messages

	rows.aligned = true;
	for (size_t r = 0, i = 0; r < runs.size(); ++r) {
		const Run &run = runs[r];

		run_of_row.push_back(i);
		for (int64_t k = run.begin; k < run.end; ++k, ++i) {
			int64_t start = starts[i];
			int64_t end = ends[i];

			if (start < 0)
				rows_error("start coordinate is negative", run, run);
			if (start >= end)
				rows_error("start coordinate exceeds or equals the end coordinate", run, run);
			if ((uint64_t)end > chromsize)
				rows_error("end coordinate exceeds chromosome's size", run, run);
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
			uint32_t &row = row_of_bin[starts[i] / m_binsize];
			if (row != UINT32_MAX)
				rows_error("rows overlap", run_of(row), run_of(i));
			row = i;
		}
	} else {
		vector<uint32_t> &order = rows.order;
		order.resize(starts.size());
		iota(order.begin(), order.end(), 0);
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
			uint64_t len = ends[i] - starts[i];
			double sum = 0;
			uint64_t covered_bases = 0;
			float v = vals[i];

			if (!std::isnan(v)) {
				sum += v * len;  // a float product, as in gtrackimportwig
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
						sum += v * overlap_len;
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

void ParquetImport::import_chroms(const vector<int> &chromids, const vector<int> &cols)
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
		prepare_rows(chromids[c], rows[c]);

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
				import.import_chroms(chroms, batch_cols[units[u].batch]);
				progress.report(next - u);
			}
			progress.report_last();
			return (char)0;
		});

		// Only now, when nothing can fail, is it safe to leave a pending message for .gcall()
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
