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
// on each chromosome (reading only the chrom columns); import_chrom() then writes the track
// files of one chromosome for any subset of the columns.
class ParquetImport {
public:
	enum { CHROM, START, END, NUM_COORDS };

	ParquetImport(const GenomeChromKey &chromkey, const vector<string> &files, const vector<string> &columns,
	              const vector<string> &dirs, unsigned binsize, double defval) :
		m_chromkey(chromkey), m_files(files), m_columns(columns), m_dirs(dirs), m_binsize(binsize), m_defval(defval) {}

	// Opens the files and finds the columns in each
	void open();

	// Fills runs, unknown chromosome names and the files they came from
	void index();

	void import_chrom(int chromid, const vector<int> &cols);

	const UnknownChroms &unknown() const { return m_unknown; }
	uint64_t num_matched_chroms() const;

	// The files that named chromosomes the database does not have, for messages
	string unknown_chroms_files() const;

private:
	const GenomeChromKey             &m_chromkey;
	vector<string>                    m_files;
	vector<string>                    m_columns;
	vector<string>                    m_dirs;
	unsigned                          m_binsize;
	double                            m_defval;
	vector<unique_ptr<ParquetReader>> m_readers;
	vector<vector<int>>               m_col_idx;  // per file: chrom, start, end, then the value columns
	vector<vector<Run>>               m_runs;     // per chromid
	UnknownChroms                     m_unknown;
	vector<bool>                      m_file_has_unknown;

	[[noreturn]] void rows_error(const char *what, const Run &run1, const Run &run2) const;
};

void ParquetImport::open()
{
	m_readers.clear();
	m_col_idx.clear();
	for (const string &file : m_files) {
		m_readers.emplace_back(new ParquetReader());
		ParquetReader &reader = *m_readers.back();
		reader.open(file);

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
}

void ParquetImport::index()
{
	vector<string> chroms;

	m_runs.assign(m_chromkey.get_num_chroms(), vector<Run>());
	m_unknown.clear();
	m_file_has_unknown.assign(m_files.size(), false);

	for (size_t f = 0; f < m_files.size(); ++f) {
		for (int rg = 0; rg < m_readers[f]->num_row_groups(); ++rg) {
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
			res += (num ? ", " : "") + m_files[f];
		++num;
	}
	if (num > (size_t)UnknownChroms::MAX_REPORTED)
		res += " (and " + to_string(num - UnknownChroms::MAX_REPORTED) + " more files)";
	return res;
}

void ParquetImport::rows_error(const char *what, const Run &run1, const Run &run2) const
{
	if (run1.file == run2.file && run1.row_group == run2.row_group)
		verror("File %s, row group %d: %s", m_files[run1.file].c_str(), run1.row_group, what);
	verror("File %s, row group %d and file %s, row group %d: %s", m_files[run1.file].c_str(), run1.row_group,
	       m_files[run2.file].c_str(), run2.row_group, what);
}

void ParquetImport::import_chrom(int chromid, const vector<int> &cols)
{
	const vector<Run> &runs = m_runs[chromid];
	uint64_t chromsize = m_chromkey.get_chrom_size(chromid);
	uint64_t num_bins = (chromsize + m_binsize - 1) / m_binsize;
	vector<int64_t> starts, ends;  // the chromosome's rows in the order of runs
	vector<size_t> run_of_row;     // index into runs of the first row of each run, for messages
	vector<int64_t> buf_starts, buf_ends;
	bool aligned = true;           // every row is exactly a bin

	// Coordinates. Runs of one row group are consecutive in runs, so each row group's start and
	// end columns are decoded once.
	for (size_t r = 0; r < runs.size(); ++r) {
		const Run &run = runs[r];
		ParquetReader &reader = *m_readers[run.file];

		if (!r || run.file != runs[r - 1].file || run.row_group != runs[r - 1].row_group) {
			reader.read_int64s(run.row_group, m_col_idx[run.file][START], buf_starts);
			reader.read_int64s(run.row_group, m_col_idx[run.file][END], buf_ends);
		}

		run_of_row.push_back(starts.size());
		for (int64_t i = run.begin; i < run.end; ++i) {
			int64_t start = buf_starts[i];
			int64_t end = buf_ends[i];

			if (start < 0)
				rows_error("start coordinate is negative", run, run);
			if (start >= end)
				rows_error("start coordinate exceeds or equals the end coordinate", run, run);
			if ((uint64_t)end > chromsize)
				rows_error("end coordinate exceeds chromosome's size", run, run);
			aligned = aligned && start % m_binsize == 0 && (uint64_t)end == min((uint64_t)start + m_binsize, chromsize);
			starts.push_back(start);
			ends.push_back(end);
		}
	}

	auto run_of = [&](size_t row) -> const Run & {
		return runs[upper_bound(run_of_row.begin(), run_of_row.end(), row) - run_of_row.begin() - 1];
	};

	// Rows that are exactly the bins (the common case) are placed by their start; any other rows
	// are sorted by start, as GenomeArraysCsv sorts them. Overlapping or duplicate rows are errors
	// either way.
	vector<uint32_t> order;
	if (aligned) {
		vector<uint32_t> row_of_bin(num_bins, UINT32_MAX);
		for (size_t i = 0; i < starts.size(); ++i) {
			uint32_t &row = row_of_bin[starts[i] / m_binsize];
			if (row != UINT32_MAX)
				rows_error("rows overlap", run_of(row), run_of(i));
			row = i;
		}
	} else {
		order.resize(starts.size());
		iota(order.begin(), order.end(), 0);
		sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) { return starts[a] < starts[b] || (starts[a] == starts[b] && ends[a] < ends[b]); });
		for (size_t k = 1; k < order.size(); ++k) {
			if (ends[order[k - 1]] > starts[order[k]])
				rows_error("rows overlap", run_of(order[k - 1]), run_of(order[k]));
		}
	}

	vector<float> vals(starts.size());
	vector<float> bins(num_bins);
	vector<double> buf_vals;
	char filename[FILENAME_MAX];

	for (int col : cols) {
		// the column's values, in the order of runs
		size_t row = 0;
		for (size_t r = 0; r < runs.size(); ++r) {
			const Run &run = runs[r];
			if (!r || run.file != runs[r - 1].file || run.row_group != runs[r - 1].row_group)
				m_readers[run.file]->read_doubles(run.row_group, m_col_idx[run.file][NUM_COORDS + col], buf_vals);
			for (int64_t i = run.begin; i < run.end; ++i)
				vals[row++] = buf_vals[i];
		}

		if (aligned) {
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
		check_interrupt();
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
// double, string, other
SEXP parquet_columns(SEXP _file, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_file) || Rf_length(_file) != 1)
			verror("File argument is not a string");

		ParquetReader reader;
		reader.open(CHAR(STRING_ELT(_file, 0)));

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
// for every chromosome of the database
SEXP gtrackimportparquet(SEXP _dirs, SEXP _files, SEXP _columns, SEXP _binsize, SEXP _defval, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_dirs) || !Rf_isString(_columns) || Rf_length(_dirs) != Rf_length(_columns) || !Rf_length(_columns))
			verror("dirs and columns must be character vectors of the same, positive, length");
		if (!Rf_isString(_files) || !Rf_length(_files))
			verror("files must be a non-empty character vector");
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

		vector<string> dirs, files, columns;
		for (int i = 0; i < Rf_length(_dirs); ++i) {
			dirs.push_back(CHAR(STRING_ELT(_dirs, i)));
			columns.push_back(CHAR(STRING_ELT(_columns, i)));
		}
		for (int i = 0; i < Rf_length(_files); ++i)
			files.push_back(CHAR(STRING_ELT(_files, i)));

		IntervUtils iu(_envir);
		const GenomeChromKey &chromkey = iu.get_chromkey();
		ParquetImport import(chromkey, files, columns, dirs, binsize, defval);

		import.open();
		import.index();

		string srcname = import.unknown_chroms_files();
		uint64_t num_matched = import.num_matched_chroms();
		verify_chroms_matched(chromkey, srcname.c_str(), import.unknown(), num_matched, _envir);

		// Units: each chromosome, largest first, times batches of columns of nearly equal size.
		// Kid i takes units i, i + num_kids, ... A batch decodes the chromosome's coordinates
		// again, so in one process all the columns are one batch.
		int num_cols = columns.size();
		int num_batches = (num_cols + MAX_UNIT_COLS - 1) / MAX_UNIT_COLS;
		int num_kids = import_num_kids(iu, chromkey.get_num_chroms() * num_batches);
		if (num_kids == 1)
			num_batches = 1;
		vector<int> chromids(chromkey.get_num_chroms());
		vector<pair<int, vector<int>>> units;

		iota(chromids.begin(), chromids.end(), 0);
		stable_sort(chromids.begin(), chromids.end(), [&](int a, int b) { return chromkey.get_chrom_size(a) > chromkey.get_chrom_size(b); });
		for (int chromid : chromids) {
			for (int batch = 0; batch < num_batches; ++batch) {
				vector<int> cols(num_cols * (batch + 1) / num_batches - num_cols * batch / num_batches);
				iota(cols.begin(), cols.end(), num_cols * batch / num_batches);
				units.emplace_back(chromid, cols);
			}
		}

		run_kids(iu, num_kids, [&](int kid) {
			// a kid reads through file descriptors of its own
			if (num_kids > 1)
				import.open();

			Progress_reporter progress;
			progress.init((units.size() - kid + num_kids - 1) / num_kids, 1);
			for (size_t u = kid; u < units.size(); u += num_kids) {
				check_interrupt();
				import.import_chrom(units[u].first, units[u].second);
				progress.report(1);
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
