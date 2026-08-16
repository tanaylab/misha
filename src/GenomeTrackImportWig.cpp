#include <cstdint>
#include <cmath>
#include <cstring>

#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

#include "GenomeArraysCsv.h"
#include "GenomeTrackFixedBin.h"
#include "GenomeTrackSparse.h"
#include "Wig.h"

using namespace std;
using namespace rdb;

// Formats up to max_shown names as "a, b, c (and N more)".
static string format_chrom_names(const vector<string> &names, uint64_t total, size_t max_shown)
{
	string res;
	size_t shown = min(names.size(), max_shown);

	for (size_t i = 0; i < shown; ++i) {
		if (i)
			res += ", ";
		res += names[i];
	}

	if (total > shown) {
		char buf[100];
		snprintf(buf, sizeof(buf), "%s(and %llu more)", shown ? " " : "", (unsigned long long)(total - shown));
		res += buf;
	}
	return res;
}

static string format_db_chroms(const GenomeChromKey &chromkey, size_t max_shown)
{
	vector<string> names;
	size_t shown = min((size_t)chromkey.get_num_chroms(), max_shown);

	for (size_t i = 0; i < shown; ++i)
		names.push_back(chromkey.id2chrom((int)i));
	return format_chrom_names(names, chromkey.get_num_chroms(), max_shown);
}

static string format_alias_hint(SEXP envir)
{
	string hint = "Chromosome aliases are resolved automatically (chr1 <-> 1, M <-> MT, plus ";

	hint += get_groot(envir);
	hint += "/chrom_aliases.tsv if present). If the file's names should map to database chromosomes, "
	        "add them to that file (tab separated: canonical<TAB>alias<TAB>source) and reload the database "
	        "with gdb.init(rescan = TRUE), or rename the chromosomes in the file.";
	return hint;
}

// Decides what an import that could not map all of the file's chromosome names should do.
//
// The importer knows two things here: how many database chromosomes the file supplied data
// for, and which names in the file matched nothing. A file that covers only part of the
// database (a chr1-only bedGraph) has no unmatched names and is silent; a file whose names
// belong to another naming convention has unmatched names and either loses everything
// (error) or loses part of its data (warning).
static void verify_chroms_matched(const GenomeChromKey &chromkey, const char *fname,
                                  const vector<string> &unknown_chroms, uint64_t num_unknown,
                                  uint64_t num_matched, SEXP envir)
{
	if (!num_unknown)
		return;

	string unknown_str = format_chrom_names(unknown_chroms, num_unknown, 5);

	if (!num_matched) {
		verror("None of the chromosome names in %s exist in the genome database, so the track would be empty (all NaN).\n"
		       "File chromosomes: %s.\n"
		       "Database chromosomes: %s.\n"
		       "%s",
		       fname, unknown_str.c_str(), format_db_chroms(chromkey, 5).c_str(), format_alias_hint(envir).c_str());
	}

	char msg[10000];
	snprintf(msg, sizeof(msg),
	         "%llu chromosome name(s) in %s do not exist in the genome database and were skipped: %s. "
	         "Data for the remaining %llu chromosome(s) was imported. %s",
	         (unsigned long long)num_unknown, fname, unknown_str.c_str(),
	         (unsigned long long)num_matched, format_alias_hint(envir).c_str());

	// Rf_warning() is unsafe here: a caller that catches the warning with an exiting
	// handler (tryCatch, or options(warn = 2)) longjmps out of the C++ frame, skipping
	// ~RdbInitializer. The message is left for .gcall() to raise once the call returns.
	SEXP rmsg = rprotect_ptr(Rf_mkString(msg));
	define_in_misha(envir, ".GPENDING.WARNING", rmsg);
	runprotect(1);
}

extern "C" {

SEXP gtrackimportwig(SEXP _track, SEXP _wig, SEXP _binsize, SEXP _defvalue, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_track) || Rf_length(_track) != 1)
			verror("Track argument is not a string");

		if (!Rf_isString(_wig) || Rf_length(_wig) != 1)
			verror("Wig argument is not a string");

		if ((!Rf_isReal(_binsize) && !Rf_isInteger(_binsize)) || Rf_length(_binsize) != 1)
			verror("Binsize argument is not a number");

		if ((!Rf_isReal(_defvalue) && !Rf_isInteger(_defvalue)) || Rf_length(_defvalue) != 1)
			verror("Defvalue argument is not a number");

		const char *track = CHAR(STRING_ELT(_track, 0));
		const char *fname = CHAR(STRING_ELT(_wig, 0));
		double dbinsize = Rf_isReal(_binsize) ? REAL(_binsize)[0] : INTEGER(_binsize)[0];
		unsigned binsize = (unsigned)dbinsize;
		double defvalue = Rf_isReal(_defvalue) ? REAL(_defvalue)[0] : INTEGER(_defvalue)[0];

		if (dbinsize < 0 || binsize != dbinsize)
			verror("Invalid value of binsize argument: %g\n", dbinsize);

		string dirname = create_track_dir(_envir, track);
		IntervUtils iu(_envir);
		Wig wig;
		GenomeArraysCsv csv;
		GIntervals data;
		GIntervals::const_iterator iinterv_begin;
		GIntervals::const_iterator iinterv_end;
		vector<float> vals;
		char filename[FILENAME_MAX];
		bool is_csv = false;
		string wig_err_msg;

		try {
			wig.init(iu.get_chromkey(), fname, true);
		} catch (TGLException &e) {
			if (e.type() != typeid(Wig) || e.code() == Wig::FILE_ERROR)
				throw e;
			wig_err_msg = e.msg();
			is_csv = true;
		}

		if (is_csv) {
			try {
				csv.init(fname, iu.get_chromkey());
				if (csv.get_colnames().size() != 1)
					verror("More than one value column appears in file %s", fname);
			} catch (TGLException &e) {
				if (e.type() != typeid(GenomeArraysCsv) || e.code() == GenomeArraysCsv::FILE_ERROR)
					throw e;
				verror("Unrecognized format of file %s.\nWIG parser error: %s\nCSV parser error: %s",
				       fname, wig_err_msg.c_str(), e.msg());
			}

			verify_chroms_matched(iu.get_chromkey(), fname, csv.get_unknown_chroms(), csv.get_num_unknown_chroms(),
			                      csv.get_num_matched_chroms(), _envir);
		} else
			verify_chroms_matched(iu.get_chromkey(), fname, wig.get_unknown_chroms(), wig.get_num_unknown_chroms(),
			                      wig.get_num_matched_chroms(), _envir);

		GIntervals all_genome_intervs;
		iu.get_all_genome_intervs(all_genome_intervs);

		Progress_reporter progress;
		progress.init(iu.get_chromkey().get_num_chroms(), 1);

		for (int chromid = 0; chromid < (int)iu.get_chromkey().get_num_chroms(); ++chromid) {
			uint64_t chromsize = iu.get_chromkey().get_chrom_size(chromid);

			check_interrupt();

			if (is_csv) {
				const GIntervals &intervals = csv.get_intervals(chromid);
				iinterv_begin = intervals.begin();
				iinterv_end = intervals.end();
			} else {
				wig.get_data(chromid, data);
				iinterv_begin = data.begin();
				iinterv_end = data.end();
			}

			snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), chromid).c_str());

			if (binsize) {  // Fixed-bin track
				GenomeTrackFixedBin gtrack;
				gtrack.init_write(filename, binsize, chromid);
				GIntervals::const_iterator iinterval = iinterv_begin;

				// Batch writing: accumulate bins and write in batches
				vector<float> bin_batch;
				bin_batch.reserve(10000);
				uint64_t bins_since_progress = 0;

				for (uint64_t start_coord = 0; start_coord < chromsize; start_coord += binsize) {
					double sum = 0;
					uint64_t covered_bases = 0;
					uint64_t end_coord = min(start_coord + binsize, chromsize);

					// Skip intervals that end before current bin
					while (iinterval != iinterv_end && (uint64_t)iinterval->end <= start_coord) {
						++iinterval;
					}

					// Process all intervals that overlap with current bin
					GIntervals::const_iterator cur_interval = iinterval;
					while (cur_interval != iinterv_end && (uint64_t)cur_interval->start < end_coord) {
						// Calculate overlap between interval and bin
						uint64_t overlap_start = max(start_coord, (uint64_t)cur_interval->start);
						uint64_t overlap_end = min(end_coord, (uint64_t)cur_interval->end);

						if (overlap_end > overlap_start) {
							float v;
							if (is_csv) {
								csv.get_sliced_vals(cur_interval, vals);
								v = vals.front();
							} else {
								memcpy(&v, &cur_interval->udata, sizeof(float));
							}

							if (!std::isnan(v)) {
								uint64_t overlap_len = overlap_end - overlap_start;
								sum += v * overlap_len;
								covered_bases += overlap_len;
							}
						}
						++cur_interval;
					}

					// Fill uncovered bases with default value
					uint64_t bin_size_actual = end_coord - start_coord;
					if (covered_bases < bin_size_actual && !std::isnan(defvalue)) {
						uint64_t uncovered = bin_size_actual - covered_bases;
						sum += defvalue * uncovered;
						covered_bases += uncovered;
					}

					// Calculate bin value
					float bin_value = covered_bases ? (sum / covered_bases) : numeric_limits<float>::quiet_NaN();
					bin_batch.push_back(bin_value);

					// Write in batches of 10000 bins
					if (bin_batch.size() >= 10000) {
						gtrack.write_next_bins(bin_batch.data(), bin_batch.size());
						bin_batch.clear();
						progress.report(0);
						check_interrupt();
						bins_since_progress = 0;
					} else {
						++bins_since_progress;
						// Report progress every 10000 bins even if batch not full
						if (bins_since_progress >= 10000) {
							progress.report(0);
							check_interrupt();
							bins_since_progress = 0;
						}
					}
				}

				// Write remaining bins
				if (!bin_batch.empty()) {
					gtrack.write_next_bins(bin_batch.data(), bin_batch.size());
				}
			} else { // Sparse track
				GenomeTrackSparse gtrack;
				gtrack.init_write(filename, chromid);

				uint64_t intervals_processed = 0;
				for (GIntervals::const_iterator iinterval = iinterv_begin; iinterval != iinterv_end; ++iinterval) {
					float v;

					if (is_csv) {
						csv.get_sliced_vals(iinterval, vals);
						v = vals.front();
					} else {
						memcpy(&v, &iinterval->udata, sizeof(float));
					}

					gtrack.write_next_interval(*iinterval, v);

					// Report progress every 10000 intervals
					++intervals_processed;
					if (intervals_processed % 10000 == 0) {
						progress.report(0);
						check_interrupt();
					}
				}
			}
			progress.report(1);
		}
		progress.report_last();
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

	return R_NilValue;
}

}
