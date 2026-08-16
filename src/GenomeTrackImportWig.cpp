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

// Formats names as "a, b, c (and N more)". "N+" when the count itself was capped.
static string format_chrom_names(const vector<string> &names, uint64_t total, bool truncated)
{
	string res;

	for (size_t i = 0; i < names.size(); ++i) {
		if (i)
			res += ", ";
		res += names[i];
	}

	if (total > names.size()) {
		char buf[100];
		snprintf(buf, sizeof(buf), "%s(and %llu%s more)", names.empty() ? "" : " ",
		         (unsigned long long)(total - names.size()), truncated ? "+" : "");
		res += buf;
	}
	return res;
}

static string format_db_chroms(const GenomeChromKey &chromkey)
{
	vector<string> names;
	size_t shown = min((size_t)chromkey.get_num_chroms(), (size_t)UnknownChroms::MAX_REPORTED);

	for (size_t i = 0; i < shown; ++i)
		names.push_back(chromkey.id2chrom((int)i));
	return format_chrom_names(names, chromkey.get_num_chroms(), false);
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

// A file that named chromosomes and matched none of them contributes nothing to the track:
// whatever the track would end up holding (all NaN, or defval everywhere) comes from nowhere
// in the file. Refuse before writing anything.
//
// A file that matched some of them is a different matter, and is left to
// report_skipped_chroms below. The distinction the importer can actually make here is
// between names in the file that resolved and names that did not - not between database
// chromosomes the file covered and ones it did not. A chr1-only bedGraph dropped nothing,
// so it is silent; it is unmatched names that prove data in the file never reached the track.
static void verify_chroms_matched(const GenomeChromKey &chromkey, const char *fname,
                                  const UnknownChroms &unknown, uint64_t num_matched, SEXP envir)
{
	if (unknown.empty() || num_matched)
		return;

	verror("None of the chromosome names in %s exist in the genome database: no data from the file would reach the track.\n"
	       "File chromosomes: %s.\n"
	       "Database chromosomes: %s.\n"
	       "%s",
	       fname, format_chrom_names(unknown.names(), unknown.num(), unknown.truncated()).c_str(),
	       format_db_chroms(chromkey).c_str(), format_alias_hint(envir).c_str());
}

// Reports the names dropped by an import that did produce data. Which channel it goes to
// depends on what was dropped, because the two cases are not equally alarming:
//
//   - a scaffold / patch / unplaced contig (chrUn_*, *_random, GL000220.1): a database built
//     from the main chromosomes deliberately does not have it, no alias can fix it, and this
//     is what every whole-genome bigWig looks like. A message: visible and in order, and it
//     cannot be aggregated away into "there were 50 or more warnings".
//   - a primary chromosome (chr7, X): the database was supposed to have that one, so the
//     naming is probably wrong. A warning, with the alias hint. Mitochondria are the first
//     kind, not this one - see UnknownChroms.
static void report_skipped_chroms(const char *fname, const UnknownChroms &unknown, uint64_t num_matched, SEXP envir)
{
	if (unknown.empty())
		return;

	char msg[10000];
	const char *slot;

	if (unknown.primary_names().empty()) {
		snprintf(msg, sizeof(msg),
		         "%llu%s chromosome name(s) in %s do not exist in the genome database and were skipped: %s. "
		         "Data for the remaining %llu chromosome(s) was imported.",
		         (unsigned long long)unknown.num(), unknown.truncated() ? "+" : "", fname,
		         format_chrom_names(unknown.names(), unknown.num(), unknown.truncated()).c_str(),
		         (unsigned long long)num_matched);
		slot = ".GPENDING.MESSAGE";
	} else {
		snprintf(msg, sizeof(msg),
		         "%llu%s chromosome name(s) in %s do not exist in the genome database and were skipped, among them "
		         "primary chromosome(s): %s. Data for the remaining %llu chromosome(s) was imported. %s",
		         (unsigned long long)unknown.num(), unknown.truncated() ? "+" : "", fname,
		         format_chrom_names(unknown.primary_names(), unknown.num_primary(), false).c_str(),
		         (unsigned long long)num_matched, format_alias_hint(envir).c_str());
		slot = ".GPENDING.WARNING";
	}

	// Neither Rf_warning() nor an R-level message() can be raised from here: a caller that
	// catches either with an exiting handler (tryCatch, or options(warn = 2)) longjmps out of
	// the C++ frame, skipping ~RdbInitializer. The text is left for .gcall() to raise once the
	// call has returned - and it is only set once the import can no longer fail, so an error
	// downstream cannot leave it behind as permanent state in .misha.
	SEXP rmsg = rprotect_ptr(Rf_mkString(msg));
	define_in_misha(envir, slot, rmsg);
	runprotect(1);
}

extern "C" {

SEXP gtrackimportwig(SEXP _track, SEXP _wig, SEXP _srcname, SEXP _binsize, SEXP _defvalue, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_track) || Rf_length(_track) != 1)
			verror("Track argument is not a string");

		if (!Rf_isString(_wig) || Rf_length(_wig) != 1)
			verror("Wig argument is not a string");

		if (!Rf_isString(_srcname) || Rf_length(_srcname) != 1)
			verror("Source name argument is not a string");

		if ((!Rf_isReal(_binsize) && !Rf_isInteger(_binsize)) || Rf_length(_binsize) != 1)
			verror("Binsize argument is not a number");

		if ((!Rf_isReal(_defvalue) && !Rf_isInteger(_defvalue)) || Rf_length(_defvalue) != 1)
			verror("Defvalue argument is not a number");

		const char *track = CHAR(STRING_ELT(_track, 0));
		const char *fname = CHAR(STRING_ELT(_wig, 0));
		// What the user asked to import. gtrack.import unzips and converts bigWig into a
		// temporary file that is deleted before the call returns, so messages that name
		// fname would point at a path that no longer exists (and, in gtrack.import_set,
		// disagree with the file name reported as failed).
		const char *srcname = CHAR(STRING_ELT(_srcname, 0));
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

		}

		const UnknownChroms &unknown_chroms = is_csv ? csv.get_unknown_chroms() : wig.get_unknown_chroms();
		uint64_t num_matched_chroms = is_csv ? csv.get_num_matched_chroms() : wig.get_num_matched_chroms();

		verify_chroms_matched(iu.get_chromkey(), srcname, unknown_chroms, num_matched_chroms, _envir);

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

		// Only now, when nothing downstream can throw, is it safe to leave a pending
		// message behind for .gcall() to raise.
		report_skipped_chroms(srcname, unknown_chroms, num_matched_chroms, _envir);
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

	return R_NilValue;
}

}
