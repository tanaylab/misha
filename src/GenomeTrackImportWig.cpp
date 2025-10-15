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

		try {
			wig.init(iu.get_chromkey(), fname, true);
		} catch (TGLException &e) {
			if (e.type() != typeid(Wig) || e.code() == Wig::FILE_ERROR)
				throw e;
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
				verror("Unrecognized format of file %s", fname);
			}
		}

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
