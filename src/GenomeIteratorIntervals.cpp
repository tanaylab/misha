#include "GIntervalsBigSet1D.h"
#include "GIntervalsBigSet2D.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"
#include "TrackExpressionScanner.h"
#include "TrackExpressionIntervalRelativeBinIterator.h"

using namespace std;
using namespace rdb;

extern "C" {

SEXP giterator_intervals(SEXP _expr, SEXP _intervals, SEXP _iterator_policy, SEXP _band, SEXP _intervals_set_out,
                         SEXP _interval_relative, SEXP _partial_bins, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_expr) || Rf_length(_expr) != 1)
			verror("Tracks expression argument must be a string");

		if (!Rf_isNull(_intervals_set_out) && (!Rf_isString(_intervals_set_out) || Rf_length(_intervals_set_out) != 1))
			verror("intervals.set.out argument is not a string");

		// Parse interval_relative parameter
		bool interval_relative = false;
		if (!Rf_isNull(_interval_relative)) {
			if (!Rf_isLogical(_interval_relative) || Rf_length(_interval_relative) != 1)
				verror("interval_relative argument must be a logical value");
			interval_relative = LOGICAL(_interval_relative)[0];
		}

		// Parse partial_bins parameter (0 = clip, 1 = exact/drop)
		TrackExpressionIntervalRelativeBinIterator::PartialBinsMode partial_bins_mode =
			TrackExpressionIntervalRelativeBinIterator::CLIP;
		if (!Rf_isNull(_partial_bins)) {
			if (!Rf_isInteger(_partial_bins) || Rf_length(_partial_bins) != 1)
				verror("partial_bins argument must be an integer");
			partial_bins_mode = (TrackExpressionIntervalRelativeBinIterator::PartialBinsMode)INTEGER(_partial_bins)[0];
		}

		string intervset_out = Rf_isNull(_intervals_set_out) ? "" : CHAR(STRING_ELT(_intervals_set_out, 0));

		IntervUtils iu(_envir);
		SEXP answer = R_NilValue;

		// Handle interval_relative mode separately
		if (interval_relative) {
			// Validate: interval_relative requires a numeric iterator (binsize)
			if (!((Rf_isReal(_iterator_policy) || Rf_isInteger(_iterator_policy)) && Rf_length(_iterator_policy) == 1))
				verror("interval_relative mode requires a numeric iterator (binsize)");

			int64_t binsize = Rf_isReal(_iterator_policy) ? (int64_t)REAL(_iterator_policy)[0] : INTEGER(_iterator_policy)[0];

			// Convert intervals to GIntervals (not using scope/unify for source intervals)
			GIntervals src_intervals;
			iu.convert_rintervs(_intervals, &src_intervals, NULL, false);
			src_intervals.sort();

			// Create and run the interval-relative iterator
			TrackExpressionIntervalRelativeBinIterator iterator;
			iterator.begin(src_intervals, binsize, partial_bins_mode, src_intervals);

			if (intervset_out.empty()) {
				// Collect results in memory with intervalID
				GIntervals res_intervs;
				vector<int64_t> interval_ids;

				for (; !iterator.isend(); iterator.next()) {
					res_intervs.push_back(iterator.last_interval());
					interval_ids.push_back(iterator.get_cur_src_interval_idx() + 1);  // 1-based for R
					iu.verify_max_data_size(res_intervs.size(), "Result");
				}

				// Convert to R data frame with extra column for intervalID
				answer = iu.convert_intervs(&res_intervs, GInterval::NUM_COLS + 1, false);

				// Add intervalID column
				SEXP ids = rprotect_ptr(RSaneAllocVector(INTSXP, interval_ids.size()));
				for (size_t i = 0; i < interval_ids.size(); ++i)
					INTEGER(ids)[i] = interval_ids[i];
				SET_VECTOR_ELT(answer, GInterval::NUM_COLS, ids);

				// Set column name
				SEXP col_names = Rf_getAttrib(answer, R_NamesSymbol);
				SET_STRING_ELT(col_names, GInterval::NUM_COLS, Rf_mkChar("intervalID"));
			} else {
				// Save to big intervals set (without intervalID for file output)
				GIntervals res_intervals;
				vector<GIntervalsBigSet1D::ChromStat> chromstats;
				char error_prefix[1000];

				GIntervalsBigSet1D::begin_save(intervset_out.c_str(), iu, chromstats);
				for (; !iterator.isend(); iterator.next()) {
					const GInterval &interval = iterator.last_interval();

					if (res_intervals.empty())
						snprintf(error_prefix, sizeof(error_prefix), "Big intervals set %s, chrom %s",
								 intervset_out.c_str(), iu.id2chrom(interval.chromid).c_str());
					else if (res_intervals.front().chromid != interval.chromid)
						GIntervalsBigSet1D::save_chrom_plain_intervals(intervset_out.c_str(), res_intervals, iu, chromstats);

					res_intervals.push_back(interval);
					iu.verify_max_data_size(res_intervals.size(), error_prefix);
				}

				GIntervalsBigSet1D::save_chrom_plain_intervals(intervset_out.c_str(), res_intervals, iu, chromstats);
				GIntervalsBigSet1D::end_save_plain_intervals(intervset_out.c_str(), iu, chromstats);
			}

			return answer;
		}

		// Original behavior: use TrackExprScanner
		GIntervalsFetcher1D *intervals1d = NULL;
		GIntervalsFetcher2D *intervals2d = NULL;
		iu.convert_rintervs(_intervals, &intervals1d, &intervals2d);
		unique_ptr<GIntervalsFetcher1D> intervals1d_guard(intervals1d);
		unique_ptr<GIntervalsFetcher2D> intervals2d_guard(intervals2d);
		intervals1d->sort();
		intervals1d->unify_overlaps();
		intervals2d->sort();
		intervals2d->verify_no_overlaps(iu.get_chromkey());

		TrackExprScanner scanner(iu);
		scanner.begin(_expr, intervals1d, intervals2d, _iterator_policy, _band);

		if (intervset_out.empty()) {
			if (scanner.get_iterator()->is_1d()) {
				GIntervals res_intervs;
				for (; !scanner.isend(); scanner.next()) {
					res_intervs.push_back(scanner.last_interval1d());
					iu.verify_max_data_size(res_intervs.size(), "Result");
				}
				answer = iu.convert_intervs(&res_intervs);
			} else {
				GIntervals2D res_intervs;
				for (; !scanner.isend(); scanner.next()) {
					res_intervs.push_back(scanner.last_interval2d());
					iu.verify_max_data_size(res_intervs.size(), "Result");
				}
				answer = iu.convert_intervs(&res_intervs);
			}
		} else {
			char error_prefix[1000];

			if (scanner.get_iterator()->is_1d()) {
				GIntervals res_intervals;
				vector<GIntervalsBigSet1D::ChromStat> chromstats;

				GIntervalsBigSet1D::begin_save(intervset_out.c_str(), iu, chromstats);
				for (; !scanner.isend(); scanner.next()) {
					const GInterval &interval = scanner.last_interval1d();

					if (res_intervals.empty())
						snprintf(error_prefix, sizeof(error_prefix), "Big intervals set %s, chrom %s",
								 intervset_out.c_str(), iu.id2chrom(interval.chromid).c_str());
					else if (res_intervals.front().chromid != interval.chromid)
						GIntervalsBigSet1D::save_chrom_plain_intervals(intervset_out.c_str(), res_intervals, iu, chromstats);

					res_intervals.push_back(interval);
					iu.verify_max_data_size(res_intervals.size(), error_prefix);
				}

				GIntervalsBigSet1D::save_chrom_plain_intervals(intervset_out.c_str(), res_intervals, iu, chromstats);
				GIntervalsBigSet1D::end_save_plain_intervals(intervset_out.c_str(), iu, chromstats);
			} else {
				GIntervals2D res_intervals;
				vector<GIntervalsBigSet2D::ChromStat> chromstats;

				GIntervalsBigSet2D::begin_save(intervset_out.c_str(), iu, chromstats);
				for (; !scanner.isend(); scanner.next()) {
					const GInterval2D &interval = scanner.last_interval2d();

					if (res_intervals.empty())
						snprintf(error_prefix, sizeof(error_prefix), "Big intervals set %s, chroms (%s, %s)",
								intervset_out.c_str(), iu.id2chrom(interval.chromid1()).c_str(), iu.id2chrom(interval.chromid2()).c_str());
					else if (!res_intervals.front().is_same_chrom(interval))
						GIntervalsBigSet2D::save_chrom_plain_intervals(intervset_out.c_str(), res_intervals, iu, chromstats);

					res_intervals.push_back(interval);
					iu.verify_max_data_size(res_intervals.size(), error_prefix);
				}

				GIntervalsBigSet2D::save_chrom_plain_intervals(intervset_out.c_str(), res_intervals, iu, chromstats);
				GIntervalsBigSet2D::end_save_plain_intervals(intervset_out.c_str(), iu, chromstats);
			}
		}

		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	return R_NilValue;
}

SEXP gcheck_iterator(SEXP _iterator_policy, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		IntervUtils iu(_envir);
		GIntervals intervals1d;
		GIntervals2D intervals2d;
		GIntervals scope1d;
		GIntervals2D scope2d;
		iu.get_all_genome_intervs(scope1d);
		iu.get_all_genome_intervs(scope2d);

		TrackExprScanner scanner(iu);
		scanner.create_expr_iterator(R_NilValue, &scope1d, &scope2d, _iterator_policy, R_NilValue);
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	return R_NilValue;
}

}
