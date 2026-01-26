/*
 * GenomeTrackExtract.cpp
 *
 *  Created on: Mar 21, 2010
 *      Author: hoichman
 */

#include <cstdint>
#include <inttypes.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "rdbinterval.h"
#include "rdbutils.h"
#include "GenomeTrack.h"
#include "GIntervalsBigSet1D.h"
#include "GIntervalsBigSet2D.h"
#include "TrackExpressionScanner.h"
#include "TrackExpressionFixedBinIterator.h"
#include "TrackExpressionFixedRectIterator.h"
#include "TrackExpressionIntervals1DIterator.h"
#include "TGLException.h"

// Optional profiling for gextract_multitask hot path.
// Enable by setting getOption("gextract.profile") to TRUE.
static bool is_gextract_profile_enabled() {
	SEXP opt = Rf_GetOption1(Rf_install("gextract.profile"));
	return !Rf_isNull(opt) && Rf_asLogical(opt) == 1;
}

static void log_gextract_timing(const char *label, double ms) {
	Rprintf("[gextract.profile] %s: %.2f ms\n", label, ms);
}

using namespace std;
using namespace rdb;

static uint64_t estimate_records_for_expr(
	IntervUtils &iu,
	SEXP exprs,
	GIntervalsFetcher1D *scope1d,
	GIntervalsFetcher2D *scope2d,
	SEXP iterator_policy,
	SEXP band)
{
	uint64_t estimated = iu.estimate_num_bins(iterator_policy, scope1d, scope2d);
	if (estimated)
		return estimated;

	// Try to infer iterator details when iterator policy is implicit.
	TrackExprScanner scanner(iu);
	TrackExpressionIteratorBase *expr_itr = scanner.create_expr_iterator(exprs, scope1d, scope2d, iterator_policy, band, true);

	if (auto *bin_itr = dynamic_cast<TrackExpressionFixedBinIterator *>(expr_itr)) {
		int64_t binsize = bin_itr->get_bin_size();
		if (binsize > 0) {
			SEXP rbinsize = PROTECT(Rf_ScalarReal((double)binsize));
			estimated = iu.estimate_num_bins(rbinsize, scope1d, scope2d);
			UNPROTECT(1);
			return estimated;
		}
	}

	if (auto *rect_itr = dynamic_cast<TrackExpressionFixedRectIterator *>(expr_itr)) {
		int64_t width = rect_itr->get_width();
		int64_t height = rect_itr->get_height();
		if (width > 0 && height > 0) {
			SEXP rrect = PROTECT(Rf_allocVector(REALSXP, 2));
			REAL(rrect)[0] = (double)width;
			REAL(rrect)[1] = (double)height;
			estimated = iu.estimate_num_bins(rrect, scope1d, scope2d);
			UNPROTECT(1);
			return estimated;
		}
	}

	if (auto *intervals_itr = dynamic_cast<TrackExpressionIntervals1DIterator *>(expr_itr)) {
		const GIntervals *intervals = intervals_itr->get_intervals();
		if (intervals)
			return intervals->size();
	}

	return 0;
}

static SEXP build_rintervals_extract(GIntervalsFetcher1D *out_intervals1d, GIntervalsFetcher2D *out_intervals2d, const vector< vector<double> > &values,
									 vector<unsigned> *interv_ids, SEXP _exprs, SEXP _colnames, IntervUtils &iu)
{
	SEXP answer;
	unsigned num_interv_cols;
	unsigned num_exprs = values.size();

	if (out_intervals1d) {
		answer = iu.convert_intervs(out_intervals1d, interv_ids ? GInterval::NUM_COLS + num_exprs + 1 : GInterval::NUM_COLS + num_exprs, false);
		num_interv_cols = GInterval::NUM_COLS;
	} else {
		answer = iu.convert_intervs(out_intervals2d, interv_ids ? GInterval2D::NUM_COLS + num_exprs + 1 : GInterval2D::NUM_COLS + num_exprs, false);
		num_interv_cols = GInterval2D::NUM_COLS;
	}

    for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
        SEXP expr_vals;
        expr_vals = rprotect_ptr(RSaneAllocVector(REALSXP, values[iexpr].size()));
		// Use memcpy for bulk copy instead of manual loop for better performance
		if (values[iexpr].size() > 0) {
			memcpy(REAL(expr_vals), values[iexpr].data(), values[iexpr].size() * sizeof(double));
		}
        SET_VECTOR_ELT(answer, num_interv_cols + iexpr, expr_vals);
	}

    SEXP col_names = rprotect_ptr(Rf_getAttrib(answer, R_NamesSymbol));
	for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
		if (Rf_isNull(_colnames))
			SET_STRING_ELT(col_names, num_interv_cols + iexpr, Rf_mkChar(get_bounded_colname(CHAR(STRING_ELT(_exprs, iexpr))).c_str()));
		else
			SET_STRING_ELT(col_names, num_interv_cols + iexpr, STRING_ELT(_colnames, iexpr));
	}

    if (interv_ids) {
        SEXP ids;
        ids = rprotect_ptr(RSaneAllocVector(INTSXP, interv_ids->size()));
		for (vector<unsigned>::const_iterator iid = interv_ids->begin(); iid != interv_ids->end(); ++iid)
			INTEGER(ids)[iid - interv_ids->begin()] = *iid;
		SET_VECTOR_ELT(answer, num_interv_cols + num_exprs, ids);

		SET_STRING_ELT(col_names, num_interv_cols + num_exprs, Rf_mkChar("intervalID"));
	}
    
    runprotect(1); // col_names
    return answer;
}


extern "C" {

SEXP C_gextract(SEXP _intervals, SEXP _exprs, SEXP _colnames, SEXP _iterator_policy, SEXP _band, SEXP _file, SEXP _intervals_set_out, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_exprs) || Rf_length(_exprs) < 1)
			verror("Track expressions argument must be a vector of strings");

		if (!Rf_isNull(_colnames)) {
			if (!Rf_isString(_colnames))
				verror("Column names argument must be a vector of strings");
			if (Rf_length(_colnames) != Rf_length(_exprs))
				verror("Number of column names must match the number of track expressions");
		}

		if (!Rf_isNull(_file) && (!Rf_isString(_file) || Rf_length(_file) != 1))
			verror("File argument must be a string or NULL");

		if (!Rf_isNull(_intervals_set_out) && (!Rf_isString(_intervals_set_out) || Rf_length(_intervals_set_out) != 1))
			verror("intervals.set.out argument is not a string");

		if (!Rf_isNull(_file) && !Rf_isNull(_intervals_set_out))
			verror("Cannot use both file and intervals.set.out arguments");

		const char *filename = Rf_isNull(_file) ? NULL : CHAR(STRING_ELT(_file, 0));
		string intervset_out = Rf_isNull(_intervals_set_out) ? "" : CHAR(STRING_ELT(_intervals_set_out, 0));
		unsigned num_exprs = (unsigned)Rf_length(_exprs);
		IntervUtils iu(_envir);

		GIntervalsFetcher1D *intervals1d = NULL;
		GIntervalsFetcher2D *intervals2d = NULL;
		iu.convert_rintervs(_intervals, &intervals1d, &intervals2d);
		unique_ptr<GIntervalsFetcher1D> intervals1d_guard(intervals1d);
		unique_ptr<GIntervalsFetcher2D> intervals2d_guard(intervals2d);
		intervals1d->sort();
		intervals2d->sort();
		intervals2d->verify_no_overlaps(iu.get_chromkey());

		TrackExprScanner scanner(iu);

		if (filename) {
			ofstream outfile;

			outfile.open(filename);
			if (outfile.fail())
				verror("Failed to open file %s for writing: %s\n", filename, strerror(errno));
			outfile << setprecision(15);

			scanner.begin(_exprs, intervals1d, intervals2d, _iterator_policy, _band);

			if (scanner.get_iterator()->is_1d()) {
				for (int i = 0; i < GInterval::NUM_COLS; ++i)
					outfile << GInterval::COL_NAMES[i] << "\t";
			} else {
				for (int i = 0; i < GInterval2D::NUM_COLS; ++i)
					outfile << GInterval2D::COL_NAMES[i] << "\t";
			}

			for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
				if (iexpr)
					outfile << "\t";
				if (Rf_isNull(_colnames))
					outfile << get_bounded_colname(CHAR(STRING_ELT(_exprs, iexpr)));
				else
					outfile << CHAR(STRING_ELT(_colnames, iexpr));
			}
			outfile << "\n";

			// Use buffered writing with snprintf for better performance than iostream
			const size_t BUFFER_SIZE = 65536;  // 64KB buffer
			char buffer[BUFFER_SIZE];
			size_t buffer_pos = 0;

			for (; !scanner.isend(); scanner.next()) {
				char line[4096];  // Temporary buffer for one line
				int len = 0;

				if (scanner.get_iterator()->is_1d()) {
					const GInterval &interval = scanner.last_interval1d();
					len = snprintf(line, sizeof(line), "%s\t%" PRId64 "\t%" PRId64,
						iu.id2chrom(interval.chromid).c_str(), interval.start, interval.end);
				} else {
					const GInterval2D &interval = scanner.last_interval2d();
					len = snprintf(line, sizeof(line), "%s\t%" PRId64 "\t%" PRId64 "\t%s\t%" PRId64 "\t%" PRId64,
						iu.id2chrom(interval.chromid1()).c_str(), interval.start1(), interval.end1(),
						iu.id2chrom(interval.chromid2()).c_str(), interval.start2(), interval.end2());
				}

				for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
					len += snprintf(line + len, sizeof(line) - len, "\t%.15g", scanner.last_real(iexpr));
				}
				len += snprintf(line + len, sizeof(line) - len, "\n");

				// Flush buffer if not enough space for this line
				if (buffer_pos + len >= BUFFER_SIZE) {
					outfile.write(buffer, buffer_pos);
					buffer_pos = 0;
				}

				memcpy(buffer + buffer_pos, line, len);
				buffer_pos += len;

				check_interrupt();
			}

			// Flush any remaining data
			if (buffer_pos > 0) {
				outfile.write(buffer, buffer_pos);
			}

			if (outfile.fail())
				verror("Failed to write to file %s: %s\n", filename, strerror(errno));

			return R_NilValue;
		}

		GIntervals out_intervals1d;
		GIntervals2D out_intervals2d;
		vector< vector<double> > values(num_exprs);

		// Pre-reserve memory for the regular (non-file) path to avoid reallocations
		uint64_t max_size = intervals1d ? intervals1d->size() : intervals2d->size();

		if (!intervset_out.empty()) {
			bool is_1d_iterator = iu.is_1d_iterator(_exprs, intervals1d, intervals2d, _iterator_policy);
			vector<GIntervalsBigSet1D::ChromStat> chromstats1d;
			vector<GIntervalsBigSet2D::ChromStat> chromstats2d;
			GInterval last_scope_interval1d;
			GInterval2D last_scope_interval2d;
			uint64_t size;
			char error_prefix[1000];

			if (is_1d_iterator)
				GIntervalsBigSet1D::begin_save(intervset_out.c_str(), iu, chromstats1d);
			else
				GIntervalsBigSet2D::begin_save(intervset_out.c_str(), iu, chromstats2d);

			scanner.begin(_exprs, intervals1d, intervals2d, _iterator_policy, _band);

			while (!scanner.isend()) {
				for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr)
					values[iexpr].push_back(scanner.last_real(iexpr));

				if (is_1d_iterator) {
					if (last_scope_interval1d.chromid != scanner.last_scope_interval1d().chromid) {
						last_scope_interval1d = scanner.last_scope_interval1d();
						snprintf(error_prefix, sizeof(error_prefix), "Big intervals set %s, chrom %s",
								intervset_out.c_str(), iu.id2chrom(last_scope_interval1d.chromid).c_str());
					}
					out_intervals1d.push_back(scanner.last_interval1d());
					size = out_intervals1d.size();
				} else {
					if (!last_scope_interval2d.is_same_chrom(scanner.last_scope_interval2d())) {
						last_scope_interval2d = scanner.last_scope_interval2d();
						snprintf(error_prefix, sizeof(error_prefix), "Big intervals set %s, chroms (%s, %s)",
								intervset_out.c_str(), iu.id2chrom(last_scope_interval2d.chromid1()).c_str(), iu.id2chrom(last_scope_interval2d.chromid2()).c_str());
					}
					out_intervals2d.push_back(scanner.last_interval2d());
					size = out_intervals2d.size();
				}

				iu.verify_max_data_size(size, error_prefix, false);

				scanner.next();

				if (is_1d_iterator) {
					if (scanner.isend() || last_scope_interval1d.chromid != scanner.last_scope_interval1d().chromid) {
						SEXP rintervals = build_rintervals_extract(&out_intervals1d, NULL, values, NULL, _exprs, _colnames, iu);
						GIntervalsBigSet1D::save_chrom(intervset_out.c_str(), &out_intervals1d, rintervals, iu, chromstats1d);
						out_intervals1d.clear();
						for (vector< vector<double> >::iterator ivalues = values.begin(); ivalues != values.end(); ++ivalues) 
							ivalues->clear();
					}
				} else {
					if (scanner.isend() || !last_scope_interval2d.is_same_chrom(scanner.last_scope_interval2d())) {
						SEXP rintervals = build_rintervals_extract(NULL, &out_intervals2d, values, NULL, _exprs, _colnames, iu);
						GIntervalsBigSet2D::save_chrom(intervset_out.c_str(), &out_intervals2d, rintervals, iu, chromstats2d);
						out_intervals2d.clear();
						for (vector< vector<double> >::iterator ivalues = values.begin(); ivalues != values.end(); ++ivalues) 
							ivalues->clear();
					}
				}
			}

			// finish saving (write meta)
			if (is_1d_iterator) {
				SEXP zeroline = build_rintervals_extract(&out_intervals1d, NULL, values, NULL, _exprs, _colnames, iu);
				GIntervalsBigSet1D::end_save(intervset_out.c_str(), zeroline, iu, chromstats1d);
			} else {
				SEXP zeroline = build_rintervals_extract(NULL, &out_intervals2d, values, NULL, _exprs, _colnames, iu);
				GIntervalsBigSet2D::end_save(intervset_out.c_str(), zeroline, iu, chromstats2d);
			}

			return R_NilValue;
		}

		vector<unsigned> interv_ids;

		// Reserve memory to avoid reallocations during the main loop
		for (auto& v : values) {
			v.reserve(max_size);
		}
		if (intervals1d)
			out_intervals1d.reserve(max_size);
		else
			out_intervals2d.reserve(max_size);
		interv_ids.reserve(max_size);

		for (scanner.begin(_exprs, intervals1d, intervals2d, _iterator_policy, _band); !scanner.isend(); scanner.next()) {
			for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr)
				values[iexpr].push_back(scanner.last_real(iexpr));

			if (scanner.get_iterator()->is_1d()) {
				out_intervals1d.push_back(scanner.last_interval1d());
				interv_ids.push_back(iu.get_orig_interv_idx(scanner.last_scope_interval1d()) + 1);
			} else {
				out_intervals2d.push_back(scanner.last_interval2d());
				interv_ids.push_back(iu.get_orig_interv_idx(scanner.last_scope_interval2d()) + 1);
			}

			iu.verify_max_data_size(values[0].size(), "Result");
			check_interrupt();
		}

		if (out_intervals1d.empty() && out_intervals2d.empty())
			return R_NilValue;

		// assemble the answer
		SEXP answer;

		if (!out_intervals1d.empty())
			answer = build_rintervals_extract(&out_intervals1d, NULL, values, &interv_ids, _exprs, _colnames, iu);
		else
			answer = build_rintervals_extract(NULL, &out_intervals2d, values, &interv_ids, _exprs, _colnames, iu);

		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	return R_NilValue;
}

SEXP gextract_multitask(SEXP _intervals, SEXP _exprs, SEXP _colnames, SEXP _iterator_policy, SEXP _band, SEXP _file, SEXP _intervals_set_out, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;
		const bool profile = is_gextract_profile_enabled();

		if (!Rf_isString(_exprs) || Rf_length(_exprs) < 1)
			verror("Tracks expressions argument must be a vector of strings");

		if (!Rf_isNull(_colnames)) {
			if (!Rf_isString(_colnames))
				verror("Column names argument must be a vector of strings");
			if (Rf_length(_colnames) != Rf_length(_exprs))
				verror("Number of column names must match the number of track expressions");
		}

		if (!Rf_isNull(_file) && (!Rf_isString(_file) || Rf_length(_file) != 1))
			verror("File argument must be a string or NULL");

		if (!Rf_isNull(_intervals_set_out) && (!Rf_isString(_intervals_set_out) || Rf_length(_intervals_set_out) != 1))
			verror("intervals.set.out argument is not a string");

		if (!Rf_isNull(_file) && !Rf_isNull(_intervals_set_out))
			verror("Cannot use both file and intervals.set.out arguments");

		const char *filename = Rf_isNull(_file) ? NULL : CHAR(STRING_ELT(_file, 0));
		unsigned num_exprs = (unsigned)Rf_length(_exprs);
		uint64_t num_intervals;
		GIntervals out_intervals1d;
		GIntervals2D out_intervals2d;
		vector<unsigned> interv_ids;
		vector< vector<double> > values(num_exprs);
		IntervUtils iu(_envir);

		GIntervalsFetcher1D *intervals1d = NULL;
		GIntervalsFetcher2D *intervals2d = NULL;
		iu.convert_rintervs(_intervals, &intervals1d, &intervals2d);
		unique_ptr<GIntervalsFetcher1D> intervals1d_guard(intervals1d);
		unique_ptr<GIntervalsFetcher2D> intervals2d_guard(intervals2d);
		intervals1d->sort();
		intervals2d->sort();
		intervals2d->verify_no_overlaps(iu.get_chromkey());

		string intervset_out = Rf_isNull(_intervals_set_out) ? "" : CHAR(STRING_ELT(_intervals_set_out, 0));

		if (filename) {
			TrackExprScanner scanner(iu);
			ofstream outfile;

			outfile.open(filename);
			if (outfile.fail())
				verror("Failed to open file %s for writing: %s\n", filename, strerror(errno));
			outfile << setprecision(15);

			scanner.begin(_exprs, intervals1d, intervals2d, _iterator_policy, _band);

			if (scanner.get_iterator()->is_1d()) {
				for (int i = 0; i < GInterval::NUM_COLS; ++i)
					outfile << GInterval::COL_NAMES[i] << "\t";
			} else {
				for (int i = 0; i < GInterval2D::NUM_COLS; ++i)
					outfile << GInterval2D::COL_NAMES[i] << "\t";
			}

			for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
				if (iexpr)
					outfile << "\t";
				if (Rf_isNull(_colnames))
					outfile << get_bounded_colname(CHAR(STRING_ELT(_exprs, iexpr)));
				else
					outfile << CHAR(STRING_ELT(_colnames, iexpr));
			}
			outfile << "\n";

			// Use buffered writing with snprintf for better performance than iostream
			const size_t BUFFER_SIZE = 65536;  // 64KB buffer
			char buffer[BUFFER_SIZE];
			size_t buffer_pos = 0;

			for (; !scanner.isend(); scanner.next()) {
				char line[4096];  // Temporary buffer for one line
				int len = 0;

				if (scanner.get_iterator()->is_1d()) {
					const GInterval &interval = scanner.last_interval1d();
					len = snprintf(line, sizeof(line), "%s\t%" PRId64 "\t%" PRId64,
						iu.id2chrom(interval.chromid).c_str(), interval.start, interval.end);
				} else {
					const GInterval2D &interval = scanner.last_interval2d();
					len = snprintf(line, sizeof(line), "%s\t%" PRId64 "\t%" PRId64 "\t%s\t%" PRId64 "\t%" PRId64,
						iu.id2chrom(interval.chromid1()).c_str(), interval.start1(), interval.end1(),
						iu.id2chrom(interval.chromid2()).c_str(), interval.start2(), interval.end2());
				}

				for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
					len += snprintf(line + len, sizeof(line) - len, "\t%.15g", scanner.last_real(iexpr));
				}
				len += snprintf(line + len, sizeof(line) - len, "\n");

				// Flush buffer if not enough space for this line
				if (buffer_pos + len >= BUFFER_SIZE) {
					outfile.write(buffer, buffer_pos);
					buffer_pos = 0;
				}

				memcpy(buffer + buffer_pos, line, len);
				buffer_pos += len;

				check_interrupt();
			}

			// Flush any remaining data
			if (buffer_pos > 0) {
				outfile.write(buffer, buffer_pos);
			}

			if (outfile.fail())
				verror("Failed to write to file %s: %s\n", filename, strerror(errno));

			return R_NilValue;
		}

		if (!iu.prepare4multitasking(_exprs, intervals1d, intervals2d, _iterator_policy, _band))
			rreturn(R_NilValue);

		bool is_1d_iterator = iu.is_1d_iterator(_exprs, intervals1d, intervals2d, _iterator_policy);
		// Estimate number of records to cap shared-memory allocation size
		uint64_t estimated_records = estimate_records_for_expr(iu, _exprs, intervals1d, intervals2d, _iterator_policy, _band);

		if (!intervset_out.empty()) {
			vector<GIntervalsBigSet1D::ChromStat> chromstats1d;
			vector<GIntervalsBigSet2D::ChromStat> chromstats2d;

			if (is_1d_iterator)
				GIntervalsBigSet1D::begin_save(intervset_out.c_str(), iu, chromstats1d);
			else
				GIntervalsBigSet2D::begin_save(intervset_out.c_str(), iu, chromstats2d);

			if (iu.distribute_task(is_1d_iterator ?
								   sizeof(GIntervalsBigSet1D::ChromStat) * chromstats1d.size() :
								   sizeof(GIntervalsBigSet2D::ChromStat) * chromstats2d.size(),
								   0))
			{ // child process
				GIntervalsFetcher1D *kid_intervals1d = iu.get_kid_intervals1d();
				GIntervalsFetcher2D *kid_intervals2d = iu.get_kid_intervals2d();
				TrackExprScanner scanner(iu);
				GInterval last_scope_interval1d;
				GInterval2D last_scope_interval2d;
				uint64_t size;
				char error_prefix[1000];

				scanner.begin(_exprs, kid_intervals1d, kid_intervals2d, _iterator_policy, _band);

				while (!scanner.isend()) {
					for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr)
						values[iexpr].push_back(scanner.last_real(iexpr));

					if (is_1d_iterator) {
						if (last_scope_interval1d.chromid != scanner.last_scope_interval1d().chromid) {
							last_scope_interval1d = scanner.last_scope_interval1d();
							snprintf(error_prefix, sizeof(error_prefix), "Big intervals set %s, chrom %s",
									intervset_out.c_str(), iu.id2chrom(last_scope_interval1d.chromid).c_str());
						}
						out_intervals1d.push_back(scanner.last_interval1d());
						size = out_intervals1d.size();
					} else {
						if (!last_scope_interval2d.is_same_chrom(scanner.last_scope_interval2d())) {
							last_scope_interval2d = scanner.last_scope_interval2d();
							snprintf(error_prefix, sizeof(error_prefix), "Big intervals set %s, chroms (%s, %s)",
									intervset_out.c_str(), iu.id2chrom(last_scope_interval2d.chromid1()).c_str(), iu.id2chrom(last_scope_interval2d.chromid2()).c_str());
						}
						out_intervals2d.push_back(scanner.last_interval2d());
						size = out_intervals2d.size();
					}

					iu.verify_max_data_size(size, error_prefix, false);

					scanner.next();

					if (is_1d_iterator) {
						if (scanner.isend() || last_scope_interval1d.chromid != scanner.last_scope_interval1d().chromid) {
							SEXP rintervals = build_rintervals_extract(&out_intervals1d, NULL, values, NULL, _exprs, _colnames, iu);
							GIntervalsBigSet1D::save_chrom(intervset_out.c_str(), &out_intervals1d, rintervals, iu, chromstats1d);
							out_intervals1d.clear();
							for (vector< vector<double> >::iterator ivalues = values.begin(); ivalues != values.end(); ++ivalues) 
								ivalues->clear();
						}
					} else {
						if (scanner.isend() || !last_scope_interval2d.is_same_chrom(scanner.last_scope_interval2d())) {
							SEXP rintervals = build_rintervals_extract(NULL, &out_intervals2d, values, NULL, _exprs, _colnames, iu);
							GIntervalsBigSet2D::save_chrom(intervset_out.c_str(), &out_intervals2d, rintervals, iu, chromstats2d);
							out_intervals2d.clear();
							for (vector< vector<double> >::iterator ivalues = values.begin(); ivalues != values.end(); ++ivalues) 
								ivalues->clear();
						}
					}
				}

				// pack the result into shared memory
				void *ptr = allocate_res(0);

				if (is_1d_iterator) 
					pack_data(ptr, chromstats1d.front(), chromstats1d.size());
				else
					pack_data(ptr, chromstats2d.front(), chromstats2d.size());
			} else { // parent process
				vector<GIntervalsBigSet1D::ChromStat> kid_chromstats1d(chromstats1d.size());
				vector<GIntervalsBigSet2D::ChromStat> kid_chromstats2d(chromstats2d.size());

				for (int i = 0; i < get_num_kids(); ++i) {
					void *ptr = get_kid_res(i);

					if (is_1d_iterator) {
						unpack_data(ptr, kid_chromstats1d.front(), kid_chromstats1d.size());
						for (vector<GIntervalsBigSet1D::ChromStat>::const_iterator istat = kid_chromstats1d.begin(); istat < kid_chromstats1d.end(); ++istat) {
							if (istat->size)
								chromstats1d[istat - kid_chromstats1d.begin()] = *istat;
						}
					} else {
						unpack_data(ptr, kid_chromstats2d.front(), kid_chromstats2d.size());
						for (vector<GIntervalsBigSet2D::ChromStat>::const_iterator istat = kid_chromstats2d.begin(); istat < kid_chromstats2d.end(); ++istat) {
							if (istat->size)
								chromstats2d[istat - kid_chromstats2d.begin()] = *istat;
						}
					}
				}

				// finish saving (write meta)
				if (is_1d_iterator) {
					SEXP zeroline = build_rintervals_extract(&out_intervals1d, NULL, values, NULL, _exprs, _colnames, iu);
					GIntervalsBigSet1D::end_save(intervset_out.c_str(), zeroline, iu, chromstats1d);
				} else {
					SEXP zeroline = build_rintervals_extract(NULL, &out_intervals2d, values, NULL, _exprs, _colnames, iu);
					GIntervalsBigSet2D::end_save(intervset_out.c_str(), zeroline, iu, chromstats2d);
				}
			}
			rreturn(R_NilValue);
		}

		bool do_child = false;
		try {
			if (estimated_records > 0) {
				do_child = iu.distribute_task(0,
											  (is_1d_iterator ? sizeof(GInterval) : sizeof(GInterval2D)) + // interval
											  sizeof(unsigned) +                                                // interval id
											  sizeof(double) * num_exprs,                                      // values
											  rdb::MT_MODE_MMAP,
											  estimated_records);
			} else {
				do_child = iu.distribute_task(0,
											  (is_1d_iterator ? sizeof(GInterval) : sizeof(GInterval2D)) + // interval
											  sizeof(unsigned) +                                                // interval id
											  sizeof(double) * num_exprs);                                     // values
			}
		} catch (TGLException &e) {
			if (string(e.msg()).find("Failed to allocate shared memory") != string::npos) {
				return C_gextract(_intervals, _exprs, _colnames, _iterator_policy, _band, R_NilValue, _intervals_set_out, _envir);
			}
			throw;
		}

		if (do_child)
		{  // child process
			GIntervalsFetcher1D *kid_intervals1d = iu.get_kid_intervals1d();
			GIntervalsFetcher2D *kid_intervals2d = iu.get_kid_intervals2d();
			TrackExprScanner scanner(iu);
			auto t_scan_start = std::chrono::steady_clock::now();

			// Pre-reserve buffers to avoid reallocations during the scan.
			uint64_t kid_size = is_1d_iterator ? kid_intervals1d->size() : kid_intervals2d->size();
			if (is_1d_iterator)
				out_intervals1d.reserve(kid_size);
			else
				out_intervals2d.reserve(kid_size);
			interv_ids.reserve(kid_size);
			for (unsigned i = 0; i < num_exprs; ++i)
				values[i].reserve(kid_size);
			
			for (scanner.begin(_exprs, kid_intervals1d, kid_intervals2d, _iterator_policy, _band); !scanner.isend(); scanner.next()) {
				for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr)
					values[iexpr].push_back(scanner.last_real(iexpr));
			
				if (is_1d_iterator) {
					out_intervals1d.push_back(scanner.last_interval1d());
					interv_ids.push_back(iu.get_orig_interv_idx(scanner.last_scope_interval1d()) + 1);
				} else {
					out_intervals2d.push_back(scanner.last_interval2d());
					interv_ids.push_back(iu.get_orig_interv_idx(scanner.last_scope_interval2d()) + 1);
				}
			
				iu.verify_max_data_size(values[0].size(), "Result");
			}
			
			// now we finally know the result size => pack the result into shared memory
			num_intervals = is_1d_iterator ? out_intervals1d.size() : out_intervals2d.size();

			void *result = allocate_res(num_intervals);
			
			if (!num_intervals)
				rreturn(R_NilValue);

			if (is_1d_iterator)
				pack_data(result, out_intervals1d.front(), num_intervals);
			else
				pack_data(result, out_intervals2d.front(), num_intervals);
			
			pack_data(result, interv_ids.front(), num_intervals);
			
			for (unsigned i = 0; i < num_exprs; ++i)
				pack_data(result, values[i].front(), num_intervals);

			if (profile) {
				auto t_end = std::chrono::steady_clock::now();
				double scan_ms = std::chrono::duration<double, std::milli>(t_end - t_scan_start).count();
				log_gextract_timing("child_scan_pack_ms", scan_ms);
			}
			
		} else {  // parent process
			auto t_gather_start = std::chrono::steady_clock::now();

			// collect results from kids
			// First pass: compute total sizes to pre-reserve and avoid repeated reallocations.
			uint64_t total_intervals = 0;
			for (int i = 0; i < get_num_kids(); ++i) {
				total_intervals += get_kid_res_size(i);
			}
			if (is_1d_iterator)
				out_intervals1d.reserve(total_intervals);
			else
				out_intervals2d.reserve(total_intervals);
			interv_ids.reserve(total_intervals);
			for (unsigned i = 0; i < num_exprs; ++i)
				values[i].reserve(total_intervals);

			for (int i = 0; i < get_num_kids(); ++i) {
				void *ptr = get_kid_res(i);
				num_intervals = get_kid_res_size(i);

				if (!num_intervals)
					continue;

				if (is_1d_iterator) {
					out_intervals1d.insert(out_intervals1d.end(), (GInterval *)ptr, (GInterval *)ptr + num_intervals);
					ptr = (GInterval *)ptr + num_intervals;
				} else {
					out_intervals2d.insert(out_intervals2d.end(), (GInterval2D *)ptr, (GInterval2D *)ptr + num_intervals);
					ptr = (GInterval2D *)ptr + num_intervals;
				}

				interv_ids.insert(interv_ids.end(), (unsigned *)ptr, (unsigned *)ptr + num_intervals);
				ptr = (unsigned *)ptr + num_intervals;

				for (unsigned i = 0; i < num_exprs; ++i) {
					values[i].insert(values[i].end(), (double *)ptr, (double *)ptr + num_intervals);
					ptr = (double *)ptr + num_intervals;
				}
			}

			if (out_intervals1d.empty() && out_intervals2d.empty()) 
				rreturn(R_NilValue);

			auto t_assemble_start = std::chrono::steady_clock::now();

			// assemble the answer
			SEXP answer;
			unsigned num_interv_cols;

			if (!out_intervals1d.empty()) {
				answer = iu.convert_intervs(&out_intervals1d, GInterval::NUM_COLS + num_exprs + 1);
				num_interv_cols = GInterval::NUM_COLS;
			} else {
				answer = iu.convert_intervs(&out_intervals2d, GInterval2D::NUM_COLS + num_exprs + 1);
				num_interv_cols = GInterval2D::NUM_COLS;
			}

            for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
                SEXP expr_vals = rprotect_ptr(RSaneAllocVector(REALSXP, values[iexpr].size()));
				// values[iexpr] is contiguous; copy in one shot
				if (!values[iexpr].empty())
					memcpy(REAL(expr_vals), &values[iexpr][0], sizeof(double) * values[iexpr].size());
                SET_VECTOR_ELT(answer, num_interv_cols + iexpr, expr_vals);
			}

            SEXP ids = rprotect_ptr(RSaneAllocVector(INTSXP, interv_ids.size()));
			if (!interv_ids.empty())
				memcpy(INTEGER(ids), &interv_ids[0], sizeof(int) * interv_ids.size());
			SET_VECTOR_ELT(answer, num_interv_cols + num_exprs, ids);

            SEXP col_names = rprotect_ptr(Rf_getAttrib(answer, R_NamesSymbol));
			for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
				if (Rf_isNull(_colnames))
					SET_STRING_ELT(col_names, num_interv_cols + iexpr, Rf_mkChar(get_bounded_colname(CHAR(STRING_ELT(_exprs, iexpr))).c_str()));
				else
					SET_STRING_ELT(col_names, num_interv_cols + iexpr, STRING_ELT(_colnames, iexpr));
			}
			SET_STRING_ELT(col_names, num_interv_cols + num_exprs, Rf_mkChar("intervalID"));

            runprotect(2); // col_names, ids

			if (profile) {
				auto t_end = std::chrono::steady_clock::now();
				double gather_ms = std::chrono::duration<double, std::milli>(t_assemble_start - t_gather_start).count();
				double assemble_ms = std::chrono::duration<double, std::milli>(t_end - t_assemble_start).count();
				log_gextract_timing("parent_gather_ms", gather_ms);
				log_gextract_timing("parent_assemble_ms", assemble_ms);
			}
            rreturn(answer);
		}
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	rreturn(R_NilValue);
}

}
