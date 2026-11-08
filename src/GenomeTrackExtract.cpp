/*
 * GenomeTrackExtract.cpp
 *
 *  Created on: Mar 21, 2010
 *      Author: hoichman
 */

#include <cctype>
#include <climits>
#include <cstdint>
#include <dirent.h>
#include <inttypes.h>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <strings.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <vector>

#include "rdbinterval.h"
#include "rdbutils.h"
#include "GenomeTrack.h"
#include "GIntervalsBigSet1D.h"
#include "GIntervalsBigSet2D.h"
#include "Parquet.h"
#include "TrackExpressionScanner.h"
#include "TrackExpressionFixedBinIterator.h"
#include "TrackExpressionFixedRectIterator.h"
#include "TrackExpressionIntervals1DIterator.h"
#include "TrackExpressionSparseIterator.h"
#include "WorkerThreads.h"
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

	if (dynamic_cast<TrackExpressionSparseIterator *>(expr_itr)) {
		return 0;
	}

	// Intervals iterators emit one record per (scope interval, iterator interval)
	// intersection, not one per iterator interval: a single long iterator interval
	// (or one block that several overlapping iterator intervals merged into) covers
	// as many scope intervals as it spans. For two sets of disjoint sorted intervals
	// the number of intersecting pairs is at most n + m - 1 per chromosome, so
	// iterator_count + scope_count is a true bound there and stays within a few
	// records of iterator_count in the common case (many small peaks, few scope
	// intervals). It is not a bound for a self-overlapping scope, which gextract
	// allows; that case is caught at runtime by the shared-memory guard, which the
	// parent handles by re-running the call serially.
	if (auto *int1d_itr = dynamic_cast<TrackExpressionIntervals1DIterator *>(expr_itr)) {
		uint64_t itr_size = int1d_itr->get_intervals() ? int1d_itr->get_intervals()->size() : 0;
		if (!itr_size)
			return 0;
		uint64_t scope_size = scope1d ? (uint64_t)scope1d->size() : 0;
		return itr_size + scope_size;
	}

	return 0;
}

static uint64_t inflate_estimated_records(IntervUtils &iu, uint64_t estimated, double factor)
{
	if (!estimated)
		return 0;

	if (!(factor > 1.0))
		return estimated;

	long double inflated = (long double)estimated * factor;
	long double max_records = (long double)iu.get_max_data_size();

	if (inflated > max_records)
		return (uint64_t)max_records;
	if (inflated > (long double)std::numeric_limits<uint64_t>::max())
		return std::numeric_limits<uint64_t>::max();

	return (uint64_t)inflated;
}

static bool allow_multichrom_range_split_for_gextract(SEXP iterator_policy, const string &intervset_out)
{
	// Keep legacy semantics for implicit iterators: splitting scope intervals across kids can
	// change row geometry for iterator modes that depend on original interval boundaries.
	return !Rf_isNull(iterator_policy) && intervset_out.empty();
}

static SEXP build_rintervals_extract(GIntervalsFetcher1D *out_intervals1d, GIntervalsFetcher2D *out_intervals2d, const vector< vector<double> > &values,
									 vector<unsigned> *interv_ids, SEXP _exprs, SEXP _colnames, IntervUtils &iu)
{
	SEXP answer;
	unsigned num_interv_cols;
	unsigned num_exprs = values.size();
	unsigned protect_mark = protect_depth();

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
			SET_STRING_ELT(col_names, num_interv_cols + iexpr, STRING_ELT(_exprs, iexpr));
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

	// Keep "answer" (the first thing convert_intervs protected) and drop everything
	// this function added on top: the value columns, col_names and ids are all
	// reachable from it now. The intervals.set.out path calls this once per contig,
	// so anything left here accumulates on a 50,000-slot stack.
	runprotect_to(protect_mark + 1);
	return answer;
}

static void set_extract_value_column_names(SEXP answer, unsigned num_interv_cols, unsigned num_exprs, SEXP _exprs, SEXP _colnames)
{
	SEXP col_names = rprotect_ptr(Rf_getAttrib(answer, R_NamesSymbol));
	for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
		if (Rf_isNull(_colnames))
			SET_STRING_ELT(col_names, num_interv_cols + iexpr, STRING_ELT(_exprs, iexpr));
		else
			SET_STRING_ELT(col_names, num_interv_cols + iexpr, STRING_ELT(_colnames, iexpr));
	}
	SET_STRING_ELT(col_names, num_interv_cols + num_exprs, Rf_mkChar("intervalID"));
	runprotect(1); // col_names
}

// intervals_join modes recognized by gextract C entrypoints.
enum class IntervalsJoinMode {
	ID,        // keep intervalID (default; current behavior)
	NONE,      // drop intervalID, attach nothing
	INTERVALS  // drop intervalID, attach all columns from input intervals data frame
};

static IntervalsJoinMode parse_intervals_join(SEXP _intervals_join) {
	if (Rf_isNull(_intervals_join))
		return IntervalsJoinMode::ID;
	if (!Rf_isString(_intervals_join) || Rf_length(_intervals_join) != 1)
		verror("intervals_join must be a single string");
	const char *s = CHAR(STRING_ELT(_intervals_join, 0));
	if (!strcmp(s, "id")) return IntervalsJoinMode::ID;
	if (!strcmp(s, "none")) return IntervalsJoinMode::NONE;
	if (!strcmp(s, "intervals")) return IntervalsJoinMode::INTERVALS;
	verror("intervals_join must be one of 'id', 'intervals', 'none' (got '%s')", s);
	return IntervalsJoinMode::ID;
}

// Pick a reasonable thread count for parallel column gathering. Capped by
// hardware concurrency and by gmax.processes (so users keep one knob to tune
// total in-process parallelism). Returns >= 1.
static unsigned pick_gather_threads(IntervUtils &iu, R_xlen_t nrows, int n_in)
{
	// No benefit from threads for tiny outputs.
	if (nrows < 200000 || n_in <= 1)
		return 1;
	unsigned hw = std::thread::hardware_concurrency();
	if (hw == 0)
		hw = 1;
	unsigned cap = (unsigned)iu.get_max_processes();
	if (cap == 0)
		cap = hw;
	unsigned n = std::min<unsigned>(hw, cap);
	// More threads than columns just wastes startup overhead.
	if (n > (unsigned)n_in)
		n = (unsigned)n_in;
	// 8 threads saturate memory bandwidth for the per-column gather; more
	// gives diminishing returns and risks contention with concurrent R work.
	if (n > 8)
		n = 8;
	return n == 0 ? 1 : n;
}

// One column of the intervals_join gather: its type and the data of the source column of
// `_intervals` and of the destination column, fetched on the calling thread - the worker
// threads must not call the R API (INTEGER() and the like may allocate, for one on an ALTREP
// column).
struct GatherColumn {
	int         type;
	const void *src;
	void       *dst;
};

// Gather a single column into its destination, using 1-based indices in `id_p`. Pure C++
// data movement - safe to run from a non-main thread. Callers must validate id_p range
// upfront.
//
// For STRSXP we bypass SET_STRING_ELT and write the SEXP* slots directly.
// This is safe because:
//   - dst_col was freshly allocated by R in the main thread (gen 0); the
//     write barrier (CHECK_OLD_TO_NEW) would be a no-op for old->new with a
//     young target.
//   - Source CHARSXPs live in `_intervals`, which is PROTECTed on the R
//     argument stack throughout the call, so they cannot be reclaimed.
//   - We do not touch any R API from worker threads, so no GC can fire mid
//     loop.
static void gather_one_column(const GatherColumn &col, const int *id_p, R_xlen_t nrows)
{
	switch (col.type) {
		case INTSXP:
		case LGLSXP: {
			const int *src = (const int *)col.src;
			int *dst = (int *)col.dst;
			for (R_xlen_t i = 0; i < nrows; ++i)
				dst[i] = src[id_p[i] - 1];
			break;
		}
		case REALSXP: {
			const double *src = (const double *)col.src;
			double *dst = (double *)col.dst;
			for (R_xlen_t i = 0; i < nrows; ++i)
				dst[i] = src[id_p[i] - 1];
			break;
		}
		case STRSXP: {
			// Direct pointer writes; see contract above.
			const SEXP *src = (const SEXP *)col.src;
			SEXP *dst = (SEXP *)col.dst;
			for (R_xlen_t i = 0; i < nrows; ++i)
				dst[i] = src[id_p[i] - 1];
			break;
		}
		default:
			// Caller must filter unsupported types before dispatch.
			break;
	}
}

// Replace the trailing intervalID column of `answer` according to `mode`.
//   ID        -> return answer unchanged
//   NONE      -> drop the intervalID column, return a new data frame without it
//   INTERVALS -> drop intervalID, attach every column of `_intervals` indexed
//                by the (1-based) intervalID values. Conflicting names get a "1" suffix.
//
// Supported input column types: INTSXP (incl. factor), REALSXP, LGLSXP, STRSXP.
// Anything else triggers a verror naming the column.
//
// In INTERVALS mode the per-column gather runs in parallel across std::thread
// workers; the bottleneck on wide-window workloads is otherwise serial
// per-cell SET_STRING_ELT on 8+ character columns.
static SEXP transform_intervals_join(SEXP answer, SEXP _intervals,
                                     IntervalsJoinMode mode, IntervUtils &iu)
{
	if (Rf_isNull(answer) || mode == IntervalsJoinMode::ID)
		return answer;

	SEXP names = Rf_getAttrib(answer, R_NamesSymbol);
	int n_old = Rf_length(answer);
	if (n_old == 0)
		return answer;

	// Locate the intervalID column (last position by construction).
	int id_idx = n_old - 1;
	if (Rf_isNull(names) || strcmp(CHAR(STRING_ELT(names, id_idx)), "intervalID") != 0)
		return answer; // unexpected layout; bail out conservatively

	SEXP id_col = VECTOR_ELT(answer, id_idx);
	if (TYPEOF(id_col) != INTSXP)
		verror("internal: intervalID column is not integer");
	const int *id_p = INTEGER(id_col);
	R_xlen_t nrows = Rf_xlength(id_col);

	if (mode == IntervalsJoinMode::NONE) {
		// Build a new VECSXP without the intervalID column.
		SEXP out = rprotect_ptr(RSaneAllocVector(VECSXP, n_old - 1));
		SEXP out_names = rprotect_ptr(RSaneAllocVector(STRSXP, n_old - 1));
		for (int i = 0; i < n_old - 1; ++i) {
			SET_VECTOR_ELT(out, i, VECTOR_ELT(answer, i));
			SET_STRING_ELT(out_names, i, STRING_ELT(names, i));
		}
		Rf_setAttrib(out, R_NamesSymbol, out_names);
		Rf_setAttrib(out, R_RowNamesSymbol, Rf_getAttrib(answer, R_RowNamesSymbol));
		Rf_setAttrib(out, R_ClassSymbol, Rf_getAttrib(answer, R_ClassSymbol));
		return out;
	}

	// INTERVALS mode: attach columns from `_intervals` indexed by id_p[i] - 1.
	if (Rf_isNull(_intervals) || !Rf_isVectorList(_intervals))
		verror("intervals_join='intervals' requires intervals to be a data frame");

	SEXP in_names = Rf_getAttrib(_intervals, R_NamesSymbol);
	int n_in = Rf_length(_intervals);
	if (n_in == 0 || Rf_isNull(in_names))
		verror("intervals_join='intervals' requires a named intervals data frame");

	R_xlen_t intervs_nrows = (n_in > 0) ? Rf_xlength(VECTOR_ELT(_intervals, 0)) : 0;

	// One-shot validation of the intervalID range. The hot gather loops then
	// run unconditional-load (no per-cell bounds checks), which is a big win
	// since they're memory-bound.
	{
		int min_id = INT_MAX;
		int max_id = INT_MIN;
		for (R_xlen_t i = 0; i < nrows; ++i) {
			int v = id_p[i];
			if (v < min_id) min_id = v;
			if (v > max_id) max_id = v;
		}
		if (nrows > 0 && (min_id < 1 || max_id > intervs_nrows))
			verror("internal: intervalID range [%d, %d] out of bounds for intervals (nrows=%lld)",
				min_id, max_id, (long long)intervs_nrows);
	}

	// Validate every input column type up front (so worker threads don't have
	// to handle errors).
	for (int j = 0; j < n_in; ++j) {
		SEXP in_col = VECTOR_ELT(_intervals, j);
		int t = TYPEOF(in_col);
		if (t != INTSXP && t != REALSXP && t != LGLSXP && t != STRSXP) {
			const char *col_nm = CHAR(STRING_ELT(in_names, j));
			verror("intervals_join='intervals': column '%s' has unsupported type (%s); supported types are integer, double, logical, character, factor",
				col_nm, Rf_type2char(t));
		}
	}

	// Build target column names: existing answer cols (without intervalID) + renamed input cols.
	int n_keep = n_old - 1;
	int n_total = n_keep + n_in;

	SEXP out = rprotect_ptr(RSaneAllocVector(VECSXP, n_total));
	SEXP out_names = rprotect_ptr(RSaneAllocVector(STRSXP, n_total));

	// Carry over answer columns + their names, except intervalID.
	for (int i = 0; i < n_keep; ++i) {
		SET_VECTOR_ELT(out, i, VECTOR_ELT(answer, i));
		SET_STRING_ELT(out_names, i, STRING_ELT(names, i));
	}

	// Build a name set for collision detection (existing answer cols only).
	// O(n_keep * n_in) is fine - typical column counts are < 50.
	auto name_clashes = [&](const char *nm) -> bool {
		for (int i = 0; i < n_keep; ++i) {
			if (!strcmp(nm, CHAR(STRING_ELT(names, i))))
				return true;
		}
		return false;
	};

	// Pre-allocate all output columns + their names in the main thread (R API
	// is not thread-safe). We park them at out[n_keep + j] so the workers can
	// just do data movement.
	for (int j = 0; j < n_in; ++j) {
		SEXP in_col = VECTOR_ELT(_intervals, j);
		const char *col_nm = CHAR(STRING_ELT(in_names, j));

		// Compute output column name (append "1" once on conflict).
		std::string out_nm = col_nm;
		if (name_clashes(col_nm))
			out_nm += "1";
		SET_STRING_ELT(out_names, n_keep + j, Rf_mkChar(out_nm.c_str()));

		int in_type = TYPEOF(in_col);
		SEXP new_col = R_NilValue;
		switch (in_type) {
			case INTSXP:
				new_col = RSaneAllocVector(INTSXP, nrows);
				break;
			case REALSXP:
				new_col = RSaneAllocVector(REALSXP, nrows);
				break;
			case LGLSXP:
				new_col = RSaneAllocVector(LGLSXP, nrows);
				break;
			case STRSXP:
				new_col = RSaneAllocVector(STRSXP, nrows);
				break;
			default:
				// Already filtered above.
				break;
		}
		// Stash into out before any further allocation so it's GC-anchored.
		SET_VECTOR_ELT(out, n_keep + j, new_col);

		// Preserve factor levels/class on integer columns.
		if (in_type == INTSXP) {
			SEXP levs = Rf_getAttrib(in_col, R_LevelsSymbol);
			if (!Rf_isNull(levs)) {
				Rf_setAttrib(new_col, R_LevelsSymbol, levs);
				Rf_setAttrib(new_col, R_ClassSymbol, Rf_getAttrib(in_col, R_ClassSymbol));
			}
		}
	}

	// Parallel column gather. Workers steal columns via an atomic counter.
	vector<GatherColumn> cols(n_in);
	for (int j = 0; j < n_in; ++j) {
		SEXP in_col = VECTOR_ELT(_intervals, j);
		SEXP dst_col = VECTOR_ELT(out, n_keep + j);

		cols[j].type = TYPEOF(in_col);
		switch (cols[j].type) {
			case INTSXP:
				cols[j].src = INTEGER(in_col);
				cols[j].dst = INTEGER(dst_col);
				break;
			case LGLSXP:
				cols[j].src = LOGICAL(in_col);
				cols[j].dst = LOGICAL(dst_col);
				break;
			case REALSXP:
				cols[j].src = REAL(in_col);
				cols[j].dst = REAL(dst_col);
				break;
			case STRSXP:
				cols[j].src = STRING_PTR_RO(in_col);
				cols[j].dst = (void *)STRING_PTR_RO(dst_col);
				break;
			default:
				cols[j].src = cols[j].dst = NULL;
				break;
		}
	}

	std::atomic<int> next{0};
	run_in_threads(pick_gather_threads(iu, nrows, n_in), [&]() {
		for (int j; (j = next.fetch_add(1, std::memory_order_relaxed)) < n_in; )
			gather_one_column(cols[j], id_p, nrows);
	});

	Rf_setAttrib(out, R_NamesSymbol, out_names);
	Rf_setAttrib(out, R_RowNamesSymbol, Rf_getAttrib(answer, R_RowNamesSymbol));
	Rf_setAttrib(out, R_ClassSymbol, Rf_getAttrib(answer, R_ClassSymbol));
	return out;
}

// Builds the gextract_multitask answer straight from the kids' result blocks in shared
// memory: the interval columns, one column per expression and intervalID, rows in kid
// order. A kid's block holds its intervals (GInterval or GInterval2D records), then their
// intervalIDs, then one run of doubles per expression (see the kid's pack_data() calls).
//
// This is the data frame convert_intervs() builds from the concatenated intervals, down
// to the attributes: factor levels from the chrom key (a separate levels vector per chrom
// column in 2D), names, class, and row names 1..N, which R stores compacted as
// c(NA, N) for N > 2. It skips the intermediate vector of every interval and the
// single-threaded walk over the rows.
//
// All R objects are allocated here, on the calling thread. The worker threads only copy
// into their memory, each taking whole kids, so the rows they write never overlap; they
// make no R API calls. As each part of a kid's block is copied its pages are handed back
// to the system, so at the peak the result is not held twice (the blocks are not read
// again; the shared memory is unmapped when the call ends). Freeing the pages takes about
// as long as the copy and the kernel does it one call at a time, so threads beyond a few
// gain little (genome-wide 16 tracks: 5.3 s on 1 thread, 2.3 s on 4 or 8 or 16).
static SEXP build_answer_from_kids(IntervUtils &iu, bool is_1d, unsigned num_exprs, uint64_t total_intervals,
								   SEXP _exprs, SEXP _colnames)
{
	const unsigned num_interv_cols = is_1d ? (unsigned)GInterval::NUM_COLS : (unsigned)GInterval2D::NUM_COLS;
	const unsigned num_cols = num_interv_cols + num_exprs + 1;
	const unsigned num_chroms = iu.get_chromkey().get_num_chroms();
	const int nkids = get_num_kids();
	const size_t interv_size = is_1d ? sizeof(GInterval) : sizeof(GInterval2D);
	unsigned protect_mark = protect_depth();

	SEXP answer = rprotect_ptr(RSaneAllocVector(VECSXP, num_cols));

	// interval columns: 1D chrom, start, end; 2D chrom1, start1, end1, chrom2, start2, end2
	vector<unsigned> chrom_col_idx;
	vector<int *> chrom_cols;
	vector<double *> coord_cols;
	for (unsigned icol = 0; icol < num_interv_cols; ++icol) {
		bool is_chrom = is_1d ? icol == GInterval::CHROM : (icol == GInterval2D::CHROM1 || icol == GInterval2D::CHROM2);
		SEXP col = RSaneAllocVector(is_chrom ? INTSXP : REALSXP, total_intervals);
		SET_VECTOR_ELT(answer, icol, col);
		if (is_chrom) {
			chrom_col_idx.push_back(icol);
			chrom_cols.push_back(INTEGER(col));
		} else
			coord_cols.push_back(REAL(col));
	}

	vector<double *> value_cols(num_exprs);
	for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
		SEXP col = RSaneAllocVector(REALSXP, total_intervals);
		SET_VECTOR_ELT(answer, num_interv_cols + iexpr, col);
		value_cols[iexpr] = REAL(col);
	}

	SEXP ids = RSaneAllocVector(INTSXP, total_intervals);
	SET_VECTOR_ELT(answer, num_interv_cols + num_exprs, ids);
	int *ids_col = INTEGER(ids);

	vector<const char *> kid_res(nkids);
	vector<uint64_t> kid_size(nkids), kid_offset(nkids);
	uint64_t offset = 0;
	for (int i = 0; i < nkids; ++i) {
		kid_res[i] = (const char *)get_kid_res(i);
		kid_size[i] = get_kid_res_size(i);
		kid_offset[i] = offset;
		offset += kid_size[i];
	}

#ifdef MADV_REMOVE
	const uintptr_t page_size = (uintptr_t)sysconf(_SC_PAGESIZE);
	std::mutex remove_mutex;
#endif

	auto copy_kid = [&](int i) {
		const uint64_t n = kid_size[i];
		const uint64_t o = kid_offset[i];
		const char *ptr = kid_res[i];

		if (!n)
			return;

		// Hands back the pages of this kid's block up to "to" once they are copied. The
		// first page of the block may hold the end of the block before it and the last
		// page the start of the block after it, which other threads may still read: the
		// start is rounded up and "to" down, so those two stay. MADV_DONTNEED unmaps the
		// pages from this process, which threads do in parallel; MADV_REMOVE then frees
		// them, which the kernel does one call at a time per mapping. The threads take
		// turns on a mutex for it: waiting in the kernel, they spin. A failure only keeps
		// the pages until the shared memory is unmapped at the end of the call.
#ifdef MADV_REMOVE
		uintptr_t released = ((uintptr_t)ptr + page_size - 1) & ~(page_size - 1);
#endif
		auto release = [&](const char *to) {
#ifdef MADV_REMOVE
			uintptr_t e = (uintptr_t)to & ~(page_size - 1);
			if (e > released) {
				madvise((void *)released, e - released, MADV_DONTNEED);
				{
					std::lock_guard<std::mutex> lock(remove_mutex);
					madvise((void *)released, e - released, MADV_REMOVE);
				}
				released = e;
			}
#endif
		};

		if (is_1d) {
			int *chroms = chrom_cols[0];
			double *starts = coord_cols[0], *ends = coord_cols[1];
			for (uint64_t j = 0; j < n; ++j) {
				// memcpy: a kid's block starts wherever the previous one ended, so the
				// records need not be aligned
				GInterval interv;
				memcpy((void *)&interv, ptr + j * interv_size, sizeof(GInterval));
				chroms[o + j] = interv.chromid + 1;
				starts[o + j] = interv.start;
				ends[o + j] = interv.end;
			}
		} else {
			int *chroms1 = chrom_cols[0], *chroms2 = chrom_cols[1];
			double *starts1 = coord_cols[0], *ends1 = coord_cols[1], *starts2 = coord_cols[2], *ends2 = coord_cols[3];
			for (uint64_t j = 0; j < n; ++j) {
				GInterval2D interv;
				memcpy((void *)&interv, ptr + j * interv_size, sizeof(GInterval2D));
				chroms1[o + j] = interv.chromid1() + 1;
				starts1[o + j] = interv.start1();
				ends1[o + j] = interv.end1();
				chroms2[o + j] = interv.chromid2() + 1;
				starts2[o + j] = interv.start2();
				ends2[o + j] = interv.end2();
			}
		}
		ptr += n * interv_size;
		release(ptr);

		static_assert(sizeof(unsigned) == sizeof(int), "intervalID is packed as unsigned");
		memcpy(ids_col + o, ptr, n * sizeof(unsigned));
		ptr += n * sizeof(unsigned);

		// column by column, so that a thread holds at most about one column of shared
		// memory that is copied but not yet released
		for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
			memcpy(value_cols[iexpr] + o, ptr, n * sizeof(double));
			ptr += n * sizeof(double);
			release(ptr);
		}
	};

	std::atomic<int> next_kid{0};
	run_in_threads(pick_gather_threads(iu, (R_xlen_t)total_intervals, nkids), [&]() {
		for (int i; (i = next_kid.fetch_add(1)) < nkids; )
			copy_kid(i);
	});

	// factor levels: as convert_intervs() sets them
	for (unsigned icol : chrom_col_idx) {
		SEXP chroms = rprotect_ptr(RSaneAllocVector(STRSXP, num_chroms));
		for (unsigned id = 0; id < num_chroms; ++id)
			SET_STRING_ELT(chroms, id, Rf_mkChar(iu.get_chromkey().id2chrom(id).c_str()));
		SEXP col = VECTOR_ELT(answer, icol);
		Rf_setAttrib(col, R_LevelsSymbol, chroms);
		Rf_setAttrib(col, R_ClassSymbol, Rf_mkString("factor"));
	}

	SEXP col_names = rprotect_ptr(RSaneAllocVector(STRSXP, num_cols));
	for (unsigned icol = 0; icol < num_interv_cols; ++icol)
		SET_STRING_ELT(col_names, icol, Rf_mkChar(is_1d ? GInterval::COL_NAMES[icol] : GInterval2D::COL_NAMES[icol]));

	// convert_intervs() fills 1..N, which R compacts to c(NA, N) when N > 2; give the
	// compact form directly and spare the N-long vector
	SEXP row_names;
	if (total_intervals > 2) {
		row_names = rprotect_ptr(RSaneAllocVector(INTSXP, 2));
		INTEGER(row_names)[0] = NA_INTEGER;
		INTEGER(row_names)[1] = (int)total_intervals;
	} else {
		row_names = rprotect_ptr(RSaneAllocVector(INTSXP, total_intervals));
		for (uint64_t i = 0; i < total_intervals; ++i)
			INTEGER(row_names)[i] = (int)(i + 1);
	}

	Rf_setAttrib(answer, R_NamesSymbol, col_names);
	Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
	Rf_setAttrib(answer, R_RowNamesSymbol, row_names);

	set_extract_value_column_names(answer, num_interv_cols, num_exprs, _exprs, _colnames);

	// everything allocated here is reachable from answer: keep only answer protected
	runprotect_to(protect_mark + 1);
	return answer;
}

static bool is_parquet_file(const char *filename)
{
	size_t len = strlen(filename);
	return len >= 8 && !strcasecmp(filename + len - 8, ".parquet");
}

// Files removed when this goes out of scope. In a multitasking child it removes nothing: the
// parent owns every temporary file of the call.
struct TmpFiles {
	vector<string> paths;
	~TmpFiles() {
		if (!RdbInitializer::is_kid()) {
			for (const string &path : paths)
				unlink(path.c_str());
		}
	}
};

// The file a write to filename replaces: filename, or where its chain of symlinks leads, whether
// that file exists or not. An ofstream writes through a link, and so does a Parquet output.
static string resolve_symlinks(const char *filename)
{
	string path(filename);
	struct stat st;

	for (int hops = 0; !lstat(path.c_str(), &st) && S_ISLNK(st.st_mode); ++hops) {
		if (hops == 40)
			verror("Failed to write %s: too many levels of symbolic links", filename);

		char text[PATH_MAX];
		ssize_t len = readlink(path.c_str(), text, sizeof(text));
		if (len <= 0 || len == (ssize_t)sizeof(text))
			verror("Failed to read symbolic link %s: %s", path.c_str(), len < 0 ? strerror(errno) : "invalid length");
		// a relative link is relative to the directory of the link (npos + 1 == 0: no directory)
		size_t slash = path.rfind('/');
		path = (text[0] == '/' ? "" : path.substr(0, slash + 1)) + string(text, len);
	}
	return path;
}

// This host as temporary names carry it, as .gdb.host_tag() writes it: the name up to the first
// dot, anything but [A-Za-z0-9_] replaced by '_'
static string host_tag()
{
	char buf[256] = "";
	gethostname(buf, sizeof(buf) - 1);
	string host(buf, strcspn(buf, "."));
	for (char &ch : host) {
		if (!isalnum((unsigned char)ch) && ch != '_')
			ch = '_';
	}
	return host.empty() ? "unknown" : host;
}

// A number as to_string() writes it, small enough for a pid
static bool is_decimal(const string &s)
{
	return !s.empty() && s.size() <= 9 && s.find_first_not_of("0123456789") == string::npos && (s[0] != '0' || s.size() == 1);
}

// Where a Parquet output goes, and the hidden file it is written to first: tmp, next to target
// (the same filesystem, so rename() works), named for this host and process, as two hosts on
// shared storage may run writers of the same pid. The parts of a multitasking write add .<kid>.
struct ParquetOutput {
	string target;  // the target, or where its symlinks lead
	string tmp;     // .<target name>.tmp.<host>.<pid>
};

// Also removes the leftovers of writers of the same target on this host that were killed
// (a SIGKILL leaves them no chance to clean up): temporary files and parts of pids that no longer
// exist. Files of another host, of a live pid, or with a name of another form are left alone.
static ParquetOutput parquet_output(const char *filename)
{
	ParquetOutput out;
	out.target = resolve_symlinks(filename);
	size_t slash = out.target.rfind('/');  // npos + 1 == 0: no directory
	string dir = out.target.substr(0, slash + 1);
	string prefix = "." + out.target.substr(slash + 1) + ".tmp." + host_tag() + ".";
	out.tmp = dir + prefix + to_string(getpid());

	DIR *dp = opendir(dir.empty() ? "." : dir.c_str());
	if (!dp)
		return out;  // writing tmp fails with the reason
	for (struct dirent *ent; (ent = readdir(dp)); ) {
		string name(ent->d_name);
		if (name.compare(0, prefix.size(), prefix))
			continue;
		string rest = name.substr(prefix.size());  // <pid> or <pid>.<kid>
		size_t dot = rest.find('.');
		string pid = rest.substr(0, dot);
		if (is_decimal(pid) && (dot == string::npos || is_decimal(rest.substr(dot + 1))) && kill((pid_t)stol(pid), 0) && errno == ESRCH)
			unlink((dir + name).c_str());
	}
	closedir(dp);
	return out;
}

// Renames the finished temporary file onto the target. A target it replaces passes on its
// permission bits, as a write in place would keep them (best effort: if chmod fails, the file
// has the permissions it was created with).
static void commit_parquet(const ParquetOutput &out)
{
	struct stat st;

	if (!stat(out.target.c_str(), &st))
		chmod(out.tmp.c_str(), st.st_mode & 0777);
	if (rename(out.tmp.c_str(), out.target.c_str()))
		verror("Failed to rename %s to %s: %s", out.tmp.c_str(), out.target.c_str(), strerror(errno));
}

// The value columns: colnames, or else the expressions
static vector<string> parquet_value_names(SEXP _exprs, SEXP _colnames)
{
	vector<string> names;
	for (int iexpr = 0; iexpr < Rf_length(_exprs); ++iexpr)
		names.push_back(CHAR(STRING_ELT(Rf_isNull(_colnames) ? _exprs : _colnames, iexpr)));
	return names;
}

// start and end are INT64 only if a chromosome of the scope does not fit in INT32
static bool parquet_int64_coords(IntervUtils &iu, GIntervalsFetcher1D *scope)
{
	const GenomeChromKey &chromkey = iu.get_chromkey();
	for (int chromid = 0; chromid < (int)chromkey.get_num_chroms(); ++chromid) {
		if (scope->size(chromid) && chromkey.get_chrom_size(chromid) > (uint64_t)numeric_limits<int32_t>::max())
			return true;
	}
	return false;
}

static string parquet_created_by()
{
	SEXP version = run_in_R("as.character(utils::packageVersion(\"misha\"))", R_BaseEnv);
	if (!Rf_isString(version) || Rf_length(version) != 1)
		verror("Failed to read the version of misha");
	return string("misha version ") + CHAR(STRING_ELT(version, 0));
}

// Writes the rows of a begun scan to a new Parquet file, path; messages call it target
static void write_parquet_rows(TrackExprScanner &scanner, IntervUtils &iu, const string &path, const char *target, const vector<string> &names,
							   bool int64_coords, const string &created_by)
{
	if (!scanner.get_iterator()->is_1d())
		verror("Parquet output does not support 2D iterators");

	ParquetWriter writer;
	vector<double> values(names.size());

	writer.open(path, names, int64_coords, created_by, target);
	for (; !scanner.isend(); scanner.next()) {
		const GInterval &interval = scanner.last_interval1d();
		for (unsigned iexpr = 0; iexpr < values.size(); ++iexpr)
			values[iexpr] = scanner.last_real(iexpr);
		writer.add_row(iu.id2chrom(interval.chromid), interval.start, interval.end, values.data());
		check_interrupt();
	}
	writer.close();
}

// Writes the scan to a Parquet file. The rows go to a temporary file renamed onto the target once
// complete, so an error or an interrupt leaves the target as it was.
static void write_parquet(TrackExprScanner &scanner, IntervUtils &iu, const char *filename, SEXP _exprs, SEXP _colnames,
						  GIntervalsFetcher1D *intervals1d, GIntervalsFetcher2D *intervals2d, SEXP _iterator_policy, SEXP _band)
{
	scanner.begin(_exprs, intervals1d, intervals2d, _iterator_policy, _band);

	ParquetOutput out = parquet_output(filename);
	TmpFiles tmp { { out.tmp } };
	write_parquet_rows(scanner, iu, out.tmp, filename, parquet_value_names(_exprs, _colnames), parquet_int64_coords(iu, intervals1d), parquet_created_by());
	commit_parquet(out);
	tmp.paths.clear();
}

extern "C" {

SEXP C_gextract(SEXP _intervals, SEXP _exprs, SEXP _colnames, SEXP _iterator_policy, SEXP _band, SEXP _file, SEXP _intervals_set_out, SEXP _intervals_join, SEXP _envir)
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

		IntervalsJoinMode join_mode = parse_intervals_join(_intervals_join);
		if (join_mode == IntervalsJoinMode::INTERVALS) {
			if (!Rf_isNull(_file))
				verror("intervals_join='intervals' cannot be combined with 'file'");
			if (!Rf_isNull(_intervals_set_out))
				verror("intervals_join='intervals' cannot be combined with 'intervals.set.out'");
		}

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

		if (filename && is_parquet_file(filename)) {
			write_parquet(scanner, iu, filename, _exprs, _colnames, intervals1d, intervals2d, _iterator_policy, _band);
			return R_NilValue;
		}

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
					outfile << CHAR(STRING_ELT(_exprs, iexpr));
				else
					outfile << CHAR(STRING_ELT(_colnames, iexpr));
			}
			outfile << "\n";

			// Use buffered writing with snprintf for better performance than iostream
			const size_t BUFFER_SIZE = 65536;  // 64KB buffer
			char buffer[BUFFER_SIZE];
			size_t buffer_pos = 0;

			for (; !scanner.isend(); scanner.next()) {
				std::string line;  // one output row; grows as needed (no fixed cap)
				char field[64];
				int n;

				if (scanner.get_iterator()->is_1d()) {
					const GInterval &interval = scanner.last_interval1d();
					line = iu.id2chrom(interval.chromid);
					n = snprintf(field, sizeof(field), "\t%" PRId64 "\t%" PRId64, interval.start, interval.end);
					line.append(field, n);
				} else {
					const GInterval2D &interval = scanner.last_interval2d();
					line = iu.id2chrom(interval.chromid1());
					n = snprintf(field, sizeof(field), "\t%" PRId64 "\t%" PRId64, interval.start1(), interval.end1());
					line.append(field, n);
					line += '\t';
					line += iu.id2chrom(interval.chromid2());
					n = snprintf(field, sizeof(field), "\t%" PRId64 "\t%" PRId64, interval.start2(), interval.end2());
					line.append(field, n);
				}

				for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
					n = snprintf(field, sizeof(field), "\t%.15g", scanner.last_real(iexpr));
					line.append(field, n);
				}
				line += '\n';

				// Flush buffer if not enough space for this line
				if (buffer_pos + line.size() >= BUFFER_SIZE) {
					outfile.write(buffer, buffer_pos);
					buffer_pos = 0;
				}

				if (line.size() >= BUFFER_SIZE) {
					outfile.write(line.data(), line.size());  // row larger than buffer: write directly
				} else {
					memcpy(buffer + buffer_pos, line.data(), line.size());
					buffer_pos += line.size();
				}

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

		// Pre-compute max_size for memory allocation
		uint64_t max_size = intervals1d ? intervals1d->size() : intervals2d->size();

		if (!intervset_out.empty()) {
			// intervals_set_out path: uses C++ vectors + build_rintervals_extract for per-chrom saving
			GIntervals out_intervals1d;
			GIntervals2D out_intervals2d;
			vector< vector<double> > values(num_exprs);
			bool is_1d_iterator = iu.is_1d_iterator(_exprs, intervals1d, intervals2d, _iterator_policy);
			vector<GIntervalsBigSet1D::ChromStat> chromstats1d;
			GIntervalsMeta2D::ChromStats2D chromstats2d;
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
						runprotect(rintervals); // one contig's worth of protect slots; the loop runs per contig
						out_intervals1d.clear();
						for (vector< vector<double> >::iterator ivalues = values.begin(); ivalues != values.end(); ++ivalues)
							ivalues->clear();
					}
				} else {
					if (scanner.isend() || !last_scope_interval2d.is_same_chrom(scanner.last_scope_interval2d())) {
						SEXP rintervals = build_rintervals_extract(NULL, &out_intervals2d, values, NULL, _exprs, _colnames, iu);
						GIntervalsBigSet2D::save_chrom(intervset_out.c_str(), &out_intervals2d, rintervals, iu, chromstats2d);
						runprotect(rintervals); // one contig's worth of protect slots; the loop runs per contig
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

		// ---- Direct-to-R accumulation path (in-memory) ----
		// When we have a good size estimate: pre-allocate R vectors and write directly.
		// When we don't (sparse/intervals iterators): accumulate in C++ vectors,
		// then copy directly to R vectors (skipping convert_intervs virtual iteration).
		uint64_t estimated = estimate_records_for_expr(iu, _exprs, intervals1d, intervals2d, _iterator_policy, _band);

		if (estimated == 0) {
			// Unknown result size — use safe C++ push_back accumulation,
			// then build R vectors directly (bypassing convert_intervs).
			GIntervals out_intervals1d;
			GIntervals2D out_intervals2d;
			vector< vector<double> > values(num_exprs);
			vector<unsigned> interv_ids;

			if (intervals1d)
				out_intervals1d.reserve(max_size);
			else
				out_intervals2d.reserve(max_size);
			for (auto& v : values)
				v.reserve(max_size);
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

			// Build R data frame directly from C++ vectors (no convert_intervs)
			bool is_1d_fb = !out_intervals1d.empty();
			uint64_t nrows = is_1d_fb ? out_intervals1d.size() : out_intervals2d.size();
			unsigned fb_interv_cols = is_1d_fb ? (unsigned)GInterval::NUM_COLS : (unsigned)GInterval2D::NUM_COLS;
			unsigned fb_total_cols = fb_interv_cols + num_exprs + 1;

			SEXP fb_answer;
			rprotect(fb_answer = RSaneAllocVector(VECSXP, fb_total_cols));

			if (is_1d_fb) {
				SEXP fb_chroms, fb_starts, fb_ends;
				rprotect(fb_chroms = RSaneAllocVector(INTSXP, nrows));
				rprotect(fb_starts = RSaneAllocVector(REALSXP, nrows));
				rprotect(fb_ends = RSaneAllocVector(REALSXP, nrows));
				int *pc = INTEGER(fb_chroms);
				double *ps = REAL(fb_starts), *pe = REAL(fb_ends);
				for (uint64_t i = 0; i < nrows; i++) {
					pc[i] = out_intervals1d[i].chromid + 1;
					ps[i] = out_intervals1d[i].start;
					pe[i] = out_intervals1d[i].end;
				}
				unsigned nc = iu.get_chromkey().get_num_chroms();
				SEXP chrom_lvls;
				rprotect(chrom_lvls = RSaneAllocVector(STRSXP, nc));
				for (unsigned id = 0; id < nc; ++id)
					SET_STRING_ELT(chrom_lvls, id, Rf_mkChar(iu.get_chromkey().id2chrom(id).c_str()));
				Rf_setAttrib(fb_chroms, R_LevelsSymbol, chrom_lvls);
				Rf_setAttrib(fb_chroms, R_ClassSymbol, Rf_mkString("factor"));
				SET_VECTOR_ELT(fb_answer, GInterval::CHROM, fb_chroms);
				SET_VECTOR_ELT(fb_answer, GInterval::START, fb_starts);
				SET_VECTOR_ELT(fb_answer, GInterval::END, fb_ends);
			} else {
				SEXP fb_c1, fb_s1, fb_e1, fb_c2, fb_s2, fb_e2;
				rprotect(fb_c1 = RSaneAllocVector(INTSXP, nrows));
				rprotect(fb_s1 = RSaneAllocVector(REALSXP, nrows));
				rprotect(fb_e1 = RSaneAllocVector(REALSXP, nrows));
				rprotect(fb_c2 = RSaneAllocVector(INTSXP, nrows));
				rprotect(fb_s2 = RSaneAllocVector(REALSXP, nrows));
				rprotect(fb_e2 = RSaneAllocVector(REALSXP, nrows));
				int *pc1 = INTEGER(fb_c1), *pc2 = INTEGER(fb_c2);
				double *ps1 = REAL(fb_s1), *pe1 = REAL(fb_e1);
				double *ps2 = REAL(fb_s2), *pe2 = REAL(fb_e2);
				for (uint64_t i = 0; i < nrows; i++) {
					pc1[i] = out_intervals2d[i].chromid1() + 1;
					ps1[i] = out_intervals2d[i].start1();
					pe1[i] = out_intervals2d[i].end1();
					pc2[i] = out_intervals2d[i].chromid2() + 1;
					ps2[i] = out_intervals2d[i].start2();
					pe2[i] = out_intervals2d[i].end2();
				}
				unsigned nc = iu.get_chromkey().get_num_chroms();
				SEXP chrom_lvls;
				rprotect(chrom_lvls = RSaneAllocVector(STRSXP, nc));
				for (unsigned id = 0; id < nc; ++id)
					SET_STRING_ELT(chrom_lvls, id, Rf_mkChar(iu.get_chromkey().id2chrom(id).c_str()));
				Rf_setAttrib(fb_c1, R_LevelsSymbol, chrom_lvls);
				Rf_setAttrib(fb_c1, R_ClassSymbol, Rf_mkString("factor"));
				Rf_setAttrib(fb_c2, R_LevelsSymbol, chrom_lvls);
				Rf_setAttrib(fb_c2, R_ClassSymbol, Rf_mkString("factor"));
				SET_VECTOR_ELT(fb_answer, GInterval2D::CHROM1, fb_c1);
				SET_VECTOR_ELT(fb_answer, GInterval2D::START1, fb_s1);
				SET_VECTOR_ELT(fb_answer, GInterval2D::END1, fb_e1);
				SET_VECTOR_ELT(fb_answer, GInterval2D::CHROM2, fb_c2);
				SET_VECTOR_ELT(fb_answer, GInterval2D::START2, fb_s2);
				SET_VECTOR_ELT(fb_answer, GInterval2D::END2, fb_e2);
			}

			// Expression value columns (memcpy from C++ vectors)
			for (unsigned iexpr = 0; iexpr < num_exprs; iexpr++) {
				SEXP rv = rprotect_ptr(RSaneAllocVector(REALSXP, nrows));
				memcpy(REAL(rv), values[iexpr].data(), nrows * sizeof(double));
				SET_VECTOR_ELT(fb_answer, fb_interv_cols + iexpr, rv);
			}

			// IntervalID column
			SEXP fb_ids;
			rprotect(fb_ids = RSaneAllocVector(INTSXP, nrows));
			for (uint64_t i = 0; i < nrows; i++)
				INTEGER(fb_ids)[i] = interv_ids[i];
			SET_VECTOR_ELT(fb_answer, fb_interv_cols + num_exprs, fb_ids);

			// Column names
			SEXP fb_colnames;
			rprotect(fb_colnames = RSaneAllocVector(STRSXP, fb_total_cols));
			if (is_1d_fb) {
				for (int i = 0; i < GInterval::NUM_COLS; i++)
					SET_STRING_ELT(fb_colnames, i, Rf_mkChar(GInterval::COL_NAMES[i]));
			} else {
				for (int i = 0; i < GInterval2D::NUM_COLS; i++)
					SET_STRING_ELT(fb_colnames, i, Rf_mkChar(GInterval2D::COL_NAMES[i]));
			}
			for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
				const char *cn = Rf_isNull(_colnames) ?
					CHAR(STRING_ELT(_exprs, iexpr)) :
					CHAR(STRING_ELT(_colnames, iexpr));
				SET_STRING_ELT(fb_colnames, fb_interv_cols + iexpr, Rf_mkChar(cn));
			}
			SET_STRING_ELT(fb_colnames, fb_interv_cols + num_exprs, Rf_mkChar("intervalID"));
			Rf_setAttrib(fb_answer, R_NamesSymbol, fb_colnames);

			// Row names and class
			SEXP fb_rownames;
			rprotect(fb_rownames = RSaneAllocVector(INTSXP, nrows));
			for (uint64_t i = 0; i < nrows; i++)
				INTEGER(fb_rownames)[i] = (int)(i + 1);
			Rf_setAttrib(fb_answer, R_RowNamesSymbol, fb_rownames);
			Rf_setAttrib(fb_answer, R_ClassSymbol, Rf_mkString("data.frame"));

			return transform_intervals_join(fb_answer, _intervals, join_mode, iu);
		}

		// We have a good estimate — use direct-to-R pre-allocation.
		// For intervals iterators this is iterator_count (exact when each
		// interval overlaps ≤1 scope interval; may underestimate with
		// fragmented scopes). The scan loop grows vectors if needed.
		uint64_t alloc_size = estimated;
		{
			uint64_t max_data = iu.get_max_data_size();
			if (alloc_size > max_data)
				alloc_size = max_data;
		}

		scanner.begin(_exprs, intervals1d, intervals2d, _iterator_policy, _band);

		if (scanner.isend())
			return R_NilValue;

		bool is_1d = scanner.get_iterator()->is_1d();
		unsigned num_interv_cols = is_1d ? (unsigned)GInterval::NUM_COLS : (unsigned)GInterval2D::NUM_COLS;
		unsigned total_cols = num_interv_cols + num_exprs + 1; // +1 for intervalID

		// Pre-allocate all R vectors at alloc_size
		SEXP answer;
		rprotect(answer = RSaneAllocVector(VECSXP, total_cols));

		SEXP r_row_names;
		rprotect(r_row_names = RSaneAllocVector(INTSXP, alloc_size));
		int *p_row_names = INTEGER(r_row_names);

		// Expression value columns
		vector<SEXP> r_expr_vals(num_exprs);
		vector<double*> p_expr_vals(num_exprs);
		for (unsigned i = 0; i < num_exprs; i++) {
			r_expr_vals[i] = rprotect_ptr(RSaneAllocVector(REALSXP, alloc_size));
			p_expr_vals[i] = REAL(r_expr_vals[i]);
		}

		// IntervalID column
		SEXP r_ids;
		rprotect(r_ids = RSaneAllocVector(INTSXP, alloc_size));
		int *p_ids = INTEGER(r_ids);

		uint64_t row = 0;

		// Helper: grow an R vector by 2× and update pointer.
		// Used when intervals iterator estimate underestimates (rare: scope splitting).
		auto grow_intsxp = [&](SEXP &vec, int* &ptr) {
			SEXP newvec = rprotect_ptr(RSaneAllocVector(INTSXP, alloc_size));
			memcpy(INTEGER(newvec), ptr, row * sizeof(int));
			vec = newvec;
			ptr = INTEGER(newvec);
		};
		auto grow_realsxp = [&](SEXP &vec, double* &ptr) {
			SEXP newvec = rprotect_ptr(RSaneAllocVector(REALSXP, alloc_size));
			memcpy(REAL(newvec), ptr, row * sizeof(double));
			vec = newvec;
			ptr = REAL(newvec);
		};

		if (is_1d) {
			// 1D interval columns
			SEXP r_chroms, r_starts, r_ends;
			rprotect(r_chroms = RSaneAllocVector(INTSXP, alloc_size));
			rprotect(r_starts = RSaneAllocVector(REALSXP, alloc_size));
			rprotect(r_ends = RSaneAllocVector(REALSXP, alloc_size));
			int *p_chroms = INTEGER(r_chroms);
			double *p_starts = REAL(r_starts);
			double *p_ends = REAL(r_ends);

			// Main scan loop for 1D
			for (; !scanner.isend(); scanner.next()) {
				// Grow all vectors if estimate was too small (rare: scope splitting)
				if (row >= alloc_size) {
					alloc_size = alloc_size * 2;
					grow_intsxp(r_chroms, p_chroms);
					grow_realsxp(r_starts, p_starts);
					grow_realsxp(r_ends, p_ends);
					grow_intsxp(r_ids, p_ids);
					grow_intsxp(r_row_names, p_row_names);
					for (unsigned i = 0; i < num_exprs; i++)
						grow_realsxp(r_expr_vals[i], p_expr_vals[i]);
				}

				const GInterval &interval = scanner.last_interval1d();
				p_chroms[row] = interval.chromid + 1;  // R factors are 1-based
				p_starts[row] = interval.start;
				p_ends[row] = interval.end;
				p_ids[row] = iu.get_orig_interv_idx(scanner.last_scope_interval1d()) + 1;

				for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr)
					p_expr_vals[iexpr][row] = scanner.last_real(iexpr);

				p_row_names[row] = (int)(row + 1);
				row++;

				iu.verify_max_data_size(row, "Result");
				check_interrupt();
			}

			// Truncate all vectors to actual size
			r_chroms = rprotect_ptr(Rf_lengthgets(r_chroms, row));
			r_starts = rprotect_ptr(Rf_lengthgets(r_starts, row));
			r_ends   = rprotect_ptr(Rf_lengthgets(r_ends, row));

			// Set chrom as factor
			unsigned num_chroms_total = iu.get_chromkey().get_num_chroms();
			SEXP chrom_levels;
			rprotect(chrom_levels = RSaneAllocVector(STRSXP, num_chroms_total));
			for (unsigned id = 0; id < num_chroms_total; ++id)
				SET_STRING_ELT(chrom_levels, id, Rf_mkChar(iu.get_chromkey().id2chrom(id).c_str()));
			Rf_setAttrib(r_chroms, R_LevelsSymbol, chrom_levels);
			Rf_setAttrib(r_chroms, R_ClassSymbol, Rf_mkString("factor"));

			SET_VECTOR_ELT(answer, GInterval::CHROM, r_chroms);
			SET_VECTOR_ELT(answer, GInterval::START, r_starts);
			SET_VECTOR_ELT(answer, GInterval::END, r_ends);
		} else {
			// 2D interval columns
			SEXP r_chroms1, r_starts1, r_ends1, r_chroms2, r_starts2, r_ends2;
			rprotect(r_chroms1 = RSaneAllocVector(INTSXP, alloc_size));
			rprotect(r_starts1 = RSaneAllocVector(REALSXP, alloc_size));
			rprotect(r_ends1 = RSaneAllocVector(REALSXP, alloc_size));
			rprotect(r_chroms2 = RSaneAllocVector(INTSXP, alloc_size));
			rprotect(r_starts2 = RSaneAllocVector(REALSXP, alloc_size));
			rprotect(r_ends2 = RSaneAllocVector(REALSXP, alloc_size));
			int *p_chroms1 = INTEGER(r_chroms1);
			double *p_starts1 = REAL(r_starts1);
			double *p_ends1 = REAL(r_ends1);
			int *p_chroms2 = INTEGER(r_chroms2);
			double *p_starts2 = REAL(r_starts2);
			double *p_ends2 = REAL(r_ends2);

			// Main scan loop for 2D
			for (; !scanner.isend(); scanner.next()) {
				if (row >= alloc_size) {
					alloc_size = alloc_size * 2;
					grow_intsxp(r_chroms1, p_chroms1);
					grow_realsxp(r_starts1, p_starts1);
					grow_realsxp(r_ends1, p_ends1);
					grow_intsxp(r_chroms2, p_chroms2);
					grow_realsxp(r_starts2, p_starts2);
					grow_realsxp(r_ends2, p_ends2);
					grow_intsxp(r_ids, p_ids);
					grow_intsxp(r_row_names, p_row_names);
					for (unsigned i = 0; i < num_exprs; i++)
						grow_realsxp(r_expr_vals[i], p_expr_vals[i]);
				}

				const GInterval2D &interval = scanner.last_interval2d();
				p_chroms1[row] = interval.chromid1() + 1;
				p_starts1[row] = interval.start1();
				p_ends1[row] = interval.end1();
				p_chroms2[row] = interval.chromid2() + 1;
				p_starts2[row] = interval.start2();
				p_ends2[row] = interval.end2();
				p_ids[row] = iu.get_orig_interv_idx(scanner.last_scope_interval2d()) + 1;

				for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr)
					p_expr_vals[iexpr][row] = scanner.last_real(iexpr);

				p_row_names[row] = (int)(row + 1);
				row++;

				iu.verify_max_data_size(row, "Result");
				check_interrupt();
			}

			// Truncate all vectors to actual size
			r_chroms1 = rprotect_ptr(Rf_lengthgets(r_chroms1, row));
			r_starts1 = rprotect_ptr(Rf_lengthgets(r_starts1, row));
			r_ends1   = rprotect_ptr(Rf_lengthgets(r_ends1, row));
			r_chroms2 = rprotect_ptr(Rf_lengthgets(r_chroms2, row));
			r_starts2 = rprotect_ptr(Rf_lengthgets(r_starts2, row));
			r_ends2   = rprotect_ptr(Rf_lengthgets(r_ends2, row));

			// Set chrom1 and chrom2 as factors
			unsigned num_chroms_total = iu.get_chromkey().get_num_chroms();
			SEXP chrom_levels1, chrom_levels2;
			rprotect(chrom_levels1 = RSaneAllocVector(STRSXP, num_chroms_total));
			rprotect(chrom_levels2 = RSaneAllocVector(STRSXP, num_chroms_total));
			for (unsigned id = 0; id < num_chroms_total; ++id) {
				SET_STRING_ELT(chrom_levels1, id, Rf_mkChar(iu.get_chromkey().id2chrom(id).c_str()));
				SET_STRING_ELT(chrom_levels2, id, Rf_mkChar(iu.get_chromkey().id2chrom(id).c_str()));
			}
			Rf_setAttrib(r_chroms1, R_LevelsSymbol, chrom_levels1);
			Rf_setAttrib(r_chroms1, R_ClassSymbol, Rf_mkString("factor"));
			Rf_setAttrib(r_chroms2, R_LevelsSymbol, chrom_levels2);
			Rf_setAttrib(r_chroms2, R_ClassSymbol, Rf_mkString("factor"));

			SET_VECTOR_ELT(answer, GInterval2D::CHROM1, r_chroms1);
			SET_VECTOR_ELT(answer, GInterval2D::START1, r_starts1);
			SET_VECTOR_ELT(answer, GInterval2D::END1, r_ends1);
			SET_VECTOR_ELT(answer, GInterval2D::CHROM2, r_chroms2);
			SET_VECTOR_ELT(answer, GInterval2D::START2, r_starts2);
			SET_VECTOR_ELT(answer, GInterval2D::END2, r_ends2);
		}

		if (row == 0)
			return R_NilValue;

		// Truncate shared columns to actual size
		for (unsigned i = 0; i < num_exprs; i++)
			r_expr_vals[i] = rprotect_ptr(Rf_lengthgets(r_expr_vals[i], row));
		r_ids = rprotect_ptr(Rf_lengthgets(r_ids, row));
		r_row_names = rprotect_ptr(Rf_lengthgets(r_row_names, row));

		// Set expression value columns
		for (unsigned iexpr = 0; iexpr < num_exprs; iexpr++)
			SET_VECTOR_ELT(answer, num_interv_cols + iexpr, r_expr_vals[iexpr]);

		// Set IntervalID column
		SET_VECTOR_ELT(answer, num_interv_cols + num_exprs, r_ids);

		// Column names
		SEXP col_names;
		rprotect(col_names = RSaneAllocVector(STRSXP, total_cols));
		if (is_1d) {
			for (int i = 0; i < GInterval::NUM_COLS; i++)
				SET_STRING_ELT(col_names, i, Rf_mkChar(GInterval::COL_NAMES[i]));
		} else {
			for (int i = 0; i < GInterval2D::NUM_COLS; i++)
				SET_STRING_ELT(col_names, i, Rf_mkChar(GInterval2D::COL_NAMES[i]));
		}
		for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
			if (Rf_isNull(_colnames))
				SET_STRING_ELT(col_names, num_interv_cols + iexpr, STRING_ELT(_exprs, iexpr));
			else
				SET_STRING_ELT(col_names, num_interv_cols + iexpr, STRING_ELT(_colnames, iexpr));
		}
		SET_STRING_ELT(col_names, num_interv_cols + num_exprs, Rf_mkChar("intervalID"));

		Rf_setAttrib(answer, R_NamesSymbol, col_names);
		Rf_setAttrib(answer, R_RowNamesSymbol, r_row_names);
		Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));

		return transform_intervals_join(answer, _intervals, join_mode, iu);
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	}
	return R_NilValue;
}

SEXP gextract_multitask(SEXP _intervals, SEXP _exprs, SEXP _colnames, SEXP _iterator_policy, SEXP _band, SEXP _file, SEXP _intervals_set_out, SEXP _intervals_join, SEXP _envir)
{
	bool single_shard = false;

	try {
		// The parts and the temporary target of a Parquet output. Declared before the
		// RdbInitializer, whose destructor kills and reaps the kids, so that they are removed
		// only once no kid can still be writing one.
		TmpFiles parquet_tmp;
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

		IntervalsJoinMode join_mode = parse_intervals_join(_intervals_join);
		if (join_mode == IntervalsJoinMode::INTERVALS) {
			if (!Rf_isNull(_file))
				verror("intervals_join='intervals' cannot be combined with 'file'");
			if (!Rf_isNull(_intervals_set_out))
				verror("intervals_join='intervals' cannot be combined with 'intervals.set.out'");
		}

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

		// Parquet: each kid writes its shard to a part file of its own, and the parent merges the
		// parts in kid order, which is genome order (prepare4multitasking() plans the shards
		// along the sorted scope)
		if (filename && is_parquet_file(filename)) {
			int num_kids = iu.prepare4multitasking(_exprs, intervals1d, intervals2d, _iterator_policy, _band,
												   allow_multichrom_range_split_for_gextract(_iterator_policy, intervset_out));

			// one shard, or an empty scope (whose file C_gextract writes too)
			if (num_kids <= 1)
				throw rdb::SingleShard();

			// decided once, on the whole scope, so that every part has the same schema
			vector<string> names = parquet_value_names(_exprs, _colnames);
			bool int64_coords = parquet_int64_coords(iu, intervals1d);
			string created_by = parquet_created_by();
			ParquetOutput out = parquet_output(filename);
			vector<string> parts;

			for (int i = 0; i < num_kids; ++i)
				parts.push_back(out.tmp + "." + to_string(i));
			parquet_tmp.paths = parts;
			parquet_tmp.paths.push_back(out.tmp);

			if (iu.distribute_task(0, 0)) {  // child process
				TrackExprScanner scanner(iu);

				scanner.begin(_exprs, iu.get_kid_intervals1d(), iu.get_kid_intervals2d(), _iterator_policy, _band);
				write_parquet_rows(scanner, iu, parts[RdbInitializer::get_kid_idx()], filename, names, int64_coords, created_by);
				rreturn(R_NilValue);
			}

			ParquetWriter::merge(out.tmp, parts, filename);
			commit_parquet(out);
			parquet_tmp.paths.pop_back();  // out.tmp is now the target
			rreturn(R_NilValue);
		}

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
					outfile << CHAR(STRING_ELT(_exprs, iexpr));
				else
					outfile << CHAR(STRING_ELT(_colnames, iexpr));
			}
			outfile << "\n";

			// Use buffered writing with snprintf for better performance than iostream
			const size_t BUFFER_SIZE = 65536;  // 64KB buffer
			char buffer[BUFFER_SIZE];
			size_t buffer_pos = 0;

			for (; !scanner.isend(); scanner.next()) {
				std::string line;  // one output row; grows as needed (no fixed cap)
				char field[64];
				int n;

				if (scanner.get_iterator()->is_1d()) {
					const GInterval &interval = scanner.last_interval1d();
					line = iu.id2chrom(interval.chromid);
					n = snprintf(field, sizeof(field), "\t%" PRId64 "\t%" PRId64, interval.start, interval.end);
					line.append(field, n);
				} else {
					const GInterval2D &interval = scanner.last_interval2d();
					line = iu.id2chrom(interval.chromid1());
					n = snprintf(field, sizeof(field), "\t%" PRId64 "\t%" PRId64, interval.start1(), interval.end1());
					line.append(field, n);
					line += '\t';
					line += iu.id2chrom(interval.chromid2());
					n = snprintf(field, sizeof(field), "\t%" PRId64 "\t%" PRId64, interval.start2(), interval.end2());
					line.append(field, n);
				}

				for (unsigned iexpr = 0; iexpr < num_exprs; ++iexpr) {
					n = snprintf(field, sizeof(field), "\t%.15g", scanner.last_real(iexpr));
					line.append(field, n);
				}
				line += '\n';

				// Flush buffer if not enough space for this line
				if (buffer_pos + line.size() >= BUFFER_SIZE) {
					outfile.write(buffer, buffer_pos);
					buffer_pos = 0;
				}

				if (line.size() >= BUFFER_SIZE) {
					outfile.write(line.data(), line.size());  // row larger than buffer: write directly
				} else {
					memcpy(buffer + buffer_pos, line.data(), line.size());
					buffer_pos += line.size();
				}

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

		// IMPORTANT: do not add estimated-size mode gating here.
		// We intentionally rely on prepare4multitasking() + runtime allocation guards,
		// because estimate-based auto-disabling previously regressed gextract throughput.
		bool allow_multichrom_1d_range_split = allow_multichrom_range_split_for_gextract(_iterator_policy, intervset_out);
		int num_kids = iu.prepare4multitasking(_exprs, intervals1d, intervals2d, _iterator_policy, _band, allow_multichrom_1d_range_split);

		if (!num_kids)
			rreturn(R_NilValue);

		if (num_kids == 1)
			throw rdb::SingleShard();

		bool is_1d_iterator = iu.is_1d_iterator(_exprs, intervals1d, intervals2d, _iterator_policy);
		// Estimate number of records to cap shared-memory allocation size
		uint64_t base_estimated_records = estimate_records_for_expr(iu, _exprs, intervals1d, intervals2d, _iterator_policy, _band);
		double max_records_factor = iu.get_multitask_max_records_factor();
		uint64_t estimated_records = inflate_estimated_records(iu, base_estimated_records, max_records_factor);

		if (!intervset_out.empty()) {
			vector<GIntervalsBigSet1D::ChromStat> chromstats1d;
			GIntervalsMeta2D::ChromStats2D chromstats2d;

			if (is_1d_iterator)
				GIntervalsBigSet1D::begin_save(intervset_out.c_str(), iu, chromstats1d);
			else {
				GIntervalsBigSet2D::begin_save(intervset_out.c_str(), iu, chromstats2d);
				chromstats2d.set_max_pairs(intervals2d ? (size_t)intervals2d->num_chrom_pairs() : 0);
			}

			if (iu.distribute_task(is_1d_iterator ?
								   sizeof(GIntervalsBigSet1D::ChromStat) * chromstats1d.size() :
								   GIntervalsMeta2D::ChromStats2D::max_pairs_bytes(chromstats2d.max_pairs()),
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
							runprotect(rintervals); // one contig's worth of protect slots; the loop runs per contig
							out_intervals1d.clear();
							for (vector< vector<double> >::iterator ivalues = values.begin(); ivalues != values.end(); ++ivalues) 
								ivalues->clear();
						}
					} else {
						if (scanner.isend() || !last_scope_interval2d.is_same_chrom(scanner.last_scope_interval2d())) {
							SEXP rintervals = build_rintervals_extract(NULL, &out_intervals2d, values, NULL, _exprs, _colnames, iu);
							GIntervalsBigSet2D::save_chrom(intervset_out.c_str(), &out_intervals2d, rintervals, iu, chromstats2d);
							runprotect(rintervals); // one contig's worth of protect slots; the loop runs per contig
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
					chromstats2d.pack(ptr);
			} else { // parent process
				vector<GIntervalsBigSet1D::ChromStat> kid_chromstats1d(chromstats1d.size());
				GIntervalsMeta2D::ChromStats2D kid_chromstats2d;

				for (int i = 0; i < get_num_kids(); ++i) {
					void *ptr = get_kid_res(i);

					if (is_1d_iterator) {
						unpack_data(ptr, kid_chromstats1d.front(), kid_chromstats1d.size());
						for (vector<GIntervalsBigSet1D::ChromStat>::const_iterator istat = kid_chromstats1d.begin(); istat < kid_chromstats1d.end(); ++istat) {
							if (istat->size)
								chromstats1d[istat - kid_chromstats1d.begin()] = *istat;
						}
					} else {
						kid_chromstats2d.clear();
						kid_chromstats2d.unpack(ptr);
						chromstats2d.merge(kid_chromstats2d);
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
		for (;;) {
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
				break;
				} catch (TGLException &e) {
					string msg(e.msg());
					if (msg.find("Failed to allocate shared memory") != string::npos) {
						// Fallback policy: keep multitasking unless shared-memory allocation fails.
						// We first back off the inflation factor; only then fall back to serial.
						if (base_estimated_records > 0 && max_records_factor > 1.0) {
							max_records_factor = max(1.0, max_records_factor / 2.0);
							estimated_records = inflate_estimated_records(iu, base_estimated_records, max_records_factor);
							continue;
						}
					// Hand the serial run to the SingleShard tail below. Calling C_gextract()
					// from here would run it while this frame still holds its RdbInitializer,
					// and an error inside it would longjmp past our destructor (see
					// rdb::SingleShard). A kid must never return to R, so it falls through to
					// the rethrow and reports through the outer catch instead.
					if (!RdbInitializer::is_kid())
						throw rdb::SingleShard();
				}
				throw;
			}
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

				uint64_t total_intervals = 0;
				for (int i = 0; i < get_num_kids(); ++i)
					total_intervals += get_kid_res_size(i);

				if (!total_intervals)
					rreturn(R_NilValue);

				SEXP answer = build_answer_from_kids(iu, is_1d_iterator, num_exprs, total_intervals, _exprs, _colnames);

				if (profile) {
					auto t_end = std::chrono::steady_clock::now();
					log_gextract_timing("parent_gather_ms", std::chrono::duration<double, std::milli>(t_end - t_gather_start).count());
				}
				rreturn(transform_intervals_join(answer, _intervals, join_mode, iu));
		}
	} catch (rdb::SingleShard &) {
		// the serial entry point below re-queues whatever it warrants
		clear_pending_diagnostics(_envir);
		single_shard = true;
	} catch (TGLException &e) {
		string msg(e.msg());
		// Only the parent may re-run the call serially. In a child process this is
		// reached when the shared-memory result arena overflows (allocate_res()), and
		// returning from here would hand control back to R inside the fork: the child
		// would run the caller's R code and then exit through R's own shutdown, whose
		// R_CleanTempDir() deletes the session tempdir the parent is still using.
		// rerror() instead publishes the message to the parent, which retries.
		if (msg.find("Result size exceeded the maximal allowed") != string::npos && !RdbInitializer::is_kid()) {
			return C_gextract(_intervals, _exprs, _colnames, _iterator_policy, _band, _file, _intervals_set_out, _intervals_join, _envir);
		}
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	}
	// Multitasking is off the table: either the scope fits in a single shard (forking one
	// kid would buy no parallelism, only a fork, a shared-memory segment and a wait) or the
	// shared-memory arena could not be allocated. Run the serial entry point instead - the
	// very path gmultitasking = FALSE takes. This must happen out here, after the try
	// block has destroyed our RdbInitializer (see rdb::SingleShard).
	if (single_shard)
		return C_gextract(_intervals, _exprs, _colnames, _iterator_policy, _band, _file, _intervals_set_out, _intervals_join, _envir);

	rreturn(R_NilValue);
}

}
