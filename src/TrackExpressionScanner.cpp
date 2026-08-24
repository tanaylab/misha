/*
 * TrackExpressionScanner.cpp
 *
 *  Created on: Apr 28, 2010
 *      Author: hoichman
 */

#include <cstdint>
#include <errno.h>
#include <sys/time.h>

#include <climits>
#include <limits>
#include <set>
#include <unordered_set>

#include "port.h"

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Parse.h>

#include "HashFunc.h"

#include "GenomeTrack.h"
#include "GenomeTrackFixedBin.h"
#include "TrackExpressionBigSet1DIterator.h"
#include "TrackExpressionBigSet2DIterator.h"
#include "TrackExpressionCartesianGridIterator.h"
#include "TrackExpressionFixedBinIterator.h"
#include "TrackExpressionFixedRectIterator.h"
#include "TrackExpressionIntervals1DIterator.h"
#include "TrackExpressionIntervals2DIterator.h"
#include "TrackExpressionIterator.h"
#include "TrackExpressionSparseIterator.h"
#include "TrackExpressionScanner.h"
#include "TrackExpressionTrackRectsIterator.h"
#include "TrackExpressionVars.h"
#include "rdbutils.h"
#include "TrackIndex.h"
#include "TrackIndex2D.h"

const int TrackExprScanner::INIT_REPORT_STEP = 10000;
// Default vectorized evaluation buffer size (the gbuf.size option's documented default).
const unsigned TrackExprScanner::DEFAULT_EVAL_BUF_SIZE = 1000;
const int TrackExprScanner::REPORT_INTERVAL = 3000;
const int TrackExprScanner::MIN_REPORT_INTERVAL = 1000;

using namespace rdb;

static uint64_t get_cur_clock()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

TrackExprScanner::TrackExprScanner(rdb::IntervUtils &iu) :
	m_iu(iu),	
	m_isend(true),
	m_expr_itr(NULL),
	m_expr_vars(iu)
{
	m_1d.cur_chromid = -1;
	m_2d.cur_chromid1 = -1;
	m_2d.cur_chromid2 = -1;
	m_do_report_progress = true;
}

TrackExprScanner::~TrackExprScanner()
{
	delete m_expr_itr;
}

void TrackExprScanner::convert_rtrack_exprs(SEXP rtrack_exprs, vector<string> &track_exprs)
{
	track_exprs.clear();

	if (!Rf_isString(rtrack_exprs) || Rf_length(rtrack_exprs) < 1)
		verror("Tracks expressions argument must be a vector of strings");

	unsigned num_track_exprs = (unsigned)Rf_length(rtrack_exprs);
	track_exprs.resize(num_track_exprs);

	for (unsigned iexpr = 0; iexpr < num_track_exprs; ++iexpr)
		track_exprs[iexpr] = CHAR(STRING_ELT(rtrack_exprs, iexpr));
}

// Size of the vectorized evaluation buffer, from the gbuf.size option.
// Unset falls back to the documented default; anything set but unusable is an error,
// because silently substituting the default hides a value the user asked for and
// casting it truncates (see the call site in begin()).
unsigned TrackExprScanner::get_eval_buf_size()
{
	SEXP gbufsize = Rf_GetOption1(Rf_install("gbuf.size"));

	if (gbufsize == R_NilValue)
		return DEFAULT_EVAL_BUF_SIZE;

	double val;

	if (Rf_isReal(gbufsize) && Rf_length(gbufsize) == 1)
		val = REAL(gbufsize)[0];
	else if (Rf_isInteger(gbufsize) && Rf_length(gbufsize) == 1)
		val = INTEGER(gbufsize)[0] == NA_INTEGER ? NA_REAL : (double)INTEGER(gbufsize)[0];
	else
		verror("gbuf.size option must be a single number");

	// INT_MAX, not UINT_MAX: past it R's own vectors become "long" ones, and both this
	// class and R's Rf_length() index the evaluation buffer with an int. A value R can
	// still allocate but not attribute used to fail with a raw Rf_error from attrib.c,
	// which is the same longjmp-past-~RdbInitializer this validation exists to prevent.
	if (!R_FINITE(val) || val < 1 || val > (double)INT_MAX)
		verror("gbuf.size option must be a finite number between 1 and %d", INT_MAX);

	return (unsigned)val;
}

void TrackExprScanner::define_r_vars(unsigned eval_buf_limit)
{
	m_eval_buf_limit = eval_buf_limit;
	m_expr_vars.define_r_vars(eval_buf_limit);

	// Prepare R-visible buffers that mirror the current chunk of iterator intervals.
	// For 1D we expose CHROM/START/END; for 2D CHROM1/START1/END1/CHROM2/START2/END2.
	// Chromosome ids are 1-based on the R side. We also preinitialize the entire
	// buffer with safe defaults so that partially filled buffers (at EOF) remain valid.

	if (m_expr_itr->is_1d()) {
		m_1d.cur_chromid = -1;
		m_1d.expr_itr_intervals.resize(m_eval_buf_limit);
		m_1d.expr_itr_scope_intervals.resize(m_eval_buf_limit);
		m_rexpr_itr_intervals = m_iu.convert_intervs(&m_1d.expr_itr_intervals);
		m_1d.expr_itr_intervals_chroms = INTEGER(VECTOR_ELT(m_rexpr_itr_intervals, GInterval::CHROM));
		m_1d.expr_itr_intervals_starts = REAL(VECTOR_ELT(m_rexpr_itr_intervals, GInterval::START));
		m_1d.expr_itr_intervals_ends = REAL(VECTOR_ELT(m_rexpr_itr_intervals, GInterval::END));
		for (unsigned i = 0; i < m_eval_buf_limit; i++)
			m_1d.expr_itr_intervals_chroms[i] = 1;
    } else {
		m_2d.cur_chromid1 = -1;
		m_2d.cur_chromid2 = -1;
		m_2d.expr_itr_intervals.resize(m_eval_buf_limit);
		m_2d.expr_itr_scope_intervals.resize(m_eval_buf_limit);
		m_rexpr_itr_intervals = m_iu.convert_intervs(&m_2d.expr_itr_intervals);
		m_2d.expr_itr_intervals_chroms1 = INTEGER(VECTOR_ELT(m_rexpr_itr_intervals, GInterval2D::CHROM1));
		m_2d.expr_itr_intervals_starts1 = REAL(VECTOR_ELT(m_rexpr_itr_intervals, GInterval2D::START1));
		m_2d.expr_itr_intervals_ends1 = REAL(VECTOR_ELT(m_rexpr_itr_intervals, GInterval2D::END1));
		m_2d.expr_itr_intervals_chroms2 = INTEGER(VECTOR_ELT(m_rexpr_itr_intervals, GInterval2D::CHROM2));
		m_2d.expr_itr_intervals_starts2 = REAL(VECTOR_ELT(m_rexpr_itr_intervals, GInterval2D::START2));
		m_2d.expr_itr_intervals_ends2 = REAL(VECTOR_ELT(m_rexpr_itr_intervals, GInterval2D::END2));
		for (unsigned i = 0; i < m_eval_buf_limit; i++)
			m_2d.expr_itr_intervals_chroms1[i] = m_2d.expr_itr_intervals_chroms2[i] = 1;
    }
	// Expose the current-batch intervals to R under GITERATOR.INTERVALS in the .misha env
    define_in_misha(m_iu.get_env(), "GITERATOR.INTERVALS", m_rexpr_itr_intervals);

    for (unsigned iexpr = 0; iexpr < m_track_exprs.size(); ++iexpr) {
        const TrackExpressionVars::Track_var *var = m_expr_vars.var(m_track_exprs[iexpr].c_str());
        if (var)  // track expression is a virtual track
            m_eval_doubles[iexpr] = REAL(var->rvar);
    }
}

void TrackExprScanner::check(string &track_expr, GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d, SEXP iterator_policy, SEXP band)
{
	check(vector<string>(1, track_expr), scope1d, scope2d, iterator_policy, band);
}

void TrackExprScanner::check(SEXP track_exprs, GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d, SEXP iterator_policy, SEXP band)
{
	vector<string> track_expr_strs;
	convert_rtrack_exprs(track_exprs, track_expr_strs);
	check(track_expr_strs, scope1d, scope2d, iterator_policy, band);
}

void TrackExprScanner::check(const vector<string> &track_exprs, GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d, SEXP iterator_policy, SEXP band)
{
	m_band = m_iu.convert_band(band);

    runprotect(m_eval_bufs);
    runprotect(m_eval_exprs);

    m_track_exprs.reserve(track_exprs.size());
	// Normalize input expressions: trim whitespace and store canonical forms.
    for (vector<string>::const_iterator iexpr = track_exprs.begin(); iexpr != track_exprs.end(); ++iexpr) {
        // trim spaces from the track expression
        string::const_iterator istr_start, istr_end;

        for (istr_start = iexpr->begin(); istr_start < iexpr->end(); ++istr_start) {
            if (!isspace(*istr_start))
                break;
        }

        for (istr_end = iexpr->end() - 1; istr_end >= iexpr->begin(); --istr_end) {
            if (!isspace(*istr_end))
                break;
        }

        m_track_exprs.push_back(iexpr->substr(istr_start - iexpr->begin(), istr_end - istr_start + 1));
    }

    m_eval_exprs.resize(m_track_exprs.size(), R_NilValue);
    m_eval_bufs.resize(m_track_exprs.size(), R_NilValue);
    m_eval_doubles.resize(m_track_exprs.size(), NULL);
    m_eval_ints.resize(m_track_exprs.size(), NULL);

	m_expr_vars.parse_exprs(m_track_exprs);

	// initiate the expression iterator
	delete m_expr_itr;
	m_expr_itr = create_expr_iterator(iterator_policy, m_expr_vars, m_track_exprs, scope1d, scope2d, m_1d.intervals, m_2d.intervals, m_band);
	m_expr_vars.init(*m_expr_itr);

	for (unsigned iexpr = 0; iexpr < m_track_exprs.size(); ++iexpr) {
        if (!m_expr_vars.var(m_track_exprs[iexpr].c_str())) {   // track expression is not a virtual track
    		SEXP expr = R_NilValue;
			SEXPCleaner expr_cleaner(expr);

			expr = rprotect_ptr(RSaneAllocVector(STRSXP, 1));
    		SET_STRING_ELT(expr, 0, Rf_mkChar(m_track_exprs[iexpr].c_str()));

			// Parse the user R expression once and keep the parsed form for repeated evaluation
    		ParseStatus status;
			SEXP parsed_expr;
			parsed_expr = rprotect_ptr(R_ParseVector(expr, -1, &status, R_NilValue));
			if (status != PARSE_OK)
    			verror("R parsing of expression \"%s\" failed", m_track_exprs[iexpr].c_str());
			m_eval_exprs[iexpr] = VECTOR_ELT(parsed_expr, 0);
        }
	}
}

bool TrackExprScanner::begin(string &track_expr, GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d, SEXP iterator_policy, SEXP band)
{
	return begin(vector<string>(1, track_expr), scope1d, scope2d, iterator_policy, band);
}

bool TrackExprScanner::begin(SEXP track_exprs, GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d, SEXP iterator_policy, SEXP band)
{
	vector<string> track_expr_strs;
	convert_rtrack_exprs(track_exprs, track_expr_strs);
	return begin(track_expr_strs, scope1d, scope2d, iterator_policy, band);
}

bool TrackExprScanner::begin(const vector<string> &track_exprs, GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d, SEXP iterator_policy, SEXP band)
{
	check(track_exprs, scope1d, scope2d, iterator_policy, band);

	// Decide evaluation buffer size. Prefer vectorized evaluation (large buffers)
	// but gracefully fall back to scalar evaluation (buffer size 1) if expressions
	// cannot produce vectors of the requested size.
    // gbuf.size is a user option, so it must be validated rather than cast. A plain
	// (unsigned)REAL(...)[0] turns 2^32 - or NA / Inf - into a buffer of 0; the scanner
	// then evaluates nothing, convert_intervs() hands back NULL and the caller's
	// VECTOR_ELT(NULL, ...) raises a raw Rf_error that longjmps past ~RdbInitializer,
	// leaving misha's SIGINT handler installed and Ctrl-C dead for the rest of the
	// session. Reject what cannot be honoured instead.
	define_r_vars(get_eval_buf_size());

    for (unsigned iexpr = 0; iexpr < m_track_exprs.size(); ++iexpr) {
        if (m_eval_exprs[iexpr] != R_NilValue) {
            SEXP res = eval_in_R(m_eval_exprs[iexpr], m_iu.get_env());

            // Name it, do not count. eval_in_R protects through rprotect(), which is
            // a no-op on R_NilValue - so an expression that evaluates to NULL
            // ("NULL", "if (FALSE) 1", any invisible(NULL)) protects nothing and a
            // runprotect(1) here popped whatever define_r_vars() had just protected:
            // the iterator-intervals data frame whose column pointers this scanner
            // caches. That is a use-after-free, and it segfaulted.
            if (Rf_length(res) != (int)m_eval_buf_limit) {
                runprotect(res);
				define_r_vars(1); // switch to scalar mode
                break;
            }
            runprotect(res);
        }
    }

	m_num_evals = 0;
	m_last_progress_reported = -1;
	m_report_step = INIT_REPORT_STEP;
	m_last_report_clock = get_cur_clock();

	m_isend = false;
	m_eval_buf_idx = m_eval_buf_limit;
	m_eval_buf_size = 0;
	m_expr_itr_scope_idx.resize(m_eval_buf_limit, -1);
	m_expr_itr_scope_chrom_idx.resize(m_eval_buf_limit, -1);

	return next();
}

bool TrackExprScanner::next()
{
	monitor_memusage();

	if (isend())
		return false;

	// Advance within the current buffer; if depleted, refill from the iterator
	if (eval_next())
		return true;

	// did we start reporting progress?
	if (m_last_progress_reported >= 0) {
		if (m_last_progress_reported != 100)
			REprintf("100%%\n");
		else
			REprintf("\n");
	}
	update_progress(100);

	runprotect(m_eval_bufs);

	return false;
}

bool TrackExprScanner::eval_next()
{
	m_eval_buf_idx++;

	// Refill the evaluation buffer when exhausted (batch processing for speed)
	if (m_eval_buf_idx >= m_eval_buf_limit) {
		m_eval_buf_idx = 0;
		if (m_expr_itr->is_1d()) {
			for (m_eval_buf_size = 0; m_eval_buf_size < m_eval_buf_limit; ++m_eval_buf_size) {
                if (m_expr_itr->isend()) {
					// Pad the remainder with sentinel values so R-side vectors remain valid
					for (unsigned i = m_eval_buf_size; i < m_eval_buf_limit; ++i) {
                        m_1d.expr_itr_intervals_chroms[i] = 1;
                        m_1d.expr_itr_intervals_starts[i] = -1.;
                        m_1d.expr_itr_intervals_ends[i] = -1.;
                    }
                    break;
                }

                const GInterval &interval = ((TrackExpression1DIterator *)m_expr_itr)->last_interval();

                m_1d.expr_itr_intervals[m_eval_buf_size] = interval;
                m_1d.expr_itr_scope_intervals[m_eval_buf_size] = ((TrackExpression1DIterator *)m_expr_itr)->last_scope_interval();
                m_1d.expr_itr_intervals_chroms[m_eval_buf_size] = interval.chromid + 1;
                m_1d.expr_itr_intervals_starts[m_eval_buf_size] = (double)interval.start;
                m_1d.expr_itr_intervals_ends[m_eval_buf_size] = (double)interval.end;
                m_1d.cur_chromid = interval.chromid;

                m_expr_itr_scope_idx[m_eval_buf_size] = ((TrackExpression1DIterator *)m_expr_itr)->get_cur_scope_idx();
                m_expr_itr_scope_chrom_idx[m_eval_buf_size] = ((TrackExpression1DIterator *)m_expr_itr)->get_cur_scope_chrom_idx();
				// Publish iterator-derived variables for this slot so expressions can use them
				m_expr_vars.set_vars(interval, m_eval_buf_size);
                m_expr_itr->next();
			}
		}
		else {
			for (m_eval_buf_size = 0; m_eval_buf_size < m_eval_buf_limit; ++m_eval_buf_size) {
				if (m_expr_itr->isend()) {
					// Pad the remainder with sentinel values in 2D as well
					for (unsigned i = m_eval_buf_size; i < m_eval_buf_limit; ++i) {
						m_2d.expr_itr_intervals_chroms1[i] = 1;
						m_2d.expr_itr_intervals_starts1[i] = -1.;
						m_2d.expr_itr_intervals_ends1[i] = -1.;
						m_2d.expr_itr_intervals_chroms2[i] = 1;
						m_2d.expr_itr_intervals_starts2[i] = -1.;
						m_2d.expr_itr_intervals_ends2[i] = -1.;
					}
					break;
				}

				const GInterval2D &interval = ((TrackExpression2DIterator *)m_expr_itr)->last_interval();

				m_2d.expr_itr_intervals[m_eval_buf_size] = interval;
				m_2d.expr_itr_scope_intervals[m_eval_buf_size] = ((TrackExpression2DIterator *)m_expr_itr)->last_scope_interval();
				m_2d.expr_itr_intervals_chroms1[m_eval_buf_size] = interval.chromid1() + 1;
				m_2d.expr_itr_intervals_starts1[m_eval_buf_size] = (double)interval.start1();
				m_2d.expr_itr_intervals_ends1[m_eval_buf_size] = (double)interval.end1();
				m_2d.expr_itr_intervals_chroms2[m_eval_buf_size] = interval.chromid2() + 1;
				m_2d.expr_itr_intervals_starts2[m_eval_buf_size] = (double)interval.start2();
				m_2d.expr_itr_intervals_ends2[m_eval_buf_size] = (double)interval.end2();
				m_2d.cur_chromid1 = interval.chromid1();
				m_2d.cur_chromid2 = interval.chromid2();

				m_expr_itr_scope_idx[m_eval_buf_size] = ((TrackExpression2DIterator *)m_expr_itr)->get_cur_scope_idx();
				m_expr_itr_scope_chrom_idx[m_eval_buf_size] = ((TrackExpression2DIterator *)m_expr_itr)->get_cur_scope_chrom_idx();
				m_expr_vars.set_vars(interval, m_band, m_eval_buf_size);
				m_expr_itr->next();
			}
		}

        check_interrupt();

		for (unsigned iexpr = 0; iexpr < m_eval_exprs.size(); ++iexpr) {
            if (m_eval_exprs[iexpr] != R_NilValue) {
                runprotect(m_eval_bufs[iexpr]);
				// Evaluate the parsed R expression over the current batch; only numeric or logical vectors are supported
				m_eval_bufs[iexpr] = eval_in_R(m_eval_exprs[iexpr], m_iu.get_env());
                if (Rf_length(m_eval_bufs[iexpr]) != (int)m_eval_buf_limit)
                    verror("Evaluation of expression \"%s\" produces a vector of size %d while expecting size %d",
                            m_track_exprs[iexpr].c_str(), Rf_length(m_eval_bufs[iexpr]), m_eval_buf_limit);
                if (Rf_isReal(m_eval_bufs[iexpr]))
                    m_eval_doubles[iexpr] = REAL(m_eval_bufs[iexpr]);
                else if (Rf_isLogical(m_eval_bufs[iexpr]))
                    m_eval_ints[iexpr] = LOGICAL(m_eval_bufs[iexpr]);
                else
                    verror("Evaluation of expression \"%s\" produces a vector of unsupported type %s",
                            m_track_exprs[iexpr].c_str(), Rf_type2char(TYPEOF(m_eval_bufs[iexpr])));
            }
        }

		report_progress();
	}

	if (m_eval_buf_idx >= m_eval_buf_size) {
		m_eval_buf_idx = m_eval_buf_limit;
		m_isend = true;
	}

	return !m_isend;
}

void TrackExprScanner::report_progress()
{
	m_num_evals += m_eval_buf_size;
	if (m_num_evals > (uint64_t)m_report_step && m_do_report_progress) {
        uint64_t curclock = get_cur_clock();
        double delta = curclock - m_last_report_clock;

		// Adapt the next reporting step to aim for roughly REPORT_INTERVAL ms between prints
		if (delta)
            m_report_step = (int)(m_report_step * (REPORT_INTERVAL / delta) + .5);
        else
            m_report_step *= 10;

        if (delta > MIN_REPORT_INTERVAL) {
			if (m_last_progress_reported < 0 && m_eval_buf_limit == 1)
                REprintf("Warning: track expression(s) cannot be evaluated as a vector. Run-times might be slow.\n");

            int progress;

            if (m_expr_itr->is_1d())
                progress = !((TrackExpression1DIterator *)m_expr_itr)->get_scope()->size() ? 0 :
                    (int)(100. * (((TrackExpression1DIterator *)m_expr_itr)->get_cur_scope_idx()) / ((TrackExpression1DIterator *)m_expr_itr)->get_scope()->size());
            else
                progress = !((TrackExpression2DIterator *)m_expr_itr)->get_scope()->size() ? 0 :
                    (int)(100. * (((TrackExpression2DIterator *)m_expr_itr)->get_cur_scope_idx()) / ((TrackExpression2DIterator *)m_expr_itr)->get_scope()->size());

			// In 2D the scope index may momentarily decrease; ensure monotonic progress display
			progress = max(progress, m_last_progress_reported);
			if (progress > 100) progress = 100;
            if (progress != 100) {
                if (progress != m_last_progress_reported) {
                    REprintf("%d%%...", progress);
                    update_progress(progress);
                } else
                    REprintf(".");
                m_last_progress_reported = progress;
            }
            m_num_evals = 0;
            m_last_report_clock = curclock;
        }
	}
}

void TrackExprScanner::verify_1d_iter(GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d) const
{
	if (!scope1d)
		verror("The function does not support 1D iterators");

	if (scope1d && scope2d && !scope1d->size() && scope2d->size())
		verror("1D iterator is used along with 2D intervals");
}

void TrackExprScanner::verify_2d_iter(GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d) const
{
	if (!scope2d)
		verror("The function does not support 2D iterators");

	if (scope1d && scope2d && scope1d->size() && !scope2d->size())
		verror("2D iterator is used along with 1D intervals");
}

TrackExpressionIteratorBase *TrackExprScanner::create_expr_iterator(SEXP rtrack_exprs, GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d,
																	SEXP iterator_policy, SEXP band, bool call_begin)
{
	m_track_exprs.resize(Rf_length(rtrack_exprs));
	for (int i = 0; i < Rf_length(rtrack_exprs); ++i) 
		m_track_exprs[i] = CHAR(STRING_ELT(rtrack_exprs, i));

	m_band = m_iu.convert_band(band);
	m_expr_vars.parse_exprs(m_track_exprs);

	// initiate the expression iterator
	delete m_expr_itr;
	m_expr_itr = create_expr_iterator(iterator_policy, m_expr_vars, m_track_exprs, scope1d, scope2d, m_1d.intervals, m_2d.intervals, m_band, call_begin);
    return m_expr_itr;
}

//---------------------------- overlapping 1D iterator intervals --------------------------------
//
// A 1D intervals iterator is unified below, so overlapping intervals are folded
// into maximal blocks and the result carries one row per block instead of one per
// interval. The 2D branch errors out on the same condition, but erroring here would
// break every script that iterates over sliding windows or unioned peak calls, so
// the merge stays and the caller gets a warning instead.
//
// The merge is not always observable: the emitted rows are (block x scope)
// intersections, so when the same intervals are also passed as the scope, every input
// interval is clipped back to itself and nothing is lost. The merge is unobservable
// exactly when every input interval the scope reaches is itself one of the emitted
// rows, and that is what is tested here - the search stops at the first interval for
// which it fails.
//
// Two shapes make it fail, and both are looked for:
//   - an emitted row that no single input interval contains: what a partial overlap
//     produces (chr1 0-300 and chr1 100-400 become chr1 0-400 under a wider scope).
//   - an input interval that is never emitted as a row of its own even though the
//     scope reaches it: what nesting produces (chr1 0-1000 and chr1 100-300 become
//     chr1 0-1000, a row the first input *does* contain, so the test above sees
//     nothing while the second peak has silently disappeared).
//
// The second test is what keeps the recommended per-interval idiom (the same
// intervals as both scope and iterator) silent: there every input interval is clipped
// by its own scope interval back to itself, nested inputs included.
//
// What still goes unreported is an exact duplicate under a wider scope: both copies
// are emitted verbatim, as one shared row, so no interval is missing from the output
// even though the row count dropped. Erroring on any overlap at all, under an opt-in
// strict mode, is the complete answer to that.

// intervals must be sorted
static bool intervals_overlap(const vector<GInterval> &intervals)
{
	for (size_t i = 1; i < intervals.size(); ++i) {
		if (intervals[i].chromid == intervals[i - 1].chromid && intervals[i].start < intervals[i - 1].end)
			return true;
	}
	return false;
}

// max_end[i] is the maximal end coordinate among the intervals of the same chromosome up to i
static bool is_row_covered(const vector<GInterval> &intervals, const vector<int64_t> &max_end, const GInterval &row)
{
	// the last interval that starts at or before the row (GInterval::operator< orders by
	// chromosome and then by start coordinate, which is how the intervals were sorted)
	vector<GInterval>::const_iterator iinterv = upper_bound(intervals.begin(), intervals.end(), row);

	if (iinterv == intervals.begin())
		return false;

	--iinterv;
	if (iinterv->chromid != row.chromid)
		return false;

	return max_end[iinterv - intervals.begin()] >= row.end;
}

// block_first[k] is the index of the first interval of unified[k], and block_first[n] is
// intervals.size(). The blocks are the maximal merged runs of the sorted intervals, so the
// intervals of a block are a contiguous range and every interval belongs to exactly one.
static void index_blocks(const vector<GInterval> &intervals, const vector<GInterval> &unified,
						 vector<size_t> &block_first)
{
	block_first.assign(unified.size() + 1, intervals.size());

	size_t iinterv = 0;
	for (size_t iblock = 0; iblock < unified.size(); ++iblock) {
		block_first[iblock] = iinterv;
		while (iinterv < intervals.size() && intervals[iinterv].chromid == unified[iblock].chromid &&
			   intervals[iinterv].start < unified[iblock].end)
			++iinterv;
	}
}

// The intervals that start before a scope interval and still reach into it: a stabbing query
// at the scope interval's start, which a prefix maximum cannot answer in less than the whole
// prefix. One interval spanning a block keeps max_end high across all of it, so the scan
// below starts at the beginning of the block and every scope interval re-walks the same
// intervals - quadratic in the size of the block, which is the shape a container interval
// with nested peaks inside it produces.
//
// This tree answers the same query in runs of consecutive indices rather than one index at
// a time, so its cost follows the number of runs the answer breaks into and not the number
// of intervals it reports. That is what keeps the dense case cheap: an overlapping window
// wider than the cap reaches back over an unbroken run of its predecessors, which is one
// descent plus the same sequential scan the walk itself would have done, while the nested
// container that motivated the tree reports a single interval from a block of millions.
//
// It is built on demand: a short scan finds the answer on its own whenever the reaching
// intervals sit near the scope interval (sliding windows, unioned peak calls), and a
// million-interval iterator should not pay for a tree it never queries.
namespace {

// a run of consecutive interval indices, all of which answer the query
struct IdxRange {
	size_t lo;
	size_t hi;
};

class EndTree {
public:
	EndTree() : m_built(false), m_size(0) {}

	void build(const vector<GInterval> &intervals) {
		if (m_built)
			return;
		m_built = true;

		m_size = 1;
		while (m_size < intervals.size())
			m_size <<= 1;

		// the two bounds of a node sit together, so testing both costs one cache line
		Node empty = { numeric_limits<int64_t>::max(), numeric_limits<int64_t>::min() };
		m_tree.assign(2 * m_size, empty);
		for (size_t i = 0; i < intervals.size(); ++i) {
			m_tree[m_size + i].min_end = intervals[i].end;
			m_tree[m_size + i].max_end = intervals[i].end;
		}
		for (size_t i = m_size - 1; i >= 1; --i) {
			m_tree[i].min_end = min(m_tree[2 * i].min_end, m_tree[2 * i + 1].min_end);
			m_tree[i].max_end = max(m_tree[2 * i].max_end, m_tree[2 * i + 1].max_end);
		}
	}

	// appends the maximal runs of consecutive indices in [lo, hi) whose interval ends after
	// `thr`, in ascending order
	void ends_after(size_t lo, size_t hi, int64_t thr, vector<IdxRange> &out) const {
		if (lo < hi)
			descend(1, 0, m_size, lo, hi, thr, out);
	}

private:
	struct Node {
		int64_t min_end;
		int64_t max_end;
	};

	void descend(size_t node, size_t node_lo, size_t node_hi, size_t lo, size_t hi, int64_t thr, vector<IdxRange> &out) const {
		if (node_hi <= lo || hi <= node_lo || m_tree[node].max_end <= thr)
			return;
		// every interval under this node answers the query, and the node is inside the queried
		// range: report the whole run without descending to its leaves. Padding leaves are
		// never inside the range, so their bounds cannot reach this test.
		if (lo <= node_lo && node_hi <= hi && m_tree[node].min_end > thr) {
			emit(node_lo, node_hi, out);
			return;
		}
		// a qualifying leaf is always covered by the test above; this is only a guard against
		// descending past the bottom of the tree
		if (node_hi - node_lo == 1) {
			emit(node_lo, node_hi, out);
			return;
		}

		size_t node_mid = node_lo + (node_hi - node_lo) / 2;
		descend(2 * node, node_lo, node_mid, lo, hi, thr, out);
		descend(2 * node + 1, node_mid, node_hi, lo, hi, thr, out);
	}

	// the descent visits leaves left to right, so a run split across several nodes arrives as
	// adjacent pieces and is glued back into one
	static void emit(size_t lo, size_t hi, vector<IdxRange> &out) {
		if (!out.empty() && out.back().hi == lo)
			out.back().hi = hi;
		else {
			IdxRange range = { lo, hi };
			out.push_back(range);
		}
	}

	bool           m_built;
	size_t         m_size;
	vector<Node>   m_tree;
};

} // anonymous namespace

// How many intervals the scan is willing to walk over before it asks the tree instead.
// Sliding windows and unioned peak calls reach back a handful of intervals - a window
// reaches back width/step of them - and walking is cheaper than a tree descent at that
// size, so the threshold sits well above any realistic reach. Past it the tree answers in
// runs, so a window that does reach back further than this is not punished for it.
static const size_t MAX_WALKED_INTERVALS = 128;

// The emitted row that swallowed intervals[iinterv], found by re-walking the scope for the
// first interval that reaches it. Only ever called once, on the way to raising the warning.
static bool find_swallowing_row(const vector<GInterval> &intervals, const vector<GInterval> &unified,
								const vector<size_t> &block_first, GIntervalsFetcher1D *scope,
								size_t iinterv, GInterval &row, GInterval &block)
{
	size_t iblock = upper_bound(block_first.begin(), block_first.end(), iinterv) - block_first.begin() - 1;
	const GInterval &interv = intervals[iinterv];

	block = unified[iblock];

	for (scope->begin_iter(); !scope->isend(); scope->next()) {
		const GInterval &scope_interv = scope->cur_interval();

		if (scope_interv.chromid != interv.chromid || scope_interv.end <= interv.start || scope_interv.start >= interv.end)
			continue;

		row = GInterval(interv.chromid, max(block.start, scope_interv.start), min(block.end, scope_interv.end), 0);
		return true;
	}
	return false;
}

// Returns the first emitted row that cost the caller something - either a row no single
// iterator interval contains, or a row that swallowed an interval the scope reached and that
// is emitted nowhere on its own - along with the block it came from.
static bool find_merged_row(const vector<GInterval> &intervals, const vector<GInterval> &unified,
							GIntervalsFetcher1D *scope, GInterval &row, GInterval &block)
{
	vector<int64_t> max_end(intervals.size());

	for (size_t i = 0; i < intervals.size(); ++i)
		max_end[i] = i && intervals[i - 1].chromid == intervals[i].chromid ? max(max_end[i - 1], intervals[i].end) : intervals[i].end;

	vector<size_t> block_first;
	index_blocks(intervals, unified, block_first);

	// An input interval is emitted verbatim when some scope interval clips it and its block to
	// the same thing; the merge then cost it nothing. Both flags are decided over the whole
	// scope, because one scope interval emitting an interval verbatim is enough even if
	// another one swallows it - which is exactly what the recommended idiom looks like.
	vector<bool> reached(intervals.size(), false);
	vector<bool> verbatim(intervals.size(), false);

	EndTree tree;
	vector<IdxRange> reaching;

	for (scope->begin_iter(); !scope->isend(); scope->next()) {
		const GInterval &scope_interv = scope->cur_interval();
		vector<GInterval>::const_iterator iblock = lower_bound(unified.begin(), unified.end(), scope_interv);

		// the preceding block might start before the scope interval and still reach into it
		if (iblock != unified.begin() && (iblock - 1)->chromid == scope_interv.chromid && (iblock - 1)->end > scope_interv.start)
			--iblock;

		for (; iblock != unified.end() && iblock->chromid == scope_interv.chromid && iblock->start < scope_interv.end; ++iblock) {
			GInterval emitted(scope_interv.chromid, max(iblock->start, scope_interv.start), min(iblock->end, scope_interv.end), 0);

			if (emitted.start >= emitted.end)
				continue;

			if (!is_row_covered(intervals, max_end, emitted)) {
				row = emitted;
				block = *iblock;
				return true;
			}

			size_t block_lo = block_first[iblock - unified.begin()];
			size_t block_hi = block_first[iblock - unified.begin() + 1];

			// the first interval of the block that begins at or after the scope interval
			// (GInterval::operator< orders by chromosome and then by start coordinate, which is
			// how the intervals were sorted)
			GInterval probe(scope_interv.chromid, scope_interv.start, scope_interv.start, 0);
			size_t istart = lower_bound(intervals.begin() + block_lo, intervals.begin() + block_hi, probe) - intervals.begin();

			// How far back the scan has to start to catch the intervals that begin before the
			// scope interval and still reach into it: max_end is a prefix maximum, so the walk
			// back stops as soon as nothing earlier can reach.
			//
			// A single interval spanning the block holds max_end high over all of it however far
			// away the scope interval is, and the walk would then reach the start of the block for
			// every scope interval that touches it - quadratic, and a container interval with peaks
			// nested inside it is exactly that shape. Being a prefix maximum is also what says in
			// advance whether that is about to happen: max_end does not decrease inside a block, so
			// if the interval MAX_WALKED_INTERVALS back still reaches, so does everything between,
			// and the walk is going to run at least that far. One look there replaces the walk with
			// a tree query, which reports the same intervals in runs.
			size_t i = istart;

			if (istart - block_lo > MAX_WALKED_INTERVALS && max_end[istart - MAX_WALKED_INTERVALS - 1] > scope_interv.start) {
				tree.build(intervals);
				reaching.clear();
				tree.ends_after(block_lo, istart, scope_interv.start, reaching);
				for (size_t k = 0; k < reaching.size(); ++k) {
					for (size_t ireach = reaching[k].lo; ireach < reaching[k].hi; ++ireach) {
						reached[ireach] = true;
						if (intervals[ireach].start <= emitted.start && intervals[ireach].end >= emitted.end)
							verbatim[ireach] = true;
					}
				}
			} else {
				// bounded by the look above: nothing MAX_WALKED_INTERVALS back reaches, or the block
				// does not extend that far
				while (i > block_lo && max_end[i - 1] > scope_interv.start)
					--i;
			}

			for (; i < block_hi && intervals[i].start < scope_interv.end; ++i) {
				if (intervals[i].end <= scope_interv.start)
					continue;

				reached[i] = true;
				// intervals[i] lies inside the block, so containing the emitted row is the
				// same as being clipped to it by this scope interval
				if (intervals[i].start <= emitted.start && intervals[i].end >= emitted.end)
					verbatim[i] = true;
			}
		}
	}

	// no row was widened, but an interval the scope reached may still have been swallowed whole
	for (size_t i = 0; i < intervals.size(); ++i) {
		if (reached[i] && !verbatim[i])
			return find_swallowing_row(intervals, unified, block_first, scope, i, row, block);
	}
	return false;
}

// The first pair of overlapping intervals, optionally restricted to those inside `block`.
static bool first_overlapping_pair(const vector<GInterval> &intervals, const GInterval *block,
								   GInterval &interv1, GInterval &interv2)
{
	int64_t max_end = -1;
	size_t imax_end = 0;
	bool started = false;

	for (size_t i = 0; i < intervals.size(); ++i) {
		if (block) {
			if (intervals[i].chromid != block->chromid || intervals[i].end <= block->start)
				continue;
			if (intervals[i].start >= block->end)
				break;
		}

		bool same_chrom = started && intervals[i].chromid == intervals[imax_end].chromid;

		if (same_chrom && intervals[i].start < max_end) {
			interv1 = intervals[imax_end];
			interv2 = intervals[i];
			return true;
		}

		if (!same_chrom || intervals[i].end > max_end) {
			max_end = intervals[i].end;
			imax_end = i;
			started = true;
		}
	}
	return false;
}

static void warn_merged_iterator(const vector<GInterval> &intervals, const vector<GInterval> &unified,
								 GIntervalsFetcher1D *scope, const IntervUtils &iu)
{
	GInterval row;
	GInterval block;

	if (!find_merged_row(intervals, unified, scope, row, block))
		return;

	// name the pair that produced the row; the fallback cannot normally happen, since a
	// row is only uncovered when its block was built out of more than one interval
	GInterval interv1;
	GInterval interv2;

	if (!first_overlapping_pair(intervals, &block, interv1, interv2) &&
		!first_overlapping_pair(intervals, NULL, interv1, interv2))
		return;

	char msg[1000];

	snprintf(msg, sizeof(msg),
			 "Overlapping intervals were used as an iterator: %llu intervals were merged into %llu non-overlapping block%s. "
			 "For example %s %lld-%lld and %s %lld-%lld are reported as a single row %s %lld-%lld, so the results do not "
			 "correspond one-to-one to the iterator intervals. To get one value per interval, call gextract() with the "
			 "same intervals as the scope (its 'intervals' argument) and map the rows back with the intervalID column.",
			 (unsigned long long)intervals.size(), (unsigned long long)unified.size(), unified.size() == 1 ? "" : "s",
			 iu.id2chrom(interv1.chromid).c_str(), (long long)interv1.start, (long long)interv1.end,
			 iu.id2chrom(interv2.chromid).c_str(), (long long)interv2.start, (long long)interv2.end,
			 iu.id2chrom(row.chromid).c_str(), (long long)row.start, (long long)row.end);

	// Rf_warning() is unsafe here (see add_pending_diagnostic for why): the text is
	// queued for .gcall() to raise once the call has returned.
	add_pending_diagnostic(iu.get_env(), "warning", msg);
}

TrackExpressionIteratorBase *TrackExprScanner::create_expr_iterator(SEXP giterator, const TrackExpressionVars &vars, const vector<string> &track_exprs,
		GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d, GIntervals &intervals1d, GIntervals2D &intervals2d, const DiagonalBand &band, bool call_begin)
{
	monitor_memusage();

	TrackExpressionIteratorBase *expr_itr = NULL;
	SEXP class_name = Rf_getAttrib(giterator, R_ClassSymbol);

	// Resolve the iterator policy from the user-provided argument when possible
	// (fixed bin/rect, cartesian grid, big sets, explicit intervals). If not provided
	// explicitly, attempt to infer a sensible default from the tracks referenced in
	// the expressions, while verifying track format consistency.

	if ((Rf_isReal(giterator) || Rf_isInteger(giterator)) && Rf_length(giterator) == 1) {            // giterator == binsize
		verify_1d_iter(scope1d, scope2d);
		expr_itr = new TrackExpressionFixedBinIterator;
		if (call_begin) 
			((TrackExpressionFixedBinIterator *)expr_itr)->begin(Rf_isReal(giterator) ? (int64_t)REAL(giterator)[0] : INTEGER(giterator)[0], *scope1d);
	} else if ((Rf_isReal(giterator) || Rf_isInteger(giterator)) && Rf_length(giterator) == 2) {     // giterator == fixed rectangle
		verify_2d_iter(scope1d, scope2d);
		expr_itr = new TrackExpressionFixedRectIterator;
		if (call_begin) {
			if (Rf_isReal(giterator))
				((TrackExpressionFixedRectIterator *)expr_itr)->begin((int64_t)REAL(giterator)[0], (int64_t)REAL(giterator)[1], *scope2d, band);
			else
				((TrackExpressionFixedRectIterator *)expr_itr)->begin(INTEGER(giterator)[0], INTEGER(giterator)[1], *scope2d, band);
		}
	} else if (Rf_isVector(giterator) && !Rf_isNull(class_name) && !strcmp(CHAR(STRING_ELT(class_name, 0)), "cartesian.grid")) {
		GIntervals intervals1;
		GIntervals intervals2;

		m_iu.convert_rintervs(VECTOR_ELT(giterator, 0), &intervals1, NULL, false);
		if (!Rf_isNull(VECTOR_ELT(giterator, 1)))
			m_iu.convert_rintervs(VECTOR_ELT(giterator, 1), &intervals2, NULL, false);

		SEXP rexpansion1 = VECTOR_ELT(giterator, 2);
		SEXP rexpansion2 = VECTOR_ELT(giterator, 3);

		if ((!Rf_isReal(rexpansion1) && !Rf_isInteger(rexpansion1)) || (!Rf_isNull(rexpansion2) && !Rf_isReal(rexpansion2) && !Rf_isInteger(rexpansion2)))
			verror("Invalid format of cartesian grid iterator");

		vector<int64_t> expansion1(Rf_length(rexpansion1));
		for (int i = 0; i < Rf_length(rexpansion1); ++i)
			expansion1[i] = Rf_isReal(rexpansion1) ? (int64_t)REAL(rexpansion1)[i] : INTEGER(rexpansion1)[i];

		vector<int64_t> expansion2;
		if (Rf_isNull(rexpansion2))
			expansion2 = expansion1;
		else {
			expansion2.resize(Rf_length(rexpansion2));
			for (int i = 0; i < Rf_length(rexpansion2); ++i)
				expansion2[i] = Rf_isReal(rexpansion2) ? (int64_t)REAL(rexpansion2)[i] : INTEGER(rexpansion2)[i];
		}

		SEXP rband_idx = VECTOR_ELT(giterator, 4);

		if ((!Rf_isReal(rband_idx) && !Rf_isInteger(rband_idx)) || Rf_length(rband_idx) != 3)
			verror("Invalid format of cartesian grid iterator");

		int64_t min_band_idx = Rf_isReal(rband_idx) ? (int64_t)REAL(rband_idx)[0] : INTEGER(rband_idx)[0];
		int64_t max_band_idx = Rf_isReal(rband_idx) ? (int64_t)REAL(rband_idx)[1] : INTEGER(rband_idx)[1];
		bool use_band_idx_limit = Rf_isReal(rband_idx) ? (bool)REAL(rband_idx)[2] : (bool)INTEGER(rband_idx)[2];

		expr_itr = new TrackExpressionCartesianGridIterator;
		if (call_begin) 
			((TrackExpressionCartesianGridIterator *)expr_itr)->begin(m_iu.get_chromkey(), &intervals1, Rf_isNull(VECTOR_ELT(giterator, 1)) ? NULL : &intervals2,
																	  expansion1, expansion2, min_band_idx, max_band_idx, use_band_idx_limit, *scope2d, band,
																	  m_iu.get_max_data_size());
	} else if (Rf_isString(giterator) && Rf_length(giterator) == 1 && GIntervalsBigSet::isbig(CHAR(STRING_ELT(giterator, 0)), m_iu)) {
		const char *intervset = CHAR(STRING_ELT(giterator, 0));
		SEXP meta = GIntervalsMeta::load_meta(interv2path(m_iu.get_env(), intervset).c_str());

		if (GIntervalsBigSet1D::is1d(meta)) {
			expr_itr = new TrackExpressionBigSet1DIterator(m_iu);
			if (call_begin) 
				((TrackExpressionBigSet1DIterator *)expr_itr)->begin(intervset, meta, *scope1d);
		} else if (GIntervalsBigSet2D::is2d(meta)) {
			expr_itr = new TrackExpressionBigSet2DIterator(m_iu);
			if (call_begin) 
				((TrackExpressionBigSet2DIterator *)expr_itr)->begin(intervset, meta, *scope2d, band, m_iu.get_max_data_size());
		}
		runprotect(meta);
	} else if ((Rf_isVector(giterator) && !Rf_isString(giterator)) ||
			   (Rf_isString(giterator) && Rf_length(giterator) == 1 && !m_iu.track_exists(CHAR(STRING_ELT(giterator, 0)))))
	{   // giterator == intervals
		unsigned intervs_type_mask;
		if (m_iu.convert_rintervs(giterator, &intervals1d, &intervals2d, true, NULL, "", &intervs_type_mask)) {
			if (intervs_type_mask == (IntervUtils::INTERVS1D | IntervUtils::INTERVS2D))
				verror("Dual (1D and 2D) intervals cannot be used as an iterator");

			if (intervs_type_mask == IntervUtils::INTERVS1D) {
				verify_1d_iter(scope1d, scope2d);
				intervals1d.sort();

				// The merge below is silent and whether it is observable depends on the scope,
				// so the check runs once per call (see rdb::once_per_call) and only pays for
				// the copy when the intervals really do overlap. The overlap test comes
				// first, and deliberately so: once_per_call consumes the key, and a call
				// whose iterator does not overlap has nothing to report - if it consumed
				// the key anyway it would silence a nested call that does overlap (an
				// inner gextract under gintervals.mapply's FUN, say). The linear scan is
				// the cheap half; the expensive half is the copy on the next line.
				vector<GInterval> orig_intervals;
				if (scope1d && intervals_overlap(intervals1d) && once_per_call("iterator1d.overlaps"))
					orig_intervals.assign(intervals1d.begin(), intervals1d.end());

				intervals1d.unify_overlaps(false);

				if (!orig_intervals.empty())
					warn_merged_iterator(orig_intervals, intervals1d, scope1d, m_iu);
				expr_itr = new TrackExpressionIntervals1DIterator();
				if (call_begin) 
					((TrackExpressionIntervals1DIterator *)expr_itr)->begin(intervals1d, *scope1d);
			} else if (intervs_type_mask == IntervUtils::INTERVS2D) {
				verify_2d_iter(scope1d, scope2d);
				intervals2d.sort();
				intervals2d.verify_no_overlaps(m_iu.get_chromkey());
				expr_itr = new TrackExpressionIntervals2DIterator();
				if (call_begin) 
					((TrackExpressionIntervals2DIterator *)expr_itr)->begin(m_iu.get_chromkey(), intervals2d, *scope2d, band, m_iu.get_max_data_size());
			}
		}
	}

	int common_track_type = -1;
	int common_binsize = -1;

	// Read-only lookup of per-chrom sizes; using the cached const ref avoids
	// a per-call O(num_contigs) copy that adds up across multitask kids on
	// large-contig genomes.
	const GIntervals &all_genome_intervs = m_iu.get_all_genome_intervs_cached1d();

	unsigned num_tracks = vars.get_num_track_vars();
	vector<string> track_names;
	vector<GenomeTrack::Type> track_types;

for (unsigned ivar = 0; ivar < vars.get_num_track_vars(); ++ivar) {
    // Skip sequence-based variables since they don't have associated tracks
    if (vars.is_seq_variable(ivar)) {
        continue;
    }
    track_names.push_back(vars.get_track_name(ivar));
    track_types.push_back(vars.get_track_type(ivar));
}

	if (Rf_isString(giterator) && !expr_itr) {
		string iter_val(CHAR(STRING_ELT(giterator, 0)));

		if (find(track_names.begin(), track_names.end(), iter_val) == track_names.end()) {
			if (Rf_isString(giterator)) {
				SEXP all_track_names;
				SEXPCleaner all_track_names_cleaner(all_track_names);

                rprotect(all_track_names = find_in_misha(m_iu.get_env(), "GTRACKS"));
				if (Rf_isString(all_track_names)) {
					int i;
					for (i = 0; i < Rf_length(all_track_names); ++i) {
						if (iter_val == CHAR(STRING_ELT(all_track_names, i)))
							break;
                }
					if (i >= Rf_length(all_track_names)) 
						verror("Invalid iterator: %s is neither a name of a track nor a name of an intervals set", iter_val.c_str());
				}
                        }
			track_names.push_back(iter_val);
			track_types.push_back(GenomeTrack::get_type(track2path(m_iu.get_env(), iter_val).c_str(), m_iu.get_chromkey()));
		}
	}

	try {
		for (vector<string>::const_iterator itrack_name = track_names.begin(); itrack_name != track_names.end(); ++itrack_name) {
			string trackpath(track2path(m_iu.get_env(), *itrack_name));
			GenomeTrack::Type track_type = track_types[itrack_name - track_names.begin()];
			vector<string> filenames;
			unsigned binsize = 0;

			// read the list of chrom files
			get_chrom_files(trackpath.c_str(), filenames);

			// Indexed 1D fast path: when track.idx is present, validate via the
			// in-memory index instead of synthesizing one filename per chromosome
			// and stat()'ing track.idx for every one. On large-contig genomes
			// (~1M+ chroms) the synthetic loop hits ~4 syscalls × N_chroms per
			// gextract setup, with stat(track.idx) dominating.
			if (filenames.empty() && GenomeTrack::is_1d(track_type)) {
				auto idx = GenomeTrack::get_track_index(trackpath);
				if (idx) {
					const auto &entries = idx->get_all_entries();
					set<int> chromids;
					unsigned bin_size = 0;

					// FIXED_BIN: bin_size is invariant across contigs; read it
					// once from the first non-empty entry, then validate each
					// non-empty contig's num_samples in-memory via entry.length.
					if (track_type == GenomeTrack::FIXED_BIN) {
						const string dat_path = trackpath + "/track.dat";
						BufferedFile bf;
						bool bf_opened = false;
						for (const auto &entry : entries) {
							chromids.insert((int)entry.chrom_id);
							if (entry.length == 0)
								continue;
							if (bin_size == 0) {
								if (bf.open(dat_path.c_str(), "rb"))
									verror("Cannot open %s: %s", dat_path.c_str(), strerror(errno));
								bf_opened = true;
								if (bf.seek((long)entry.offset, SEEK_SET))
									verror("Failed to seek in %s", dat_path.c_str());
								if (bf.read(&bin_size, sizeof(bin_size)) != sizeof(bin_size) || bin_size == 0)
									verror("Invalid fixed-bin header in %s", dat_path.c_str());
							}
							if (entry.length < sizeof(bin_size))
								verror("Invalid fixed-bin format in %s", dat_path.c_str());
							const uint64_t data_bytes = entry.length - sizeof(bin_size);
							if ((data_bytes % sizeof(float)) != 0)
								verror("Invalid fixed-bin format in %s", dat_path.c_str());
							const int64_t num_samples = (int64_t)(data_bytes / sizeof(float));
							const int chromid = (int)entry.chrom_id;
							// chrom_id comes straight from the on-disk index; validate it
							// against the current genome before indexing (a stale index
							// built under a different chrom set could point past the end).
							if (chromid < 0 || chromid >= (int)all_genome_intervs.size())
								verror("Track %s index references chromosome id %d outside the current genome (stale index?)",
								       itrack_name->c_str(), chromid);
							const int64_t expected_num_bins = (int64_t)ceil(all_genome_intervs[chromid].end / (double)bin_size);
							if (num_samples != expected_num_bins)
								verror("Number of bins in track %s, chrom %s do not match the chromosome size (expecting: %ld, reading: %ld)",
								       itrack_name->c_str(),
								       m_iu.get_chromkey().id2chrom(chromid).c_str(),
								       expected_num_bins, num_samples);
						}
						if (bf_opened)
							bf.close();
					} else {
						// SPARSE / ARRAYS: no bin_size invariant to enforce; just
						// collect chromids for the missing-chrom check.
						for (const auto &entry : entries)
							chromids.insert((int)entry.chrom_id);
					}

					// Update common_track_type / common_binsize using the same
					// semantics as the legacy per-chrom loop's first-iteration
					// branch (line ~669 below).
					if (Rf_isString(giterator)) {
						if (*itrack_name == CHAR(STRING_ELT(giterator, 0)))
							common_track_type = track_type;
					} else {
						if (itrack_name == track_names.begin())
							common_track_type = track_type;
						else if (common_track_type != track_type)
							common_track_type = -1;
					}
					if (track_type == GenomeTrack::FIXED_BIN) {
						if (Rf_isString(giterator)) {
							if (*itrack_name == CHAR(STRING_ELT(giterator, 0)))
								common_binsize = (int)bin_size;
						} else {
							if (itrack_name == track_names.begin())
								common_binsize = (int)bin_size;
							else if (common_binsize != (int)bin_size)
								common_binsize = -1;
						}
					}

					// Equivalent of the "missing-chrom" check that ran after the
					// legacy per-chrom loop. The index entry table covers the
					// full chromkey by construction, but verify in case of a
					// stale/partial index.
					for (GIntervals::const_iterator iinterv = all_genome_intervs.begin(); iinterv != all_genome_intervs.end(); ++iinterv) {
						if (chromids.find(iinterv->chromid) == chromids.end())
							verror("Chrom %s presented in the global chrom list is missing in track %s",
							       m_iu.id2chrom(iinterv->chromid).c_str(), itrack_name->c_str());
					}

					continue;  // Skip the legacy per-chrom loop below.
				}
				// No index found: fall through to legacy per-chrom synthesis.
				for (unsigned i = 0; i < m_iu.get_chromkey().get_num_chroms(); ++i) {
					filenames.push_back(m_iu.get_chromkey().id2chrom(i));
				}
			} else if (filenames.empty() && GenomeTrack::is_2d(track_type)) {
				// For indexed 2D tracks, populate filenames from the index entries
				try {
					auto idx2d = TrackIndex2D::get_track_index_2d(trackpath);
					if (idx2d && idx2d->is_loaded()) {
						for (const auto &entry : idx2d->entries()) {
							filenames.push_back(GenomeTrack::get_2d_filename(m_iu.get_chromkey(), entry.chrom1_id, entry.chrom2_id));
						}
					}
				} catch (...) {
					// Fall through - filenames stays empty
				}
			}

			sort(filenames.begin(), filenames.end());

			if (GenomeTrack::is_1d(track_type)) {
				set<int> chromids;
				// Declare outside the loop so init_read can reuse mmap for indexed
				// tracks (single track.dat) instead of re-mmapping per chromosome.
				GenomeTrackFixedBin gtrack_fbin;

				for (vector<string>::const_iterator ifilename = filenames.begin(); ifilename != filenames.end(); ++ifilename) {
					int chromid = -1;

					try {
						 chromid = GenomeTrack::get_chromid_1d(m_iu.get_chromkey(), *ifilename);
					} catch (TGLException &e) {
						verror("Track %s: %s\n", itrack_name->c_str(), e.msg());
					}

					chromids.insert(chromid);
					if (track_type == GenomeTrack::FIXED_BIN) {
                        // Metadata-only: this loop only validates bin_size and
                        // num_samples, so skip the mmap setup that init_read
                        // otherwise does. With many tracks × many chromosomes
                        // the per-call open+mmap+madvise+close was dominating
                        // gextract setup time.
                        gtrack_fbin.init_read_metadata((trackpath + "/" + *ifilename).c_str(), chromid);
                        int64_t expected_num_bins = (int64_t)ceil(all_genome_intervs[chromid].end / (double)gtrack_fbin.get_bin_size());

                        if (gtrack_fbin.get_num_samples() != expected_num_bins){
                            verror("Number of bins in track %s, chrom %s do not match the chromosome size (expecting: %ld, reading: %ld)",
                                    itrack_name->c_str(), m_iu.get_chromkey().id2chrom(chromid).c_str(), expected_num_bins, gtrack_fbin.get_num_samples());
						}
					}

					if (ifilename == filenames.begin()) {
						if (Rf_isString(giterator)) {
							if (*itrack_name == CHAR(STRING_ELT(giterator, 0)))
								common_track_type = track_type;
						} else {
							if (itrack_name == track_names.begin())
								common_track_type = track_type;
							else if (common_track_type != track_type)
								common_track_type = -1;
						}

						if (track_type == GenomeTrack::FIXED_BIN) {
							binsize = gtrack_fbin.get_bin_size();

							if (Rf_isString(giterator)) {
								if (*itrack_name == CHAR(STRING_ELT(giterator, 0)))
									common_binsize = binsize;
							} else {
								if (itrack_name == track_names.begin())
									common_binsize = binsize;
								else if (common_binsize != (int)binsize)
									common_binsize = -1;
							}
						}
					} else if (track_type == GenomeTrack::FIXED_BIN && binsize != gtrack_fbin.get_bin_size())
						verror("Track %s: bin size of chroms %s and %s differ (%d and %d respectively)",
								itrack_name->c_str(), m_iu.id2chrom(GenomeTrack::get_chromid_1d(m_iu.get_chromkey(), *(ifilename - 1))).c_str(),
								m_iu.id2chrom(chromid).c_str(), binsize, gtrack_fbin.get_bin_size());
				}

				for (GIntervals::const_iterator iinterv = all_genome_intervs.begin(); iinterv != all_genome_intervs.end(); ++iinterv) {
					if (chromids.find(iinterv->chromid) == chromids.end())
						verror("Chrom %s presented in the global chrom list is missing in track %s", m_iu.id2chrom(iinterv->chromid).c_str(), itrack_name->c_str());
				}
			} else if (GenomeTrack::is_2d(track_type)) {
				for (vector<string>::const_iterator ifilename = filenames.begin(); ifilename != filenames.end(); ++ifilename) {
					try {
						 GenomeTrack::get_chromid_2d(m_iu.get_chromkey(), *ifilename);
					} catch (TGLException &e) {
						verror("Track %s: %s\n", itrack_name->c_str(), e.msg());
					}

					if (ifilename == filenames.begin()) {
						if (Rf_isString(giterator)) {
							if (*itrack_name == CHAR(STRING_ELT(giterator, 0)))
								common_track_type = track_type;
						} else {
							if (itrack_name == track_names.begin())
								common_track_type = track_type;
							else if (common_track_type != track_type)
								common_track_type = -1;
						}
					}
				}
			} else
				verror("Tracks of type %s are not supported in track expression", GenomeTrack::TYPE_NAMES[track_type]);
		}
	} catch (TGLException &e) {
		verror("%s\n", e.msg());
	}

	// if the iterator type is track or NULL (set implicitly), the iterator should be initialized now
	if (!expr_itr) {
		if (common_track_type < 0) {
			if (num_tracks) {
				if (track_exprs.size() == 1)
					verror("Cannot implicitly determine iterator policy:\ntrack expression \"%s\" contains a pwm virtual track or tracks in different formats.\n", track_exprs.front().c_str());
				verror("Cannot implicitly determine iterator policy: track expressions contain tracks in different formats.\n");
			}
			if (track_exprs.size() == 1)
				verror("Cannot implicitly determine iterator policy:\ntrack expression \"%s\" does not contain any tracks.\n", track_exprs.front().c_str());
			verror("Cannot implicitly determine iterator policy: track expressions do not contain any tracks.\n");
		}

		if (common_track_type == GenomeTrack::FIXED_BIN && common_binsize < 0) {
			if (track_exprs.size() == 1)
				verror("Cannot implicitly determine iterator policy:\ntrack expression \"%s\" contains tracks with different bin sizes.\n", track_exprs.front().c_str());
			verror("Cannot implicitly determine iterator policy: track expressions contain tracks with different bin sizes.\n");
		}

		if (!Rf_isString(giterator) &&
			(common_track_type == GenomeTrack::SPARSE || common_track_type == GenomeTrack::ARRAYS ||
			 common_track_type == GenomeTrack::POINTS || common_track_type == GenomeTrack::RECTS ||
			 common_track_type == GenomeTrack::COMPUTED) &&
			num_tracks > 1)
		{
			for (unsigned ivar = 1; ivar < vars.get_num_track_vars(); ++ivar) {
				if (vars.get_track_name(ivar) != vars.get_track_name(ivar - 1)) {
					if (track_exprs.size() == 1)
						verror("Cannot implicitly determine iterator policy: track expression \"%s\" contains more than one %s track.\n",
								track_exprs.front().c_str(), GenomeTrack::TYPE_NAMES[common_track_type]);
					verror("Cannot implicitly determine iterator policy: track expressions contain more than one %s track.\n", GenomeTrack::TYPE_NAMES[common_track_type]);
				}
			}
		}

		if (common_track_type == GenomeTrack::FIXED_BIN) {
			verify_1d_iter(scope1d, scope2d);
			expr_itr = new TrackExpressionFixedBinIterator();
			if (call_begin) 
				((TrackExpressionFixedBinIterator *)expr_itr)->begin(common_binsize, *scope1d);
		} else if (common_track_type == GenomeTrack::SPARSE || common_track_type == GenomeTrack::ARRAYS) {
			verify_1d_iter(scope1d, scope2d);
			expr_itr = new TrackExpressionSparseIterator(m_iu, (GenomeTrack::Type)common_track_type);
			if (call_begin) {
				if (Rf_isString(giterator))
					((TrackExpressionSparseIterator *)expr_itr)->begin(track2path(m_iu.get_env(), CHAR(STRING_ELT(giterator, 0))), *scope1d);
				else
					((TrackExpressionSparseIterator *)expr_itr)->begin(track2path(m_iu.get_env(), vars.get_track_name(0)), *scope1d);
			}
		} else if (common_track_type == GenomeTrack::RECTS || common_track_type == GenomeTrack::POINTS || common_track_type == GenomeTrack::COMPUTED) {
			verify_2d_iter(scope1d, scope2d);
			expr_itr = new TrackExpressionTrackRectsIterator(m_iu);
			if (call_begin) {
				if (Rf_isString(giterator))
					((TrackExpressionTrackRectsIterator *)expr_itr)->begin(track2path(m_iu.get_env(), CHAR(STRING_ELT(giterator, 0))), (GenomeTrack::Type)common_track_type, *scope2d, band, m_iu.get_max_data_size());
				else
					((TrackExpressionTrackRectsIterator *)expr_itr)->begin(track2path(m_iu.get_env(), vars.get_track_name(0)), (GenomeTrack::Type)common_track_type, *scope2d, band, m_iu.get_max_data_size());
			}
		} else
		verror("Unrecognized type of iterator");
	}

	return expr_itr;
}
