/*
 * GenomeTrackCor.cpp
 *
 *  Created on: Sep 20, 2025
 *      Author: codex
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "rdbinterval.h"
#include "rdbutils.h"
#include "TrackExpressionScanner.h"

using namespace std;
using namespace rdb;

struct CorSummary {
    double num_bins;
    double num_non_nan_bins;
    double sum_x;
    double sum_y;
    double sum_x2;
    double sum_y2;
    double sum_xy;

    void update(double x, double y) {
        ++num_bins;
        if (!std::isnan(x) && !std::isnan(y)) {
            ++num_non_nan_bins;
            sum_x += x;
            sum_y += y;
            sum_x2 += x * x;
            sum_y2 += y * y;
            sum_xy += x * y;
        }
    }

    void merge(const CorSummary &obj) {
        num_bins += obj.num_bins;
        num_non_nan_bins += obj.num_non_nan_bins;
        sum_x += obj.sum_x;
        sum_y += obj.sum_y;
        sum_x2 += obj.sum_x2;
        sum_y2 += obj.sum_y2;
        sum_xy += obj.sum_xy;
    }

    double mean_x() const { return sum_x / num_non_nan_bins; }
    double mean_y() const { return sum_y / num_non_nan_bins; }

    double cov() const {
        double n = num_non_nan_bins;
        return (sum_xy - (sum_x * sum_y) / n) / (n - 1);
    }

    double var_x() const {
        double n = num_non_nan_bins;
        double var = (sum_x2 - (sum_x * sum_x) / n) / (n - 1);
        return var;
    }

    double var_y() const {
        double n = num_non_nan_bins;
        double var = (sum_y2 - (sum_y * sum_y) / n) / (n - 1);
        return var;
    }

    CorSummary() :
        num_bins(0),
        num_non_nan_bins(0),
        sum_x(0),
        sum_y(0),
        sum_x2(0),
        sum_y2(0),
        sum_xy(0) {}
};

enum CorSummaryCols { TOTAL_BINS, TOTAL_NAN_BINS, MEAN1, MEAN2, SD1, SD2, COV, CORR, NUM_COLS };

static const char *CorSummaryColNames[NUM_COLS] = { "n", "n.na", "mean1", "mean2", "sd1", "sd2", "cov", "cor" };

static void fill_values(const CorSummary &summary, double *out)
{
    out[TOTAL_BINS] = summary.num_bins;
    out[TOTAL_NAN_BINS] = summary.num_bins - summary.num_non_nan_bins;

    if (summary.num_non_nan_bins == 0) {
        out[MEAN1] = numeric_limits<double>::quiet_NaN();
        out[MEAN2] = numeric_limits<double>::quiet_NaN();
        out[SD1] = numeric_limits<double>::quiet_NaN();
        out[SD2] = numeric_limits<double>::quiet_NaN();
        out[COV] = numeric_limits<double>::quiet_NaN();
        out[CORR] = numeric_limits<double>::quiet_NaN();
        return;
    }

    double mean1 = summary.mean_x();
    double mean2 = summary.mean_y();
    out[MEAN1] = mean1;
    out[MEAN2] = mean2;

    if (summary.num_non_nan_bins < 2) {
        out[SD1] = numeric_limits<double>::quiet_NaN();
        out[SD2] = numeric_limits<double>::quiet_NaN();
        out[COV] = numeric_limits<double>::quiet_NaN();
        out[CORR] = numeric_limits<double>::quiet_NaN();
        return;
    }

    double var1 = summary.var_x();
    double var2 = summary.var_y();

    if (var1 < 0 && var1 > -1e-12) {
        var1 = 0;
    }
    if (var2 < 0 && var2 > -1e-12) {
        var2 = 0;
    }

    double sd1 = var1 >= 0 ? sqrt(var1) : numeric_limits<double>::quiet_NaN();
    double sd2 = var2 >= 0 ? sqrt(var2) : numeric_limits<double>::quiet_NaN();
    out[SD1] = sd1;
    out[SD2] = sd2;

    double cov = summary.cov();
    out[COV] = cov;

    if (sd1 > 0 && sd2 > 0) {
        out[CORR] = cov / (sd1 * sd2);
    } else {
        out[CORR] = numeric_limits<double>::quiet_NaN();
    }
}

static void set_colnames(SEXP colnames)
{
    for (int i = 0; i < NUM_COLS; i++)
        SET_STRING_ELT(colnames, i, Rf_mkChar(CorSummaryColNames[i]));
}

static SEXP build_result_vector(const CorSummary &summary)
{
    SEXP answer = rprotect_ptr(RSaneAllocVector(REALSXP, NUM_COLS));
    SEXP colnames = rprotect_ptr(RSaneAllocVector(STRSXP, NUM_COLS));

    fill_values(summary, REAL(answer));
    set_colnames(colnames);
    Rf_setAttrib(answer, R_NamesSymbol, colnames);

    runprotect(2);
    return answer;
}

static SEXP build_result_matrix(const vector<CorSummary> &summaries)
{
    uint64_t num_pairs = summaries.size();
    SEXP answer = rprotect_ptr(RSaneAllocVector(REALSXP, num_pairs * NUM_COLS));
    double *panswer = REAL(answer);

    for (uint64_t i = 0; i < num_pairs; ++i) {
        double vals[NUM_COLS];
        fill_values(summaries[i], vals);
        for (int col = 0; col < NUM_COLS; ++col) {
            panswer[num_pairs * col + i] = vals[col];
        }
    }

    SEXP dim = rprotect_ptr(RSaneAllocVector(INTSXP, 2));
    SEXP dimnames = rprotect_ptr(RSaneAllocVector(VECSXP, 2));
    SEXP colnames = rprotect_ptr(RSaneAllocVector(STRSXP, NUM_COLS));

    INTEGER(dim)[0] = num_pairs;
    INTEGER(dim)[1] = NUM_COLS;
    SET_VECTOR_ELT(dimnames, 0, R_NilValue);
    set_colnames(colnames);
    SET_VECTOR_ELT(dimnames, 1, colnames);

    Rf_setAttrib(answer, R_DimSymbol, dim);
    Rf_setAttrib(answer, R_DimNamesSymbol, dimnames);

    runprotect(4);
    return answer;
}

extern "C" {

SEXP gtrackcor(SEXP _track_exprs, SEXP _intervals, SEXP _iterator_policy, SEXP _band, SEXP _envir)
{
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_track_exprs) || Rf_length(_track_exprs) < 2 || (Rf_length(_track_exprs) % 2) != 0)
            verror("Track expression argument must be a string vector of even length (pairs)");
        unsigned num_pairs = Rf_length(_track_exprs) / 2;

        IntervUtils iu(_envir);
        TrackExprScanner scanner(iu);
        GIntervalsFetcher1D *intervals1d = NULL;
        GIntervalsFetcher2D *intervals2d = NULL;
        iu.convert_rintervs(_intervals, &intervals1d, &intervals2d);
        unique_ptr<GIntervalsFetcher1D> intervals1d_guard(intervals1d);
        unique_ptr<GIntervalsFetcher2D> intervals2d_guard(intervals2d);
        intervals1d->sort();
        intervals1d->unify_overlaps();
        intervals2d->sort();
        intervals2d->verify_no_overlaps(iu.get_chromkey());

        vector<CorSummary> summaries(num_pairs);

        for (scanner.begin(_track_exprs, intervals1d, intervals2d, _iterator_policy, _band); !scanner.isend(); scanner.next()) {
            for (unsigned i = 0; i < num_pairs; ++i) {
                summaries[i].update(scanner.last_real(2 * i), scanner.last_real(2 * i + 1));
            }
        }

        if (num_pairs == 1)
            return build_result_vector(summaries[0]);

        return build_result_matrix(summaries);
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    return R_NilValue;
}

SEXP gtrackcor_multitask(SEXP _track_exprs, SEXP _intervals, SEXP _iterator_policy, SEXP _band, SEXP _envir)
{
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_track_exprs) || Rf_length(_track_exprs) < 2 || (Rf_length(_track_exprs) % 2) != 0)
            verror("Track expression argument must be a string vector of even length (pairs)");
        unsigned num_pairs = Rf_length(_track_exprs) / 2;

        IntervUtils iu(_envir);
        GIntervalsFetcher1D *intervals1d = NULL;
        GIntervalsFetcher2D *intervals2d = NULL;
        iu.convert_rintervs(_intervals, &intervals1d, &intervals2d);
        unique_ptr<GIntervalsFetcher1D> intervals1d_guard(intervals1d);
        unique_ptr<GIntervalsFetcher2D> intervals2d_guard(intervals2d);
        intervals1d->sort();
        intervals1d->unify_overlaps();
        intervals2d->sort();
        intervals2d->verify_no_overlaps(iu.get_chromkey());

        if (!iu.prepare4multitasking(_track_exprs, intervals1d, intervals2d, _iterator_policy, _band))
            rreturn(R_NilValue);

        if (iu.distribute_task(num_pairs * sizeof(CorSummary), 0)) { // child process
            CorSummary *summaries = (CorSummary *)allocate_res(0);
            TrackExprScanner scanner(iu);

            for (unsigned i = 0; i < num_pairs; ++i)
                summaries[i] = CorSummary();

            for (scanner.begin(_track_exprs, iu.get_kid_intervals1d(), iu.get_kid_intervals2d(), _iterator_policy, _band); !scanner.isend(); scanner.next()) {
                for (unsigned i = 0; i < num_pairs; ++i) {
                    summaries[i].update(scanner.last_real(2 * i), scanner.last_real(2 * i + 1));
                }
            }
        } else { // parent process
            vector<CorSummary> summaries(num_pairs);

            for (int i = 0; i < get_num_kids(); ++i) {
                CorSummary *kid_summaries = (CorSummary *)get_kid_res(i);
                for (unsigned j = 0; j < num_pairs; ++j)
                    summaries[j].merge(kid_summaries[j]);
            }

            if (num_pairs == 1)
                rreturn(build_result_vector(summaries[0]));

            rreturn(build_result_matrix(summaries));
        }
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(R_NilValue);
}

}
