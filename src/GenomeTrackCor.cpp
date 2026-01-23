/*
 * GenomeTrackCor.cpp
 *
 *  Created on: Sep 20, 2025
 *      Author: codex
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "rdbinterval.h"
#include "rdbutils.h"
#include "TrackExpressionScanner.h"
#include "StreamSampler.h"

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

// ============================================================================
// Spearman correlation support
// ============================================================================

// Spearman result has fewer columns (no mean/sd/cov for ranks)
enum SpearmanCols { SPEARMAN_TOTAL_BINS, SPEARMAN_TOTAL_NAN_BINS, SPEARMAN_CORR, SPEARMAN_NUM_COLS };
static const char *SpearmanColNames[SPEARMAN_NUM_COLS] = { "n", "n.na", "cor" };

static void set_spearman_colnames(SEXP colnames)
{
    for (int i = 0; i < SPEARMAN_NUM_COLS; i++)
        SET_STRING_ELT(colnames, i, Rf_mkChar(SpearmanColNames[i]));
}

// Structure to hold pair of values for exact Spearman
struct ValuePair {
    double x;
    double y;
};

// Compute ranks with average rank for ties (same as R's rank with ties.method="average")
static void compute_ranks(vector<double> &values, vector<double> &ranks)
{
    uint64_t n = values.size();
    ranks.resize(n);

    // Create index array
    vector<uint64_t> indices(n);
    for (uint64_t i = 0; i < n; ++i)
        indices[i] = i;

    // Sort indices by values
    sort(indices.begin(), indices.end(), [&values](uint64_t a, uint64_t b) {
        return values[a] < values[b];
    });

    // Assign ranks, handling ties with average
    uint64_t i = 0;
    while (i < n) {
        uint64_t j = i;
        // Find all elements with the same value (ties)
        while (j < n && values[indices[j]] == values[indices[i]])
            ++j;

        // Average rank for this group
        double avg_rank = (i + 1 + j) / 2.0;  // ranks are 1-based

        // Assign average rank to all ties
        for (uint64_t k = i; k < j; ++k)
            ranks[indices[k]] = avg_rank;

        i = j;
    }
}

// Compute Pearson correlation from ranks
static double pearson_from_values(const vector<double> &x, const vector<double> &y)
{
    uint64_t n = x.size();
    if (n < 2)
        return numeric_limits<double>::quiet_NaN();

    double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0, sum_xy = 0;
    for (uint64_t i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_x2 += x[i] * x[i];
        sum_y2 += y[i] * y[i];
        sum_xy += x[i] * y[i];
    }

    double var_x = (sum_x2 - (sum_x * sum_x) / n) / (n - 1);
    double var_y = (sum_y2 - (sum_y * sum_y) / n) / (n - 1);

    if (var_x <= 0 || var_y <= 0)
        return numeric_limits<double>::quiet_NaN();

    double cov = (sum_xy - (sum_x * sum_y) / n) / (n - 1);
    return cov / (sqrt(var_x) * sqrt(var_y));
}

// Build result vector for Spearman (single pair)
static SEXP build_spearman_result_vector(double num_bins, double num_nan_bins, double cor)
{
    SEXP answer = rprotect_ptr(RSaneAllocVector(REALSXP, SPEARMAN_NUM_COLS));
    SEXP colnames = rprotect_ptr(RSaneAllocVector(STRSXP, SPEARMAN_NUM_COLS));

    REAL(answer)[SPEARMAN_TOTAL_BINS] = num_bins;
    REAL(answer)[SPEARMAN_TOTAL_NAN_BINS] = num_nan_bins;
    REAL(answer)[SPEARMAN_CORR] = cor;

    set_spearman_colnames(colnames);
    Rf_setAttrib(answer, R_NamesSymbol, colnames);

    runprotect(2);
    return answer;
}

// Build result matrix for Spearman (multiple pairs)
static SEXP build_spearman_result_matrix(const vector<double> &num_bins,
                                          const vector<double> &num_nan_bins,
                                          const vector<double> &cors)
{
    uint64_t num_pairs = cors.size();
    SEXP answer = rprotect_ptr(RSaneAllocVector(REALSXP, num_pairs * SPEARMAN_NUM_COLS));
    double *panswer = REAL(answer);

    for (uint64_t i = 0; i < num_pairs; ++i) {
        panswer[num_pairs * SPEARMAN_TOTAL_BINS + i] = num_bins[i];
        panswer[num_pairs * SPEARMAN_TOTAL_NAN_BINS + i] = num_nan_bins[i];
        panswer[num_pairs * SPEARMAN_CORR + i] = cors[i];
    }

    SEXP dim = rprotect_ptr(RSaneAllocVector(INTSXP, 2));
    SEXP dimnames = rprotect_ptr(RSaneAllocVector(VECSXP, 2));
    SEXP colnames = rprotect_ptr(RSaneAllocVector(STRSXP, SPEARMAN_NUM_COLS));

    INTEGER(dim)[0] = num_pairs;
    INTEGER(dim)[1] = SPEARMAN_NUM_COLS;
    SET_VECTOR_ELT(dimnames, 0, R_NilValue);
    set_spearman_colnames(colnames);
    SET_VECTOR_ELT(dimnames, 1, colnames);

    Rf_setAttrib(answer, R_DimSymbol, dim);
    Rf_setAttrib(answer, R_DimNamesSymbol, dimnames);

    runprotect(4);
    return answer;
}

// Estimate rank of a value using binary search in sorted sample
static double estimate_rank(double value, const vector<double> &sorted_sample, uint64_t total_count)
{
    if (sorted_sample.empty())
        return numeric_limits<double>::quiet_NaN();

    // Binary search to find position
    auto it = lower_bound(sorted_sample.begin(), sorted_sample.end(), value);
    uint64_t pos = it - sorted_sample.begin();

    // Count equal values for tie handling
    auto upper = upper_bound(sorted_sample.begin(), sorted_sample.end(), value);
    uint64_t count_equal = upper - it;

    // Estimate rank: scale position to total count
    double scale = (double)total_count / sorted_sample.size();

    if (count_equal > 0) {
        // Average rank for ties
        double low_rank = pos * scale + 1;
        double high_rank = (pos + count_equal) * scale;
        return (low_rank + high_rank) / 2.0;
    } else {
        return pos * scale + 0.5;  // Interpolate
    }
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

// ============================================================================
// Exact Spearman correlation (two-pass, collects all values)
// ============================================================================

SEXP gtrackcor_spearman_exact(SEXP _track_exprs, SEXP _intervals, SEXP _iterator_policy, SEXP _band, SEXP _envir)
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

        // Collect all (x, y) pairs for each track pair
        vector<vector<ValuePair>> all_pairs(num_pairs);
        vector<double> total_bins(num_pairs, 0);

        for (scanner.begin(_track_exprs, intervals1d, intervals2d, _iterator_policy, _band); !scanner.isend(); scanner.next()) {
            for (unsigned i = 0; i < num_pairs; ++i) {
                double x = scanner.last_real(2 * i);
                double y = scanner.last_real(2 * i + 1);
                ++total_bins[i];
                if (!std::isnan(x) && !std::isnan(y)) {
                    all_pairs[i].push_back({x, y});
                    iu.verify_max_data_size(all_pairs[i].size(), "Result");
                }
            }
        }

        // Compute Spearman correlation for each pair
        vector<double> num_nan_bins(num_pairs);
        vector<double> cors(num_pairs);

        for (unsigned i = 0; i < num_pairs; ++i) {
            num_nan_bins[i] = total_bins[i] - all_pairs[i].size();

            if (all_pairs[i].size() < 2) {
                cors[i] = numeric_limits<double>::quiet_NaN();
                continue;
            }

            // Extract x and y values
            vector<double> x_vals(all_pairs[i].size());
            vector<double> y_vals(all_pairs[i].size());
            for (uint64_t j = 0; j < all_pairs[i].size(); ++j) {
                x_vals[j] = all_pairs[i][j].x;
                y_vals[j] = all_pairs[i][j].y;
            }

            // Compute ranks
            vector<double> x_ranks, y_ranks;
            compute_ranks(x_vals, x_ranks);
            compute_ranks(y_vals, y_ranks);

            // Compute Pearson correlation on ranks
            cors[i] = pearson_from_values(x_ranks, y_ranks);
        }

        if (num_pairs == 1)
            return build_spearman_result_vector(total_bins[0], num_nan_bins[0], cors[0]);

        return build_spearman_result_matrix(total_bins, num_nan_bins, cors);
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    return R_NilValue;
}

SEXP gtrackcor_spearman_exact_multitask(SEXP _track_exprs, SEXP _intervals, SEXP _iterator_policy, SEXP _band, SEXP _envir)
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

        // For exact Spearman with multitasking, children collect pairs and send to parent
        // Parent merges all pairs and computes final correlation

        // Fixed header size per pair: total_bins (uint64_t) + num_values (uint64_t)
        uint64_t header_size = sizeof(uint64_t) * 2 * num_pairs;
        // Each record is a ValuePair (two doubles)
        uint64_t record_size = sizeof(ValuePair);

        if (iu.distribute_task(header_size, record_size)) {
            // Child process: collect all pairs
            TrackExprScanner scanner(iu);
            vector<vector<ValuePair>> all_pairs(num_pairs);
            vector<double> total_bins(num_pairs, 0);

            for (scanner.begin(_track_exprs, iu.get_kid_intervals1d(), iu.get_kid_intervals2d(), _iterator_policy, _band); !scanner.isend(); scanner.next()) {
                for (unsigned i = 0; i < num_pairs; ++i) {
                    double x = scanner.last_real(2 * i);
                double y = scanner.last_real(2 * i + 1);
                ++total_bins[i];
                if (!std::isnan(x) && !std::isnan(y)) {
                    all_pairs[i].push_back({x, y});
                    iu.verify_max_data_size(all_pairs[i].size(), "Result");
                }
            }
        }

            // Pack results: for each pair, send total_bins, num_values, then values
            uint64_t total_values = 0;
            for (unsigned i = 0; i < num_pairs; ++i)
                total_values += all_pairs[i].size();

            void *ptr = allocate_res(total_values);  // Allocates header_size + total_values * record_size

            for (unsigned i = 0; i < num_pairs; ++i) {
                uint64_t tb = (uint64_t)total_bins[i];
                uint64_t nv = all_pairs[i].size();
                pack_data(ptr, tb, 1);
                pack_data(ptr, nv, 1);
                if (nv > 0)
                    pack_data(ptr, all_pairs[i][0], nv);
            }
        } else {
            // Parent process: collect from all children, merge, and compute
            vector<vector<ValuePair>> all_pairs(num_pairs);
            vector<double> total_bins(num_pairs, 0);

            for (int k = 0; k < get_num_kids(); ++k) {
                void *ptr = get_kid_res(k);

                for (unsigned i = 0; i < num_pairs; ++i) {
                    uint64_t tb, nv;
                    unpack_data(ptr, tb, 1);
                    unpack_data(ptr, nv, 1);
                    total_bins[i] += tb;

                    if (nv > 0) {
                        uint64_t old_size = all_pairs[i].size();
                        uint64_t new_size = old_size + nv;
                        iu.verify_max_data_size(new_size, "Result");
                        all_pairs[i].resize(new_size);
                        unpack_data(ptr, all_pairs[i][old_size], nv);
                    }
                }
            }

            // Compute Spearman correlation for each pair
            vector<double> num_nan_bins(num_pairs);
            vector<double> cors(num_pairs);

            for (unsigned i = 0; i < num_pairs; ++i) {
                num_nan_bins[i] = total_bins[i] - all_pairs[i].size();

                if (all_pairs[i].size() < 2) {
                    cors[i] = numeric_limits<double>::quiet_NaN();
                    continue;
                }

                // Extract x and y values
                vector<double> x_vals(all_pairs[i].size());
                vector<double> y_vals(all_pairs[i].size());
                for (uint64_t j = 0; j < all_pairs[i].size(); ++j) {
                    x_vals[j] = all_pairs[i][j].x;
                    y_vals[j] = all_pairs[i][j].y;
                }

                // Compute ranks
                vector<double> x_ranks, y_ranks;
                compute_ranks(x_vals, x_ranks);
                compute_ranks(y_vals, y_ranks);

                // Compute Pearson correlation on ranks
                cors[i] = pearson_from_values(x_ranks, y_ranks);
            }

            if (num_pairs == 1)
                rreturn(build_spearman_result_vector(total_bins[0], num_nan_bins[0], cors[0]));

            rreturn(build_spearman_result_matrix(total_bins, num_nan_bins, cors));
        }
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(R_NilValue);
}

// ============================================================================
// Approximate Spearman correlation (streaming with sampling)
// ============================================================================

SEXP gtrackcor_spearman(SEXP _track_exprs, SEXP _intervals, SEXP _iterator_policy, SEXP _band, SEXP _envir)
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

        uint64_t sample_size = iu.get_max_data_size();

        // For each pair: sample (x,y) pairs, and separately sample x and y for rank estimation
        vector<StreamSampler<ValuePair>> pair_samplers(num_pairs);
        vector<StreamSampler<double>> x_samplers(num_pairs);
        vector<StreamSampler<double>> y_samplers(num_pairs);
        vector<double> total_bins(num_pairs, 0);
        vector<uint64_t> non_nan_counts(num_pairs, 0);

        for (unsigned i = 0; i < num_pairs; ++i) {
            pair_samplers[i].init(sample_size);
            x_samplers[i].init(sample_size);
            y_samplers[i].init(sample_size);
        }

        for (scanner.begin(_track_exprs, intervals1d, intervals2d, _iterator_policy, _band); !scanner.isend(); scanner.next()) {
            for (unsigned i = 0; i < num_pairs; ++i) {
                double x = scanner.last_real(2 * i);
                double y = scanner.last_real(2 * i + 1);
                ++total_bins[i];
                if (!std::isnan(x) && !std::isnan(y)) {
                    ++non_nan_counts[i];
                    pair_samplers[i].add({x, y}, unif_rand);
                    x_samplers[i].add(x, unif_rand);
                    y_samplers[i].add(y, unif_rand);
                }
            }
        }

        // Compute approximate Spearman correlation for each pair
        vector<double> num_nan_bins(num_pairs);
        vector<double> cors(num_pairs);

        for (unsigned i = 0; i < num_pairs; ++i) {
            num_nan_bins[i] = total_bins[i] - non_nan_counts[i];

            if (non_nan_counts[i] < 2) {
                cors[i] = numeric_limits<double>::quiet_NaN();
                continue;
            }

            // Sort samples for rank estimation
            vector<double> &x_samples = x_samplers[i].samples();
            vector<double> &y_samples = y_samplers[i].samples();
            sort(x_samples.begin(), x_samples.end());
            sort(y_samples.begin(), y_samples.end());

            // Estimate ranks for sampled pairs
            const vector<ValuePair> &pairs = pair_samplers[i].samples();
            vector<double> x_ranks(pairs.size());
            vector<double> y_ranks(pairs.size());

            for (uint64_t j = 0; j < pairs.size(); ++j) {
                x_ranks[j] = estimate_rank(pairs[j].x, x_samples, non_nan_counts[i]);
                y_ranks[j] = estimate_rank(pairs[j].y, y_samples, non_nan_counts[i]);
            }

            // Compute Pearson correlation on estimated ranks
            cors[i] = pearson_from_values(x_ranks, y_ranks);
        }

        if (num_pairs == 1)
            return build_spearman_result_vector(total_bins[0], num_nan_bins[0], cors[0]);

        return build_spearman_result_matrix(total_bins, num_nan_bins, cors);
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    return R_NilValue;
}

SEXP gtrackcor_spearman_multitask(SEXP _track_exprs, SEXP _intervals, SEXP _iterator_policy, SEXP _band, SEXP _envir)
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

        int num_kids = iu.prepare4multitasking(_track_exprs, intervals1d, intervals2d, _iterator_policy, _band);
        if (!num_kids)
            rreturn(R_NilValue);

        uint64_t sample_size = (uint64_t)ceil(iu.get_max_data_size() / (double)num_kids);

        // Estimate result size: for each pair, we send:
        // - total_bins (uint64_t) + non_nan_counts (uint64_t)
        // - pair_samples size (uint64_t) + pair_samples data (sample_size * sizeof(ValuePair))
        // - x_samples size (uint64_t) + x_samples data (sample_size * sizeof(double))
        // - y_samples size (uint64_t) + y_samples data (sample_size * sizeof(double))
        uint64_t per_pair_size = 2 * sizeof(uint64_t) +                           // total_bins, non_nan_counts
                                 sizeof(uint64_t) + sample_size * sizeof(ValuePair) +  // pair_samples
                                 sizeof(uint64_t) + sample_size * sizeof(double) +     // x_samples
                                 sizeof(uint64_t) + sample_size * sizeof(double);      // y_samples
        uint64_t result_size = num_pairs * per_pair_size;

        if (iu.distribute_task(result_size, 0)) {
            // Child process: collect samples
            TrackExprScanner scanner(iu);

            vector<StreamSampler<ValuePair>> pair_samplers(num_pairs);
            vector<StreamSampler<double>> x_samplers(num_pairs);
            vector<StreamSampler<double>> y_samplers(num_pairs);
            vector<double> total_bins(num_pairs, 0);
            vector<uint64_t> non_nan_counts(num_pairs, 0);

            for (unsigned i = 0; i < num_pairs; ++i) {
                pair_samplers[i].init(sample_size);
                x_samplers[i].init(sample_size);
                y_samplers[i].init(sample_size);
            }

            for (scanner.begin(_track_exprs, iu.get_kid_intervals1d(), iu.get_kid_intervals2d(), _iterator_policy, _band); !scanner.isend(); scanner.next()) {
                for (unsigned i = 0; i < num_pairs; ++i) {
                    double x = scanner.last_real(2 * i);
                    double y = scanner.last_real(2 * i + 1);
                    ++total_bins[i];
                    if (!std::isnan(x) && !std::isnan(y)) {
                        ++non_nan_counts[i];
                        pair_samplers[i].add({x, y}, unif_rand);
                        x_samplers[i].add(x, unif_rand);
                        y_samplers[i].add(y, unif_rand);
                    }
                }
            }

            // Pack results
            void *ptr = allocate_res(0);

            for (unsigned i = 0; i < num_pairs; ++i) {
                uint64_t tb = (uint64_t)total_bins[i];
                uint64_t nnc = non_nan_counts[i];
                pack_data(ptr, tb, 1);
                pack_data(ptr, nnc, 1);

                // Pack pair samples
                uint64_t ps = pair_samplers[i].samples().size();
                pack_data(ptr, ps, 1);
                if (ps > 0)
                    pack_data(ptr, pair_samplers[i].samples()[0], ps);

                // Pack x samples
                uint64_t xs = x_samplers[i].samples().size();
                pack_data(ptr, xs, 1);
                if (xs > 0)
                    pack_data(ptr, x_samplers[i].samples()[0], xs);

                // Pack y samples
                uint64_t ys = y_samplers[i].samples().size();
                pack_data(ptr, ys, 1);
                if (ys > 0)
                    pack_data(ptr, y_samplers[i].samples()[0], ys);
            }
        } else {
            // Parent process: merge samples from all children
            vector<double> total_bins(num_pairs, 0);
            vector<uint64_t> non_nan_counts(num_pairs, 0);
            vector<vector<ValuePair>> merged_pairs(num_pairs);
            vector<vector<double>> merged_x(num_pairs);
            vector<vector<double>> merged_y(num_pairs);
            vector<double> min_sampling_rate_pairs(num_pairs, 1.0);
            vector<double> min_sampling_rate_x(num_pairs, 1.0);
            vector<double> min_sampling_rate_y(num_pairs, 1.0);

            // First pass: collect all data and find minimum sampling rates
            struct KidData {
                vector<uint64_t> tb, nnc;
                vector<vector<ValuePair>> pairs;
                vector<vector<double>> x_samp, y_samp;
            };
            vector<KidData> kid_data(get_num_kids());

            for (int k = 0; k < get_num_kids(); ++k) {
                void *ptr = get_kid_res(k);
                kid_data[k].tb.resize(num_pairs);
                kid_data[k].nnc.resize(num_pairs);
                kid_data[k].pairs.resize(num_pairs);
                kid_data[k].x_samp.resize(num_pairs);
                kid_data[k].y_samp.resize(num_pairs);

                for (unsigned i = 0; i < num_pairs; ++i) {
                    unpack_data(ptr, kid_data[k].tb[i], 1);
                    unpack_data(ptr, kid_data[k].nnc[i], 1);

                    total_bins[i] += kid_data[k].tb[i];
                    non_nan_counts[i] += kid_data[k].nnc[i];

                    uint64_t ps, xs, ys;
                    unpack_data(ptr, ps, 1);
                    if (ps > 0) {
                        kid_data[k].pairs[i].resize(ps);
                        unpack_data(ptr, kid_data[k].pairs[i][0], ps);
                        if (kid_data[k].nnc[i] > 0) {
                            double rate = (double)ps / kid_data[k].nnc[i];
                            min_sampling_rate_pairs[i] = min(min_sampling_rate_pairs[i], rate);
                        }
                    }

                    unpack_data(ptr, xs, 1);
                    if (xs > 0) {
                        kid_data[k].x_samp[i].resize(xs);
                        unpack_data(ptr, kid_data[k].x_samp[i][0], xs);
                        if (kid_data[k].nnc[i] > 0) {
                            double rate = (double)xs / kid_data[k].nnc[i];
                            min_sampling_rate_x[i] = min(min_sampling_rate_x[i], rate);
                        }
                    }

                    unpack_data(ptr, ys, 1);
                    if (ys > 0) {
                        kid_data[k].y_samp[i].resize(ys);
                        unpack_data(ptr, kid_data[k].y_samp[i][0], ys);
                        if (kid_data[k].nnc[i] > 0) {
                            double rate = (double)ys / kid_data[k].nnc[i];
                            min_sampling_rate_y[i] = min(min_sampling_rate_y[i], rate);
                        }
                    }
                }
            }

            // Second pass: merge samples with uniform sampling rate
            for (int k = 0; k < get_num_kids(); ++k) {
                for (unsigned i = 0; i < num_pairs; ++i) {
                    if (kid_data[k].nnc[i] == 0)
                        continue;

                    // Merge pairs
                    double kid_rate = kid_data[k].pairs[i].empty() ? 0 : (double)kid_data[k].pairs[i].size() / kid_data[k].nnc[i];
                    if (kid_rate > 0) {
                        double keep_prob = min_sampling_rate_pairs[i] / kid_rate;
                        for (const auto &p : kid_data[k].pairs[i]) {
                            if (keep_prob >= 1.0 || unif_rand() < keep_prob)
                                merged_pairs[i].push_back(p);
                        }
                    }

                    // Merge x samples
                    kid_rate = kid_data[k].x_samp[i].empty() ? 0 : (double)kid_data[k].x_samp[i].size() / kid_data[k].nnc[i];
                    if (kid_rate > 0) {
                        double keep_prob = min_sampling_rate_x[i] / kid_rate;
                        for (double v : kid_data[k].x_samp[i]) {
                            if (keep_prob >= 1.0 || unif_rand() < keep_prob)
                                merged_x[i].push_back(v);
                        }
                    }

                    // Merge y samples
                    kid_rate = kid_data[k].y_samp[i].empty() ? 0 : (double)kid_data[k].y_samp[i].size() / kid_data[k].nnc[i];
                    if (kid_rate > 0) {
                        double keep_prob = min_sampling_rate_y[i] / kid_rate;
                        for (double v : kid_data[k].y_samp[i]) {
                            if (keep_prob >= 1.0 || unif_rand() < keep_prob)
                                merged_y[i].push_back(v);
                        }
                    }
                }
            }

            // Compute approximate Spearman correlation for each pair
            vector<double> num_nan_bins(num_pairs);
            vector<double> cors(num_pairs);

            for (unsigned i = 0; i < num_pairs; ++i) {
                num_nan_bins[i] = total_bins[i] - non_nan_counts[i];

                if (non_nan_counts[i] < 2 || merged_pairs[i].size() < 2) {
                    cors[i] = numeric_limits<double>::quiet_NaN();
                    continue;
                }

                // Sort samples for rank estimation
                sort(merged_x[i].begin(), merged_x[i].end());
                sort(merged_y[i].begin(), merged_y[i].end());

                // Estimate ranks for sampled pairs
                vector<double> x_ranks(merged_pairs[i].size());
                vector<double> y_ranks(merged_pairs[i].size());

                for (uint64_t j = 0; j < merged_pairs[i].size(); ++j) {
                    x_ranks[j] = estimate_rank(merged_pairs[i][j].x, merged_x[i], non_nan_counts[i]);
                    y_ranks[j] = estimate_rank(merged_pairs[i][j].y, merged_y[i], non_nan_counts[i]);
                }

                // Compute Pearson correlation on estimated ranks
                cors[i] = pearson_from_values(x_ranks, y_ranks);
            }

            if (num_pairs == 1)
                rreturn(build_spearman_result_vector(total_bins[0], num_nan_bins[0], cors[0]));

            rreturn(build_spearman_result_matrix(total_bins, num_nan_bins, cors));
        }
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(R_NilValue);
}

}
