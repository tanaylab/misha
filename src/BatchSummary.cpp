// BatchSummary.cpp — batched multi-track genome-wide summary statistics.
//
// Uses BatchTrackScan<Summary>. Each (track, chrom) worker accumulates
// per-position counts/sums/min/max; per-track merging combines them;
// the main thread emits one row per track with columns
// {n, n_nan, min, max, sum, mean, sd}.

#include "BatchTrackScan.h"
#include "rdbinterval.h"
#include "rdbutils.h"
#include "GenomeTrack.h"

#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <vector>

using namespace rdb;
using namespace batchscan;

// ---------------------------------------------------------------------------
// Summary reducer. needs_pruning = false → driver skips sliding-deque
// maintenance (Summary needs every value).
// ---------------------------------------------------------------------------
struct Summary {
    struct Config {};

    struct State {
        uint64_t n = 0, n_nan = 0;
        double sum = 0.0, sum_sq = 0.0;
        double min_val = std::numeric_limits<double>::infinity();
        double max_val = -std::numeric_limits<double>::infinity();

        void init(const Config &, int, int) {}

        // Called only for non-NaN aggregate values (driver filters NaN
        // before calling accept). Bumps total count and updates aggregates.
        void accept(float v, int64_t) {
            ++n;
            double d = (double)v;
            sum += d;
            sum_sq += d * d;
            if (d < min_val) min_val = d;
            if (d > max_val) max_val = d;
        }

        // Called by the driver when the aggregate is NaN (all-NaN window).
        // Matches legacy gsummary semantics: every evaluated scan position
        // contributes to num_bins, and NaN-aggregate positions also
        // contribute to num_nan.
        void nan_seen() { ++n; ++n_nan; }

        // Count a pruned valid position — required interface (no-op here
        // because needs_pruning is false, this never gets called).
        void count_pruned() {}

        void boundary() {}

        bool prune(float, float) const { return false; }

        void merge(const State &o) {
            n += o.n;
            n_nan += o.n_nan;
            sum += o.sum;
            sum_sq += o.sum_sq;
            if (o.min_val < min_val) min_val = o.min_val;
            if (o.max_val > max_val) max_val = o.max_val;
        }
    };

    struct Result {
        double n, n_nan, min, max, sum, mean, sd;
    };

    static constexpr bool needs_pruning = false;
    static constexpr bool needs_lower_bound = false;
};

static Summary::Result summary_finalize(const Summary::State &s)
{
    Summary::Result r;
    r.n = (double)s.n;
    r.n_nan = (double)s.n_nan;
    uint64_t n_valid = s.n - s.n_nan;
    if (n_valid == 0) {
        r.min = r.max = r.sum = r.mean = r.sd =
            std::numeric_limits<double>::quiet_NaN();
        return r;
    }
    r.min = s.min_val;
    r.max = s.max_val;
    r.sum = s.sum;
    double mean = s.sum / (double)n_valid;
    r.mean = mean;
    if (n_valid > 1) {
        // Bessel-corrected stdev, matching IntervalSummary::get_stdev().
        double var =
            (s.sum_sq - (double)n_valid * mean * mean) / (double)(n_valid - 1);
        if (var < 0) var = 0;  // float-precision guard
        r.sd = std::sqrt(var);
    } else {
        r.sd = std::numeric_limits<double>::quiet_NaN();
    }
    return r;
}

static WindowAggFunc parse_func(SEXP _func)
{
    if (!Rf_isString(_func) || Rf_length(_func) != 1)
        verror("func must be a character scalar");
    const char *s = CHAR(STRING_ELT(_func, 0));
    if (!std::strcmp(s, "lse")) return WindowAggFunc::LSE;
    if (!std::strcmp(s, "avg")) return WindowAggFunc::AVG;
    if (!std::strcmp(s, "sum")) return WindowAggFunc::SUM;
    if (!std::strcmp(s, "max")) return WindowAggFunc::MAX;
    if (!std::strcmp(s, "min")) return WindowAggFunc::MIN;
    verror("func must be one of: lse, avg, sum, max, min");
    return WindowAggFunc::LSE;  // unreachable
}

// ---------------------------------------------------------------------------
// .Call entry: C_gsummary_multi
//
// Args (8):
//   1. _track_names   character vector
//   2. _iterator      integer scalar (bp step)
//   3. _sshift        integer scalar
//   4. _eshift        integer scalar
//   5. _n_threads     integer scalar (0 = auto)
//   6. _func          character scalar: "lse"|"avg"|"sum"|"max"|"min"
//   7. _intervals     data.frame or NULL (whole-genome)
//   8. _envir         R environment
//
// Returns: REALSXP matrix with dimensions n_tracks x 7, rownames = track
// names, colnames = {"n","n_nan","min","max","sum","mean","sd"}.
// ---------------------------------------------------------------------------
extern "C" SEXP C_gsummary_multi(
    SEXP _track_names,
    SEXP _iterator,
    SEXP _sshift,
    SEXP _eshift,
    SEXP _n_threads,
    SEXP _func,
    SEXP _intervals,
    SEXP _envir)
{
    try {
        RdbInitializer rdb_init;
        IntervUtils iu(_envir);

        int n_tracks = Rf_length(_track_names);
        std::vector<std::string> track_names(n_tracks);
        for (int i = 0; i < n_tracks; ++i)
            track_names[i] = CHAR(STRING_ELT(_track_names, i));

        int iterator_step = INTEGER(_iterator)[0];
        int sshift = INTEGER(_sshift)[0];
        int eshift = INTEGER(_eshift)[0];
        int n_threads_req = INTEGER(_n_threads)[0];
        if (iterator_step <= 0)
            verror("iterator must be a positive integer");

        WindowAggFunc func = parse_func(_func);

        const GenomeChromKey &chromkey = iu.get_chromkey();
        SEXP envir = iu.get_env();

        std::vector<std::string> track_dirs(n_tracks);
        std::vector<GenomeTrack::Type> track_types(n_tracks);
        for (int m = 0; m < n_tracks; ++m) {
            track_dirs[m] = track2path(envir, track_names[m]);
            track_types[m] = GenomeTrack::get_type(
                track_dirs[m].c_str(), chromkey, false);
        }

        const int n_chroms = (int)chromkey.get_num_chroms();
        std::vector<std::vector<GInterval>> per_chrom(n_chroms);
        bool use_intervals = !Rf_isNull(_intervals);
        if (use_intervals) {
            GIntervalsFetcher1D *i1d = nullptr;
            GIntervalsFetcher2D *i2d = nullptr;
            iu.convert_rintervs(_intervals, &i1d, &i2d);
            std::unique_ptr<GIntervalsFetcher1D> g1(i1d);
            std::unique_ptr<GIntervalsFetcher2D> g2(i2d);
            if (i2d && i2d->size() > 0)
                verror("intervals must be 1D");
            i1d->sort();
            i1d->unify_overlaps();
            for (int c = 0; c < n_chroms; ++c) {
                if (i1d->size(c) == 0) continue;
                i1d->begin_chrom_iter(c);
                for (auto it = i1d->get_chrom_begin();
                     it != i1d->get_chrom_end(); ++it) {
                    per_chrom[c].push_back(*it);
                }
            }
        }

        std::vector<Summary::Config> configs(n_tracks);

        unsigned hw = std::thread::hardware_concurrency();
        if (hw == 0) hw = 4;
        int n_threads;
        if (n_threads_req > 0) {
            n_threads = std::min(n_threads_req, n_tracks);
        } else {
            n_threads = (int)std::min((unsigned)n_tracks, std::min(hw, 40u));
        }
        if (n_threads < 1) n_threads = 1;

        ScanConfig scan;
        scan.func = func;
        scan.iterator_step = iterator_step;
        scan.sshift = sshift;
        scan.eshift = eshift;
        scan.per_chrom_intervals = use_intervals ? &per_chrom : nullptr;
        scan.n_threads = n_threads;

        BatchTrackScanResult<Summary> scan_result;
        run_batch_scan<Summary>(track_names, track_dirs, track_types,
                                configs, scan, chromkey, scan_result);

        for (int m = 0; m < n_tracks; ++m) {
            if (!scan_result.error_messages[m].empty())
                verror("Error processing track %s: %s",
                       track_names[m].c_str(),
                       scan_result.error_messages[m].c_str());
        }

        // n_tracks x 7 matrix with columns n, n_nan, min, max, sum, mean, sd.
        SEXP result;
        rprotect(result = Rf_allocMatrix(REALSXP, n_tracks, 7));
        double *out = REAL(result);
        for (int m = 0; m < n_tracks; ++m) {
            auto r = summary_finalize(scan_result.per_track_states[m]);
            out[m + 0 * n_tracks] = r.n;
            out[m + 1 * n_tracks] = r.n_nan;
            out[m + 2 * n_tracks] = r.min;
            out[m + 3 * n_tracks] = r.max;
            out[m + 4 * n_tracks] = r.sum;
            out[m + 5 * n_tracks] = r.mean;
            out[m + 6 * n_tracks] = r.sd;
        }

        SEXP rownames;
        rprotect(rownames = Rf_allocVector(STRSXP, n_tracks));
        for (int i = 0; i < n_tracks; ++i)
            SET_STRING_ELT(rownames, i, STRING_ELT(_track_names, i));

        SEXP colnames;
        rprotect(colnames = Rf_allocVector(STRSXP, 7));
        const char *cn[] = {"n", "n_nan", "min", "max", "sum", "mean", "sd"};
        for (int i = 0; i < 7; ++i) SET_STRING_ELT(colnames, i, Rf_mkChar(cn[i]));

        SEXP dimnames;
        rprotect(dimnames = Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, rownames);
        SET_VECTOR_ELT(dimnames, 1, colnames);
        Rf_setAttrib(result, R_DimNamesSymbol, dimnames);

        return result;

    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const std::bad_alloc &) {
        rerror("Out of memory");
    } catch (const std::exception &e) {
        rerror("%s", e.what());
    }
    return R_NilValue;
}
