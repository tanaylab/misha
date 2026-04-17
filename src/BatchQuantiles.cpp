// BatchQuantiles.cpp — batched multi-track genome-wide quantile scan.
//
// Phase 2: TopKQuantile reducer with top-K / bottom-K heap-backed vector
// plus full-vector fallback, sliding-max upper-bound pruning, and
// per-aggregator (`lse`/`avg`/`sum`/`max`/`min`) templating. Supports
// whole-genome scans and caller-provided intervals.

#include "BatchTrackScan.h"
#include "rdbinterval.h"
#include "rdbutils.h"
#include "GenomeTrack.h"

#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <vector>

using namespace rdb;
using namespace batchscan;

// ---------------------------------------------------------------------------
// TopKQuantile reducer.
//
// Mode (set at Config build time):
//   use_fallback = true  -> State::buf holds every accepted value; final
//                           nth_element runs over the full buffer.
//   use_fallback = false -> State::buf is a heap-backed vector of size K.
//                           top_side = true  -> min-heap of largest K so far
//                           top_side = false -> max-heap of smallest K so far
//
// Heap build is lazy: we append until buf.size() == K, then make_heap once,
// and thereafter swap-replace the extremum. This avoids heap overhead while
// the buffer is filling and amortizes the log(K) cost across pushes.
// ---------------------------------------------------------------------------
struct TopKQuantile {
    struct Config {
        std::vector<double> percentiles;
        uint32_t K = 0;
        bool use_fallback = true;
        bool top_side = true;
    };

    struct State {
        const Config *cfg = nullptr;
        uint64_t n_total = 0;
        std::vector<float> buf;
        bool heap_built = false;

        void init(const Config &c, int /*chromid*/, int /*iterator_step*/) {
            cfg = &c;
            buf.reserve(c.use_fallback ? 1024
                                       : std::min<size_t>(c.K, 1u << 20));
        }

        // Called by the scan driver for every non-NaN position, whether it
        // gets pushed into the heap or not. This keeps n_total consistent
        // between fallback and heap modes — without it, heap-mode would
        // undercount pruned-but-valid positions, shifting the rank
        // calculation in topk_finalize.
        void count_pruned() { ++n_total; }

        // NaN-aggregate positions have no rank. Ignored for quantile math;
        // they neither contribute to n_total nor the heap.
        void nan_seen() {}

        void accept(float v, int64_t /*pos*/) {
            ++n_total;
            if (!cfg || cfg->use_fallback) {
                buf.push_back(v);
                return;
            }
            uint32_t K = cfg->K;
            if (!heap_built) {
                buf.push_back(v);
                // Use >= K so that a merged accumulator (which may already
                // hold K elements before the first accept) heapifies on its
                // very first accept rather than silently growing past K.
                if (buf.size() >= (size_t)K) {
                    if (cfg->top_side)
                        std::make_heap(buf.begin(), buf.end(),
                                       std::greater<float>{});
                    else
                        std::make_heap(buf.begin(), buf.end(),
                                       std::less<float>{});
                    heap_built = true;
                }
                return;
            }
            // Heap full; swap-replace if v improves on the extremum.
            if (cfg->top_side) {
                if (v > buf.front()) {
                    std::pop_heap(buf.begin(), buf.end(),
                                  std::greater<float>{});
                    buf.back() = v;
                    std::push_heap(buf.begin(), buf.end(),
                                   std::greater<float>{});
                }
            } else {
                if (v < buf.front()) {
                    std::pop_heap(buf.begin(), buf.end(),
                                  std::less<float>{});
                    buf.back() = v;
                    std::push_heap(buf.begin(), buf.end(),
                                   std::less<float>{});
                }
            }
        }

        void boundary() {}

        // Runtime hint to the driver: is pruning ever possible for this
        // State? Returning false lets the driver skip sliding-deque
        // maintenance entirely (each position saves 1-2 deque pushes +
        // advance_front calls — measurable on large fallback-mode scans).
        // For TopKQuantile, fallback mode never prunes, so maintaining
        // the max/min deques is pure overhead.
        bool pruning_active() const {
            return cfg && !cfg->use_fallback;
        }

        // Pruning: only meaningful in heap mode, and only once the heap is
        // full. For top-K (all percentiles >= 0.5), skip positions whose
        // upper bound is below the current K-th largest. For bottom-K,
        // skip positions whose lower bound is above the current K-th smallest.
        bool prune(float upper, float lower) const {
            if (!cfg || cfg->use_fallback || !heap_built) return false;
            return cfg->top_side ? (upper < buf.front())
                                 : (lower > buf.front());
        }

        // Merge: concat, then if in heap mode and over capacity, trim to K
        // via nth_element (cheaper than N·K heap pushes).
        void merge(const State &o) {
            n_total += o.n_total;
            if (!cfg && o.cfg) cfg = o.cfg;  // pick up Config on first merge

            if (cfg && !cfg->use_fallback) {
                buf.insert(buf.end(), o.buf.begin(), o.buf.end());
                if (buf.size() > (size_t)cfg->K) {
                    uint32_t K = cfg->K;
                    if (cfg->top_side) {
                        // Keep largest K: order buf so that the largest K
                        // occupy positions [N-K, N); then discard prefix.
                        std::nth_element(buf.begin(), buf.begin() + (buf.size() - K),
                                         buf.end());
                        buf.erase(buf.begin(), buf.begin() + (buf.size() - K));
                    } else {
                        std::nth_element(buf.begin(), buf.begin() + K - 1,
                                         buf.end());
                        buf.resize(K);
                    }
                }
                heap_built = false;  // re-heapify lazily on next accept
            } else {
                buf.insert(buf.end(), o.buf.begin(), o.buf.end());
            }
        }
    };

    struct Result { std::vector<double> quantile_vals; };

    // Pruning is meaningful only in heap mode; setting needs_pruning=true
    // unconditionally is safe (State::prune returns false in fallback mode
    // because `heap_built` stays false and the early-return fires). In
    // fallback mode the scan driver still maintains the sliding deques —
    // measurable overhead flagged as future work (addendum, Phase 2 review).
    static constexpr bool needs_pruning = true;
    // needs_lower_bound = true is required for bottom-K pruning
    // (top_side=false, all percentiles < 0.5), which consults `lower` in
    // State::prune. Spec originally had false; implementation flipped to
    // true to make bottom-K pruning actually fire. See addendum D-phase-2.
    static constexpr bool needs_lower_bound = true;
    // Quantiles need the actual value, not just a pass-flag, so the
    // certain-pass shortcut doesn't apply.
    static constexpr bool supports_certain_pass = false;
};

// ---------------------------------------------------------------------------
// finalize: compute the quantile value for each requested percentile, using
// the merged State::buf. In fallback mode buf holds all N values; in heap
// mode buf holds the K extreme values.
// ---------------------------------------------------------------------------
static std::vector<double> topk_finalize(TopKQuantile::State &s)
{
    const auto &pctiles = s.cfg->percentiles;
    std::vector<double> out(pctiles.size(),
                            std::numeric_limits<double>::quiet_NaN());
    int64_t N_total = (int64_t)s.n_total;     // all accepted (non-NaN) values
    int64_t N_buf   = (int64_t)s.buf.size();
    if (N_buf == 0) return out;
    for (size_t i = 0; i < pctiles.size(); ++i) {
        double p = pctiles[i];
        int64_t rank = (int64_t)std::floor(p * (double)(N_total - 1));
        if (rank < 0) rank = 0;
        if (rank >= N_total) rank = N_total - 1;

        int64_t buf_rank;
        if (s.cfg->use_fallback) {
            buf_rank = rank;
        } else {
            // Heap mode: buf holds K extreme values.
            //   top_side = true  -> buf holds the largest K values;
            //                       rank r among all N maps to buf index
            //                       r - (N - K).
            //   top_side = false -> buf holds the smallest K values;
            //                       rank r maps to buf index r.
            if (s.cfg->top_side)
                buf_rank = rank - (N_total - N_buf);
            else
                buf_rank = rank;
        }
        if (buf_rank < 0) buf_rank = 0;
        if (buf_rank >= N_buf) buf_rank = N_buf - 1;
        std::nth_element(s.buf.begin(), s.buf.begin() + buf_rank,
                         s.buf.end());
        out[i] = (double)s.buf[buf_rank];
    }
    return out;
}

// ---------------------------------------------------------------------------
// Parse `_func` SEXP into WindowAggFunc.
// ---------------------------------------------------------------------------
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
// .Call entry: C_gquantiles_multi
//
// Args (9):
//   1. _track_names   character vector
//   2. _percentiles   numeric vector
//   3. _iterator      integer scalar (bp step)
//   4. _sshift        integer scalar
//   5. _eshift        integer scalar
//   6. _n_threads     integer scalar (0 = auto)
//   7. _func          character scalar: "lse"|"avg"|"sum"|"max"|"min"
//   8. _intervals     data.frame (or NULL for whole-genome)
//   9. _envir         R environment
// ---------------------------------------------------------------------------
extern "C" SEXP C_gquantiles_multi(
    SEXP _track_names,
    SEXP _percentiles,
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

        int n_pctiles = Rf_length(_percentiles);
        std::vector<double> pctiles(REAL(_percentiles),
                                    REAL(_percentiles) + n_pctiles);
        for (double p : pctiles)
            if (p < 0.0 || p > 1.0) verror("percentile %g not in [0, 1]", p);

        int iterator_step = INTEGER(_iterator)[0];
        int sshift = INTEGER(_sshift)[0];
        int eshift = INTEGER(_eshift)[0];
        int n_threads_req = INTEGER(_n_threads)[0];
        if (iterator_step <= 0)
            verror("iterator must be a positive integer");

        WindowAggFunc func = parse_func(_func);

        const GenomeChromKey &chromkey = iu.get_chromkey();
        SEXP envir = iu.get_env();

        // Main-thread resolution of track paths and types.
        std::vector<std::string> track_dirs(n_tracks);
        std::vector<GenomeTrack::Type> track_types(n_tracks);
        for (int m = 0; m < n_tracks; ++m) {
            track_dirs[m] = track2path(envir, track_names[m]);
            track_types[m] = GenomeTrack::get_type(
                track_dirs[m].c_str(), chromkey, false);
        }

        // Build per-chrom intervals list (optional restriction).
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

        // Estimate N (positions per track) for adaptive K.
        uint64_t total_bp = 0;
        if (use_intervals) {
            for (auto &v : per_chrom)
                for (auto &g : v)
                    total_bp += (uint64_t)(g.end - g.start);
        } else {
            for (int c = 0; c < n_chroms; ++c)
                total_bp += chromkey.get_chrom_size(c);
        }
        uint64_t N_est = total_bp / (uint64_t)iterator_step;
        if (N_est == 0) N_est = 1;

        constexpr uint32_t K_MAX = 10'000'000u;

        double min_p = *std::min_element(pctiles.begin(), pctiles.end());
        double max_p = *std::max_element(pctiles.begin(), pctiles.end());
        bool all_top = (min_p >= 0.5);
        bool all_bot = (max_p <  0.5);
        bool mixed_tail = !(all_top || all_bot);

        double tail = 0.0;
        if (all_top)      tail = 1.0 - min_p;
        else if (all_bot) tail = max_p;
        double K_needed_f = std::ceil(tail * (double)N_est * 1.2);
        if (K_needed_f < 1.0) K_needed_f = 1.0;
        bool clamp = K_needed_f > (double)K_MAX;
        uint32_t K = clamp ? K_MAX
                           : (uint32_t)K_needed_f;

        bool use_fallback = mixed_tail || clamp;

        std::vector<TopKQuantile::Config> configs(n_tracks);
        for (int m = 0; m < n_tracks; ++m) {
            configs[m].percentiles = pctiles;
            configs[m].K = K;
            configs[m].use_fallback = use_fallback;
            configs[m].top_side = all_top;
        }

        if (mixed_tail)
            Rf_warning("percentiles span both tails (<0.5 and >=0.5); "
                       "top-K pruning disabled, falling back to full storage");
        if (clamp)
            Rf_warning("quantile K=%g exceeds K_MAX=%u; "
                       "falling back to full storage "
                       "(memory ~= 4 B * N_positions * n_tracks)",
                       K_needed_f, K_MAX);

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

        BatchTrackScanResult<TopKQuantile> scan_result;
        run_batch_scan<TopKQuantile>(track_names, track_dirs, track_types,
                                     configs, scan, chromkey, scan_result);

        for (int m = 0; m < n_tracks; ++m) {
            if (!scan_result.error_messages[m].empty())
                verror("Error processing track %s: %s",
                       track_names[m].c_str(),
                       scan_result.error_messages[m].c_str());
        }
        auto &per_track = scan_result.per_track_states;

        SEXP result;
        if (n_pctiles == 1)
            rprotect(result = Rf_allocVector(REALSXP, n_tracks));
        else
            rprotect(result = Rf_allocMatrix(REALSXP, n_tracks, n_pctiles));
        double *out = REAL(result);

        for (int m = 0; m < n_tracks; ++m) {
            auto qs = topk_finalize(per_track[m]);
            for (int p = 0; p < n_pctiles; ++p) {
                if (n_pctiles == 1)
                    out[m] = qs[p];
                else
                    out[m + (int64_t)p * n_tracks] = qs[p];
            }
        }

        if (n_pctiles == 1) {
            SEXP names;
            rprotect(names = Rf_allocVector(STRSXP, n_tracks));
            for (int i = 0; i < n_tracks; ++i)
                SET_STRING_ELT(names, i, STRING_ELT(_track_names, i));
            Rf_setAttrib(result, R_NamesSymbol, names);
        } else {
            SEXP rownames;
            rprotect(rownames = Rf_allocVector(STRSXP, n_tracks));
            for (int i = 0; i < n_tracks; ++i)
                SET_STRING_ELT(rownames, i, STRING_ELT(_track_names, i));

            SEXP colnames;
            rprotect(colnames = Rf_allocVector(STRSXP, n_pctiles));
            for (int i = 0; i < n_pctiles; ++i) {
                char buf[64];
                snprintf(buf, sizeof(buf), "%.6g", pctiles[i]);
                SET_STRING_ELT(colnames, i, Rf_mkChar(buf));
            }

            SEXP dimnames;
            rprotect(dimnames = Rf_allocVector(VECSXP, 2));
            SET_VECTOR_ELT(dimnames, 0, rownames);
            SET_VECTOR_ELT(dimnames, 1, colnames);
            Rf_setAttrib(result, R_DimNamesSymbol, dimnames);
        }

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
