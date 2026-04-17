#ifndef BATCHTRACKSCAN_TPP_
#define BATCHTRACKSCAN_TPP_

// Included from BatchTrackScan.h. Contains templated scan driver and
// per-track inner loops. Reducer methods are called via task.state.*.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>

#include "GenomeTrack1D.h"
#include "GenomeTrackFixedBin.h"
#include "GenomeTrackSparse.h"

namespace batchscan {

// Array-backed monotonic deque. CAP must be >= number of bins in the window.
// Entries are indexed by "bin index"; caller is responsible for pushing bins
// in increasing index order and advancing the front as the window slides.
template <int CAP, bool IsMax>
struct SlidingExtremum {
    int idx[CAP];
    float val[CAP];
    int head = 0;
    int tail = 0;

    void reset() { head = tail = 0; }

    void advance_front(int front_bin) {
        while (head < tail && idx[head % CAP] < front_bin) ++head;
    }

    void push_bin(int b, float v) {
        // Evict worse entries from the back: for max-deque, anything <= v is
        // dominated; for min-deque, anything >= v.
        while (head < tail) {
            float back = val[(tail - 1) % CAP];
            bool dominated = IsMax ? (back <= v) : (back >= v);
            if (!dominated) break;
            --tail;
        }
        // Caller is required to hold the live window to <= CAP bins (enforced
        // by the scan_fixedbin_inner guard at window setup time).
        // assert((size_t)(tail - head) < (size_t)CAP);
        idx[tail % CAP] = b;
        val[tail % CAP] = v;
        ++tail;
    }

    float current_extremum() const {
        if (head >= tail) {
            return IsMax ? -std::numeric_limits<float>::infinity()
                         :  std::numeric_limits<float>::infinity();
        }
        return val[head % CAP];
    }
};

// Window capacity: supports windows up to 256 bins (e.g. 256*20bp = 5120bp).
// No realistic vtrack call exceeds this.
static constexpr int WINDOW_CAP = 256;

using SlidingMax = SlidingExtremum<WINDOW_CAP, true>;
using SlidingMin = SlidingExtremum<WINDOW_CAP, false>;

// Open per-chrom handle. Called from worker thread — the main thread has
// already resolved track_dir and track_type (which are R-env operations).
inline void open_chrom_track(const std::string &track_dir,
                             GenomeTrack::Type type, int chromid,
                             const GenomeChromKey &chromkey,
                             std::shared_ptr<GenomeTrack> &out_owner,
                             GenomeTrackFixedBin **out_fb,
                             GenomeTrackSparse **out_sp)
{
    std::string resolved =
        GenomeTrack::find_existing_1d_filename(chromkey, track_dir, chromid);
    std::string filename = track_dir + "/" + resolved;
    *out_fb = nullptr;
    *out_sp = nullptr;
    if (type == GenomeTrack::FIXED_BIN) {
        auto t = std::make_shared<GenomeTrackFixedBin>();
        t->init_read(filename.c_str(), chromid);
        out_owner = t;
        *out_fb = t.get();
    } else if (type == GenomeTrack::SPARSE) {
        auto t = std::make_shared<GenomeTrackSparse>();
        t->init_read(filename.c_str(), chromid);
        out_owner = t;
        *out_sp = t.get();
    } else {
        throw std::runtime_error(std::string("track ") + track_dir +
                                 " has unsupported type (expected dense or sparse)");
    }
}

// -----------------------------------------------------------------------------
// FixedBin inner scan — sliding max/min + pruning gated on Reducer traits.
// -----------------------------------------------------------------------------
template <typename Reducer, WindowAggFunc F>
void scan_fixedbin_inner(GenomeTrackFixedBin *fb, unsigned bin_size,
                         int64_t chrom_size, int iterator_step,
                         int sshift, int eshift,
                         const std::vector<GInterval> *allowed_intervals,
                         typename Reducer::State &state)
{
    int64_t total_bins = (int64_t)((chrom_size + bin_size - 1) / bin_size);
    if (total_bins <= 0) return;

    int64_t out_count = 0;
    const float *all_bins = fb->get_mmap_bins_ptr(0, total_bins, out_count);
    if (!all_bins || out_count <= 0) return;

    const int n_bins_per_window = std::max(
        1, (int)(((eshift - sshift) + (int)bin_size - 1) / (int)bin_size));
    if (n_bins_per_window > WINDOW_CAP) {
        throw std::runtime_error(
            "BatchTrackScan: window exceeds WINDOW_CAP bins");
    }
    const float pre_const = aggregate_precomputed_const(F, n_bins_per_window);

    SlidingMax smax;
    SlidingMin smin;
    int next_bin_to_push = 0;   // next bin index to enter the deque(s)

    size_t interval_cursor = 0;
    bool prev_in_mask = false;

    for (int64_t c = 0; c < chrom_size; c += iterator_step) {
        // Interval-mask check with monotonic cursor.
        if (allowed_intervals) {
            while (interval_cursor < allowed_intervals->size() &&
                   (*allowed_intervals)[interval_cursor].end <= c)
                ++interval_cursor;
            bool in = (interval_cursor < allowed_intervals->size() &&
                       c >= (*allowed_intervals)[interval_cursor].start);
            if (!in) {
                if (prev_in_mask) state.boundary();
                prev_in_mask = false;
                if (Reducer::needs_pruning) {
                    smax.reset();
                    smin.reset();
                    next_bin_to_push = 0;  // rebuild on next in-mask position
                }
                continue;
            }
            prev_in_mask = true;
        }

        int64_t win_s = c + sshift;
        int64_t win_e = c + eshift;
        if (win_s < 0) win_s = 0;
        if (win_e <= win_s) continue;
        int64_t sbin = win_s / (int64_t)bin_size;
        int64_t ebin = (win_e + (int64_t)bin_size - 1) / (int64_t)bin_size;
        if (sbin < 0) sbin = 0;
        if (ebin > out_count) ebin = out_count;
        if (sbin >= ebin) continue;

        if constexpr (Reducer::needs_pruning) {
            // Forward-jump past bins in a mask gap: after the deque reset
            // on a gap, next_bin_to_push is 0, and we don't want to walk
            // through the skipped region — skip directly to sbin.
            if (next_bin_to_push < (int)sbin) next_bin_to_push = (int)sbin;
            // Push all bins up to ebin that we haven't seen yet.
            while (next_bin_to_push < (int)ebin) {
                int b = next_bin_to_push;
                float v = all_bins[b];
                smax.push_bin(b, std::isnan(v) ? -std::numeric_limits<float>::infinity() : v);
                if constexpr (Reducer::needs_lower_bound) {
                    smin.push_bin(b, std::isnan(v) ? std::numeric_limits<float>::infinity() : v);
                }
                ++next_bin_to_push;
            }
            // Advance fronts to current window start.
            smax.advance_front((int)sbin);
            if constexpr (Reducer::needs_lower_bound) {
                smin.advance_front((int)sbin);
            }

            float wmax = smax.current_extremum();
            float wmin = Reducer::needs_lower_bound
                         ? smin.current_extremum()
                         : std::numeric_limits<float>::quiet_NaN();
            float upper = aggregate_upper_bound(F, wmax, pre_const);
            float lower = Reducer::needs_lower_bound
                          ? aggregate_lower_bound(F, wmin, pre_const)
                          : std::numeric_limits<float>::quiet_NaN();
            if (state.prune(upper, lower)) {
                // Pruned position still counts as a valid sample if at
                // least one bin in the window is non-NaN (wmax > -inf).
                // Without this, the effective N used for rank
                // calculation in heap-mode reducers drifts off from the
                // fallback N.
                if (wmax > -std::numeric_limits<float>::infinity())
                    state.count_pruned();
                continue;
            }
        }

        float val = aggregate_window(F, all_bins, sbin, ebin);
        if (std::isnan(val)) state.nan_seen();
        else state.accept(val, c);
    }
}

// -----------------------------------------------------------------------------
// Sparse inner scan — no sliding-window optimization (irregular positions).
// Pruning is never triggered (bounds reported as +/-inf).
// Pattern mirrors bq_aggregate_lse_sparse from the original GlmBatchQuantiles.
// -----------------------------------------------------------------------------
template <typename Reducer, WindowAggFunc F>
void scan_sparse_inner(GenomeTrackSparse *sp, int64_t chrom_size,
                       int iterator_step, int sshift, int eshift,
                       const std::vector<GInterval> *allowed_intervals,
                       typename Reducer::State &state)
{
    const GIntervals &intervals = sp->get_intervals();
    const std::vector<float> &vals = sp->get_vals();

    size_t interval_cursor = 0;
    bool prev_in_mask = false;

    for (int64_t c = 0; c < chrom_size; c += iterator_step) {
        if (allowed_intervals) {
            while (interval_cursor < allowed_intervals->size() &&
                   (*allowed_intervals)[interval_cursor].end <= c)
                ++interval_cursor;
            bool in = (interval_cursor < allowed_intervals->size() &&
                       c >= (*allowed_intervals)[interval_cursor].start);
            if (!in) {
                if (prev_in_mask) state.boundary();
                prev_in_mask = false;
                continue;
            }
            prev_in_mask = true;
        }

        int64_t win_s = c + sshift;
        int64_t win_e = c + eshift;
        if (win_s < 0) win_s = 0;
        if (win_e <= win_s) continue;

        // Binary-search first interval whose end > win_s.
        size_t lo = 0, hi = intervals.size();
        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if ((int64_t)intervals[mid].end <= win_s) lo = mid + 1;
            else hi = mid;
        }

        // Scan overlapping sparse intervals and compute the aggregator inline.
        // We don't build a dense bins array here (sparse positions are
        // irregular); aggregate state is kept directly. SUM/AVG accumulate
        // in double for numerical stability, mirroring aggregate_window's
        // fixedbin path. No pruning on sparse — we always call accept().
        float acc_lse = -std::numeric_limits<float>::infinity();
        double acc_sum = 0.0;
        float acc_max = -std::numeric_limits<float>::infinity();
        float acc_min = std::numeric_limits<float>::infinity();
        uint64_t n = 0;
        for (size_t i = lo; i < intervals.size(); ++i) {
            if ((int64_t)intervals[i].start >= win_e) break;
            float v = vals[i];
            if (std::isnan(v)) continue;
            ++n;
            if constexpr (F == WindowAggFunc::LSE) lse_accumulate(acc_lse, v);
            else if constexpr (F == WindowAggFunc::SUM) acc_sum += v;
            else if constexpr (F == WindowAggFunc::AVG) acc_sum += v;
            else if constexpr (F == WindowAggFunc::MAX) { if (v > acc_max) acc_max = v; }
            else if constexpr (F == WindowAggFunc::MIN) { if (v < acc_min) acc_min = v; }
        }
        if (n == 0) {
            // Zero non-NaN bins in window → NaN aggregate. Still a valid
            // scan position; reducers that count bins (Summary) need to
            // see it via nan_seen(). Reducers that only care about
            // non-NaN values (TopKQuantile, ThresholdScreen behavior on
            // gap) ignore.
            state.nan_seen();
            continue;
        }

        float val;
        if constexpr (F == WindowAggFunc::LSE) val = acc_lse;
        else if constexpr (F == WindowAggFunc::SUM) val = (float)acc_sum;
        else if constexpr (F == WindowAggFunc::AVG) val = (float)(acc_sum / (double)n);
        else if constexpr (F == WindowAggFunc::MAX) val = acc_max;
        else if constexpr (F == WindowAggFunc::MIN) val = acc_min;

        // No pruning for sparse path.
        if (std::isnan(val)) state.nan_seen();
        else state.accept(val, c);
    }
}

// -----------------------------------------------------------------------------
// Dispatcher: picks the right aggregator F template based on ScanConfig.func.
// -----------------------------------------------------------------------------
template <typename Reducer>
void scan_task_dispatch(BatchTrackScanTask<Reducer> &task,
                        const ScanConfig &scan,
                        const GenomeChromKey &chromkey,
                        const std::vector<GInterval> *allowed_intervals)
{
    std::shared_ptr<GenomeTrack> owner;
    GenomeTrackFixedBin *fb = nullptr;
    GenomeTrackSparse *sp = nullptr;
    open_chrom_track(task.track_dir, task.track_type, task.chromid,
                     chromkey, owner, &fb, &sp);
    int64_t chrom_size = chromkey.get_chrom_size(task.chromid);

#define DISPATCH_F(FENUM)                                                      \
    case WindowAggFunc::FENUM:                                                  \
        if (fb) scan_fixedbin_inner<Reducer, WindowAggFunc::FENUM>(             \
                    fb, fb->get_bin_size(), chrom_size, scan.iterator_step,    \
                    scan.sshift, scan.eshift, allowed_intervals, task.state);   \
        else if (sp) scan_sparse_inner<Reducer, WindowAggFunc::FENUM>(          \
                    sp, chrom_size, scan.iterator_step, scan.sshift,            \
                    scan.eshift, allowed_intervals, task.state);                \
        break;

    switch (scan.func) {
        DISPATCH_F(LSE)
        DISPATCH_F(AVG)
        DISPATCH_F(SUM)
        DISPATCH_F(MAX)
        DISPATCH_F(MIN)
    }
#undef DISPATCH_F

    task.state.boundary();  // flush any pending run at chrom end
}

// -----------------------------------------------------------------------------
// Top-level driver. Tasks = all (track, chrom) pairs with non-zero chrom_size.
// Work queue via std::atomic<size_t>. Workers catch all exceptions locally.
//
// Memory discipline: each worker, after completing a (track, chrom) task,
// merges its state into the corresponding per-track accumulator (under a
// per-track mutex) and frees the task's state. This keeps peak memory
// bounded by (per_track_accumulator * n_tracks + concurrent_task_buffers *
// n_threads), not by (all_task_buffers summed). Critical for the
// Phase 1 fallback mode where task state holds the full value vector.
// -----------------------------------------------------------------------------
template <typename Reducer>
struct BatchTrackScanResult {
    std::vector<typename Reducer::State> per_track_states;   // length n_tracks
    std::vector<std::string> error_messages;                 // one per track, empty => OK
};

template <typename Reducer>
void run_batch_scan(
    const std::vector<std::string> &track_names,
    const std::vector<std::string> &track_dirs,
    const std::vector<GenomeTrack::Type> &track_types,
    const std::vector<typename Reducer::Config> &per_track_configs,
    const ScanConfig &scan,
    const GenomeChromKey &chromkey,
    BatchTrackScanResult<Reducer> &out)
{
    const int n_tracks = (int)track_names.size();
    const int n_chroms = (int)chromkey.get_num_chroms();

    // Per-track accumulators + per-track mutexes. Accumulators are
    // initialized with the track's Config so their cfg pointer outlives
    // any merge call.
    out.per_track_states.clear();
    out.per_track_states.resize(n_tracks);
    out.error_messages.assign(n_tracks, std::string());
    for (int t = 0; t < n_tracks; ++t)
        out.per_track_states[t].init(per_track_configs[t], /*chromid=*/0,
                                     scan.iterator_step);
    std::vector<std::mutex> per_track_mutex(n_tracks);

    // Build (track, chrom) task list.
    std::vector<BatchTrackScanTask<Reducer>> tasks;
    tasks.reserve((size_t)n_tracks * n_chroms);
    for (int t = 0; t < n_tracks; ++t) {
        for (int c = 0; c < n_chroms; ++c) {
            if (chromkey.get_chrom_size(c) == 0) continue;
            std::string resolved;
            try {
                resolved = GenomeTrack::find_existing_1d_filename(
                    chromkey, track_dirs[t], c);
            } catch (...) {
                continue;  // no file for this (track, chrom) → skip task
            }
            if (resolved.empty()) continue;

            BatchTrackScanTask<Reducer> task;
            task.track_name = track_names[t];
            task.track_dir = track_dirs[t];
            task.track_type = track_types[t];
            task.track_idx = t;
            task.chromid = c;
            task.state.init(per_track_configs[t], c, scan.iterator_step);
            tasks.push_back(std::move(task));
        }
    }

    std::atomic<size_t> next_task{0};
    const size_t n_tasks = tasks.size();

    auto worker = [&]() {
        for (;;) {
            size_t i = next_task.fetch_add(1);
            if (i >= n_tasks) return;
            auto &task = tasks[i];
            try {
                const std::vector<GInterval> *allowed =
                    scan.per_chrom_intervals
                        ? &(*scan.per_chrom_intervals)[task.chromid]
                        : nullptr;
                scan_task_dispatch<Reducer>(task, scan, chromkey, allowed);
                // Merge into the per-track accumulator and release task state.
                {
                    std::lock_guard<std::mutex> lk(per_track_mutex[task.track_idx]);
                    out.per_track_states[task.track_idx].merge(task.state);
                }
                task.state = typename Reducer::State{};  // free buffers
            } catch (const std::exception &e) {
                task.error_msg = e.what();
            } catch (...) {
                task.error_msg = "unknown error in worker thread";
            }
        }
    };

    int n_threads = scan.n_threads;
    if (n_threads < 1) n_threads = 1;
    if ((size_t)n_threads > n_tasks) n_threads = (int)n_tasks;
    if (n_threads > 0) {
        std::vector<std::thread> threads;
        threads.reserve(n_threads);
        for (int i = 0; i < n_threads; ++i) threads.emplace_back(worker);
        for (auto &t : threads) t.join();
    }

    // Surface the first error per track (if any).
    for (auto &task : tasks) {
        if (!task.error_msg.empty() && out.error_messages[task.track_idx].empty())
            out.error_messages[task.track_idx] = task.error_msg;
    }
}

}  // namespace batchscan

#endif  // BATCHTRACKSCAN_TPP_
