#ifndef BATCHTRACKSCAN_H_
#define BATCHTRACKSCAN_H_

// Shared skeleton for batched multi-track scans used by:
//   - BatchQuantiles.cpp (TopKQuantile reducer)
//   - BatchSummary.cpp   (Summary reducer, Phase 3)
//   - BatchScreen.cpp    (ThresholdScreen reducer, Phase 4)
//
// The scan driver is templated over a Reducer concept. Reducers must expose:
//   - struct Config, State, Result
//   - static constexpr bool needs_pruning       (enables sliding-max + prune())
//   - static constexpr bool needs_lower_bound   (also enables sliding-min)
//   - State::init(const Config&, int chromid, int iterator_step)
//   - State::accept(float val, int64_t pos)
//   - State::boundary()                          (flush at mask gap / chrom end)
//   - State::prune(float upper, float lower)     (true => skip position)
//   - State::merge(const State& other)
//
// Thread-safety: workers MUST NOT call into the R C API. Everything that
// touches R (track name resolution, intervals conversion, SEXP allocation)
// happens on the main thread before/after run_batch_scan.

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "GenomeChromKey.h"
#include "GenomeTrack.h"
#include "GInterval.h"

namespace batchscan {

enum class WindowAggFunc { LSE, AVG, SUM, MAX, MIN };

struct ScanConfig {
    WindowAggFunc func;
    int iterator_step;
    int sshift;
    int eshift;
    // nullptr => whole genome. Otherwise indexed by chromid; each inner vector
    // is the sorted, non-overlapping intervals restricting the scan on that
    // chromosome. Positions whose scan center falls outside all intervals are
    // skipped (and reducer.boundary() is called at the gap).
    const std::vector<std::vector<GInterval>> *per_chrom_intervals;
    int n_threads;
};

// One task per (track, chrom). Workers write into State; main thread merges
// per-track states after join.
template <typename Reducer>
struct BatchTrackScanTask {
    std::string track_name;
    std::string track_dir;
    GenomeTrack::Type track_type;
    int track_idx;              // index into the original track_names vector
    int chromid;
    typename Reducer::State state;
    std::string error_msg;      // set by worker on exception; checked on main
};

// Aggregator math helpers — safe to call from any thread (pure, no R).
float aggregate_window(WindowAggFunc func, const float *bins, int64_t sbin,
                       int64_t ebin);

float aggregate_upper_bound(WindowAggFunc func, float window_max,
                            float precomputed_const);

float aggregate_lower_bound(WindowAggFunc func, float window_min,
                            float precomputed_const);

// log(n_bins) for LSE, (float)n_bins for SUM, 0 otherwise.
float aggregate_precomputed_const(WindowAggFunc func, int n_bins);

// Top-level driver. `track_dirs` and `track_types` must be resolved on the
// main thread before calling; workers never touch R to resolve names.
template <typename Reducer>
void run_batch_scan(
    const std::vector<std::string> &track_names,
    const std::vector<std::string> &track_dirs,
    const std::vector<GenomeTrack::Type> &track_types,
    const std::vector<typename Reducer::Config> &per_track_configs,
    const ScanConfig &scan,
    const GenomeChromKey &chromkey,
    std::vector<BatchTrackScanTask<Reducer>> &out_tasks);

}  // namespace batchscan

#include "BatchTrackScan.tpp"

#endif  // BATCHTRACKSCAN_H_
