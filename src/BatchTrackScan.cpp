#include "BatchTrackScan.h"
#include "GenomeTrack1D.h"  // lse_accumulate

#include <cmath>
#include <limits>

namespace batchscan {

float aggregate_window(WindowAggFunc func, const float *bins, int64_t sbin,
                       int64_t ebin)
{
    // The legacy read path (GenomeTrackFixedBin::read_next_bin) converts
    // +/-inf bins to NaN, so a non-finite bin is "missing" and excluded from
    // every reducer. We read raw mmap bins here, so apply the same rule via
    // std::isfinite (false for NaN and inf) to stay bit-identical to the slow
    // path on tracks that contain inf bins.
    switch (func) {
        case WindowAggFunc::LSE: {
            float lse = -std::numeric_limits<float>::infinity();
            uint64_t n = 0;
            for (int64_t b = sbin; b < ebin; ++b) {
                if (std::isfinite(bins[b])) { lse_accumulate(lse, bins[b]); ++n; }
            }
            return n ? lse : std::numeric_limits<float>::quiet_NaN();
        }
        case WindowAggFunc::AVG: {
            double s = 0; uint64_t n = 0;
            for (int64_t b = sbin; b < ebin; ++b)
                if (std::isfinite(bins[b])) { s += bins[b]; ++n; }
            return n ? (float)(s / (double)n)
                     : std::numeric_limits<float>::quiet_NaN();
        }
        case WindowAggFunc::SUM: {
            double s = 0; uint64_t n = 0;
            for (int64_t b = sbin; b < ebin; ++b)
                if (std::isfinite(bins[b])) { s += bins[b]; ++n; }
            return n ? (float)s : std::numeric_limits<float>::quiet_NaN();
        }
        case WindowAggFunc::MAX: {
            float m = -std::numeric_limits<float>::infinity();
            uint64_t n = 0;
            for (int64_t b = sbin; b < ebin; ++b)
                if (std::isfinite(bins[b])) { if (bins[b] > m) m = bins[b]; ++n; }
            return n ? m : std::numeric_limits<float>::quiet_NaN();
        }
        case WindowAggFunc::MIN: {
            float m = std::numeric_limits<float>::infinity();
            uint64_t n = 0;
            for (int64_t b = sbin; b < ebin; ++b)
                if (std::isfinite(bins[b])) { if (bins[b] < m) m = bins[b]; ++n; }
            return n ? m : std::numeric_limits<float>::quiet_NaN();
        }
    }
    return std::numeric_limits<float>::quiet_NaN();
}

float aggregate_precomputed_const(WindowAggFunc func, int n_bins)
{
    switch (func) {
        case WindowAggFunc::LSE: return std::log((float)n_bins);
        case WindowAggFunc::SUM: return (float)n_bins;
        default: return 0.0f;
    }
}

// Per-aggregator upper bound derived from the window max.
// Used by reducers that prune extreme quantiles/threshold checks; tight but
// cheap. For MIN we return +inf because max-based bound is uninformative.
float aggregate_upper_bound(WindowAggFunc func, float window_max,
                            float precomputed_const)
{
    switch (func) {
        case WindowAggFunc::LSE: return window_max + precomputed_const;
        // SUM upper bound. precomputed_const is an OVER-estimate of the bin
        // count (ceil(W/B)+1, and the true non-NaN count can be lower still
        // for partial/NaN bins). window_max * count is only a valid upper
        // bound when window_max >= 0 (more bins -> larger sum); when
        // window_max < 0 the largest possible sum is achieved with the FEWEST
        // bins, so sum <= 0. The universally valid bound is
        // max(window_max, 0) * count. (Was window_max * count, which is not a
        // bound for negative-valued tracks.)
        case WindowAggFunc::SUM:
            return (window_max > 0.0f ? window_max : 0.0f) * precomputed_const;
        case WindowAggFunc::AVG: return window_max;
        case WindowAggFunc::MAX: return window_max;
        case WindowAggFunc::MIN: return std::numeric_limits<float>::infinity();
    }
    return std::numeric_limits<float>::infinity();
}

// Per-aggregator lower bound derived from the window min.
// Used for bottom-percentile / LT / LE pruning. For LSE, window_min is a
// valid lower bound because LSE >= max >= min of the non-NaN bins. For MAX
// we return -inf (max-based aggregator has no useful min lower bound).
float aggregate_lower_bound(WindowAggFunc func, float window_min,
                            float precomputed_const)
{
    switch (func) {
        case WindowAggFunc::LSE: return window_min;
        // SUM lower bound. Mirror of the upper bound: window_min * count is a
        // valid lower bound only when window_min <= 0; when window_min > 0 the
        // smallest possible sum uses the FEWEST bins (0) -> sum >= 0. The
        // universally valid bound is min(window_min, 0) * count. (Was
        // window_min * count, which over-bounds even an all-positive track
        // whenever the count is over-estimated, wrongly pruning LT/LE passes.)
        case WindowAggFunc::SUM:
            return (window_min < 0.0f ? window_min : 0.0f) * precomputed_const;
        case WindowAggFunc::AVG: return window_min;
        case WindowAggFunc::MIN: return window_min;
        case WindowAggFunc::MAX: return -std::numeric_limits<float>::infinity();
    }
    return -std::numeric_limits<float>::infinity();
}

}  // namespace batchscan
