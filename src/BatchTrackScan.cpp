#include "BatchTrackScan.h"
#include "GenomeTrack1D.h"  // lse_accumulate

#include <cmath>
#include <limits>

namespace batchscan {

float aggregate_window(WindowAggFunc func, const float *bins, int64_t sbin,
                       int64_t ebin)
{
    switch (func) {
        case WindowAggFunc::LSE: {
            float lse = -std::numeric_limits<float>::infinity();
            uint64_t n = 0;
            for (int64_t b = sbin; b < ebin; ++b) {
                if (!std::isnan(bins[b])) { lse_accumulate(lse, bins[b]); ++n; }
            }
            return n ? lse : std::numeric_limits<float>::quiet_NaN();
        }
        case WindowAggFunc::AVG: {
            double s = 0; uint64_t n = 0;
            for (int64_t b = sbin; b < ebin; ++b)
                if (!std::isnan(bins[b])) { s += bins[b]; ++n; }
            return n ? (float)(s / (double)n)
                     : std::numeric_limits<float>::quiet_NaN();
        }
        case WindowAggFunc::SUM: {
            double s = 0; uint64_t n = 0;
            for (int64_t b = sbin; b < ebin; ++b)
                if (!std::isnan(bins[b])) { s += bins[b]; ++n; }
            return n ? (float)s : std::numeric_limits<float>::quiet_NaN();
        }
        case WindowAggFunc::MAX: {
            float m = -std::numeric_limits<float>::infinity();
            uint64_t n = 0;
            for (int64_t b = sbin; b < ebin; ++b)
                if (!std::isnan(bins[b])) { if (bins[b] > m) m = bins[b]; ++n; }
            return n ? m : std::numeric_limits<float>::quiet_NaN();
        }
        case WindowAggFunc::MIN: {
            float m = std::numeric_limits<float>::infinity();
            uint64_t n = 0;
            for (int64_t b = sbin; b < ebin; ++b)
                if (!std::isnan(bins[b])) { if (bins[b] < m) m = bins[b]; ++n; }
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
        case WindowAggFunc::SUM: return window_max * precomputed_const;
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
        case WindowAggFunc::SUM: return window_min * precomputed_const;
        case WindowAggFunc::AVG: return window_min;
        case WindowAggFunc::MIN: return window_min;
        case WindowAggFunc::MAX: return -std::numeric_limits<float>::infinity();
    }
    return -std::numeric_limits<float>::infinity();
}

}  // namespace batchscan
