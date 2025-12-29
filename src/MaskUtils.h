/*
 * MaskUtils.h
 *
 * Utility functions for working with interval masks during
 * genome scanning operations (training, sampling, k-mer counting).
 */

#ifndef MASK_UTILS_H_
#define MASK_UTILS_H_

#include <cstddef>
#include <cstdint>
#include <vector>

#include "GInterval.h"

/**
 * Check if a position is within any masked interval.
 *
 * Uses a cursor for efficient sequential access when scanning positions
 * in order. The cursor is advanced past intervals that end before the
 * current position.
 *
 * Note on overlapping intervals: If mask intervals overlap, the cursor
 * will still work correctly but may not find all overlapping intervals.
 * For best results, sort and merge overlapping mask intervals before use.
 *
 * @param pos Genomic position to check
 * @param mask_intervals Sorted vector of mask intervals (by start position)
 * @param cursor [in/out] Current position in the mask_intervals vector.
 *               Updated to skip past intervals that end before pos.
 * @return true if pos is within any mask interval, false otherwise
 */
inline bool is_position_masked(int64_t pos,
                               const std::vector<GInterval>& mask_intervals,
                               size_t& cursor) {
    // Advance cursor past intervals that end before this position
    while (cursor < mask_intervals.size() &&
           mask_intervals[cursor].end <= pos) {
        ++cursor;
    }

    // Check if position is within current interval
    if (cursor < mask_intervals.size() &&
        mask_intervals[cursor].start <= pos &&
        pos < mask_intervals[cursor].end) {
        return true;
    }

    return false;
}

#endif // MASK_UTILS_H_
