#ifndef FILTER_H_
#define FILTER_H_

#include <cstdint>
#include <vector>
#include "GInterval.h"

/**
 * Genome1DFilter: Represents a genomic mask stored as sorted, non-overlapping intervals per chromosome.
 * Used to exclude masked regions from virtual track evaluation.
 */
struct Genome1DFilter {
    // Per-chromosome sorted, merged mask intervals (half-open [start, end))
    std::vector<std::vector<GInterval>> per_chrom;
    bool empty = true;

    /**
     * Subtract mask from interval I, returning the unmasked parts.
     * @param I Input interval (must be valid with chromid < per_chrom.size())
     * @param out Output vector to store resulting unmasked intervals (same chrom as I)
     */
    void subtract(const GInterval& I, std::vector<GInterval>& out) const;

    /**
     * Cursor for efficient sequential subtraction on sorted intervals.
     * Maintains position in mask for faster lookup when intervals are processed in order.
     */
    struct Cursor {
        size_t idx = 0;
        int chromid = -1;
    };

    /**
     * Subtract mask from interval I using a cursor for efficiency.
     * @param I Input interval
     * @param cursor Maintains state across calls for sequential intervals
     * @param out Output vector to store resulting unmasked intervals
     */
    void subtract_with_cursor(const GInterval& I, Cursor& cursor, std::vector<GInterval>& out) const;

    /**
     * Get total number of masked bases across all chromosomes.
     */
    int64_t total_masked_bases() const;

    /**
     * Get number of chromosomes with masks.
     */
    size_t num_masked_chroms() const;
};

#endif // FILTER_H_
