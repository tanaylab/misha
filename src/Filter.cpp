#include "Filter.h"
#include <algorithm>

void Genome1DFilter::subtract(const GInterval& I, std::vector<GInterval>& out) const {
    out.clear();

    // Validate chromid
    if (I.chromid < 0 || (size_t)I.chromid >= per_chrom.size()) {
        out.push_back(I);
        return;
    }

    const auto& mask_vec = per_chrom[I.chromid];
    if (mask_vec.empty()) {
        out.push_back(I);
        return;
    }

    // Find first mask interval that could overlap with I
    // lower_bound finds first mask where mask.end > I.start
    auto it = std::lower_bound(mask_vec.begin(), mask_vec.end(), I.start,
                               [](const GInterval& m, int64_t pos) {
                                   return m.end <= pos;
                               });

    int64_t cur = I.start;

    // Process overlapping mask intervals
    for (; it != mask_vec.end() && it->start < I.end; ++it) {
        // Emit unmasked region before this mask
        if (it->start > cur) {
            out.emplace_back(I.chromid, cur, std::min<int64_t>(it->start, I.end), I.strand);
        }
        // Advance cursor past the mask
        if (it->end > cur) {
            cur = it->end;
        }
        if (cur >= I.end) break;
    }

    // Emit final unmasked region
    if (cur < I.end) {
        out.emplace_back(I.chromid, cur, I.end, I.strand);
    }
}

void Genome1DFilter::subtract_with_cursor(const GInterval& I, Cursor& cursor, std::vector<GInterval>& out) const {
    out.clear();

    // Validate chromid
    if (I.chromid < 0 || (size_t)I.chromid >= per_chrom.size()) {
        out.push_back(I);
        return;
    }

    const auto& mask_vec = per_chrom[I.chromid];
    if (mask_vec.empty()) {
        out.push_back(I);
        return;
    }

    // Reset cursor if we changed chromosomes or moved backwards
    if (cursor.chromid != I.chromid || cursor.idx >= mask_vec.size() ||
        (cursor.idx > 0 && mask_vec[cursor.idx].start > I.start)) {
        cursor.chromid = I.chromid;
        cursor.idx = 0;
    }

    // Find first mask interval that could overlap with I
    while (cursor.idx < mask_vec.size() && mask_vec[cursor.idx].end <= I.start) {
        cursor.idx++;
    }

    int64_t cur = I.start;
    size_t idx = cursor.idx;

    // Process overlapping mask intervals
    for (; idx < mask_vec.size() && mask_vec[idx].start < I.end; ++idx) {
        // Emit unmasked region before this mask
        if (mask_vec[idx].start > cur) {
            out.emplace_back(I.chromid, cur, std::min<int64_t>(mask_vec[idx].start, I.end), I.strand);
        }
        // Advance cursor past the mask
        if (mask_vec[idx].end > cur) {
            cur = mask_vec[idx].end;
        }
        if (cur >= I.end) break;
    }

    // Emit final unmasked region
    if (cur < I.end) {
        out.emplace_back(I.chromid, cur, I.end, I.strand);
    }
}

int64_t Genome1DFilter::total_masked_bases() const {
    int64_t total = 0;
    for (const auto& chrom_masks : per_chrom) {
        for (const auto& mask : chrom_masks) {
            total += mask.range();
        }
    }
    return total;
}

size_t Genome1DFilter::num_masked_chroms() const {
    size_t count = 0;
    for (const auto& chrom_masks : per_chrom) {
        if (!chrom_masks.empty()) {
            count++;
        }
    }
    return count;
}
