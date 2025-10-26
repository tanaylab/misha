/*
 * IntervalsIndex2D.h
 *
 * Index for single-file 2D interval set storage
 * Maps chromosome pair (chrom1_id, chrom2_id) to offset/length in intervals2d.dat
 *
 * Format specification: intervals2d.idx v1.0
 * - Global header (40 bytes):
 *   - Magic: "MISHAI2D" (8 bytes)
 *   - Version: uint32_t (4 bytes)
 *   - NumEntries: uint32_t (4 bytes) - number of chromosome pairs (may include empty pairs)
 *                                       Actual number of pairs found in per-chromosome format
 *   - Flags: uint64_t (8 bytes) - bit 0: IS_LITTLE_ENDIAN
 *   - Checksum: uint64_t (8 bytes) - CRC64-ECMA of entry table
 *   - Reserved: uint64_t (8 bytes)
 *
 * - Per-Pair Entry Table (28 bytes per entry):
 *   - chrom1_id: uint32_t (4 bytes)
 *   - chrom2_id: uint32_t (4 bytes)
 *   - offset: uint64_t (8 bytes)
 *   - length: uint64_t (8 bytes)
 *   - reserved: uint32_t (4 bytes)
 */

#ifndef INTERVALSINDEX2D_H_
#define INTERVALSINDEX2D_H_

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// Index entry for a single chromosome pair
struct IntervalsPairEntry {
    uint32_t chrom1_id;    // First chromosome ID from GenomeChromKey
    uint32_t chrom2_id;    // Second chromosome ID from GenomeChromKey
    uint64_t offset;       // Byte offset in intervals2d.dat
    uint64_t length;       // Byte length of pair data (0 if empty)
    uint32_t reserved;     // Reserved for future use

    IntervalsPairEntry() : chrom1_id(0), chrom2_id(0), offset(0), length(0), reserved(0) {}
    IntervalsPairEntry(uint32_t id1, uint32_t id2, uint64_t off, uint64_t len)
        : chrom1_id(id1), chrom2_id(id2), offset(off), length(len), reserved(0) {}
};

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!
class IntervalsIndex2D {
public:
    enum Errors {
        FILE_READ_FAILED,
        INVALID_FORMAT,
        VERSION_MISMATCH,
        CHECKSUM_FAILED,
        ENDIAN_MISMATCH
    };

    IntervalsIndex2D();
    ~IntervalsIndex2D();

    // Load index from file
    // Returns true on success, false if file doesn't exist
    // Throws TGLException on format errors
    bool load(const string &index_path);

    // Get index entry for a chromosome pair
    // Returns nullptr if pair not found
    const IntervalsPairEntry* get_entry(uint32_t chrom1_id, uint32_t chrom2_id) const;

    // Check if index is loaded
    bool is_loaded() const { return m_loaded; }

    // Get all entries (for iteration)
    const vector<IntervalsPairEntry>& get_all_entries() const { return m_entries; }

    // Get number of pairs
    size_t get_num_pairs() const { return m_entries.size(); }

private:
    bool m_loaded;
    vector<IntervalsPairEntry> m_entries;
    unordered_map<uint64_t, size_t> m_pair_to_index;  // key = (chrom1_id << 32) | chrom2_id

    // Magic header and version
    static const char MAGIC_HEADER[8];
    static const uint32_t INDEX_VERSION = 1;

    // Flags
    static const uint64_t FLAG_LITTLE_ENDIAN = 0x01;

    // Helper: create lookup key from pair
    static uint64_t make_pair_key(uint32_t chrom1_id, uint32_t chrom2_id) {
        return (static_cast<uint64_t>(chrom1_id) << 32) | chrom2_id;
    }

    // Helper: check if system is little-endian
    static bool is_little_endian();

    // Helper: compute checksum of index entries
    uint64_t compute_checksum(const vector<IntervalsPairEntry> &entries);
};

#endif /* INTERVALSINDEX2D_H_ */
