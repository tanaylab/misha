/*
 * TrackIndex2D.h
 *
 * Index for single-file 2D track storage
 * Maps chromosome pair (chrom1_id, chrom2_id) to offset/length in track.dat
 *
 * Format specification: track.idx v1.0 (2D)
 * - Global header (44 bytes):
 *   - Magic: "MISHT2D\0" (8 bytes)
 *   - Version: uint32_t (4 bytes) -- value 1
 *   - TrackType: uint32_t (4 bytes) -- 0=RECTS, 1=POINTS
 *   - NumPairs: uint32_t (4 bytes) -- number of entries in entry table
 *   - Flags: uint64_t (8 bytes) -- bit 0: IS_LITTLE_ENDIAN
 *   - Checksum: uint64_t (8 bytes) -- CRC64-ECMA of entry table
 *   - Reserved: uint64_t (8 bytes) -- for future use
 *
 * - Per-Pair Entry Table (28 bytes per entry):
 *   - chrom1_id: uint32_t (4 bytes)
 *   - chrom2_id: uint32_t (4 bytes)
 *   - offset: uint64_t (8 bytes)
 *   - length: uint64_t (8 bytes)
 *   - reserved: uint32_t (4 bytes)
 */

#ifndef TRACKINDEX2D_H_
#define TRACKINDEX2D_H_

#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// Track type enumeration for 2D tracks
enum class MishaTrack2DType : uint32_t {
    RECTS  = 0,   // Rectangle track
    POINTS = 1,   // Points track
};

// Index entry for a single chromosome pair
struct Track2DPairEntry {
    uint32_t chrom1_id;    // First chromosome ID from GenomeChromKey
    uint32_t chrom2_id;    // Second chromosome ID from GenomeChromKey
    uint64_t offset;       // Byte offset in track.dat
    uint64_t length;       // Byte length of pair data (0 if empty)
    uint32_t reserved;     // Reserved for future use

    Track2DPairEntry() : chrom1_id(0), chrom2_id(0), offset(0), length(0), reserved(0) {}
    Track2DPairEntry(uint32_t id1, uint32_t id2, uint64_t off, uint64_t len)
        : chrom1_id(id1), chrom2_id(id2), offset(off), length(len), reserved(0) {}
};

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!
class TrackIndex2D {
public:
    enum Errors {
        FILE_READ_FAILED,
        INVALID_FORMAT,
        VERSION_MISMATCH,
        CHECKSUM_FAILED,
        ENDIAN_MISMATCH
    };

    TrackIndex2D();
    ~TrackIndex2D();

    // Load index from file
    // Returns true on success, false if file doesn't exist
    // Throws TGLException on format errors
    bool load(const string &index_path);

    // Get index entry for a chromosome pair
    // Returns nullptr if pair not found
    const Track2DPairEntry* get_entry(uint32_t chrom1_id, uint32_t chrom2_id) const;

    // Check if a chromosome pair exists in the index
    bool has_entry(uint32_t chrom1_id, uint32_t chrom2_id) const;

    // Get track type from index
    MishaTrack2DType get_track_type() const { return m_track_type; }

    // Check if index is loaded
    bool is_loaded() const { return m_loaded; }

    // Get all entries (for iteration)
    const vector<Track2DPairEntry>& entries() const { return m_entries; }

    // Get number of pairs
    size_t num_entries() const { return m_entries.size(); }

    // --- Static cache management ---

    // Get-or-load the 2D track index for a track directory (thread-safe)
    static std::shared_ptr<TrackIndex2D> get_track_index_2d(const std::string &track_dir);

    // Clear the index cache
    static void clear_cache();

    // Write a 2D index file
    // entries must be populated with chrom1_id, chrom2_id, offset, length
    static void write_index(const string &index_path,
                            MishaTrack2DType track_type,
                            const vector<Track2DPairEntry> &entries);

private:
    bool m_loaded;
    MishaTrack2DType m_track_type;
    vector<Track2DPairEntry> m_entries;
    unordered_map<uint64_t, size_t> m_pair_to_index;  // key = (chrom1_id << 32) | chrom2_id

    // Magic header and version
    static const char MAGIC_HEADER[8];  // "MISHT2D\0"
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
    static uint64_t compute_checksum(const vector<Track2DPairEntry> &entries);

    // Static index cache (thread-safe)
    static std::map<std::string, std::shared_ptr<TrackIndex2D>> s_index_cache;
    static std::mutex s_cache_mutex;
};

#endif /* TRACKINDEX2D_H_ */
