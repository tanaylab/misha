/*
 * TrackIndex.h
 *
 * Index for single-file track storage
 * Maps contig ID to offset/length in track.dat
 *
 * Format specification: track.idx v2.0
 * - Global header (36 bytes):
 *   - Magic: "MISHATDX" (8 bytes)
 *   - Version: uint32_t (4 bytes)
 *   - TrackType: uint32_t (4 bytes) - 0=dense, 1=sparse, 2=array
 *   - NumContigs: uint32_t (4 bytes) - total chromosome count (may include empty chromosomes)
 *                                       Number of entries in entry table matches genome size
 *   - Flags: uint64_t (8 bytes) - bit 0: IS_LITTLE_ENDIAN
 *   - Checksum: uint64_t (8 bytes) - CRC64-ECMA of entry table
 *
 * - Per-Contig Entry Table (24 bytes per entry):
 *   - chrom_id: uint32_t (4 bytes)
 *   - offset: uint64_t (8 bytes)
 *   - length: uint64_t (8 bytes)
 *   - reserved: uint32_t (4 bytes)
 */

#ifndef TRACKINDEX_H_
#define TRACKINDEX_H_

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// Track type enumeration (matches GenomeTrack::Type order for 1D tracks)
enum class MishaTrackType : uint32_t {
    DENSE = 0,   // Fixed-bin / dense track
    SPARSE = 1,  // Sparse track
    ARRAY = 2,   // Array track
    // Values 3-255 reserved for future 2D types
};

// Index entry for a single contig
struct TrackContigEntry {
    uint32_t chrom_id;     // Chromosome ID from GenomeChromKey
    uint64_t offset;       // Byte offset in track.dat
    uint64_t length;       // Byte length of contig data (0 if empty)
    uint32_t reserved;     // Reserved for future use

    TrackContigEntry() : chrom_id(0), offset(0), length(0), reserved(0) {}
    TrackContigEntry(uint32_t id, uint64_t off, uint64_t len)
        : chrom_id(id), offset(off), length(len), reserved(0) {}
};

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!
class TrackIndex {
public:
    enum Errors {
        FILE_READ_FAILED,
        INVALID_FORMAT,
        VERSION_MISMATCH,
        CHECKSUM_FAILED,
        ENDIAN_MISMATCH
    };

    TrackIndex();
    ~TrackIndex();

    // Load index from file
    // Returns true on success, false if file doesn't exist
    // Throws TGLException on format errors
    bool load(const string &index_path);

    // Get index entry for a chromosome
    // Returns nullptr if chromosome not found
    const TrackContigEntry* get_entry(uint32_t chromid) const;

    // Get track type from index
    MishaTrackType get_track_type() const { return m_track_type; }

    // Check if index is loaded
    bool is_loaded() const { return m_loaded; }

    // Get all entries (for iteration)
    const vector<TrackContigEntry>& get_all_entries() const { return m_entries; }

    // Get number of contigs
    size_t get_num_contigs() const { return m_entries.size(); }

private:
    bool m_loaded;
    MishaTrackType m_track_type;
    vector<TrackContigEntry> m_entries;
    unordered_map<uint32_t, size_t> m_chromid_to_index;

    // Magic header and version
    static const char MAGIC_HEADER[8];
    static const uint32_t INDEX_VERSION = 1;

    // Flags
    static const uint64_t FLAG_LITTLE_ENDIAN = 0x01;

    // Helper: check if system is little-endian
    static bool is_little_endian();

    // Helper: compute checksum of index entries
    uint64_t compute_checksum(const vector<TrackContigEntry> &entries);
};

#endif /* TRACKINDEX_H_ */
