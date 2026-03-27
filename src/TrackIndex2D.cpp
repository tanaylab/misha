/*
 * TrackIndex2D.cpp
 *
 * Implementation of 2D track index loader, writer, and cache
 */

#include <errno.h>
#include <cstring>
#include <sys/stat.h>
#include <unistd.h>

#include "TrackIndex2D.h"
#include "CRC64.h"
#include "TGLException.h"

// Magic header: "MISHT2D\0" (8 bytes)
const char TrackIndex2D::MAGIC_HEADER[8] = {'M','I','S','H','T','2','D','\0'};

// Static cache members
std::map<std::string, std::shared_ptr<TrackIndex2D>> TrackIndex2D::s_index_cache;
std::mutex TrackIndex2D::s_cache_mutex;

TrackIndex2D::TrackIndex2D() : m_loaded(false), m_track_type(MishaTrack2DType::RECTS) {}

TrackIndex2D::~TrackIndex2D() {}

bool TrackIndex2D::is_little_endian() {
    uint32_t test = 1;
    return *reinterpret_cast<uint8_t*>(&test) == 1;
}

bool TrackIndex2D::load(const string &index_path) {
    // Check if file exists
    struct stat st;
    if (stat(index_path.c_str(), &st) != 0) {
        // File doesn't exist - this is not an error, just return false
        if (errno == ENOENT) {
            return false;
        }
        // Other stat errors are real errors
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to stat index file %s: %s", index_path.c_str(), strerror(errno));
    }

    FILE *fp = fopen(index_path.c_str(), "rb");
    if (!fp) {
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to open index file %s: %s", index_path.c_str(), strerror(errno));
    }

    // Read entire fixed header in one call (44 bytes)
#pragma pack(push, 1)
    struct Index2DHeader {
        char magic[8];
        uint32_t version;
        uint32_t track_type_raw;
        uint32_t num_pairs;
        uint64_t flags;
        uint64_t stored_checksum;
        uint64_t reserved;
    };
#pragma pack(pop)

    Index2DHeader hdr;
    if (fread(&hdr, sizeof(hdr), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to read index header from %s", index_path.c_str());
    }

    // Validate magic header
    if (memcmp(hdr.magic, MAGIC_HEADER, 8) != 0) {
        fclose(fp);
        TGLError<TrackIndex2D>(INVALID_FORMAT,
            "Invalid index file header in %s (expected MISHT2D magic)", index_path.c_str());
    }

    // Validate version
    if (hdr.version != INDEX_VERSION) {
        fclose(fp);
        TGLError<TrackIndex2D>(VERSION_MISMATCH,
            "Index version %u not supported (expected %u) in %s",
            hdr.version, INDEX_VERSION, index_path.c_str());
    }

    // Validate track type
    if (hdr.track_type_raw > 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(INVALID_FORMAT,
            "Invalid 2D track type %u in %s (expected 0=RECTS or 1=POINTS)",
            hdr.track_type_raw, index_path.c_str());
    }
    m_track_type = static_cast<MishaTrack2DType>(hdr.track_type_raw);

    uint32_t num_pairs = hdr.num_pairs;
    uint64_t stored_checksum = hdr.stored_checksum;

    // Sanity check: for 2D tracks the max is N*(N+1)/2 pairs.
    // For a genome with 50K contigs that is ~1.25 billion, but indexed
    // format is sparse -- only non-empty pairs are stored.
    // 20 million is a generous upper bound.
    if (num_pairs > 20000000) {
        fclose(fp);
        TGLError<TrackIndex2D>(INVALID_FORMAT,
            "Number of pairs %u exceeds maximum (20000000) in %s",
            num_pairs, index_path.c_str());
    }

    // Check endianness
    bool index_is_little_endian = (hdr.flags & FLAG_LITTLE_ENDIAN) != 0;
    if (index_is_little_endian != is_little_endian()) {
        fclose(fp);
        TGLError<TrackIndex2D>(ENDIAN_MISMATCH,
            "Index file %s has incompatible endianness", index_path.c_str());
    }

    // Read pair entries
    m_entries.clear();
    m_entries.reserve(num_pairs);
    m_pair_to_index.clear();

    // Packed struct matching the on-disk per-entry layout (28 bytes)
#pragma pack(push, 1)
    struct DiskPairEntry {
        uint32_t chrom1_id;
        uint32_t chrom2_id;
        uint64_t offset;
        uint64_t length;
        uint32_t reserved;
    };
#pragma pack(pop)

    for (uint32_t i = 0; i < num_pairs; ++i) {
        // Read all entry fields in one call (28 bytes)
        DiskPairEntry disk_entry;
        if (fread(&disk_entry, sizeof(disk_entry), 1, fp) != 1) {
            fclose(fp);
            TGLError<TrackIndex2D>(FILE_READ_FAILED,
                "Failed to read entry %u in %s", i, index_path.c_str());
        }

        Track2DPairEntry entry;
        entry.chrom1_id = disk_entry.chrom1_id;
        entry.chrom2_id = disk_entry.chrom2_id;
        entry.offset = disk_entry.offset;
        entry.length = disk_entry.length;
        entry.reserved = disk_entry.reserved;

        // Validate offset+length for overflow
        if (entry.offset + entry.length < entry.offset) {
            fclose(fp);
            TGLError<TrackIndex2D>(INVALID_FORMAT,
                "Offset+length overflow for pair (%u,%u) in %s",
                entry.chrom1_id, entry.chrom2_id, index_path.c_str());
        }

        // Check for duplicate pair
        uint64_t pair_key = make_pair_key(entry.chrom1_id, entry.chrom2_id);
        auto insert_result = m_pair_to_index.insert({pair_key, i});
        if (!insert_result.second) {
            fclose(fp);
            TGLError<TrackIndex2D>(INVALID_FORMAT,
                "Duplicate chromosome pair (%u,%u) in 2D track index %s",
                entry.chrom1_id, entry.chrom2_id, index_path.c_str());
        }

        m_entries.push_back(entry);
    }

    fclose(fp);

    // Validate checksum
    uint64_t computed_checksum = compute_checksum(m_entries);
    if (computed_checksum != stored_checksum) {
        TGLError<TrackIndex2D>(CHECKSUM_FAILED,
            "Index file checksum mismatch in %s (expected %016llX, got %016llX). "
            "Index may be corrupt.",
            index_path.c_str(),
            (unsigned long long)stored_checksum,
            (unsigned long long)computed_checksum);
    }

    m_loaded = true;
    return true;
}

const Track2DPairEntry* TrackIndex2D::get_entry(uint32_t chrom1_id, uint32_t chrom2_id) const {
    uint64_t key = make_pair_key(chrom1_id, chrom2_id);
    auto it = m_pair_to_index.find(key);
    if (it == m_pair_to_index.end()) {
        return nullptr;
    }
    return &m_entries[it->second];
}

bool TrackIndex2D::has_entry(uint32_t chrom1_id, uint32_t chrom2_id) const {
    uint64_t key = make_pair_key(chrom1_id, chrom2_id);
    return m_pair_to_index.find(key) != m_pair_to_index.end();
}

uint64_t TrackIndex2D::compute_checksum(const vector<Track2DPairEntry> &entries) {
    // Use CRC64-ECMA for checksum (same as TrackIndex and IntervalsIndex2D)
    misha::CRC64 crc64;
    uint64_t checksum = crc64.init_incremental();

    for (const auto &entry : entries) {
        // Hash all fields in order (excluding reserved for future compatibility)
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.chrom1_id, sizeof(entry.chrom1_id));
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.chrom2_id, sizeof(entry.chrom2_id));
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.offset, sizeof(entry.offset));
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.length, sizeof(entry.length));
        // Note: reserved field is intentionally NOT included in checksum
        // to allow future use without breaking compatibility
    }

    return crc64.finalize_incremental(checksum);
}

// --- Static cache management ---

std::shared_ptr<TrackIndex2D> TrackIndex2D::get_track_index_2d(const std::string &track_dir) {
    // Fast path: check cache under lock
    {
        std::lock_guard<std::mutex> lock(s_cache_mutex);
        auto it = s_index_cache.find(track_dir);
        if (it != s_index_cache.end()) {
            return it->second;
        }
    }

    // Load outside lock to avoid blocking other threads during I/O
    auto idx = std::make_shared<TrackIndex2D>();
    std::string idx_path = track_dir + "/track.idx";

    if (!idx->load(idx_path)) {
        // Index file doesn't exist - return nullptr
        return nullptr;
    }

    // Re-acquire lock and insert; if another thread inserted first, use its entry
    {
        std::lock_guard<std::mutex> lock(s_cache_mutex);
        auto [it, inserted] = s_index_cache.emplace(track_dir, idx);
        return it->second;  // return existing if another thread inserted first
    }
}

void TrackIndex2D::clear_cache() {
    std::lock_guard<std::mutex> lock(s_cache_mutex);
    s_index_cache.clear();
}

// --- Static write method ---

void TrackIndex2D::write_index(const string &index_path,
                               MishaTrack2DType track_type,
                               const vector<Track2DPairEntry> &entries) {
    FILE *fp = fopen(index_path.c_str(), "wb");
    if (!fp) {
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to open index file %s for writing: %s",
            index_path.c_str(), strerror(errno));
    }

    // --- Write header (44 bytes) in one call ---
#pragma pack(push, 1)
    struct WriteIndex2DHeader {
        char magic[8];
        uint32_t version;
        uint32_t track_type_raw;
        uint32_t num_pairs;
        uint64_t flags;
        uint64_t checksum;
        uint64_t reserved;
    };
#pragma pack(pop)

    WriteIndex2DHeader hdr;
    memcpy(hdr.magic, MAGIC_HEADER, 8);
    hdr.version = INDEX_VERSION;
    hdr.track_type_raw = static_cast<uint32_t>(track_type);
    hdr.num_pairs = static_cast<uint32_t>(entries.size());
    hdr.flags = is_little_endian() ? FLAG_LITTLE_ENDIAN : 0;
    hdr.checksum = compute_checksum(entries);
    hdr.reserved = 0;

    if (fwrite(&hdr, sizeof(hdr), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to write index header to %s", index_path.c_str());
    }

    // --- Write entry table ---
#pragma pack(push, 1)
    struct WriteDiskPairEntry {
        uint32_t chrom1_id;
        uint32_t chrom2_id;
        uint64_t offset;
        uint64_t length;
        uint32_t reserved;
    };
#pragma pack(pop)

    for (const auto &entry : entries) {
        WriteDiskPairEntry disk_entry;
        disk_entry.chrom1_id = entry.chrom1_id;
        disk_entry.chrom2_id = entry.chrom2_id;
        disk_entry.offset = entry.offset;
        disk_entry.length = entry.length;
        disk_entry.reserved = entry.reserved;

        if (fwrite(&disk_entry, sizeof(disk_entry), 1, fp) != 1) {
            fclose(fp);
            TGLError<TrackIndex2D>(FILE_READ_FAILED,
                "Failed to write entry for pair (%u,%u) to %s",
                entry.chrom1_id, entry.chrom2_id, index_path.c_str());
        }
    }

    // Flush and sync
    if (fflush(fp) != 0) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to flush index file %s: %s", index_path.c_str(), strerror(errno));
    }

    int fd = fileno(fp);
    if (fd >= 0) {
        fsync(fd);
    }

    fclose(fp);
}
