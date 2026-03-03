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

    // Read and validate magic header
    char header[8];
    if (fread(header, 1, 8, fp) != 8 || memcmp(header, MAGIC_HEADER, 8) != 0) {
        fclose(fp);
        TGLError<TrackIndex2D>(INVALID_FORMAT,
            "Invalid index file header in %s (expected MISHT2D magic)", index_path.c_str());
    }

    // Read version
    uint32_t version;
    if (fread(&version, sizeof(version), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to read index version from %s", index_path.c_str());
    }
    if (version != INDEX_VERSION) {
        fclose(fp);
        TGLError<TrackIndex2D>(VERSION_MISMATCH,
            "Index version %u not supported (expected %u) in %s",
            version, INDEX_VERSION, index_path.c_str());
    }

    // Read track type
    uint32_t track_type_raw;
    if (fread(&track_type_raw, sizeof(track_type_raw), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to read track type from %s", index_path.c_str());
    }
    if (track_type_raw > 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(INVALID_FORMAT,
            "Invalid 2D track type %u in %s (expected 0=RECTS or 1=POINTS)",
            track_type_raw, index_path.c_str());
    }
    m_track_type = static_cast<MishaTrack2DType>(track_type_raw);

    // Read number of pairs
    uint32_t num_pairs;
    if (fread(&num_pairs, sizeof(num_pairs), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to read pair count from %s", index_path.c_str());
    }

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

    // Read flags
    uint64_t flags;
    if (fread(&flags, sizeof(flags), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to read flags from %s", index_path.c_str());
    }

    // Check endianness
    bool index_is_little_endian = (flags & FLAG_LITTLE_ENDIAN) != 0;
    if (index_is_little_endian != is_little_endian()) {
        fclose(fp);
        TGLError<TrackIndex2D>(ENDIAN_MISMATCH,
            "Index file %s has incompatible endianness", index_path.c_str());
    }

    // Read stored checksum (will validate later)
    uint64_t stored_checksum;
    if (fread(&stored_checksum, sizeof(stored_checksum), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to read checksum from %s", index_path.c_str());
    }

    // Read reserved field
    uint64_t reserved;
    if (fread(&reserved, sizeof(reserved), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to read reserved field from %s", index_path.c_str());
    }

    // Read pair entries
    m_entries.clear();
    m_entries.reserve(num_pairs);
    m_pair_to_index.clear();

    for (uint32_t i = 0; i < num_pairs; ++i) {
        Track2DPairEntry entry;

        // Read chrom1_id, chrom2_id, offset, length, reserved (28 bytes total)
        if (fread(&entry.chrom1_id, sizeof(entry.chrom1_id), 1, fp) != 1 ||
            fread(&entry.chrom2_id, sizeof(entry.chrom2_id), 1, fp) != 1 ||
            fread(&entry.offset, sizeof(entry.offset), 1, fp) != 1 ||
            fread(&entry.length, sizeof(entry.length), 1, fp) != 1 ||
            fread(&entry.reserved, sizeof(entry.reserved), 1, fp) != 1) {
            fclose(fp);
            TGLError<TrackIndex2D>(FILE_READ_FAILED,
                "Failed to read entry %u in %s", i, index_path.c_str());
        }

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
    std::lock_guard<std::mutex> lock(s_cache_mutex);

    // Check if already in cache
    auto it = s_index_cache.find(track_dir);
    if (it != s_index_cache.end()) {
        return it->second;
    }

    // Create new index and try to load
    auto idx = std::make_shared<TrackIndex2D>();
    std::string idx_path = track_dir + "/track.idx";

    if (!idx->load(idx_path)) {
        // Index file doesn't exist - return nullptr
        return nullptr;
    }

    // Cache and return
    s_index_cache[track_dir] = idx;
    return idx;
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

    // --- Write header (44 bytes) ---

    // Magic header
    if (fwrite(MAGIC_HEADER, 1, 8, fp) != 8) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to write magic header to %s", index_path.c_str());
    }

    // Version
    uint32_t version = INDEX_VERSION;
    if (fwrite(&version, sizeof(version), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to write version to %s", index_path.c_str());
    }

    // Track type
    uint32_t track_type_raw = static_cast<uint32_t>(track_type);
    if (fwrite(&track_type_raw, sizeof(track_type_raw), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to write track type to %s", index_path.c_str());
    }

    // Number of pairs
    uint32_t num_pairs = static_cast<uint32_t>(entries.size());
    if (fwrite(&num_pairs, sizeof(num_pairs), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to write pair count to %s", index_path.c_str());
    }

    // Flags (little-endian flag)
    uint64_t flags = is_little_endian() ? FLAG_LITTLE_ENDIAN : 0;
    if (fwrite(&flags, sizeof(flags), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to write flags to %s", index_path.c_str());
    }

    // Compute checksum
    uint64_t checksum = compute_checksum(entries);
    if (fwrite(&checksum, sizeof(checksum), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to write checksum to %s", index_path.c_str());
    }

    // Reserved
    uint64_t reserved = 0;
    if (fwrite(&reserved, sizeof(reserved), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex2D>(FILE_READ_FAILED,
            "Failed to write reserved field to %s", index_path.c_str());
    }

    // --- Write entry table ---

    for (const auto &entry : entries) {
        if (fwrite(&entry.chrom1_id, sizeof(entry.chrom1_id), 1, fp) != 1 ||
            fwrite(&entry.chrom2_id, sizeof(entry.chrom2_id), 1, fp) != 1 ||
            fwrite(&entry.offset, sizeof(entry.offset), 1, fp) != 1 ||
            fwrite(&entry.length, sizeof(entry.length), 1, fp) != 1 ||
            fwrite(&entry.reserved, sizeof(entry.reserved), 1, fp) != 1) {
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
