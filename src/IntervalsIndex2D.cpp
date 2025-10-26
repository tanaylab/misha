/*
 * IntervalsIndex2D.cpp
 *
 * Implementation of 2D intervals index loader
 */

#include <errno.h>
#include <cstring>
#include <sys/stat.h>

#include "IntervalsIndex2D.h"
#include "CRC64.h"
#include "TGLException.h"

const char IntervalsIndex2D::MAGIC_HEADER[8] = {'M','I','S','H','A','I','2','D'};

IntervalsIndex2D::IntervalsIndex2D() : m_loaded(false) {}

IntervalsIndex2D::~IntervalsIndex2D() {}

bool IntervalsIndex2D::is_little_endian() {
    uint32_t test = 1;
    return *reinterpret_cast<uint8_t*>(&test) == 1;
}

bool IntervalsIndex2D::load(const string &index_path) {
    // Check if file exists
    struct stat st;
    if (stat(index_path.c_str(), &st) != 0) {
        // File doesn't exist - this is not an error, just return false
        if (errno == ENOENT) {
            return false;
        }
        // Other stat errors are real errors
        TGLError<IntervalsIndex2D>(FILE_READ_FAILED,
            "Failed to stat index file %s: %s", index_path.c_str(), strerror(errno));
    }

    FILE *fp = fopen(index_path.c_str(), "rb");
    if (!fp) {
        TGLError<IntervalsIndex2D>(FILE_READ_FAILED,
            "Failed to open index file %s: %s", index_path.c_str(), strerror(errno));
    }

    // Read and validate magic header
    char header[8];
    if (fread(header, 1, 8, fp) != 8 || memcmp(header, MAGIC_HEADER, 8) != 0) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(INVALID_FORMAT,
            "Invalid index file header in %s (expected MISHAI2D magic)", index_path.c_str());
    }

    // Read version
    uint32_t version;
    if (fread(&version, sizeof(version), 1, fp) != 1) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(FILE_READ_FAILED,
            "Failed to read index version from %s", index_path.c_str());
    }
    if (version != INDEX_VERSION) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(VERSION_MISMATCH,
            "Index version %u not supported (expected %u) in %s",
            version, INDEX_VERSION, index_path.c_str());
    }

    // Read number of entries
    uint32_t num_entries;
    if (fread(&num_entries, sizeof(num_entries), 1, fp) != 1) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(FILE_READ_FAILED,
            "Failed to read entry count from %s", index_path.c_str());
    }

    // Read flags
    uint64_t flags;
    if (fread(&flags, sizeof(flags), 1, fp) != 1) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(FILE_READ_FAILED,
            "Failed to read flags from %s", index_path.c_str());
    }

    // Check endianness
    bool index_is_little_endian = (flags & FLAG_LITTLE_ENDIAN) != 0;
    if (index_is_little_endian != is_little_endian()) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(ENDIAN_MISMATCH,
            "Index file %s has incompatible endianness", index_path.c_str());
    }

    // Read stored checksum (will validate later)
    uint64_t stored_checksum;
    if (fread(&stored_checksum, sizeof(stored_checksum), 1, fp) != 1) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(FILE_READ_FAILED,
            "Failed to read checksum from %s", index_path.c_str());
    }

    // Read reserved field
    uint64_t reserved;
    if (fread(&reserved, sizeof(reserved), 1, fp) != 1) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(FILE_READ_FAILED,
            "Failed to read reserved field from %s", index_path.c_str());
    }

    // Read pair entries
    m_entries.clear();
    m_entries.reserve(num_entries);
    m_pair_to_index.clear();

    for (uint32_t i = 0; i < num_entries; ++i) {
        IntervalsPairEntry entry;

        // Read chrom1_id, chrom2_id, offset, length, reserved (28 bytes total)
        if (fread(&entry.chrom1_id, sizeof(entry.chrom1_id), 1, fp) != 1 ||
            fread(&entry.chrom2_id, sizeof(entry.chrom2_id), 1, fp) != 1 ||
            fread(&entry.offset, sizeof(entry.offset), 1, fp) != 1 ||
            fread(&entry.length, sizeof(entry.length), 1, fp) != 1 ||
            fread(&entry.reserved, sizeof(entry.reserved), 1, fp) != 1) {
            fclose(fp);
            TGLError<IntervalsIndex2D>(FILE_READ_FAILED,
                "Failed to read entry %u in %s", i, index_path.c_str());
        }

        // Validate offset+length for overflow
        if (entry.offset + entry.length < entry.offset) {
            fclose(fp);
            TGLError<IntervalsIndex2D>(INVALID_FORMAT,
                "Offset+length overflow for pair (%u,%u) in %s",
                entry.chrom1_id, entry.chrom2_id, index_path.c_str());
        }

        // Check for duplicate pair
        uint64_t pair_key = make_pair_key(entry.chrom1_id, entry.chrom2_id);
        auto insert_result = m_pair_to_index.insert({pair_key, i});
        if (!insert_result.second) {
            fclose(fp);
            TGLError<IntervalsIndex2D>(INVALID_FORMAT,
                "Duplicate chromosome pair (%u,%u) in intervals2d index %s",
                entry.chrom1_id, entry.chrom2_id, index_path.c_str());
        }

        m_entries.push_back(entry);
    }

    fclose(fp);

    // Validate checksum
    uint64_t computed_checksum = compute_checksum(m_entries);
    if (computed_checksum != stored_checksum) {
        TGLError<IntervalsIndex2D>(CHECKSUM_FAILED,
            "Index file checksum mismatch in %s (expected %016llX, got %016llX). "
            "Index may be corrupt.",
            index_path.c_str(),
            (unsigned long long)stored_checksum,
            (unsigned long long)computed_checksum);
    }

    m_loaded = true;
    return true;
}

const IntervalsPairEntry* IntervalsIndex2D::get_entry(uint32_t chrom1_id, uint32_t chrom2_id) const {
    uint64_t key = make_pair_key(chrom1_id, chrom2_id);
    auto it = m_pair_to_index.find(key);
    if (it == m_pair_to_index.end()) {
        return nullptr;
    }
    return &m_entries[it->second];
}

uint64_t IntervalsIndex2D::compute_checksum(const vector<IntervalsPairEntry> &entries) {
    // Use CRC64-ECMA for checksum (same as TrackIndex)
    misha::CRC64 crc64;
    uint64_t checksum = crc64.init_incremental();

    for (const auto &entry : entries) {
        // Hash all fields in order
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
