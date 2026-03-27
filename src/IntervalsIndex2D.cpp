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

    // Read entire fixed header in one call (40 bytes)
#pragma pack(push, 1)
    struct Index2DHeader {
        char magic[8];
        uint32_t version;
        uint32_t num_entries;
        uint64_t flags;
        uint64_t stored_checksum;
        uint64_t reserved;
    };
#pragma pack(pop)

    Index2DHeader hdr;
    if (fread(&hdr, sizeof(hdr), 1, fp) != 1) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(FILE_READ_FAILED,
            "Failed to read index header from %s", index_path.c_str());
    }

    // Validate magic header
    if (memcmp(hdr.magic, MAGIC_HEADER, 8) != 0) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(INVALID_FORMAT,
            "Invalid index file header in %s (expected MISHAI2D magic)", index_path.c_str());
    }

    // Validate version
    if (hdr.version != INDEX_VERSION) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(VERSION_MISMATCH,
            "Index version %u not supported (expected %u) in %s",
            hdr.version, INDEX_VERSION, index_path.c_str());
    }

    uint32_t num_entries = hdr.num_entries;
    uint64_t stored_checksum = hdr.stored_checksum;

    // Check endianness
    bool index_is_little_endian = (hdr.flags & FLAG_LITTLE_ENDIAN) != 0;
    if (index_is_little_endian != is_little_endian()) {
        fclose(fp);
        TGLError<IntervalsIndex2D>(ENDIAN_MISMATCH,
            "Index file %s has incompatible endianness", index_path.c_str());
    }

    // Read pair entries
    m_entries.clear();
    m_entries.reserve(num_entries);
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

    for (uint32_t i = 0; i < num_entries; ++i) {
        // Read all entry fields in one call (28 bytes)
        DiskPairEntry disk_entry;
        if (fread(&disk_entry, sizeof(disk_entry), 1, fp) != 1) {
            fclose(fp);
            TGLError<IntervalsIndex2D>(FILE_READ_FAILED,
                "Failed to read entry %u in %s", i, index_path.c_str());
        }

        IntervalsPairEntry entry;
        entry.chrom1_id = disk_entry.chrom1_id;
        entry.chrom2_id = disk_entry.chrom2_id;
        entry.offset = disk_entry.offset;
        entry.length = disk_entry.length;
        entry.reserved = disk_entry.reserved;

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
