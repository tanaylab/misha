/*
 * TrackIndex.cpp
 *
 * Implementation of track index loader
 */

#include <errno.h>
#include <cstring>
#include <sys/stat.h>

#include "TrackIndex.h"
#include "CRC64.h"
#include "TGLException.h"

const char TrackIndex::MAGIC_HEADER[8] = {'M','I','S','H','A','T','D','X'};

TrackIndex::TrackIndex() : m_loaded(false), m_track_type(MishaTrackType::DENSE) {}

TrackIndex::~TrackIndex() {}

bool TrackIndex::is_little_endian() {
    uint32_t test = 1;
    return *reinterpret_cast<uint8_t*>(&test) == 1;
}

bool TrackIndex::load(const string &index_path) {
    // Check if file exists
    struct stat st;
    if (stat(index_path.c_str(), &st) != 0) {
        // File doesn't exist - this is not an error, just return false
        if (errno == ENOENT) {
            return false;
        }
        // Other stat errors are real errors
        TGLError<TrackIndex>(FILE_READ_FAILED,
            "Failed to stat index file %s: %s", index_path.c_str(), strerror(errno));
    }

    FILE *fp = fopen(index_path.c_str(), "rb");
    if (!fp) {
        TGLError<TrackIndex>(FILE_READ_FAILED,
            "Failed to open index file %s: %s", index_path.c_str(), strerror(errno));
    }

    // Read and validate magic header
    char header[8];
    if (fread(header, 1, 8, fp) != 8 || memcmp(header, MAGIC_HEADER, 8) != 0) {
        fclose(fp);
        TGLError<TrackIndex>(INVALID_FORMAT,
            "Invalid index file header in %s (expected MISHATDX magic)", index_path.c_str());
    }

    // Read version
    uint32_t version;
    if (fread(&version, sizeof(version), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex>(FILE_READ_FAILED,
            "Failed to read index version from %s", index_path.c_str());
    }
    if (version != INDEX_VERSION) {
        fclose(fp);
        TGLError<TrackIndex>(VERSION_MISMATCH,
            "Index version %u not supported (expected %u) in %s",
            version, INDEX_VERSION, index_path.c_str());
    }

    // Read track type
    uint32_t track_type_raw;
    if (fread(&track_type_raw, sizeof(track_type_raw), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex>(FILE_READ_FAILED,
            "Failed to read track type from %s", index_path.c_str());
    }
    if (track_type_raw > 2) {
        fclose(fp);
        TGLError<TrackIndex>(INVALID_FORMAT,
            "Invalid track type %u in %s (expected 0-2)", track_type_raw, index_path.c_str());
    }
    m_track_type = static_cast<MishaTrackType>(track_type_raw);

    // Read number of contigs
    uint32_t num_contigs;
    if (fread(&num_contigs, sizeof(num_contigs), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex>(FILE_READ_FAILED,
            "Failed to read contig count from %s", index_path.c_str());
    }

    // Sanity check: 20 million contigs should be more than enough for any genome
    if (num_contigs > 20000000) {
        fclose(fp);
        TGLError<TrackIndex>(INVALID_FORMAT,
            "Number of contigs %u exceeds maximum (20000000) in %s", num_contigs, index_path.c_str());
    }

    // Read flags
    uint64_t flags;
    if (fread(&flags, sizeof(flags), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex>(FILE_READ_FAILED,
            "Failed to read flags from %s", index_path.c_str());
    }

    // Check endianness
    bool index_is_little_endian = (flags & FLAG_LITTLE_ENDIAN) != 0;
    if (index_is_little_endian != is_little_endian()) {
        fclose(fp);
        TGLError<TrackIndex>(ENDIAN_MISMATCH,
            "Index file %s has incompatible endianness", index_path.c_str());
    }

    // Read stored checksum (will validate later)
    uint64_t stored_checksum;
    if (fread(&stored_checksum, sizeof(stored_checksum), 1, fp) != 1) {
        fclose(fp);
        TGLError<TrackIndex>(FILE_READ_FAILED,
            "Failed to read checksum from %s", index_path.c_str());
    }

    // Read contig entries
    m_entries.clear();
    m_entries.reserve(num_contigs);
    m_chromid_to_index.clear();

    for (uint32_t i = 0; i < num_contigs; ++i) {
        TrackContigEntry entry;

        // Read chrom_id, offset, length, reserved (24 bytes total)
        if (fread(&entry.chrom_id, sizeof(entry.chrom_id), 1, fp) != 1 ||
            fread(&entry.offset, sizeof(entry.offset), 1, fp) != 1 ||
            fread(&entry.length, sizeof(entry.length), 1, fp) != 1 ||
            fread(&entry.reserved, sizeof(entry.reserved), 1, fp) != 1) {
            fclose(fp);
            TGLError<TrackIndex>(FILE_READ_FAILED,
                "Failed to read entry %u in %s", i, index_path.c_str());
        }

        // Validate offset+length for overflow
        if (entry.offset + entry.length < entry.offset) {
            fclose(fp);
            TGLError<TrackIndex>(INVALID_FORMAT,
                "Offset+length overflow for chromid %u in %s",
                entry.chrom_id, index_path.c_str());
        }

        // Check for duplicate chromid
        auto insert_result = m_chromid_to_index.insert({entry.chrom_id, i});
        if (!insert_result.second) {
            fclose(fp);
            TGLError<TrackIndex>(INVALID_FORMAT,
                "Duplicate chromosome ID %u in track index %s",
                entry.chrom_id, index_path.c_str());
        }

        m_entries.push_back(entry);
    }

    fclose(fp);

    // Validate checksum
    uint64_t computed_checksum = compute_checksum(m_entries);
    if (computed_checksum != stored_checksum) {
        TGLError<TrackIndex>(CHECKSUM_FAILED,
            "Index file checksum mismatch in %s (expected %016llX, got %016llX). "
            "Index may be corrupt.",
            index_path.c_str(),
            (unsigned long long)stored_checksum,
            (unsigned long long)computed_checksum);
    }

    m_loaded = true;
    return true;
}

const TrackContigEntry* TrackIndex::get_entry(uint32_t chromid) const {
    auto it = m_chromid_to_index.find(chromid);
    if (it == m_chromid_to_index.end()) {
        return nullptr;
    }
    return &m_entries[it->second];
}

uint64_t TrackIndex::compute_checksum(const vector<TrackContigEntry> &entries) {
    // Use CRC64-ECMA for checksum (same as GenomeIndex)
    misha::CRC64 crc64;
    uint64_t checksum = crc64.init_incremental();

    for (const auto &entry : entries) {
        // Hash all fields in order
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.chrom_id, sizeof(entry.chrom_id));
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.offset, sizeof(entry.offset));
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.length, sizeof(entry.length));
        // Note: reserved field is intentionally NOT included in checksum
        // to allow future use without breaking compatibility
    }

    return crc64.finalize_incremental(checksum);
}
