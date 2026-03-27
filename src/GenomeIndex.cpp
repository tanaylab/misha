/*
 * GenomeIndex.cpp
 *
 * Implementation of genome index loader
 */

#include <errno.h>
#include <cstring>
#include <sys/stat.h>

#include "GenomeIndex.h"
#include "CRC64.h"
#include "TGLException.h"

const char GenomeIndex::MAGIC_HEADER[8] = {'M','I','S','H','A','I','D','X'};

GenomeIndex::GenomeIndex() : m_loaded(false) {}

GenomeIndex::~GenomeIndex() {}

void GenomeIndex::load(const string &index_path) {
    FILE *fp = fopen(index_path.c_str(), "rb");
    if (!fp) {
        TGLError<GenomeIndex>(FILE_READ_FAILED,
            "Failed to open index file %s: %s", index_path.c_str(), strerror(errno));
    }

    // Read entire fixed header in one call (24 bytes)
#pragma pack(push, 1)
    struct IndexHeader {
        char magic[8];
        uint32_t version;
        uint32_t num_contigs;
        uint64_t stored_checksum;
    };
#pragma pack(pop)

    IndexHeader hdr;
    if (fread(&hdr, sizeof(hdr), 1, fp) != 1) {
        fclose(fp);
        TGLError<GenomeIndex>(FILE_READ_FAILED,
            "Failed to read index header from %s", index_path.c_str());
    }

    // Validate magic header
    if (memcmp(hdr.magic, MAGIC_HEADER, 8) != 0) {
        fclose(fp);
        TGLError<GenomeIndex>(INVALID_FORMAT,
            "Invalid index file header in %s", index_path.c_str());
    }

    // Validate version
    if (hdr.version != INDEX_VERSION) {
        fclose(fp);
        TGLError<GenomeIndex>(VERSION_MISMATCH,
            "Index version %u not supported (expected %u) in %s",
            hdr.version, INDEX_VERSION, index_path.c_str());
    }

    uint32_t num_contigs = hdr.num_contigs;
    uint64_t stored_checksum = hdr.stored_checksum;

    // Sanity check: 20 million contigs should be more than enough for any genome
    if (num_contigs > 20000000) {
        fclose(fp);
        TGLError<GenomeIndex>(INVALID_FORMAT,
            "Number of contigs %u exceeds maximum (20000000) in %s", num_contigs, index_path.c_str());
    }

    // Read contig entries
    m_entries.clear();
    m_entries.reserve(num_contigs);
    m_chromid_to_index.clear();

    for (uint32_t i = 0; i < num_contigs; ++i) {
        ContigIndexEntry entry;

        // Read chromid
        if (fread(&entry.chromid, sizeof(entry.chromid), 1, fp) != 1) {
            fclose(fp);
            TGLError<GenomeIndex>(FILE_READ_FAILED,
                "Failed to read chromid at entry %u in %s", i, index_path.c_str());
        }

        // Read name length
        uint16_t name_length;
        if (fread(&name_length, sizeof(name_length), 1, fp) != 1) {
            fclose(fp);
            TGLError<GenomeIndex>(FILE_READ_FAILED,
                "Failed to read name length at entry %u in %s", i, index_path.c_str());
        }

        // Validate name length (reasonable limit for contig names)
        if (name_length > 1024) {
            fclose(fp);
            TGLError<GenomeIndex>(INVALID_FORMAT,
                "Contig name length %u exceeds maximum (1024) at entry %u in %s",
                name_length, i, index_path.c_str());
        }

        // Read name
        if (name_length > 0) {
            vector<char> name_buf(name_length + 1);
            if (fread(name_buf.data(), 1, name_length, fp) != name_length) {
                fclose(fp);
                TGLError<GenomeIndex>(FILE_READ_FAILED,
                    "Failed to read contig name at entry %u in %s", i, index_path.c_str());
            }
            name_buf[name_length] = '\0';
            entry.name = string(name_buf.data());
        }

        // Read offset, length, and reserved in one call
#pragma pack(push, 1)
        struct EntryTail {
            uint64_t offset;
            uint64_t length;
            uint64_t reserved;
        };
#pragma pack(pop)
        EntryTail tail;
        if (fread(&tail, sizeof(tail), 1, fp) != 1) {
            fclose(fp);
            TGLError<GenomeIndex>(FILE_READ_FAILED,
                "Failed to read offset/length/reserved at entry %u in %s", i, index_path.c_str());
        }
        entry.offset = tail.offset;
        entry.length = tail.length;
        entry.reserved = tail.reserved;

        // Validate offset+length for overflow
        if (entry.offset + entry.length < entry.offset) {
            fclose(fp);
            TGLError<GenomeIndex>(INVALID_FORMAT,
                "Offset+length overflow for contig '%s' (chromid %u) in %s",
                entry.name.c_str(), entry.chromid, index_path.c_str());
        }

        // Check for duplicate chromid
        auto insert_result = m_chromid_to_index.insert({entry.chromid, i});
        if (!insert_result.second) {
            fclose(fp);
            TGLError<GenomeIndex>(INVALID_FORMAT,
                "Duplicate chromosome ID %u (contig '%s') in index %s",
                entry.chromid, entry.name.c_str(), index_path.c_str());
        }

        m_entries.push_back(entry);
    }

    fclose(fp);

    // Validate checksum
    uint64_t computed_checksum = compute_checksum(m_entries);
    if (computed_checksum != stored_checksum) {
        TGLError<GenomeIndex>(CHECKSUM_FAILED,
            "Index file checksum mismatch in %s (expected %016llX, got %016llX). "
            "Index may be corrupt.",
            index_path.c_str(),
            (unsigned long long)stored_checksum,
            (unsigned long long)computed_checksum);
    }

    m_loaded = true;
}

const ContigIndexEntry* GenomeIndex::get_entry(uint32_t chromid) const {
    auto it = m_chromid_to_index.find(chromid);
    if (it == m_chromid_to_index.end()) {
        return nullptr;
    }
    return &m_entries[it->second];
}

uint64_t GenomeIndex::compute_checksum(const vector<ContigIndexEntry> &entries) {
    // Use CRC64-ECMA for checksum
    misha::CRC64 crc64;
    uint64_t checksum = crc64.init_incremental();

    for (const auto &entry : entries) {
        // Hash chromid
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.chromid, sizeof(entry.chromid));

        // Hash name
        if (!entry.name.empty()) {
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)entry.name.c_str(), entry.name.size());
        }

        // Hash offset and length
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.offset, sizeof(entry.offset));
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.length, sizeof(entry.length));
    }

    return crc64.finalize_incremental(checksum);
}
