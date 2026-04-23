/*
 * GenomeSeqRenameIndex.cpp
 *
 * C++ entry point used by gdb.rename_chroms() to rewrite seq/genome.idx
 * with renamed chromosomes. Offsets/lengths are preserved; only names
 * change. The CRC64 checksum is recomputed from scratch using the same
 * formula as GenomeIndex::compute_checksum (chromid, name-if-nonempty,
 * offset, length -- in that order).
 */

#include <errno.h>
#include <unistd.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#include "CRC64.h"
#include "GenomeIndex.h"
#include "TGLException.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

extern "C" {

SEXP C_gdb_rewrite_genome_idx(SEXP _idx_path, SEXP _old_names, SEXP _new_names, SEXP _envir)
{
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_idx_path) || Rf_length(_idx_path) != 1)
            verror("idx_path must be a character(1)");
        if (!Rf_isString(_old_names) || !Rf_isString(_new_names))
            verror("old_names and new_names must be character vectors");
        if (Rf_length(_old_names) != Rf_length(_new_names))
            verror("old_names and new_names must have the same length");

        string idx_path = CHAR(STRING_ELT(_idx_path, 0));

        unordered_map<string, string> rename_map;
        int n_rename = Rf_length(_old_names);
        for (int i = 0; i < n_rename; ++i) {
            rename_map[CHAR(STRING_ELT(_old_names, i))] = CHAR(STRING_ELT(_new_names, i));
        }

        GenomeIndex index;
        index.load(idx_path);

        const vector<ContigIndexEntry> &old_entries = index.get_all_entries();
        vector<ContigIndexEntry> new_entries;
        new_entries.reserve(old_entries.size());

        for (const auto &e : old_entries) {
            ContigIndexEntry ne = e;
            auto it = rename_map.find(e.name);
            if (it != rename_map.end()) {
                ne.name = it->second;
            }
            new_entries.push_back(ne);
        }

        // Recompute CRC64 over renamed entries. Must match
        // GenomeIndex::compute_checksum exactly: chromid, name (if non-empty),
        // offset, length -- all hashed incrementally in that order.
        misha::CRC64 crc64;
        uint64_t checksum = crc64.init_incremental();
        for (const auto &e : new_entries) {
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&e.chromid, sizeof(e.chromid));
            if (!e.name.empty()) {
                checksum = crc64.compute_incremental(checksum,
                    (const unsigned char*)e.name.c_str(), e.name.size());
            }
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&e.offset, sizeof(e.offset));
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&e.length, sizeof(e.length));
        }
        uint64_t stored_checksum = crc64.finalize_incremental(checksum);

        string tmp_path = idx_path + ".tmp." + to_string((long long)getpid());
        FILE *fp = fopen(tmp_path.c_str(), "wb");
        if (!fp) {
            verror("failed to open %s: %s", tmp_path.c_str(), strerror(errno));
        }

        #pragma pack(push, 1)
        struct IndexHeader {
            char     magic[8];
            uint32_t version;
            uint32_t num_contigs;
            uint64_t stored_checksum;
        };
        #pragma pack(pop)

        IndexHeader hdr;
        memcpy(hdr.magic, "MISHAIDX", 8);
        hdr.version = 1;
        hdr.num_contigs = (uint32_t)new_entries.size();
        hdr.stored_checksum = stored_checksum;
        if (fwrite(&hdr, sizeof(hdr), 1, fp) != 1) {
            fclose(fp); unlink(tmp_path.c_str());
            verror("failed to write header to %s", tmp_path.c_str());
        }

        #pragma pack(push, 1)
        struct EntryTail {
            uint64_t offset;
            uint64_t length;
            uint64_t reserved;
        };
        #pragma pack(pop)

        for (const auto &e : new_entries) {
            if (fwrite(&e.chromid, sizeof(e.chromid), 1, fp) != 1) goto write_err;
            {
                uint16_t name_len = (uint16_t)e.name.size();
                if (fwrite(&name_len, sizeof(name_len), 1, fp) != 1) goto write_err;
                if (name_len > 0) {
                    if (fwrite(e.name.data(), 1, name_len, fp) != name_len) goto write_err;
                }
            }
            EntryTail tail;
            tail.offset   = e.offset;
            tail.length   = e.length;
            tail.reserved = e.reserved;
            if (fwrite(&tail, sizeof(tail), 1, fp) != 1) goto write_err;
        }

        fflush(fp);
        fsync(fileno(fp));
        fclose(fp);

        if (rename(tmp_path.c_str(), idx_path.c_str()) != 0) {
            unlink(tmp_path.c_str());
            verror("failed to rename %s to %s: %s",
                   tmp_path.c_str(), idx_path.c_str(), strerror(errno));
        }

        return R_NilValue;

    write_err:
        fclose(fp);
        unlink(tmp_path.c_str());
        verror("failed to write entry to %s", tmp_path.c_str());
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

} // extern "C"
