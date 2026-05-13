// src/GenomeTrackSplitIndexed.cpp
//
// Splits a 1D indexed-format track (track.dat + track.idx) back into per-chromosome
// files in the same directory, named by the supplied chrom names.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <vector>

// Note: do NOT #include <Rinternals.h> directly. It defines a length(x) macro
// that collides with TrackContigEntry::length in TrackIndex.h. R's headers are
// pulled in transitively via rdbutils.h.
#include "TrackIndex.h"
#include "TrackIndexWriter.h"
#include "TGLException.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

extern "C" {

SEXP gtrack_split_indexed_to_per_chrom(SEXP _track_dir, SEXP _chrom_names, SEXP _remove_indexed) {
    vector<string> tmp_files_to_cleanup;
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_track_dir) || Rf_length(_track_dir) != 1)
            verror("track_dir must be a single string");
        if (!Rf_isString(_chrom_names))
            verror("chrom_names must be a character vector");

        const string track_dir = CHAR(STRING_ELT(_track_dir, 0));
        const int n_chroms = Rf_length(_chrom_names);
        vector<string> chrom_names(n_chroms);
        for (int i = 0; i < n_chroms; ++i)
            chrom_names[i] = CHAR(STRING_ELT(_chrom_names, i));
        const bool remove_indexed = Rf_asLogical(_remove_indexed) == TRUE;

        const string idx_path = track_dir + "/track.idx";
        const string dat_path = track_dir + "/track.dat";

        TrackIndex idx;
        if (!idx.load(idx_path))
            verror("track.idx not found in %s", track_dir.c_str());

        FILE *dat_fp = fopen(dat_path.c_str(), "rb");
        if (!dat_fp)
            verror("Failed to open %s: %s", dat_path.c_str(), strerror(errno));

        const size_t BUF = 1 << 20; // 1 MiB
        vector<char> buffer(BUF);

        for (const TrackContigEntry &entry : idx.get_all_entries()) {
            if (entry.chrom_id >= (uint32_t)n_chroms) {
                fclose(dat_fp);
                verror("track.idx references chrom_id %u but only %d chrom names supplied "
                       "(internal mismatch or corrupt index)",
                       entry.chrom_id, n_chroms);
            }

            const string out_path     = track_dir + "/" + chrom_names[entry.chrom_id];
            const string out_path_tmp = out_path + ".tmp";
            tmp_files_to_cleanup.push_back(out_path_tmp);

            FILE *out_fp = fopen(out_path_tmp.c_str(), "wb");
            if (!out_fp) {
                fclose(dat_fp);
                verror("Failed to create %s: %s", out_path_tmp.c_str(), strerror(errno));
            }

            // Length=0 entries arise when the source had no per-chrom file at convert time.
            // We still touch an output file (atomic via tmp+rename above) so that downstream
            // per-chrom invariants hold, but we write zero bytes. In practice this is
            // unreachable for tracks created via gtrack.create_*, which always writes a
            // 4-byte format-signature header for every chromosome.
            if (entry.length > 0) {
                if (fseeko(dat_fp, (off_t)entry.offset, SEEK_SET) != 0) {
                    fclose(out_fp); fclose(dat_fp);
                    verror("Failed to seek to offset %llu in %s",
                           (unsigned long long)entry.offset, dat_path.c_str());
                }

                uint64_t remaining = entry.length;
                while (remaining > 0) {
                    size_t to_read = (size_t)min((uint64_t)BUF, remaining);
                    size_t got = fread(buffer.data(), 1, to_read, dat_fp);
                    if (got != to_read) {
                        fclose(out_fp); fclose(dat_fp);
                        verror("Short read from %s at offset %llu",
                               dat_path.c_str(), (unsigned long long)entry.offset);
                    }
                    if (fwrite(buffer.data(), 1, got, out_fp) != got) {
                        fclose(out_fp); fclose(dat_fp);
                        verror("Failed to write %s: %s", out_path_tmp.c_str(), strerror(errno));
                    }
                    remaining -= got;
                }
            }

            // Per-file fsync + atomic rename: each per-chrom file is canonical db state;
            // we accept N fsyncs (one per contig) for crash-safe individual files.
            fflush(out_fp);
            fsync(fileno(out_fp));
            fclose(out_fp);

            if (rename(out_path_tmp.c_str(), out_path.c_str()) != 0) {
                fclose(dat_fp);
                verror("Failed to rename %s to %s: %s",
                       out_path_tmp.c_str(), out_path.c_str(), strerror(errno));
            }
            tmp_files_to_cleanup.pop_back(); // succeeded
        }

        fclose(dat_fp);

        if (remove_indexed) {
            unlink(dat_path.c_str());
            unlink(idx_path.c_str());
        }

        return R_NilValue;
    } catch (TGLException &e) {
        for (const string &p : tmp_files_to_cleanup) unlink(p.c_str());
        verror("%s", e.msg());
    } catch (const bad_alloc &) {
        for (const string &p : tmp_files_to_cleanup) unlink(p.c_str());
        verror("Out of memory");
    }
    return R_NilValue;
}

// Pack per-chromosome files in track_dir into track.dat + track.idx.
// Mirrors gtrack_convert_to_indexed_format but takes explicit args (track_dir,
// chrom_names, track_type) so it works without GROOT/ALLGENOME context.
SEXP gtrack_pack_per_chrom_to_indexed(SEXP _track_dir, SEXP _chrom_names, SEXP _track_type) {
    string dat_path_tmp;
    string idx_path_tmp;
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_track_dir) || Rf_length(_track_dir) != 1)
            verror("track_dir must be a single string");
        if (!Rf_isString(_chrom_names))
            verror("chrom_names must be a character vector");
        if (!Rf_isString(_track_type) || Rf_length(_track_type) != 1)
            verror("track_type must be a single string ('dense', 'sparse', or 'array')");

        const string track_dir = CHAR(STRING_ELT(_track_dir, 0));
        const int n_chroms = Rf_length(_chrom_names);
        vector<string> chrom_names(n_chroms);
        for (int i = 0; i < n_chroms; ++i)
            chrom_names[i] = CHAR(STRING_ELT(_chrom_names, i));

        const string type_str = CHAR(STRING_ELT(_track_type, 0));
        MishaTrackType track_type = MishaTrackType::DENSE;
        if (type_str == "dense")       track_type = MishaTrackType::DENSE;
        else if (type_str == "sparse") track_type = MishaTrackType::SPARSE;
        else if (type_str == "array")  track_type = MishaTrackType::ARRAY;
        else verror("Unsupported track_type '%s'; expected dense/sparse/array", type_str.c_str());

        dat_path_tmp = track_dir + "/track.dat.tmp";
        idx_path_tmp = track_dir + "/track.idx.tmp";
        const string dat_path = track_dir + "/track.dat";
        const string idx_path = track_dir + "/track.idx";

        FILE *dat_fp = fopen(dat_path_tmp.c_str(), "wb");
        if (!dat_fp)
            verror("Failed to create %s: %s", dat_path_tmp.c_str(), strerror(errno));
        FILE *idx_fp = fopen(idx_path_tmp.c_str(), "wb");
        if (!idx_fp) {
            fclose(dat_fp);
            verror("Failed to create %s: %s", idx_path_tmp.c_str(), strerror(errno));
        }

        // Header (checksum=0 for now; updated at end). Shared format
        // definition lives in TrackIndexWriter / TrackIndex.h.
        try {
            TrackIndexWriter::write_header(idx_fp, track_type, (uint32_t)n_chroms);
        } catch (TGLException &) {
            fclose(dat_fp); fclose(idx_fp);
            throw;
        }

        vector<TrackContigEntry> entries;
        vector<string> chr_files_to_remove;
        uint64_t current_offset = 0;

        const size_t BUF = 1 << 20;
        vector<char> buffer(BUF);

        for (int chromid = 0; chromid < n_chroms; ++chromid) {
            const string chr_file = track_dir + "/" + chrom_names[chromid];

            TrackContigEntry entry;
            entry.chrom_id = (uint32_t)chromid;
            entry.offset = current_offset;
            entry.length = 0;
            entry.reserved = 0;

            FILE *src_fp = fopen(chr_file.c_str(), "rb");
            if (src_fp) {
                if (fseeko(src_fp, 0, SEEK_END) != 0) {
                    fclose(src_fp); fclose(dat_fp); fclose(idx_fp);
                    verror("Failed to size %s", chr_file.c_str());
                }
                const uint64_t file_size = (uint64_t)ftello(src_fp);
                rewind(src_fp);

                uint64_t remaining = file_size;
                while (remaining > 0) {
                    size_t to_read = (size_t)min((uint64_t)BUF, remaining);
                    size_t got = fread(buffer.data(), 1, to_read, src_fp);
                    if (got != to_read) {
                        fclose(src_fp); fclose(dat_fp); fclose(idx_fp);
                        verror("Short read from %s", chr_file.c_str());
                    }
                    if (fwrite(buffer.data(), 1, got, dat_fp) != got) {
                        fclose(src_fp); fclose(dat_fp); fclose(idx_fp);
                        verror("Failed to write track.dat");
                    }
                    remaining -= got;
                }
                fclose(src_fp);
                entry.length = file_size;
                current_offset += file_size;
                chr_files_to_remove.push_back(chr_file);
            }
            // else: entry stays length=0, no per-chrom file present.

            try {
                TrackIndexWriter::write_entry(idx_fp, entry);
            } catch (TGLException &) {
                fclose(dat_fp); fclose(idx_fp);
                verror("Failed to write index entry for %s", chrom_names[chromid].c_str());
            }
            entries.push_back(entry);
        }

        // Compute and patch checksum.
        try {
            TrackIndexWriter::finalize_checksum(idx_fp, entries);
        } catch (TGLException &) {
            fclose(dat_fp); fclose(idx_fp);
            throw;
        }

        fflush(dat_fp); fflush(idx_fp);
        fsync(fileno(dat_fp)); fsync(fileno(idx_fp));
        fclose(dat_fp); fclose(idx_fp);

        if (rename(dat_path_tmp.c_str(), dat_path.c_str()) != 0)
            verror("Failed to rename %s to %s: %s", dat_path_tmp.c_str(), dat_path.c_str(), strerror(errno));
        if (rename(idx_path_tmp.c_str(), idx_path.c_str()) != 0)
            verror("Failed to rename %s to %s: %s", idx_path_tmp.c_str(), idx_path.c_str(), strerror(errno));

        // Validate track.dat size matches what we wrote, before destroying source files.
        struct stat dat_stat;
        if (stat(dat_path.c_str(), &dat_stat) != 0)
            verror("Failed to stat %s after pack: %s", dat_path.c_str(), strerror(errno));
        if ((uint64_t)dat_stat.st_size != current_offset)
            verror("track.dat size mismatch after pack: expected %llu bytes, got %llu bytes",
                   (unsigned long long)current_offset,
                   (unsigned long long)dat_stat.st_size);

        // Remove old per-chrom files (always remove; this is a destructive pack)
        for (const string &p : chr_files_to_remove) unlink(p.c_str());

        return R_NilValue;
    } catch (TGLException &e) {
        unlink(dat_path_tmp.c_str());
        unlink(idx_path_tmp.c_str());
        verror("%s", e.msg());
    } catch (const bad_alloc &) {
        unlink(dat_path_tmp.c_str());
        unlink(idx_path_tmp.c_str());
        verror("Out of memory");
    }
    return R_NilValue;
}

} // extern "C"
