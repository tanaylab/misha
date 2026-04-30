// src/GenomeTrackSplitIndexed.cpp
//
// Splits a 1D indexed-format track (track.dat + track.idx) back into per-chromosome
// files in the same directory, named by the supplied chrom names.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <errno.h>
#include <unistd.h>
#include <string>
#include <vector>

// Note: do NOT #include <Rinternals.h> directly. It defines a length(x) macro
// that collides with TrackContigEntry::length in TrackIndex.h. R's headers are
// pulled in transitively via rdbutils.h.
#include "TrackIndex.h"
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
            if (entry.length == 0) continue;
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

            // Per-file fsync + atomic rename: each per-chrom file is canonical db state;
            // we accept N fsyncs (one per non-empty contig) for crash-safe individual files.
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

} // extern "C"
