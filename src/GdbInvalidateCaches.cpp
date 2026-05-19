/*
 * Drops the per-directory entries from every process-static index cache
 * that misha maintains (1D track index, 2D track index, 1D bigset index,
 * 2D bigset index). Each cache is keyed by the absolute directory path of
 * the on-disk track or interval set; if that directory is removed,
 * recreated, or converted between formats, the cached entry becomes
 * stale and downstream readers will route to the wrong file layout
 * (e.g., open a non-existent track.dat on what is now a per-chrom track).
 *
 * Exposed to R as `.gdb.invalidate_dir_cache(dir)`. R-level mutation
 * hooks (gtrack.rm, .gtrack.create_atomic post-rename,
 * gtrack.convert_to_indexed, gintervals.rm, gintervals.convert_to_indexed,
 * the 2D variants of each) call this whenever the contents of a track or
 * interval set directory change.
 */

#include "GenomeTrack.h"
#include "TrackIndex2D.h"
#include "GIntervalsBigSet1D.h"
#include "GIntervalsBigSet2D.h"
#include "rdbutils.h"

extern "C" {

SEXP gdb_invalidate_dir_cache(SEXP _dir, SEXP _envir)
{
    try {
        if (!Rf_isString(_dir) || Rf_length(_dir) < 1)
            verror("'dir' must be a character vector");

        const int n = Rf_length(_dir);
        for (int i = 0; i < n; ++i) {
            const char *path = CHAR(STRING_ELT(_dir, i));
            if (path == nullptr || path[0] == '\0')
                continue;
            const std::string dir(path);
            GenomeTrack::invalidate_index_cache(dir);
            TrackIndex2D::invalidate_cache(dir);
            GIntervalsBigSet1D::invalidate_index_cache(dir);
            GIntervalsBigSet2D::invalidate_index_cache(dir);
        }
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    }
    return R_NilValue;
}

// Wipe every entry from every process-static index cache. Used by
// gdb.reload(rescan = TRUE) so out-of-process mutations (a sibling R
// session, a manual rm/cp, an external rebuild) cannot leave the
// caller routing reads through a stale track-type / layout entry.
SEXP gdb_clear_all_dir_caches(SEXP _envir)
{
    try {
        GenomeTrack::clear_index_cache();
        TrackIndex2D::clear_cache();
        GIntervalsBigSet1D::clear_index_cache();
        GIntervalsBigSet2D::clear_index_cache();
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    }
    return R_NilValue;
}

}
