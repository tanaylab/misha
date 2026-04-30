// src/GenomeTrackSplitIndexed.cpp
//
// Splits a 1D indexed-format track (track.dat + track.idx) back into per-chromosome
// files in the same directory, named by the supplied chrom names.

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
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
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_track_dir) || Rf_length(_track_dir) != 1)
            verror("track_dir must be a single string");
        if (!Rf_isString(_chrom_names))
            verror("chrom_names must be a character vector");

        // Stub: do nothing yet.
        return R_NilValue;
    } catch (TGLException &e) {
        verror("%s", e.msg());
    } catch (const bad_alloc &) {
        verror("Out of memory");
    }
    return R_NilValue;
}

} // extern "C"
