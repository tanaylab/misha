# Drop the process-static index cache entries that misha maintains for
# `dir` (a track or interval-set directory). Call this from R-level
# mutation points (rm, create-atomic rename, convert) so subsequent
# readers don't route through a stale `dir -> shared_ptr<TrackIndex>`
# entry. Without this, "rm <track>; create_*<track>" cycles on indexed
# DBs can route reads to a non-existent track.dat.
#
# Best-effort: failures swallowed - cache miss is not fatal.
.gdb.invalidate_dir_cache <- function(dir) {
    if (is.null(dir) || length(dir) == 0L) {
        return(invisible(NULL))
    }
    dir <- as.character(dir)
    dir <- dir[nzchar(dir)]
    if (!length(dir)) {
        return(invisible(NULL))
    }
    tryCatch(
        .gcall("gdb_invalidate_dir_cache", dir, .misha_env()),
        error = function(e) NULL
    )
    invisible(NULL)
}
