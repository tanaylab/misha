# Atomic-rename-then-async-unlink helper.
#
# Why: on databases with millions of contigs, a single .track directory can
# contain >1M files. unlink(recursive=TRUE) blocks for minutes, which makes
# gtrack.rm and gtrack.create's failure cleanup feel hung. file.rename to a
# hidden sibling completes in microseconds; the actual unlink runs in a
# detached child so it doesn't tie up the R session.

.gdb.trash <- function(path, async = TRUE) {
    if (!file.exists(path)) {
        return(invisible(FALSE))
    }
    parent <- dirname(path)
    base <- basename(path)
    rand <- basename(tempfile(""))
    trash <- file.path(parent, sprintf(".trash.%s.%d.%s", base, Sys.getpid(), rand))

    # suppressWarnings: cross-filesystem rename returns FALSE and warns; the
    # FALSE return is the documented signal we rely on, the warning is noise.
    if (!suppressWarnings(file.rename(path, trash))) {
        unlink(path, recursive = TRUE, force = TRUE)
        return(invisible(TRUE))
    }

    if (async) {
        cmd <- sprintf("nohup rm -rf -- %s >/dev/null 2>&1 &", shQuote(trash))
        rc <- system(cmd, wait = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
        if (rc != 0) {
            # Fork or shell launch failed - unlink synchronously so the trash
            # dir doesn't linger until the next sweep_old.
            unlink(trash, recursive = TRUE, force = TRUE)
        }
    } else {
        unlink(trash, recursive = TRUE, force = TRUE)
    }
    invisible(TRUE)
}

# Sweep stale .trash.* siblings from a parent directory. Called by gdb.init
# (Phase 1 Task 1.3) to clean up leftovers from prior sessions.
.gdb.trash_sweep_old <- function(parent, max_age_hours = 24) {
    if (!dir.exists(parent)) {
        return(invisible(0L))
    }
    entries <- list.files(parent, pattern = "^\\.trash\\.", all.files = TRUE, full.names = TRUE)
    if (!length(entries)) {
        return(invisible(0L))
    }
    cutoff <- Sys.time() - as.difftime(max_age_hours, units = "hours")
    info <- file.info(entries)
    # NA mtimes (typically permission-denied stat) are intentionally treated as
    # "not stale" so we don't loop on entries we can't even inspect.
    stale <- entries[!is.na(info$mtime) & info$mtime < cutoff]
    for (s in stale) {
        try(unlink(s, recursive = TRUE, force = TRUE), silent = TRUE)
    }
    invisible(length(stale))
}
