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
        # Fallback: rename failed (cross-filesystem, or a concurrent mover
        # beat us to it). Try a direct unlink. If even that doesn't fully
        # clear path (e.g. permission denied on a sub-file), report failure
        # so callers don't scrub their caches prematurely.
        unlink(path, recursive = TRUE, force = TRUE)
        if (file.exists(path)) {
            return(invisible(FALSE))
        }
        .gdb.invalidate_dir_cache(path)
        return(invisible(TRUE))
    }

    # Path's contents are gone; drop any cached index entry keyed by it
    # before recreate-at-same-path resurrects a stale entry.
    .gdb.invalidate_dir_cache(path)

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

# The host tag embedded in a staging directory name. Only entries carrying the
# current host's tag may be judged by whether their pid is still alive - a pid
# from another machine says nothing about a process here, and the databases are
# on shared NFS.
.gdb.host_tag <- function() {
    node <- tryCatch(Sys.info()[["nodename"]], error = function(e) "")
    if (!length(node) || is.na(node) || !nzchar(node)) {
        node <- "unknown"
    }
    node <- sub("\\..*$", "", node)
    node <- gsub("[^A-Za-z0-9_]", "_", node)
    if (!nzchar(node)) "unknown" else node
}

# Name of the hidden staging directory a writer builds a track in before
# renaming it into place. Encodes host and pid so that a leftover from a
# process that is provably gone can be removed at once (see .gdb.trash_sweep_old).
.gdb.staging_name <- function(base) {
    sprintf(".%s.tmp.%s-%d.%s", base, .gdb.host_tag(), Sys.getpid(), basename(tempfile("")))
}

# TRUE when a staging directory's owning process is provably gone: same host,
# and no process with that pid. Anything else (another host, a live pid, an
# unparseable or pre-5.11.22 name without a host tag) is left to the age rule.
# The host tag is what keeps this safe on a shared database: a pid is only
# meaningful next to the machine - and the pid namespace - that issued it, and
# a container reports its own nodename, so its leftovers fall to the age rule.
.gdb.staging_owner_gone <- function(name) {
    m <- regmatches(name, regexec("\\.tmp\\.([A-Za-z0-9_]+)-([0-9]+)\\.[^.]+$", name))[[1]]
    if (length(m) != 3L) {
        return(FALSE)
    }
    if (!identical(m[2], .gdb.host_tag())) {
        return(FALSE)
    }
    if (!dir.exists("/proc")) {
        return(FALSE)
    }
    pid <- suppressWarnings(as.integer(m[3]))
    if (is.na(pid) || pid == Sys.getpid()) {
        return(FALSE)
    }
    !file.exists(file.path("/proc", pid))
}

# Sweep leftovers from a parent directory: .trash.* handoffs and the hidden
# staging dirs the atomic writers build tracks in. Called by gdb.init and by
# every writer before it stages, so a track killed mid-write does not leave a
# partial copy sitting in the database forever.
#
# A staging dir is removed when its owning process is provably gone, or, failing
# that, once it is older than max_age_hours - a live writer on another host may
# still be filling it.
.gdb.trash_sweep_old <- function(parent, max_age_hours = 24) {
    if (!dir.exists(parent)) {
        return(invisible(0L))
    }
    trash <- list.files(parent, pattern = "^\\.trash\\.", all.files = TRUE, full.names = TRUE)
    staging <- list.files(parent, pattern = "^\\..*\\.tmp\\.[A-Za-z0-9_-]+\\.[^.]+$", all.files = TRUE, full.names = TRUE)
    staging <- setdiff(staging, trash)
    entries <- c(trash, staging)
    if (!length(entries)) {
        return(invisible(0L))
    }
    cutoff <- Sys.time() - as.difftime(max_age_hours, units = "hours")
    info <- file.info(entries)
    # NA mtimes (typically permission-denied stat) are intentionally treated as
    # "not stale" so we don't loop on entries we can't even inspect.
    stale <- !is.na(info$mtime) & info$mtime < cutoff
    orphan <- entries %in% staging & vapply(basename(entries), .gdb.staging_owner_gone, logical(1))
    stale <- entries[stale | orphan]
    for (s in stale) {
        try(unlink(s, recursive = TRUE, force = TRUE), silent = TRUE)
    }
    invisible(length(stale))
}
