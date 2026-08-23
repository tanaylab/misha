#' Create a minimal test database for dataset/multi-db tests
#'
#' Creates a minimal misha database with the specified chromosome sizes.
#' Useful for tests that need isolated databases without external dependencies.
#'
#' @param path Path where the database should be created
#' @param chrom_sizes Data frame with chrom and size columns
#' @return Invisible path to the created database
create_test_db <- function(path, chrom_sizes = data.frame(chrom = c("chr1", "chr2"), size = c(10000, 10000))) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(path, "tracks"), showWarnings = FALSE)
    dir.create(file.path(path, "seq"), showWarnings = FALSE)

    write.table(chrom_sizes, file.path(path, "chrom_sizes.txt"),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )

    # Create dummy sequence files
    for (chr in chrom_sizes$chrom) {
        writeLines(
            paste0(rep("A", chrom_sizes$size[chrom_sizes$chrom == chr]), collapse = ""),
            file.path(path, "seq", paste0(chr, ".seq"))
        )
    }

    invisible(path)
}

#' Path to the shared, read-only test database
#'
#' This is the *source* the per-run overlays are built from. Nothing in the
#' test suite should root a session here directly - use `test_db_root()`.
shared_test_db_path <- function(indexed = getOption("gmulticontig.indexed_format", FALSE)) {
    if (indexed) {
        "/net/mraid20/ifs/wisdom/tanay_lab/tgdata/db/tgdb/misha_test_db_indexed/"
    } else {
        "/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"
    }
}

# ---------------------------------------------------------------------------
# Per-run overlays of the shared test database
# ---------------------------------------------------------------------------
#
# The shared database is 5.9 TB in 106k files, so it cannot be copied. It does
# not have to be. Everything a test writes lands in a *namespace* directory
# under tracks/ (`test`, `temp`, `global`, ...), and misha tells a namespace
# from a leaf purely by the name: src/FindTracksNIntervals.cpp descends into
# any directory whose name has no dot in it and FTS_SKIPs anything ending in
# `.track` / `.interv`. The scan runs with FTS_LOGICAL, so it follows symlinks.
#
# So an overlay only has to mirror the namespace directories - a handful of
# them, plus the interval-set directories (see .mirror_db_dir) - and symlink
# the leaves. It reads exactly like the shared database, costs a few tenths of
# a second and about a megabyte of symlinks, and every write a test makes lands
# in the overlay instead of in lab data.

# Recursively mirror `src` into `dst`: directories misha would descend into
# become real (writable) directories, everything else becomes a symlink.
.mirror_db_dir <- function(src, dst, skip = character(0)) {
    if (!dir.create(dst, recursive = TRUE, showWarnings = FALSE) && !dir.exists(dst)) {
        stop(sprintf("test db overlay: cannot create %s", dst), call. = FALSE)
    }
    # Hidden entries matter (a big interval set is unloadable without its
    # `.meta`), except `.trash.*` siblings - those are the debris .gdb.trash()'s
    # async unlink left behind, which misha's scanner ignores anyway.
    entries <- setdiff(list.files(src, all.files = TRUE, no.. = TRUE), skip)
    entries <- entries[!startsWith(entries, ".trash.")]
    if (!length(entries)) {
        return(invisible(0L))
    }
    paths <- file.path(src, entries)
    # Namespace directories, plus interval-set directories: those are only a
    # couple of dozen chromosome files each, and gintervals.convert_to_indexed()
    # writes an index into an existing set in place - through the symlink, into
    # lab data, if the set were linked rather than mirrored. Track directories
    # stay symlinked; they hold the 106k files, and nothing in the suite
    # rewrites a fixture track in place.
    descend <- dir.exists(paths) &
        (!grepl(".", entries, fixed = TRUE) | grepl("\\.interv$", entries))

    n <- 0L
    if (any(!descend)) {
        ok <- file.symlink(paths[!descend], file.path(dst, entries[!descend]))
        if (!all(ok)) {
            stop(sprintf(
                "test db overlay: failed to link %s into %s",
                paste(entries[!descend][!ok], collapse = ", "), dst
            ), call. = FALSE)
        }
        n <- sum(!descend)
    }
    for (i in which(descend)) {
        n <- n + .mirror_db_dir(paths[i], file.path(dst, entries[i]))
    }
    invisible(n)
}

#' Build a private, writable overlay of the shared test database
#'
#' Reads identically to `source_db`; every write goes to `path`. Nothing under
#' `source_db` is created, modified or removed.
#'
#' @param path Directory to build the overlay in (created if absent)
#' @param source_db The shared database to overlay
#' @param top_level Restrict the overlay to these top-level entries of
#'   `source_db`. NULL (the default) brings over everything.
#' @param tracks_skip Top-level entries of `source_db/tracks` to leave out.
#' @return Invisible path to the overlay
build_test_db_overlay <- function(path, source_db = shared_test_db_path(), top_level = NULL,
                                  tracks_skip = character(0)) {
    if (!dir.create(path, recursive = TRUE, showWarnings = FALSE) && !dir.exists(path)) {
        stop(sprintf(
            "build_test_db_overlay(): cannot create %s. Is tempdir() (%s) out of space?",
            path, tempdir()
        ), call. = FALSE)
    }

    # Directories misha writes into, or that hold per-database state a test may
    # rewrite. Mirrored so the writes stay local.
    mirrored <- c("tracks", "intervs", "pssms", "tmp")
    # Small per-database files a test may rewrite (gdb.set_readonly_attrs,
    # gdb.create). Copied, not linked, so the rewrite stays local.
    copied <- c("chrom_sizes.txt", ".ro_attributes")
    # The overlay builds its own.
    ignored <- c(".db.cache", ".db.cache.dirty", ".db.cache.tmp")

    present <- list.files(source_db, all.files = TRUE, no.. = TRUE)
    if (!is.null(top_level)) {
        present <- intersect(present, top_level)
        mirrored <- intersect(mirrored, top_level)
    }
    entries <- setdiff(present, c(mirrored, ignored))
    to_copy <- intersect(entries, copied)
    # Everything else at the top level is bulk read-only data (seq/, the SAM
    # and export fixtures, the vtracks file): link it.
    to_link <- setdiff(entries, to_copy)

    if (length(to_copy)) {
        if (!all(file.copy(file.path(source_db, to_copy), file.path(path, to_copy)))) {
            stop(sprintf(
                "build_test_db_overlay(): failed to copy %s into %s (out of space?)",
                paste(to_copy, collapse = ", "), path
            ), call. = FALSE)
        }
        # file.copy preserves the mode, and the indexed database ships
        # chrom_sizes.txt read-only. The overlay's copy has to be writable.
        Sys.chmod(file.path(path, to_copy), "0644", use_umask = FALSE)
    }
    if (length(to_link) && !all(file.symlink(file.path(source_db, to_link), file.path(path, to_link)))) {
        stop(sprintf("build_test_db_overlay(): failed to link %s into %s", paste(to_link, collapse = ", "), path),
            call. = FALSE
        )
    }

    for (d in mirrored) {
        if (dir.exists(file.path(source_db, d))) {
            # tracks/temp is scratch: start empty rather than linking whatever a
            # previous era of the suite left in the shared database.
            .mirror_db_dir(file.path(source_db, d), file.path(path, d),
                skip = if (d == "tracks") c("temp", tracks_skip) else character(0)
            )
        }
    }
    dir.create(file.path(path, "tracks", "temp"), recursive = TRUE, showWarnings = FALSE)

    invisible(path)
}

# One overlay per R process, keyed by the source database. Lives under
# tempdir(), which R removes on exit; a killed run leaves behind ~700 KB of
# symlinks under TMPDIR and nothing at all in the shared database.
.test_db_overlays <- new.env(parent = emptyenv())

#' Path to this process's writable view of the shared test database
#'
#' Root sessions here, never at `shared_test_db_path()`. Returns NULL when the
#' shared database is not mounted, so callers can skip.
#'
#' @param source_db The shared database to overlay
#' @return Path to the overlay, or NULL
test_db_root <- function(source_db = shared_test_db_path()) {
    if (!dir.exists(source_db)) {
        return(NULL)
    }
    key <- normalizePath(source_db, mustWork = FALSE)
    cached <- .test_db_overlays[[key]]
    if (!is.null(cached) && dir.exists(cached)) {
        return(cached)
    }
    path <- build_test_db_overlay(
        tempfile(pattern = "misha_dbroot_", tmpdir = tempdir()),
        source_db = source_db
    )
    assign(key, path, envir = .test_db_overlays)
    path
}

#' Point the session at this process's overlay of the shared test database
#'
#' Historically this rooted the session at the lab database in place and wiped
#' `tracks/temp` inside it, which meant only one test suite could run at a
#' time. It now roots at a private overlay (see `build_test_db_overlay()`), so
#' concurrent runs cannot see each other's tracks, interval sets or attributes.
#'
#' The overlay is per process and shared by every file in that process, so
#' `tracks/temp` is still emptied on each call - but it is the overlay's temp
#' directory, not the shared database's.
load_test_db <- function() {
    root <- test_db_root()
    if (is.null(root)) {
        return(invisible(NULL))
    }

    if (getOption("misha.test.verbose", FALSE)) {
        if (getOption("gmulticontig.indexed_format", FALSE)) {
            message("Loading indexed test database overlay at ", root)
        } else {
            message("Loading per-chromosome test database overlay at ", root)
        }
    }
    gsetroot(root)

    temp_dir <- file.path(root, "tracks", "temp")
    if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE)
    }
    gdb.reload()
    gdir.create("temp", showWarnings = FALSE)
    invisible(root)
}

#' Save and restore the current database state
#'
#' This is useful for tests that temporarily switch to a different database.
#' Uses withr-style automatic cleanup.
#'
#' @param env The environment to use for defer (defaults to parent frame)
local_db_state <- function(env = parent.frame()) {
    # Save current state
    original_groot <- if (exists("GROOT", envir = .misha, inherits = FALSE)) {
        get("GROOT", envir = .misha)
    } else {
        NULL
    }

    # Register cleanup. Guard against the saved GROOT having been deleted
    # between save and restore -- happens when a prior test file used
    # create_isolated_test_db(), whose deferred cleanup unlinks the tmpdir
    # but leaves .misha$GROOT pointing at it. Tests that don't actually use
    # GROOT (e.g. ggenome.implant with an explicit fasta) shouldn't fail on
    # this stale pointer.
    withr::defer(
        {
            if (!is.null(original_groot) && dir.exists(original_groot)) {
                suppressMessages(gdb.init(original_groot))
            }
        },
        envir = env
    )
}

#' Guarantee this test file leaves a usable GROOT behind
#'
#' Files that build their own temporary database and re-root into it leave
#' .misha$GROOT dangling as soon as that directory is removed (withr deletes
#' it at the end of the test or the file). Under TESTTHAT_PARALLEL the next
#' file in the same worker process inherits the dangling root and fails with
#' something unrelated to its own subject matter - "Database directory does
#' not exist", "Chromosome chr1 does not exist ... Known chromosomes: chrA",
#' or "Cannot delete track from read-only database" (a deleted directory is
#' not writable, so it reads as read-only).
#'
#' Call this once at the top of any file that re-roots. It is idempotent and
#' does nothing when the file leaves a valid root behind.
ensure_valid_groot <- function() {
    groot <- if (exists("GROOT", envir = .misha, inherits = FALSE)) {
        get("GROOT", envir = .misha)
    } else {
        NULL
    }
    if (is.null(groot) || !dir.exists(groot)) {
        root <- test_db_root()
        if (!is.null(root)) {
            suppressMessages(gdb.init(root))
        }
    }
    invisible(NULL)
}

restore_groot_on_exit <- function(envir = parent.frame()) {
    withr::defer(ensure_valid_groot(), envir = envir)
}

#' Unconditionally re-root at this process's test-database overlay
#'
#' For files that finish rooted somewhere private (a `tempfile()` database of
#' their own) and must hand the next file in the parallel worker a sane root
#' back, whether or not that private database still exists.
reset_test_db_root <- function() {
    root <- test_db_root()
    if (!is.null(root)) {
        suppressMessages(gdb.init(root))
    }
    invisible(NULL)
}

#' Create an isolated test database for this test file
#'
#' Like `test_db_root()`, but built fresh for the calling file and removed when
#' that file finishes, so a track one file leaves behind is invisible to the
#' next. Namespace directories under tracks/ (`test`, `global`, ...) are real
#' directories whose contents are symlinked, so a track created as `test.foo`
#' lands here and not in the shared database.
#'
#' @return Path to the isolated test database
create_isolated_test_db <- function() {
    source_db <- shared_test_db_path()

    # Create unique temp dir for this test file/process
    testdb_dir <- tempfile(pattern = "misha_testdb_", tmpdir = tempdir())

    # Same top-level content this helper has always exposed - deliberately less
    # than the full database, so switching it to an overlay does not change what
    # the 100+ files using it can see. Failures are raised, not swallowed: a
    # silently half-built DB (a full tempdir() is the usual cause) surfaces much
    # later as a baffling "Interval test.fixedbin does not exist" elsewhere.
    # ... including the historical quirk that only *directories* under tracks/
    # were brought over, so the single-file `testintervs` set stays invisible
    # here (test-gintervals1.R asserts the exact gintervals.ls() listing).
    tracks_src <- file.path(source_db, "tracks")
    tracks_files <- setdiff(
        list.files(tracks_src),
        list.dirs(tracks_src, full.names = FALSE, recursive = FALSE)
    )
    build_test_db_overlay(testdb_dir,
        source_db = source_db,
        top_level = c("chrom_sizes.txt", "intervs", "pssms", "seq", "tracks"),
        tracks_skip = tracks_files
    )

    # Capture the previous GROOT so we can restore it on teardown. Without
    # this, the next test file in the same parallel-test worker inherits a
    # missing/stale GROOT and any file-scope misha calls (or setup_db-style
    # helpers that cache GWD before re-rooting) blow up.
    prev_groot <- if (exists("GROOT", envir = .misha, inherits = FALSE)) {
        get("GROOT", envir = .misha)
    } else {
        NULL
    }

    # Set as root and reload to recognize symlinked tracks
    gsetroot(testdb_dir)
    gdb.reload()

    withr::defer(
        {
            unlink(testdb_dir, recursive = TRUE)
            current_groot <- if (exists("GROOT", envir = .misha, inherits = FALSE)) {
                get("GROOT", envir = .misha)
            } else {
                NULL
            }
            # Only restore if GROOT still points at the dir we just deleted;
            # the test may have switched to another db in the meantime.
            if (!is.null(current_groot) && identical(current_groot, testdb_dir)) {
                fallback <- test_db_root()
                if (!is.null(prev_groot) && dir.exists(prev_groot)) {
                    suppressMessages(gdb.init(prev_groot))
                } else if (!is.null(fallback)) {
                    # prev_groot was itself an isolated db that has already been
                    # torn down. Fall back to this process's overlay rather than
                    # unsetting GROOT: in a parallel run the next file in this
                    # worker may be one that relies on the ambient root, and it
                    # would otherwise fail with a confusing "does not exist" or
                    # "Chromosome chr1 does not exist" error.
                    suppressMessages(gdb.init(fallback))
                } else {
                    rm("GROOT", envir = .misha)
                }
            }
        },
        envir = parent.frame()
    )

    invisible(testdb_dir)
}
