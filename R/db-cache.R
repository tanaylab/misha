# Cache management functions
#' Reloads database from the disk
#'
#' Reloads database from disk: list of tracks, intervals, etc.
#'
#' Reloads Genomic Database from disk: list of tracks, intervals, etc. Use this
#' function if you manually add tracks or if for any reason the database
#' becomes corrupted. If 'rescan' is 'TRUE', the list of tracks and intervals
#' is achieved by rescanning directory structure under the current current
#' working directory. Otherwise 'gdb.reload' attempts to use the cached list
#' that resides in 'GROOT/.db.cache' file.
#'
#' @param rescan indicates whether the file structure should be rescanned
#' @seealso \code{\link{gdb.init}}, \code{\link{gdb.create}},
#' \code{\link{gdir.cd}},
#' @keywords ~db
#' @return No return value, called for side effects.
#' @export gdb.reload
gdb.reload <- function(rescan = TRUE) {
    if (!exists("GROOT", envir = .misha)) {
        stop("gdb.init() must be called beforehand.", call. = FALSE)
    }

    # Clear session scan cache when rescanning
    if (rescan && exists(".gdb.clear_scan_cache", mode = "function")) {
        .gdb.clear_scan_cache()
    }

    # Drop every entry from the process-static index caches the C++ side
    # maintains (track type / layout, bigset layout). Without this, an
    # out-of-process mutation (sibling R session, manual rm/cp, external
    # rebuild) leaves the caller routing reads through stale entries -
    # e.g. a track converted from sparse to dense reads back as the cached
    # sparse type until R restarts. gdb.reload(rescan = TRUE) is the
    # documented "rescan from disk" hook; the C++ caches must follow.
    if (rescan && exists(".gdb.clear_all_dir_caches", mode = "function")) {
        .gdb.clear_all_dir_caches()
    }

    assign("GTRACKS", NULL, envir = .misha)
    assign("GINTERVS", NULL, envir = .misha)
    assign("GTRACKS_SRC", NULL, envir = .misha)
    assign("GTRACK_DATASET", NULL, envir = .misha)
    assign("GINTERVALS_DATASET", NULL, envir = .misha)

    gwd <- get("GWD", envir = .misha)

    # Get working database and loaded datasets
    groot <- get("GROOT", envir = .misha)
    gdatasets <- get("GDATASETS", envir = .misha)
    if (is.null(gdatasets)) gdatasets <- character(0)
    groots <- c(groot, gdatasets)

    # Check if GWD is a root tracks directory of any database
    # If GWD is a subdirectory within a database, we scan only from GWD
    # (original single-db behavior - no multi-db aggregation)
    root_tracks_dirs <- file.path(groots, "tracks")
    is_root_tracks <- any(gwd == root_tracks_dirs)

    all_tracks <- character(0)
    all_intervals <- character(0)
    track_db <- character(0) # track_name -> db_path
    track_dbs <- list() # track_name -> db_paths (connection order)
    intervals_db <- character(0) # intervals_name -> db_path

    if (!is_root_tracks) {
        # GWD is a subdirectory - scan only from GWD (original behavior)
        # Force rescan when in a subdirectory
        rescan <- TRUE
        res <- .gcall("gfind_tracks_n_intervals", gwd, .misha_env())

        if (!is.null(res)) {
            all_tracks <- res[[1]]
            all_intervals <- res[[2]]
            # Determine which database this GWD belongs to
            gwd_db <- .gdb.gwd_db()
            if (!is.na(gwd_db)) {
                groot <- gwd_db
                if (length(all_tracks)) {
                    track_db[all_tracks] <- rep(groot, length(all_tracks))
                    track_dbs[all_tracks] <- rep(list(groot), length(all_tracks))
                }
                if (length(all_intervals)) {
                    intervals_db[all_intervals] <- rep(groot, length(all_intervals))
                }
            }
        }
    } else if (length(gdatasets) == 0) {
        # Single working database fast path (avoid per-track dataset mapping)
        tracks_dir <- file.path(groot, "tracks")
        db.filename <- file.path(groot, ".db.cache")
        res <- NULL

        use_cache <- !rescan && !.gdb.cache_is_dirty(groot)
        suppressWarnings({
            if (use_cache) {
                retv <- try(
                    {
                        f <- file(db.filename, "rb")
                        res <- unserialize(f)
                        close(f)
                    },
                    silent = TRUE
                )
                if (inherits(retv, "try-error")) {
                    res <- NULL
                }
            }

            if (is.null(res)) {
                res <- .gdb.cache_rescan(groot)
            }
        })

        if (!is.null(res)) {
            all_tracks <- res[[1]]
            all_intervals <- res[[2]]
        }
    } else {
        # GWD is a root tracks directory - do multi-db aggregation
        # Scan each database in order (working db first, then datasets)
        # Working db always wins for collision resolution
        working_db <- groot
        for (g in groots) {
            tracks_dir <- file.path(g, "tracks")
            if (!dir.exists(tracks_dir)) {
                next
            }

            db.filename <- file.path(g, ".db.cache")
            res <- NULL

            # Check if we can use cache for this database
            use_cache <- !rescan && !.gdb.cache_is_dirty(g)

            suppressWarnings({
                if (use_cache) {
                    retv <- try(
                        {
                            f <- file(db.filename, "rb")
                            res <- unserialize(f)
                            close(f)
                        },
                        silent = TRUE
                    )

                    if (inherits(retv, "try-error")) {
                        res <- NULL
                    }
                }

                if (is.null(res)) {
                    # Scan this database, and write its cache
                    res <- .gdb.cache_rescan(g)
                } else {
                    # Ensure cached lists are sorted/unique; rewrite if needed
                    norm_tracks <- .gdb.normalize_cache_list(res[[1]])
                    norm_intervals <- .gdb.normalize_cache_list(res[[2]])
                    if (!identical(norm_tracks, res[[1]]) || !identical(norm_intervals, res[[2]])) {
                        res[[1]] <- norm_tracks
                        res[[2]] <- norm_intervals
                        try(
                            {
                                f <- .gwith_umask(file(db.filename, "wb"))
                                serialize(res, f)
                                close(f)
                            },
                            silent = TRUE
                        )
                    }
                }
            })

            if (!is.null(res)) {
                db_tracks <- res[[1]]
                db_intervals <- res[[2]]

                # Working db always wins; for datasets, later ones override earlier
                # Use vectorized operations to avoid O(n^2) performance
                if (length(db_tracks)) {
                    # Get existing track names efficiently
                    existing_track_names <- names(track_db)
                    # Find tracks that need to be added/updated
                    # Add if: not in track_db OR not from working_db
                    if (length(existing_track_names) == 0) {
                        # No existing tracks - add all
                        new_entries <- setNames(rep(g, length(db_tracks)), db_tracks)
                        track_db <- c(track_db, new_entries)
                    } else {
                        # Check which tracks are new or need override
                        is_new <- !(db_tracks %in% existing_track_names)
                        # For existing tracks, check if they're NOT from working_db
                        existing_tracks <- db_tracks[!is_new]
                        if (length(existing_tracks) > 0) {
                            existing_dbs <- track_db[existing_tracks]
                            should_override <- existing_dbs != working_db
                            tracks_to_update <- existing_tracks[should_override]
                        } else {
                            tracks_to_update <- character(0)
                        }
                        tracks_to_add <- c(db_tracks[is_new], tracks_to_update)
                        if (length(tracks_to_add) > 0) {
                            # Remove tracks that will be overridden before adding new entries
                            if (length(tracks_to_update) > 0) {
                                track_db <- track_db[!(names(track_db) %in% tracks_to_update)]
                            }
                            new_entries <- setNames(rep(g, length(tracks_to_add)), tracks_to_add)
                            track_db <- c(track_db, new_entries)
                        }
                    }
                    all_tracks <- c(all_tracks, db_tracks)
                }
                if (length(db_intervals)) {
                    # Same vectorized approach for intervals
                    existing_interval_names <- names(intervals_db)
                    if (length(existing_interval_names) == 0) {
                        new_entries <- setNames(rep(g, length(db_intervals)), db_intervals)
                        intervals_db <- c(intervals_db, new_entries)
                    } else {
                        is_new <- !(db_intervals %in% existing_interval_names)
                        existing_intervals <- db_intervals[!is_new]
                        if (length(existing_intervals) > 0) {
                            existing_dbs <- intervals_db[existing_intervals]
                            should_override <- existing_dbs != working_db
                            intervals_to_update <- existing_intervals[should_override]
                        } else {
                            intervals_to_update <- character(0)
                        }
                        intervals_to_add <- c(db_intervals[is_new], intervals_to_update)
                        if (length(intervals_to_add) > 0) {
                            # Remove intervals that will be overridden before adding new entries
                            if (length(intervals_to_update) > 0) {
                                intervals_db <- intervals_db[!(names(intervals_db) %in% intervals_to_update)]
                            }
                            new_entries <- setNames(rep(g, length(intervals_to_add)), intervals_to_add)
                            intervals_db <- c(intervals_db, new_entries)
                        }
                    }
                    all_intervals <- c(all_intervals, db_intervals)
                }
            }
        }
    }

    tracks <- all_tracks
    if (length(tracks) > 1) {
        if (is.unsorted(tracks, strictly = FALSE) || anyDuplicated(tracks)) {
            tracks <- sort(unique(tracks))
        }
    }

    intervals <- all_intervals
    if (length(intervals) > 1) {
        if (is.unsorted(intervals, strictly = FALSE) || anyDuplicated(intervals)) {
            intervals <- sort(unique(intervals))
        }
    }

    res <- intersect(tracks, intervals)
    if (length(res) > 0) {
        stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))
    }

    assign("GTRACKS", tracks, envir = .misha)
    assign("GINTERVS", intervals, envir = .misha)
    # Provenance stamp: which database set this listing was actually scanned
    # from. .gdb.cache_update_lists() refuses to persist the in-memory listing
    # to <groot>/.db.cache unless the stamp still matches, so a listing that
    # belongs to another database can never be published to this one.
    assign("GTRACKS_SRC", groots, envir = .misha)
    if (length(gdatasets) > 0) {
        assign("GTRACK_DATASET", track_db, envir = .misha)
        assign("GINTERVALS_DATASET", intervals_db, envir = .misha)
    } else {
        assign("GTRACK_DATASET", NULL, envir = .misha)
        assign("GINTERVALS_DATASET", NULL, envir = .misha)
    }
}

# Multi-database context switching helper.
# Temporarily switches GWD to the correct database for a track operation,
# restoring the original GWD on exit. This is necessary because C++ functions
# use GWD to resolve track paths.
.with_track_context <- function(trackname, fn) {
    # in a subdirectory (gdir.cd()) every listed name is under GWD
    if (!is.null(.gdb.gwd_db())) {
        return(fn())
    }
    .with_db_context(.gtrack_db_path(trackname), fn)
}

# Get the database path for a track (returns NULL if not found)
.gtrack_db_path <- function(trackname) {
    track_db <- get("GTRACK_DATASET", envir = .misha)
    if (is.null(track_db) || !(trackname %in% names(track_db))) {
        return(NULL)
    }
    track_db[[trackname]]
}

# Get the database path for an intervals set (returns NULL if not found)
.gintervals_db_path <- function(intervalsname) {
    intervals_db <- get("GINTERVALS_DATASET", envir = .misha)
    if (is.null(intervals_db) || !(intervalsname %in% names(intervals_db))) {
        return(NULL)
    }
    intervals_db[[intervalsname]]
}

# Check write permission for a path and provide clear error message.
# Called before write operations to give user-friendly errors when trying
# to write to read-only databases.
.gcheck_write_permission <- function(path, operation = "write to") {
    dir_path <- if (dir.exists(path)) path else dirname(path)

    if (file.access(dir_path, 2) == 0) {
        return(invisible(TRUE))
    }

    # Find which database this path belongs to for a better error message
    groot <- get("GROOT", envir = .misha)
    gdatasets <- get("GDATASETS", envir = .misha)
    if (is.null(gdatasets)) gdatasets <- character(0)
    groots <- c(groot, gdatasets)

    db_path <- NULL
    for (g in groots) {
        if (.gpath_is_within(path, g)) {
            db_path <- g
            break
        }
    }

    if (!is.null(db_path)) {
        stop(sprintf("Cannot %s read-only database: %s", operation, db_path), call. = FALSE)
    }
    stop(sprintf("Cannot %s read-only path: %s", operation, dir_path), call. = FALSE)
}


.gdb.convert_attrs <- function() {
    .gcheckroot()

    ro_attrs <- c("created.by", "created.date", "created.user")
    .gcall_noninteractive(gdb.set_readonly_attrs, ro_attrs)

    for (track in .misha$GTRACKS) {
        for (attr in ro_attrs) {
            try(
                {
                    if (.gcall_noninteractive(.gtrack.var.exists, track, attr)) {
                        .gcall_noninteractive(.gtrack.attr.set, track, attr, as.character(.gtrack.var.get(track, attr))[1], TRUE)
                        .gcall_noninteractive(gtrack.var.rm, track, attr)
                    }
                },
                silent = TRUE
            )
        }
        message(track)
    }
}

.gdb.convert_tracks <- function() {
    .gcheckroot()

    for (track in .misha$GTRACKS) {
        try(
            {
                retv <- try(.gcall_noninteractive(gtrack.info, track))
                if (inherits(retv, "try-error") && length(grep("obsolete", retv)) > 0) {
                    message(sprintf("Converting track %s", track))
                    .gcall_noninteractive(gtrack.convert, track)
                }
            },
            silent = TRUE
        )
    }
}


.gconfirmtrackcreate <- function(track) {
    tracks <- get("GTRACKS", envir = .misha)
    if (!is.null(tracks) && track %in% tracks) {
        gwd <- get("GWD", envir = .misha)
        # in a subdirectory (gdir.cd()), the database the session lists the names in
        target_db <- .gdb.gwd_db() %||% .gdb.resolve_db_for_path(gwd)
        if (is.null(target_db) || is.na(target_db)) {
            target_db <- get("GROOT", envir = .misha)
        }

        track_db <- get("GTRACK_DATASET", envir = .misha)
        if (is.null(track_db) || !(track %in% names(track_db))) {
            stop(sprintf("Track %s already exists", track), call. = FALSE)
        }

        existing_db <- track_db[[track]]
        if (!is.null(existing_db) && identical(existing_db, target_db)) {
            stop(sprintf("Track %s already exists", track), call. = FALSE)
        }

        # The name is taken off the dataset for the creation to write into the working database
        # (.track_dir()); a successful one maps it there (.gdb.add_track()). If the creation fails, the
        # dataset's track is mapped back when the creating function exits.
        track_db <- track_db[names(track_db) != track]
        assign("GTRACK_DATASET", track_db, envir = .misha)
        restore <- function() {
            track_db <- get("GTRACK_DATASET", envir = .misha)
            if (!is.null(track_db) && !(track %in% names(track_db)) && track %in% get("GTRACKS", envir = .misha)) {
                track_db[[track]] <- existing_db
                assign("GTRACK_DATASET", track_db, envir = .misha)
            }
        }
        # (the function itself in the call, not a name looked up in the caller's frame)
        do.call(on.exit, list(as.call(list(restore)), add = TRUE), envir = parent.frame())
    }

    # In a subdirectory (gdir.cd()), a track of an indexed working database is written in the session's
    # chromosome ids, so only into a loaded database that uses them (.gdb.same_order()): not through a
    # link or into a database kept below GWD (.gdb.path_db()), nor into a dataset of unknown order
    if (!is.null(.gdb.gwd_db()) && .gdb.is_indexed() && !.gdb.same_order(.gdb.path_db(.track_dir(track)))) {
        stop(sprintf("Track %s would be written into a database that may number the chromosomes differently from this session; create it in a session of that database.", track), call. = FALSE)
    }

    path <- gsub(".", "/", track, fixed = TRUE)
    dir <- dirname(path)
    fulldir <- paste(get("GWD", envir = .misha), dir, sep = "/")
    fullpath <- sprintf("%s.track", paste(get("GWD", envir = .misha), path, sep = "/"))

    if (!file.exists(fulldir)) {
        stop(sprintf("Directory %s does not exist", dir), call. = FALSE)
    }

    if (file.exists(fullpath)) {
        stop(sprintf("File %s already exists", path), call. = FALSE)
    }

    if (!is.na(match(track, get("GINTERVS", envir = .misha)))) {
        stop(sprintf("Interval %s already exists", track), call. = FALSE)
    }

    if (!is.na(match(track, gvtrack.ls()))) {
        stop(sprintf("Virtual track %s already exists", track), call. = FALSE)
    }

    if (.ggetOption(".gautocompletion", FALSE) && exists(track)) {
        stop(sprintf("Variable \"%s\" shadows the name of the new track.\nPlease remove this variable from the environment or switch off autocompletion mode.", track), call. = FALSE)
    }
}

# Get the directory path for a track, in database `db` (default: the one it is listed in). Names are
# relative to the database's tracks/, or to GWD when that is a subdirectory of it (gdir.cd()), where
# the session lists only the names under GWD (.gdb.gwd_db()). Falls back to GWD for new tracks not
# yet in any database.
.track_dir <- function(trackname, db = NULL) {
    if (is.null(db) && is.null(.gdb.gwd_db())) db <- .gtrack_db_path(trackname)
    file.path(.gdb.names_dir(db), paste0(gsub("\\.", "/", trackname), ".track"))
}

# Get the directory path for an intervals set, as .track_dir()
.intervals_dir <- function(intervalsname, db = NULL) {
    if (is.null(db) && is.null(.gdb.gwd_db())) db <- .gintervals_db_path(intervalsname)
    file.path(.gdb.names_dir(db), paste0(gsub("\\.", "/", intervalsname), ".interv"))
}

.gdb.cache_path <- function(groot = NULL) {
    if (is.null(groot)) {
        if (!exists("GROOT", envir = .misha, inherits = FALSE)) {
            return(NULL)
        }
        groot <- get("GROOT", envir = .misha)
    }
    if (is.null(groot) || groot == "") {
        return(NULL)
    }
    file.path(groot, ".db.cache")
}

# The dirty flags of a database's .db.cache: one file per mark, <groot>/.db.cache.dirty.<token>, so
# that no mark overwrites another. A flag is removed only by the session that made it
# (.gdb.cache_clear_dirty) or by a rescan that started after it existed (.gdb.cache_rescan), and the
# cache is dirty while any exists. <groot>/.db.cache.dirty, the single flag of earlier versions, counts
# as one too (and is removed the same way). No sweep of leftovers touches them: a flag left by a
# killed session costs one rescan, which removes it.
.gdb.cache_dirty_files <- function(groot = NULL) {
    cache_path <- .gdb.cache_path(groot)
    if (is.null(cache_path)) {
        return(character(0))
    }
    list.files(dirname(cache_path), pattern = "^\\.db\\.cache\\.dirty(\\..+)?$", all.files = TRUE, full.names = TRUE, no.. = TRUE)
}

.gdb.cache_is_dirty <- function(groot = NULL) {
    length(.gdb.cache_dirty_files(groot)) > 0
}

# Marks the cache dirty with a flag file of this mark's own. Returns the mark's token, or NULL if
# it was not made: quietly for a database root that cannot be written (nothing writes its cache
# then either).
.gdb.cache_mark_dirty <- function(groot = NULL) {
    cache_path <- .gdb.cache_path(groot)
    if (is.null(cache_path) || file.access(dirname(cache_path), 2) != 0) {
        return(invisible(NULL))
    }
    token <- paste(.gdb.host_tag(), Sys.getpid(), basename(tempfile("")), sep = "-")
    success <- tryCatch(
        {
            .gwith_umask(writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%OS6"), paste0(cache_path, ".dirty.", token)))
            TRUE
        },
        warning = function(w) {
            warning(sprintf("Failed to mark cache dirty: %s", w$message), call. = FALSE)
            FALSE
        },
        error = function(e) {
            warning(sprintf("Failed to mark cache dirty: %s", e$message), call. = FALSE)
            FALSE
        }
    )
    invisible(if (success) token)
}

# Removes the dirty flag this session made with `token` (.gdb.cache_mark_dirty()); the others stay.
# NULL removes nothing.
.gdb.cache_clear_dirty <- function(groot = NULL, token = NULL) {
    cache_path <- .gdb.cache_path(groot)
    if (is.null(cache_path) || is.null(token)) {
        return(invisible(FALSE))
    }
    path <- paste0(cache_path, ".dirty.", token)
    if (!file.exists(path)) {
        return(invisible(FALSE))
    }
    if (unlink(path) != 0 || file.exists(path)) {
        warning(sprintf("Failed to clear cache dirty flag at %s", path), call. = FALSE)
        return(invisible(FALSE))
    }
    invisible(TRUE)
}

# Normalize cached track/interval lists to sorted unique vectors.
.gdb.normalize_cache_list <- function(x) {
    if (length(x) <= 1) {
        return(x)
    }
    if (!is.unsorted(x, strictly = FALSE) && !anyDuplicated(x)) {
        return(x)
    }
    sort(unique(x))
}

# The cache is updated one name at a time: a write that created or removed a track (kind "tracks")
# or an interval set ("intervals") reads the database's .db.cache as it is on disk, adds or removes
# that one name, and writes it back (a temporary file renamed over it). The session's own lists
# are not written: another session may have created or removed names since this one loaded them.
# A cache that is dirty (marked by anyone: some change is missing from it), missing or unreadable
# is left dirty for the next gsetroot()/gdb.reload() to rescan. Two sessions writing at the same
# instant can still lose one change (there is no lock); the flag each sets is cleared only by its
# own writer, so the later write leaves it set when the other got in between.
.gdb.cache_update_lists <- function(groot = NULL, name = NULL, kind = c("tracks", "intervals"), add = TRUE) {
    kind <- match.arg(kind)
    if (!exists("GROOT", envir = .misha, inherits = FALSE) ||
        !exists("GTRACKS", envir = .misha, inherits = FALSE) ||
        !exists("GINTERVS", envir = .misha, inherits = FALSE)) {
        return(invisible(FALSE))
    }

    working_db <- get("GROOT", envir = .misha)
    gdatasets <- get("GDATASETS", envir = .misha)
    if (is.null(gdatasets)) gdatasets <- character(0)
    groots <- c(working_db, gdatasets)

    if (is.null(groot)) {
        groot <- working_db
    }
    if (is.null(groot) || groot == "") {
        return(invisible(FALSE))
    }

    # GTRACKS/GINTERVS may describe a *different* database than the one we are
    # about to write to. gdb.reload() is the only thing that scans a database and
    # it stamps GTRACKS_SRC with the groots it scanned; if GROOT moved on without
    # a completed reload - gsetroot() interrupted mid-reload, or a caller that
    # restores a memoized session snapshot (misha.ext::gset_genome) - the stamp
    # no longer matches, and neither may the names. Mark dirty instead and let
    # the next gsetroot()/gdb.reload() rescan from disk.
    src <- if (exists("GTRACKS_SRC", envir = .misha, inherits = FALSE)) {
        get("GTRACKS_SRC", envir = .misha)
    } else {
        NULL
    }
    if (!identical(src, groots)) {
        .gdb.cache_mark_dirty(groot)
        return(invisible(FALSE))
    }

    # If the working directory is a subdirectory of the database (not a root
    # "tracks" dir), names are relative to it (see .gdb.reload, which scans only
    # from GWD in that case), not to tracks/ as the cache's are; and a call that
    # names no change cannot be merged. Mark the cache dirty (forcing a rescan).
    gwd <- get("GWD", envir = .misha)
    if (!any(gwd == file.path(groots, "tracks")) || is.null(name)) {
        .gdb.cache_mark_dirty(groot)
        return(invisible(FALSE))
    }

    cache_path <- .gdb.cache_path(groot)
    if (is.null(cache_path)) {
        return(invisible(FALSE))
    }

    # read, merged and written back under the lock (.gdb.cache_lock); without it, left dirty
    lock <- .gdb.cache_lock(groot)
    if (is.null(lock)) {
        .gdb.cache_mark_dirty(groot)
        return(invisible(FALSE))
    }
    # a lock broken meanwhile (held past the stale age) may have let another write in
    on.exit(if (!.gdb.cache_unlock(groot, lock)) .gdb.cache_mark_dirty(groot), add = TRUE)

    # A dirty cache misses some change: marked again, with a flag of this writer's own, which only a
    # rescan that starts after it removes, as this name is missing too
    if (.gdb.cache_is_dirty(groot)) {
        .gdb.cache_mark_dirty(groot)
        return(invisible(FALSE))
    }

    # dirty while it is rewritten: a session killed meanwhile leaves it for a rescan
    token <- .gdb.cache_mark_dirty(groot)
    if (is.null(token)) {
        return(invisible(FALSE))
    }
    read_cache <- function() {
        f <- file(cache_path, "rb")
        on.exit(close(f))
        unserialize(f)
    }
    lists <- tryCatch(read_cache(), error = function(e) NULL, warning = function(w) NULL)
    if (!is.list(lists) || length(lists) < 2) {
        return(invisible(FALSE))
    }
    i <- if (kind == "tracks") 1L else 2L
    names_now <- as.character(lists[[i]])
    lists[[i]] <- if (add) .gdb.normalize_cache_list(c(names_now, name)) else setdiff(names_now, name)

    success <- tryCatch(
        .gdb.cache_write_file(cache_path, lists),
        error = function(e) {
            warning(sprintf("Failed to write database cache: %s", e$message), call. = FALSE)
            FALSE
        }
    )
    if (success) {
        .gdb.cache_clear_dirty(groot, token)
    }
    invisible(success)
}

# Writes the lists of a .db.cache (a temporary file of this writer's renamed over it)
.gdb.cache_write_file <- function(cache_path, lists) {
    tmp_path <- paste0(cache_path, ".tmp.", Sys.getpid(), ".", basename(tempfile("")))
    on.exit(if (file.exists(tmp_path)) try(unlink(tmp_path), silent = TRUE), add = TRUE)
    # The rename below carries this file's mode to cache_path, so the umask has to be right here,
    # at the open.
    .gwith_umask(.gserialize_file(lists, tmp_path))
    if (!suppressWarnings(file.rename(tmp_path, cache_path))) {
        stop("Failed to replace .db.cache file", call. = FALSE)
    }
    TRUE
}

# The lock of a database's .db.cache: <groot>/.db.cache.lock, created exclusively (open mode "wx",
# O_EXCL) and holding its holder's token. It is held only around a writer's read-merge-rename and a
# rescan's marking and final write. A lock older than `stale` seconds belongs to a holder that was
# killed, and is broken (.gdb.cache_break_lock). Returns the token, or NULL when the lock could not
# be taken within `wait` seconds (or cannot be created at all: a read-only database).
.gdb.cache_lock <- function(groot, wait = 5, stale = 30) {
    cache_path <- .gdb.cache_path(groot)
    if (is.null(cache_path)) {
        return(NULL)
    }
    path <- paste0(cache_path, ".lock")
    token <- paste(.gdb.host_tag(), Sys.getpid(), basename(tempfile("")))
    deadline <- Sys.time() + wait
    repeat {
        con <- tryCatch(suppressWarnings(.gwith_umask(file(path, "wx"))), error = function(e) NULL)
        if (!is.null(con)) {
            written <- tryCatch(
                {
                    writeLines(token, con)
                    TRUE
                },
                error = function(e) FALSE
            )
            closed <- tryCatch(
                {
                    close(con)
                    TRUE
                },
                error = function(e) FALSE
            )
            if (written && closed) {
                return(token)
            }
            unlink(path)
            return(NULL)
        }
        if (!file.exists(path) && file.access(dirname(path), 2) != 0) {
            return(NULL)
        }
        # a stale lock is broken and the lock taken again at once; on every other path (held, or a
        # break that failed: an unreadable lock, a root that cannot be written) the deadline is
        # checked and the next attempt waits
        mtime <- file.info(path)$mtime
        broken <- !is.na(mtime) && difftime(Sys.time(), mtime, units = "secs") > stale && .gdb.cache_break_lock(path)
        if (Sys.time() > deadline) {
            return(NULL)
        }
        if (!broken) {
            Sys.sleep(stats::runif(1, 0.005, 0.03))
        }
    }
}

# Breaks a stale lock: it is renamed to a name of its own (of two breakers, one wins the rename),
# then removed if it is still the lock that was judged stale. A fresh lock taken in between is put
# back, or, if the name is taken again, removed: its holder then finds the lock not its own when it
# unlocks, and marks the cache dirty (.gdb.cache_unlock).
.gdb.cache_break_lock <- function(path) {
    read_lock <- function(p) tryCatch(suppressWarnings(readLines(p, warn = FALSE)), error = function(e) NULL)
    held <- read_lock(path)
    if (is.null(held)) {
        return(invisible(FALSE))
    }
    moved <- paste0(path, ".stale.", basename(tempfile("")))
    if (!suppressWarnings(file.rename(path, moved))) {
        return(invisible(FALSE))
    }
    if (!identical(read_lock(moved), held) && !file.exists(path) && suppressWarnings(file.rename(moved, path))) {
        return(invisible(FALSE))
    }
    unlink(moved)
    invisible(TRUE)
}

# Releases the lock if it is still this holder's (`token`). FALSE if it is not: it was broken
# meanwhile, and another writer may have written the cache too.
.gdb.cache_unlock <- function(groot, token) {
    path <- paste0(.gdb.cache_path(groot), ".lock")
    held <- tryCatch(suppressWarnings(readLines(path, warn = FALSE)), error = function(e) NULL)
    if (!identical(held, token)) {
        return(invisible(FALSE))
    }
    unlink(path)
    invisible(TRUE)
}

# A rescan of a database's tracks and interval sets (gdb.reload, .gdb.scan_db_fast), written to its
# .db.cache. Under the lock, before the scan, it notes the dirty flags there are (their changes were
# made before the scan, which sees them) and adds a flag of its own, so that no writer merges into
# the cache during the scan (.gdb.cache_update_lists leaves a dirty cache alone). Under the lock
# again, it writes the scan only if its own flag is still there: a rescan that started after this one
# removed it, and wrote a newer scan. Having written, it removes its own flag and the ones it noted;
# a flag made during the scan stays. Nothing is written to a database that cannot be written.
# Returns the scan.
.gdb.cache_rescan <- function(groot) {
    tracks_dir <- file.path(groot, "tracks")
    cache_path <- .gdb.cache_path(groot)
    writable <- file.access(groot, 2) == 0
    lock <- NULL
    # (an interrupt does not leave the lock behind for the stale age)
    on.exit(if (!is.null(lock)) .gdb.cache_unlock(groot, lock), add = TRUE)
    token <- NULL
    seen <- character(0)
    if (writable) {
        lock <- .gdb.cache_lock(groot)
        seen <- .gdb.cache_dirty_files(groot)
        token <- .gdb.cache_mark_dirty(groot)
        if (!is.null(lock)) {
            held <- .gdb.cache_unlock(groot, lock)
            lock <- NULL
            if (!held) {
                .gdb.cache_mark_dirty(groot)
            }
        }
    }
    res <- .gcall("gfind_tracks_n_intervals", tracks_dir, .misha_env())
    if (!is.null(res)) {
        res[[1]] <- .gdb.normalize_cache_list(res[[1]])
        res[[2]] <- .gdb.normalize_cache_list(res[[2]])
    }
    if (!is.null(token) && !is.null(res)) {
        lock <- .gdb.cache_lock(groot)
        if (!is.null(lock)) {
            if (file.exists(paste0(cache_path, ".dirty.", token))) {
                written <- tryCatch(.gdb.cache_write_file(cache_path, res), error = function(e) FALSE)
                if (isTRUE(written)) {
                    unlink(seen)
                    .gdb.cache_clear_dirty(groot, token)
                }
            }
            held <- .gdb.cache_unlock(groot, lock)
            lock <- NULL
            if (!held) {
                .gdb.cache_mark_dirty(groot)
            }
        }
    }
    res
}

#' Mark cached track list as dirty
#'
#' When tracks or interval sets are modified outside of misha (e.g. files copied
#' manually), the cached inventory may become out of date. Calling this helper
#' marks the cache as dirty so the next \code{gsetroot()} forces a rescan.
#'
#' @return Invisible \code{TRUE} if the dirty flag was written, \code{FALSE} otherwise.
#' @seealso \code{\link{gdb.reload}}, \code{\link{gsetroot}}
#' @export gdb.mark_cache_dirty
gdb.mark_cache_dirty <- function() {
    .gcheckroot()
    invisible(!is.null(.gdb.cache_mark_dirty()))
}

# Names of files that live inside a track directory but are not per-chromosome
# data: the indexed-format pair, attribute/var sidecars, and the meta file used
# by big interval sets. Keep in sync with src/ when adding new sidecar files.
.TRACK_INTERNAL_FILES <- c("track.idx", "track.dat", ".attributes", "vars", ".meta")

# The per-chromosome files of a track or intervals set directory, as the readers find them (C++
# GenomeTrack::ChromFiles): a row per file named by a chromosome (by a pair if is_2d), named by its
# canonical name or by an alias, with the canonical chrom (chrom1, chrom2), the file, and used:
# whether it is the file that the readers open
.gchrom_files <- function(dir, is_2d) {
    as.data.frame(.gcall("gchrom_files", dir, is_2d, .misha_env()), stringsAsFactors = FALSE)
}

# Check whether a database at the given path is in indexed format.
# Unlike .gdb.is_indexed(), does not depend on the loaded GROOT.
.gdb.is_indexed_at <- function(groot) {
    if (is.null(groot) || !nzchar(groot)) {
        return(FALSE)
    }
    # A loaded dataset with no sequence of its own (its seq/ link gone, or an empty folder) is in
    # the loaded database's format; a dataset that is a database in its own right keeps its own
    loaded <- c(get0("GROOT", envir = .misha, ifnotfound = NULL), get0("GDATASETS", envir = .misha, ifnotfound = NULL))
    if (length(loaded) && normalizePath(groot, mustWork = FALSE) %in% normalizePath(loaded, mustWork = FALSE)) {
        first_chrom <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom[1])
        # sequence as .gdb.chrom_names_at() probes it
        own_seq <- all(file.exists(file.path(groot, "seq", c("genome.idx", "genome.seq")))) ||
            any(file.exists(file.path(groot, "seq", paste0(c(first_chrom, paste0("chr", first_chrom), sub("^chr", "", first_chrom)), ".seq"))))
        if (!own_seq) {
            groot <- get("GROOT", envir = .misha)
        }
    }
    seq_dir <- file.path(groot, "seq")
    file.exists(file.path(seq_dir, "genome.idx")) && file.exists(file.path(seq_dir, "genome.seq"))
}

# The chromosome names of the database at groot in chrom id order, as gsetroot() gives
# them (.gdb.chrom_order): chrom id i (0-based) is element i + 1. An indexed track of the
# database is keyed by these ids. In a per-chromosome database the names get the "chr" prefix
# and follow the sorted names, not the chrom_sizes.txt order. The loaded database and its
# loaded datasets (gdataset.load() requires the same chrom_sizes.txt) take the order gsetroot()
# gave, held in ALLGENOME, rather than probe seq/: a dataset's seq/ link may be gone or empty.
# Any other database is read from disk, without loading it, and so is a loaded one with from_disk =
# TRUE (the order its own files use); one whose seq/ holds no sequence cannot show which order its
# names take, and is an error rather than a guess.
.gdb.chrom_names_at <- function(groot, from_disk = FALSE) {
    loaded <- c(get0("GROOT", envir = .misha, ifnotfound = NULL), get0("GDATASETS", envir = .misha, ifnotfound = NULL))
    if (!from_disk && length(loaded) && normalizePath(groot, mustWork = FALSE) %in% normalizePath(loaded, mustWork = FALSE)) {
        # a dataset loaded with an order of its own (gdataset.load()) keeps it
        own <- get0("GDATASET_CHROMS", envir = .misha, ifnotfound = NULL)[[normalizePath(groot, mustWork = FALSE)]]
        return(if (is.null(own)) as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom) else own)
    }
    cs <- file.path(groot, "chrom_sizes.txt")
    if (!file.exists(cs)) {
        stop(sprintf("chrom_sizes.txt missing in %s", groot), call. = FALSE)
    }
    # read as gsetroot() reads it
    chromsizes <- utils::read.csv(
        cs,
        sep = "\t", header = FALSE, col.names = c("chrom", "size"), colClasses = c("character", "numeric")
    )
    chrom_order <- .gdb.chrom_order(groot, chromsizes)
    # Names mostly without the "chr" prefix are sorted in a per-chromosome database and kept in
    # chrom_sizes.txt order in an indexed one (.is_per_chromosome_db); a seq/ with neither
    # genome.idx and genome.seq nor the first chromosome's .seq file cannot tell which
    seq_dir <- file.path(groot, "seq")
    first_chrom <- chromsizes$chrom[1]
    # sequence as .gdb.is_indexed_at() probes it
    no_seq <- !(all(file.exists(file.path(seq_dir, c("genome.idx", "genome.seq")))) ||
        any(file.exists(file.path(seq_dir, paste0(c(first_chrom, paste0("chr", first_chrom), sub("^chr", "", first_chrom)), ".seq")))))
    if (!chrom_order$per_chromosome && nrow(chromsizes) && no_seq &&
        mean(!startsWith(chromsizes$chrom, "chr")) >= 0.8) {
        stop(sprintf(
            "The chromosome order of %s cannot be told: its chrom_sizes.txt names lack the \"chr\" prefix, so its indexed tracks are keyed by the sorted names if it is a per-chromosome database and by the chrom_sizes.txt order if it is indexed, and its seq/ has no sequence to show which (as in a copy of a database's tracks without its seq/). Put back its seq/ (a link to the original database's seq/ will do), or load it with gdataset.load() as a dataset of its database.",
            groot
        ), call. = FALSE)
    }
    # as gsetroot() checks it: an indexed seq/ linked to a database converted since gives the order
    # of chrom_sizes.txt, not of the index its tracks would be read with
    .gdb.check_genome_idx(groot, chrom_order$names[chrom_order$id_order], chromsizes$size[chrom_order$id_order])
    chrom_order$names[chrom_order$id_order]
}
