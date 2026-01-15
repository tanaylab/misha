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

    assign("GTRACKS", NULL, envir = .misha)
    assign("GINTERVS", NULL, envir = .misha)
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
            groot_idx <- which(vapply(root_tracks_dirs, function(root) .gpath_is_within(gwd, root), logical(1)))[1]
            if (!is.na(groot_idx)) {
                groot <- groots[groot_idx]
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
                res <- .gcall("gfind_tracks_n_intervals", tracks_dir, .misha_env())
                if (!is.null(res)) {
                    res[[1]] <- .gdb.normalize_cache_list(res[[1]])
                    res[[2]] <- .gdb.normalize_cache_list(res[[2]])
                }
                try(
                    {
                        f <- file(db.filename, "wb")
                        serialize(res, f)
                        close(f)
                        .gdb.cache_clear_dirty(groot)
                    },
                    silent = TRUE
                )
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
                    # Scan this database
                    res <- .gcall("gfind_tracks_n_intervals", tracks_dir, .misha_env())

                    # Normalize cached lists for faster subsequent loads
                    if (!is.null(res)) {
                        res[[1]] <- .gdb.normalize_cache_list(res[[1]])
                        res[[2]] <- .gdb.normalize_cache_list(res[[2]])
                    }

                    # Write cache for this db
                    try(
                        {
                            f <- file(db.filename, "wb")
                            serialize(res, f)
                            close(f)
                            .gdb.cache_clear_dirty(g)
                        },
                        silent = TRUE
                    )
                } else {
                    # Ensure cached lists are sorted/unique; rewrite if needed
                    norm_tracks <- .gdb.normalize_cache_list(res[[1]])
                    norm_intervals <- .gdb.normalize_cache_list(res[[2]])
                    if (!identical(norm_tracks, res[[1]]) || !identical(norm_intervals, res[[2]])) {
                        res[[1]] <- norm_tracks
                        res[[2]] <- norm_intervals
                        try(
                            {
                                f <- file(db.filename, "wb")
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
    db_path <- .gtrack_db_path(trackname)
    if (is.null(db_path)) {
        return(fn())
    }

    correct_gwd <- file.path(db_path, "tracks")
    current_gwd <- get("GWD", envir = .misha)
    if (correct_gwd == current_gwd) {
        return(fn())
    }

    assign("GWD", correct_gwd, envir = .misha)
    on.exit(assign("GWD", current_gwd, envir = .misha), add = TRUE)
    fn()
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
                if (inherits(retv, "try-error") & length(grep("obsolete", retv)) > 0) {
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
        target_db <- .gdb.resolve_db_for_path(gwd)
        if (is.null(target_db)) {
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

        track_db <- track_db[names(track_db) != track]
        assign("GTRACK_DATASET", track_db, envir = .misha)
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

# Get the directory path for a track, resolving to the correct database.
# Falls back to current GWD for new tracks not yet in any database.
.track_dir <- function(trackname) {
    db_path <- .gtrack_db_path(trackname)
    base <- if (!is.null(db_path)) file.path(db_path, "tracks") else get("GWD", envir = .misha)
    file.path(base, paste0(gsub("\\.", "/", trackname), ".track"))
}

# Get the directory path for an intervals set, resolving to the correct database.
# Falls back to current GWD for new intervals not yet in any database.
.intervals_dir <- function(intervalsname) {
    db_path <- .gintervals_db_path(intervalsname)
    base <- if (!is.null(db_path)) file.path(db_path, "tracks") else get("GWD", envir = .misha)
    file.path(base, paste0(gsub("\\.", "/", intervalsname), ".interv"))
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

.gdb.cache_dirty_path <- function(groot = NULL) {
    cache_path <- .gdb.cache_path(groot)
    if (is.null(cache_path)) {
        return(NULL)
    }
    paste0(cache_path, ".dirty")
}

.gdb.cache_is_dirty <- function(groot = NULL) {
    dirty_path <- .gdb.cache_dirty_path(groot)
    !is.null(dirty_path) && file.exists(dirty_path)
}

.gdb.cache_mark_dirty <- function(groot = NULL) {
    dirty_path <- .gdb.cache_dirty_path(groot)
    if (is.null(dirty_path)) {
        return(invisible(FALSE))
    }
    dir.create(dirname(dirty_path), recursive = TRUE, showWarnings = FALSE)

    # Write timestamp to dirty file for debugging and to reduce race conditions
    success <- tryCatch(
        {
            writeLines(as.character(Sys.time()), dirty_path)
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
    invisible(success)
}

.gdb.cache_clear_dirty <- function(groot = NULL) {
    dirty_path <- .gdb.cache_dirty_path(groot)
    if (is.null(dirty_path)) {
        return(invisible(FALSE))
    }
    if (file.exists(dirty_path)) {
        result <- try(unlink(dirty_path), silent = TRUE)
        if (inherits(result, "try-error") || (is.integer(result) && result != 0)) {
            warning(sprintf("Failed to clear cache dirty flag at %s", dirty_path), call. = FALSE)
            return(invisible(FALSE))
        }
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

# Cache write uses atomic file rename to minimize corruption risk.
# Concurrency note: Multiple processes writing to the same database may race,
# but the dirty flag and atomic rename provide basic safety. In concurrent
# scenarios, the last writer wins, which is acceptable since both should scan
# the same filesystem state. For true multi-process safety, external locking
# would be needed.
.gdb.cache_write_lists <- function(tracks, intervals, groot = NULL) {
    cache_path <- .gdb.cache_path(groot)
    if (is.null(cache_path)) {
        return(invisible(FALSE))
    }

    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)

    tmp_path <- paste0(cache_path, ".tmp")
    .gdb.cache_mark_dirty(groot)

    success <- FALSE
    error_msg <- NULL
    tryCatch(
        {
            con <- file(tmp_path, "wb")
            on.exit(close(con), add = TRUE)
            serialize(list(tracks, intervals), con)
            close(con)
            on.exit(NULL, add = FALSE)

            if (!suppressWarnings(file.rename(tmp_path, cache_path))) {
                unlink(cache_path)
                if (!suppressWarnings(file.rename(tmp_path, cache_path))) {
                    stop("Failed to replace .db.cache file")
                }
            }
            success <- TRUE
        },
        error = function(e) {
            error_msg <<- e$message
            success <<- FALSE
        },
        finally = {
            if (file.exists(tmp_path)) {
                try(unlink(tmp_path), silent = TRUE)
            }
        }
    )

    if (!success) {
        warning(sprintf("Failed to write database cache: %s", error_msg), call. = FALSE)
    } else {
        .gdb.cache_clear_dirty(groot)
    }
    invisible(success)
}

.gdb.cache_update_lists <- function(groot = NULL) {
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

    tracks <- get("GTRACKS", envir = .misha)
    intervals <- get("GINTERVS", envir = .misha)

    if (length(groots) > 1) {
        track_db <- get("GTRACK_DATASET", envir = .misha)
        if (!is.null(track_db) && length(tracks)) {
            tracks <- tracks[!is.na(track_db[tracks]) & track_db[tracks] == groot]
        }

        intervals_db <- get("GINTERVALS_DATASET", envir = .misha)
        if (!is.null(intervals_db) && length(intervals)) {
            intervals <- intervals[!is.na(intervals_db[intervals]) & intervals_db[intervals] == groot]
        }
    }

    .gdb.cache_write_lists(tracks, intervals, groot)
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
    invisible(.gdb.cache_mark_dirty())
}
