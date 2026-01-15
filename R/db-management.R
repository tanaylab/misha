# Track and intervals set management functions

.gpath_is_within <- function(path, root) {
    if (is.null(path) || is.null(root)) {
        return(FALSE)
    }

    path <- normalizePath(path, mustWork = FALSE)
    root <- normalizePath(root, mustWork = FALSE)

    if (path == root) {
        return(TRUE)
    }

    startsWith(path, paste0(root, .Platform$file.sep))
}

.gdb.resolve_db_for_path <- function(path) {
    groot <- get("GROOT", envir = .misha)
    gdatasets <- get("GDATASETS", envir = .misha)
    if (is.null(gdatasets)) gdatasets <- character(0)
    groots <- c(groot, gdatasets)

    if (length(groots) == 0 || is.null(groot)) {
        return(NULL)
    }

    root_tracks_dirs <- file.path(groots, "tracks")
    idx <- which(vapply(root_tracks_dirs, function(root) .gpath_is_within(path, root), logical(1)))[1]
    if (!is.na(idx)) {
        groots[idx]
    } else {
        NULL
    }
}

.gdb.resolve_track_db <- function(track) {
    groot <- get("GROOT", envir = .misha)
    gdatasets <- get("GDATASETS", envir = .misha)
    if (is.null(gdatasets)) gdatasets <- character(0)
    groots <- c(groot, gdatasets)

    if (length(groots) == 0 || is.null(groot)) {
        return(NULL)
    }

    rel_path <- paste0(gsub("\\.", "/", track), ".track")
    # Working db first (wins), then datasets in load order
    for (g in groots) {
        if (file.exists(file.path(g, "tracks", rel_path))) {
            return(g)
        }
    }
    NULL
}

.gdb.ensure_dataset_maps <- function() {
    .gcheckroot()

    groot <- get("GROOT", envir = .misha)
    tracks <- get("GTRACKS", envir = .misha)
    intervals <- get("GINTERVS", envir = .misha)

    track_db <- get("GTRACK_DATASET", envir = .misha)
    intervals_db <- get("GINTERVALS_DATASET", envir = .misha)

    if (is.null(track_db)) {
        track_db <- character(0)
        if (length(tracks)) {
            track_db[tracks] <- groot
        }
        assign("GTRACK_DATASET", track_db, envir = .misha)
    }

    if (is.null(intervals_db)) {
        intervals_db <- character(0)
        if (length(intervals)) {
            intervals_db[intervals] <- groot
        }
        assign("GINTERVALS_DATASET", intervals_db, envir = .misha)
    }

    invisible(TRUE)
}

.gdb.add_track <- function(track, db = NULL) {
    .gcheckroot()

    # Get the track directory, using db if specified
    if (!is.null(db)) {
        trackdir <- file.path(db, "tracks", paste0(gsub("\\.", "/", track), ".track"))
    } else {
        trackdir <- .track_dir(track)
    }

    if (file.exists(trackdir)) {
        tracks <- sort(c(get("GTRACKS", envir = .misha), track))
        intervals <- sort(get("GINTERVS", envir = .misha))

        res <- intersect(tracks, intervals)
        if (length(res) > 0) {
            stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))
        }

        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(track, envir = .misha)) {
                stop(sprintf("Variable \"%s\" shadows the name of identically named track.\nPlease remove this variable from the environment or switch off autocompletion mode.", track), call. = FALSE)
            }

            if (.ggetOption(".ginteractive", FALSE)) { # set track to NULL otherwise evaluation of track expression pmin(track, 2) will produce a string "2"
                assign(track, NULL, envir = .misha)
            } else {
                assign(track, track, envir = .misha)
            }
        }

        assign("GTRACKS", tracks, envir = .misha)

        if (is.null(db)) {
            db <- .gdb.resolve_db_for_path(trackdir)
            if (is.null(db)) {
                db <- get("GROOT", envir = .misha)
            }
        }

        track_db <- get("GTRACK_DATASET", envir = .misha)
        if (is.null(track_db)) {
            track_db <- character(0)
        }
        track_db[[track]] <- db
        assign("GTRACK_DATASET", track_db, envir = .misha)

        .gdb.cache_update_lists(db)
    }
}

# Clean up empty directories after track operations
.cleanup_empty_dirs <- function(dir) {
    db_root <- .gdb.resolve_db_for_path(dir)
    if (!is.null(db_root)) {
        tracks_dir <- file.path(db_root, "tracks")
    } else {
        tracks_dir <- get("GWD", envir = .misha)
    }

    # Walk up from dir, removing empty directories until we hit tracks_dir
    current <- dir
    while (current != tracks_dir && nchar(current) > nchar(tracks_dir)) {
        if (dir.exists(current) && length(list.files(current, all.files = TRUE, no.. = TRUE)) == 0) {
            unlink(current, recursive = TRUE)
            current <- dirname(current)
        } else {
            break
        }
    }
}

.rm_track_dir <- function(trackname) {
    dirname <- .track_dir(trackname)
    unlink(dirname, recursive = TRUE)
    if (dir.exists(dirname)) {
        message(sprintf("Failed to remove track directory %s", dirname))
        invisible(FALSE)
    }
    invisible(TRUE)
}

.gdb.rm_track <- function(track, trackdir = NULL, db = NULL) {
    .gcheckroot()

    if (is.null(trackdir)) {
        trackdir <- .track_dir(track)
    }
    if (is.null(db)) {
        db <- .gtrack_db_path(track)
        if (is.null(db)) {
            db <- .gdb.resolve_db_for_path(trackdir)
        }
        if (is.null(db)) {
            db <- get("GROOT", envir = .misha)
        }
    }

    if (!file.exists(trackdir)) {
        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(track, envir = .misha)) {
                remove(list = track, envir = .misha)
            }
        }

        new_db <- .gdb.resolve_track_db(track)

        tracks <- get("GTRACKS", envir = .misha)
        track_db <- get("GTRACK_DATASET", envir = .misha)

        if (is.null(new_db)) {
            tracks <- tracks[tracks != track]
            if (!is.null(track_db) && track %in% names(track_db)) {
                track_db <- track_db[names(track_db) != track]
            }
        } else {
            if (!(track %in% tracks)) {
                tracks <- sort(c(tracks, track))
            }
            if (is.null(track_db)) {
                track_db <- character(0)
            }
            track_db[[track]] <- new_db
        }

        assign("GTRACKS", tracks, envir = .misha)
        assign("GTRACK_DATASET", track_db, envir = .misha)

        .gdb.cache_update_lists(db)
        if (!is.null(new_db) && new_db != db) {
            .gdb.cache_update_lists(new_db)
        }
    }
}

.gdb.add_intervals.set <- function(intervals.set) {
    .gcheckroot()

    fname <- sprintf("%s.interv", paste(get("GWD", envir = .misha), gsub("\\.", "/", intervals.set), sep = "/"))
    if (file.exists(fname)) {
        tracks <- get("GTRACKS", envir = .misha)
        intervals <- sort(c(get("GINTERVS", envir = .misha), intervals.set))

        res <- intersect(tracks, intervals)
        if (length(res) > 0) {
            stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))
        }

        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(intervals.set, envir = .misha)) {
                stop(sprintf("Variable \"%s\" shadows the name of identically named intervals set.\nPlease remove this variable from the environment or switch off autocompletion mode.", intervals.set), call. = FALSE)
            }

            assign(intervals.set, intervals.set, envir = .misha)
        }

        assign("GINTERVS", intervals, envir = .misha)

        db <- .gdb.resolve_db_for_path(fname)
        if (is.null(db)) {
            db <- get("GROOT", envir = .misha)
        }

        intervals_db <- get("GINTERVALS_DATASET", envir = .misha)
        if (is.null(intervals_db)) {
            intervals_db <- character(0)
        }
        intervals_db[[intervals.set]] <- db
        assign("GINTERVALS_DATASET", intervals_db, envir = .misha)

        .gdb.cache_update_lists(db)
    }
}

.gdb.rm_intervals.set <- function(intervals.set) {
    .gcheckroot()

    fname <- sprintf("%s.interv", paste(get("GWD", envir = .misha), gsub("\\.", "/", intervals.set), sep = "/"))
    db <- .gintervals_db_path(intervals.set)
    if (is.null(db)) {
        db <- .gdb.resolve_db_for_path(fname)
        if (is.null(db)) {
            db <- get("GROOT", envir = .misha)
        }
    }
    if (!file.exists(fname)) {
        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(intervals.set, envir = .misha)) {
                remove(list = intervals.set, envir = .misha)
            }
        }

        intervals <- get("GINTERVS", envir = .misha)
        intervals <- intervals[intervals != intervals.set]
        assign("GINTERVS", intervals, envir = .misha)

        intervals_db <- get("GINTERVALS_DATASET", envir = .misha)
        if (!is.null(intervals_db) && intervals.set %in% names(intervals_db)) {
            intervals_db <- intervals_db[names(intervals_db) != intervals.set]
            assign("GINTERVALS_DATASET", intervals_db, envir = .misha)
        }

        .gdb.cache_update_lists(db)
    }
}
