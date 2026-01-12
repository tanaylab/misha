# Track and intervals set management functions

.gdb.resolve_db_for_path <- function(path) {
    groots <- get("GROOTS", envir = .misha)
    if (is.null(groots)) {
        groots <- get("GROOT", envir = .misha)
    }
    if (is.null(groots)) {
        return(NULL)
    }

    root_tracks_dirs <- file.path(groots, "tracks")
    normalized_path <- normalizePath(path, mustWork = FALSE)
    idx <- which(startsWith(normalized_path, root_tracks_dirs))[1]
    if (!is.na(idx)) {
        groots[idx]
    } else {
        NULL
    }
}

.gdb.resolve_track_db <- function(track) {
    groots <- get("GROOTS", envir = .misha)
    if (is.null(groots)) {
        groots <- get("GROOT", envir = .misha)
    }
    if (is.null(groots)) {
        return(NULL)
    }

    rel_path <- paste0(gsub("\\.", "/", track), ".track")
    for (groot in rev(groots)) {
        if (file.exists(file.path(groot, "tracks", rel_path))) {
            return(groot)
        }
    }
    NULL
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

        track_db <- get("GTRACK_DB", envir = .misha)
        if (is.null(track_db)) {
            track_db <- list()
        }
        track_db[[track]] <- db
        assign("GTRACK_DB", track_db, envir = .misha)

        .gdb.cache_update_lists(db)
    }
}

# Clean up empty directories after track operations
.cleanup_empty_dirs <- function(dir) {
    gwd <- get("GWD", envir = .misha)
    tracks_dir <- gwd

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
        track_db <- get("GTRACK_DB", envir = .misha)

        if (is.null(new_db)) {
            tracks <- tracks[tracks != track]
            if (!is.null(track_db) && track %in% names(track_db)) {
                track_db[[track]] <- NULL
            }
        } else {
            if (!(track %in% tracks)) {
                tracks <- sort(c(tracks, track))
            }
            if (is.null(track_db)) {
                track_db <- list()
            }
            track_db[[track]] <- new_db
        }

        assign("GTRACKS", tracks, envir = .misha)
        assign("GTRACK_DB", track_db, envir = .misha)

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

        intervals_db <- get("GINTERVALS_DB", envir = .misha)
        if (is.null(intervals_db)) {
            intervals_db <- list()
        }
        intervals_db[[intervals.set]] <- db
        assign("GINTERVALS_DB", intervals_db, envir = .misha)

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

        intervals_db <- get("GINTERVALS_DB", envir = .misha)
        if (!is.null(intervals_db) && intervals.set %in% names(intervals_db)) {
            intervals_db[[intervals.set]] <- NULL
            assign("GINTERVALS_DB", intervals_db, envir = .misha)
        }

        .gdb.cache_update_lists(db)
    }
}
