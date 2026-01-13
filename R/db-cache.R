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
    assign("GTRACK_DB", NULL, envir = .misha)
    assign("GINTERVALS_DB", NULL, envir = .misha)
    assign("GTRACK_DBS", NULL, envir = .misha)

    gwd <- get("GWD", envir = .misha)

    # Get all connected databases
    groots <- get("GROOTS", envir = .misha)
    if (is.null(groots)) {
        # Fallback for backward compatibility
        groots <- get("GROOT", envir = .misha)
    }

    # Check if GWD is a root tracks directory of any database
    # If GWD is a subdirectory within a database, we scan only from GWD
    # (original single-db behavior - no multi-db aggregation)
    root_tracks_dirs <- file.path(groots, "tracks")
    is_root_tracks <- any(gwd == root_tracks_dirs)

    all_tracks <- character(0)
    all_intervals <- character(0)
    track_db <- list() # track_name -> db_path
    track_dbs <- list() # track_name -> db_paths (connection order)
    intervals_db <- list() # intervals_name -> db_path

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
                # Check if this database has a prefix
                db_prefix <- .gdb_to_prefix(groot)

                if (length(all_tracks)) {
                    # For prefixed databases, qualify track names
                    if (!is.null(db_prefix)) {
                        qualified_tracks <- paste0(db_prefix, .GPREFIX_SEP, all_tracks)
                    } else {
                        qualified_tracks <- all_tracks
                    }
                    track_db[qualified_tracks] <- rep(list(groot), length(qualified_tracks))
                    track_dbs[qualified_tracks] <- rep(list(groot), length(qualified_tracks))
                    all_tracks <- qualified_tracks
                }
                if (length(all_intervals)) {
                    # For prefixed databases, qualify interval names
                    if (!is.null(db_prefix)) {
                        qualified_intervals <- paste0(db_prefix, .GPREFIX_SEP, all_intervals)
                    } else {
                        qualified_intervals <- all_intervals
                    }
                    intervals_db[qualified_intervals] <- rep(list(groot), length(qualified_intervals))
                    all_intervals <- qualified_intervals
                }
            }
        }
    } else {
        # GWD is a root tracks directory - do multi-db aggregation
        # Scan each database in order (later dbs override earlier ones - "last wins")
        for (groot in groots) {
            tracks_dir <- file.path(groot, "tracks")
            if (!dir.exists(tracks_dir)) {
                next
            }

            db.filename <- file.path(groot, ".db.cache")
            res <- NULL

            # Check if we can use cache for this database
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
                    # Scan this database
                    res <- .gcall("gfind_tracks_n_intervals", tracks_dir, .misha_env())

                    # Write cache for this db
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
                db_tracks <- res[[1]]
                db_intervals <- res[[2]]

                # Check if this database has a prefix
                db_prefix <- .gdb_to_prefix(groot)

                # "Last wins": later dbs override earlier ones
                if (length(db_tracks)) {
                    # For prefixed databases, qualify track names
                    if (!is.null(db_prefix)) {
                        qualified_tracks <- paste0(db_prefix, .GPREFIX_SEP, db_tracks)
                    } else {
                        qualified_tracks <- db_tracks
                    }

                    track_db[qualified_tracks] <- rep(list(groot), length(qualified_tracks))
                    for (i in seq_along(db_tracks)) {
                        track_dbs[[qualified_tracks[i]]] <- c(track_dbs[[qualified_tracks[i]]], groot)
                    }
                    all_tracks <- c(all_tracks, qualified_tracks)
                }
                if (length(db_intervals)) {
                    # For prefixed databases, qualify interval names
                    if (!is.null(db_prefix)) {
                        qualified_intervals <- paste0(db_prefix, .GPREFIX_SEP, db_intervals)
                    } else {
                        qualified_intervals <- db_intervals
                    }

                    intervals_db[qualified_intervals] <- rep(list(groot), length(qualified_intervals))
                    all_intervals <- c(all_intervals, qualified_intervals)
                }
            }
        }
    }

    tracks <- sort(unique(all_tracks))
    intervals <- sort(unique(all_intervals))

    res <- intersect(tracks, intervals)
    if (length(res) > 0) {
        stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))
    }

    assign("GTRACKS", tracks, envir = .misha)
    assign("GINTERVS", intervals, envir = .misha)
    assign("GTRACK_DB", track_db, envir = .misha)
    assign("GINTERVALS_DB", intervals_db, envir = .misha)
    assign("GTRACK_DBS", track_dbs, envir = .misha)
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
    track_db <- get("GTRACK_DB", envir = .misha)
    if (is.null(track_db) || !(trackname %in% names(track_db))) {
        return(NULL)
    }
    track_db[[trackname]]
}

# Get the database path for an intervals set (returns NULL if not found)
.gintervals_db_path <- function(intervalsname) {
    intervals_db <- get("GINTERVALS_DB", envir = .misha)
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
    groots <- get("GROOTS", envir = .misha)
    db_path <- NULL
    if (!is.null(groots)) {
        for (g in groots) {
            if (.gpath_is_within(path, g)) {
                db_path <- g
                break
            }
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


# Resolve track creation context, handling prefixed names
#
# @param track Track name (possibly prefixed)
# @return List with:
#   - base_name: track name without prefix
#   - qualified_name: full track name with prefix if applicable
#   - db_path: target database path
#   - gwd: GWD path for creation
.gresolve_track_creation_context <- function(track) {
    parsed <- .gparse_prefixed_name(track)

    if (!is.null(parsed$prefix)) {
        db_path <- .gprefix_to_db(parsed$prefix)
        if (is.null(db_path)) {
            .gstop_unknown_prefix(parsed$prefix)
        }
        gwd <- file.path(db_path, "tracks")
        qualified_name <- track
    } else {
        # Unprefixed name - create in current GWD database
        gwd <- get("GWD", envir = .misha)
        db_path <- .gdb.resolve_db_for_path(gwd)
        if (is.null(db_path)) {
            db_path <- get("GROOT", envir = .misha)
        }
        # Qualify name if database has prefix
        qualified_name <- .gqualify_name(track, db_path)
    }

    list(
        base_name = parsed$name,
        qualified_name = qualified_name,
        db_path = db_path,
        gwd = gwd
    )
}

.gconfirmtrackcreate <- function(track) {
    # Handle prefixed track names
    ctx <- .gresolve_track_creation_context(track)

    tracks <- get("GTRACKS", envir = .misha)
    if (!is.null(tracks) && ctx$qualified_name %in% tracks) {
        track_db <- get("GTRACK_DB", envir = .misha)
        if (is.null(track_db) || !(ctx$qualified_name %in% names(track_db))) {
            stop(sprintf("Track %s already exists", ctx$qualified_name), call. = FALSE)
        }

        existing_db <- track_db[[ctx$qualified_name]]
        if (!is.null(existing_db) && identical(existing_db, ctx$db_path)) {
            stop(sprintf("Track %s already exists", ctx$qualified_name), call. = FALSE)
        }

        # Remove from cache so we can create override
        track_db[[ctx$qualified_name]] <- NULL
        assign("GTRACK_DB", track_db, envir = .misha)
    }

    # Use base name for path construction
    path <- gsub(".", "/", ctx$base_name, fixed = TRUE)
    dir <- dirname(path)
    fulldir <- paste(ctx$gwd, dir, sep = "/")
    fullpath <- sprintf("%s.track", paste(ctx$gwd, path, sep = "/"))

    if (!file.exists(fulldir)) {
        stop(sprintf("Directory %s does not exist", dir), call. = FALSE)
    }

    if (file.exists(fullpath)) {
        stop(sprintf("File %s already exists", path), call. = FALSE)
    }

    if (!is.na(match(ctx$qualified_name, get("GINTERVS", envir = .misha)))) {
        stop(sprintf("Interval %s already exists", ctx$qualified_name), call. = FALSE)
    }

    if (!is.na(match(ctx$qualified_name, gvtrack.ls()))) {
        stop(sprintf("Virtual track %s already exists", ctx$qualified_name), call. = FALSE)
    }

    if (.ggetOption(".gautocompletion", FALSE) && exists(ctx$base_name)) {
        stop(sprintf("Variable \"%s\" shadows the name of the new track.\nPlease remove this variable from the environment or switch off autocompletion mode.", ctx$base_name), call. = FALSE)
    }

    # Return context for use by caller
    ctx
}

# Internal helper to resolve item path, handling prefixed names.
# Falls back to current GWD for new items not yet in any database.
#
# @param name Track or intervals name (possibly prefixed)
# @param extension File extension (".track" or ".interv")
# @param db_lookup_fn Function to look up database for unprefixed names
# @return Full path to the item
.gitem_path <- function(name, extension, db_lookup_fn) {
    parsed <- .gparse_prefixed_name(name)
    base_name <- parsed$name

    if (!is.null(parsed$prefix)) {
        db_path <- .gprefix_to_db(parsed$prefix)
        if (is.null(db_path)) {
            .gstop_unknown_prefix(parsed$prefix)
        }
        base <- file.path(db_path, "tracks")
    } else {
        db_path <- db_lookup_fn(name)
        base <- if (!is.null(db_path)) file.path(db_path, "tracks") else get("GWD", envir = .misha)
    }
    file.path(base, paste0(gsub("\\.", "/", base_name), extension))
}


# Get the directory path for a track, resolving to the correct database.
.track_dir <- function(trackname) {
    .gitem_path(trackname, ".track", .gtrack_db_path)
}


# Get the directory path for an intervals set, resolving to the correct database.
.intervals_dir <- function(intervalsname) {
    .gitem_path(intervalsname, ".interv", .gintervals_db_path)
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

    groots <- get("GROOTS", envir = .misha)
    if (is.null(groot)) {
        groot <- get("GROOT", envir = .misha)
    }
    if (is.null(groot) || groot == "") {
        return(invisible(FALSE))
    }

    tracks <- get("GTRACKS", envir = .misha)
    intervals <- get("GINTERVS", envir = .misha)

    if (!is.null(groots) && length(groots) > 1) {
        track_db <- get("GTRACK_DB", envir = .misha)
        if (!is.null(track_db) && length(tracks)) {
            track_db_vec <- unlist(track_db, use.names = TRUE)
            if (length(track_db_vec)) {
                tracks <- tracks[!is.na(track_db_vec[tracks]) & track_db_vec[tracks] == groot]
            }
        }

        intervals_db <- get("GINTERVALS_DB", envir = .misha)
        if (!is.null(intervals_db) && length(intervals)) {
            intervals_db_vec <- unlist(intervals_db, use.names = TRUE)
            if (length(intervals_db_vec)) {
                intervals <- intervals[!is.na(intervals_db_vec[intervals]) & intervals_db_vec[intervals] == groot]
            }
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
