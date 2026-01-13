# Track management functions (exists, path, ls, info, rm)

#' Tests for a track existence
#'
#' Tests for a track existence.
#'
#' This function returns 'TRUE' if a track exists in Genomic Database.
#'
#' @param track track name
#' @return 'TRUE' if a track exists. Otherwise 'FALSE'.
#' @seealso \code{\link{gtrack.ls}}, \code{\link{gtrack.info}},
#' \code{\link{gtrack.create}}, \code{\link{gtrack.rm}}
#' @keywords ~track
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gtrack.exists("dense_track")
#'
#' @export gtrack.exists
gtrack.exists <- function(track = NULL) {
    if (is.null(substitute(track))) {
        stop("Usage: gtrack.exists(track)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    !is.na(match(trackstr, get("GTRACKS", envir = .misha)))
}


#' Returns the path on disk of a track
#'
#' Returns the path on disk of a track.
#'
#' This function returns the actual file system path where a track is stored.
#' The function works with a single track name or a vector of track names.
#'
#' @param track track name or a vector of track names
#' @return A character vector containing the full paths to the tracks on disk.
#' @seealso \code{\link{gtrack.exists}}, \code{\link{gtrack.ls}},
#' \code{\link{gintervals.path}}
#' @keywords ~track ~path
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gtrack.path("dense_track")
#' gtrack.path(c("dense_track", "sparse_track"))
#'
#' @export gtrack.path
gtrack.path <- function(track = NULL) {
    if (is.null(substitute(track))) {
        stop("Usage: gtrack.path(track)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())

    # Handle vectorized input
    if (length(trackstr) == 0) {
        return(character(0))
    }

    # Use .track_dir function for each track
    paths <- vapply(trackstr, .track_dir, character(1), USE.NAMES = FALSE)

    paths
}

#' Returns the database path(s) for track(s)
#'
#' Returns the database path where a track is stored.
#'
#' This function returns the database path(s) for one or more tracks. When
#' multiple databases are connected, each track is resolved to its source
#' database using "last wins" semantics - if the same track exists in multiple
#' databases, the last connected database takes precedence.
#'
#' @param track track name or a vector of track names
#' @return A character vector containing the database paths for each track.
#' Returns NA for tracks that don't exist in any connected database.
#' @seealso \code{\link{gtrack.exists}}, \code{\link{gtrack.path}},
#' \code{\link{gtrack.ls}}, \code{\link{gdb.ls}}
#' @keywords ~track ~path ~database
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gtrack.db("dense_track")
#' gtrack.db(c("dense_track", "sparse_track"))
#'
#' @export gtrack.db
gtrack.db <- function(track = NULL) {
    if (is.null(substitute(track))) {
        stop("Usage: gtrack.db(track)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    if (length(trackstr) == 0) {
        return(character(0))
    }

    track_db <- get("GTRACK_DB", envir = .misha)
    if (is.null(track_db) || length(track_db) == 0) {
        return(rep(NA_character_, length(trackstr)))
    }

    track_db_vec <- unlist(track_db, use.names = TRUE)
    if (is.null(track_db_vec) || length(track_db_vec) == 0) {
        return(rep(NA_character_, length(trackstr)))
    }

    unname(track_db_vec[trackstr])
}

#' Returns the database paths that contain track(s)
#'
#' Returns all database paths that contain a version of a track.
#'
#' When multiple databases are connected, a track can exist in more than one
#' database. This function returns all such database paths in connection order.
#'
#' @param track track name or a vector of track names
#' @param dataframe return a data frame with columns \code{track} and \code{db}
#' instead of a named character vector.
#' @return A named character vector of database paths for each track. If
#' \code{dataframe} is TRUE, returns a data frame with columns \code{track} and
#' \code{db}, with multiple rows per track when it appears in multiple databases.
#' @seealso \code{\link{gtrack.db}}, \code{\link{gtrack.exists}},
#' \code{\link{gtrack.ls}}, \code{\link{gdb.ls}}
#' @keywords ~track ~database ~path
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gtrack.dbs("dense_track")
#' gtrack.dbs(gtrack.ls(), dataframe = TRUE)
#'
#' @export gtrack.dbs
gtrack.dbs <- function(track = NULL, dataframe = FALSE) {
    if (is.null(substitute(track))) {
        stop("Usage: gtrack.dbs(track)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    if (length(trackstr) == 0) {
        if (dataframe) {
            return(data.frame(track = character(0), db = character(0), stringsAsFactors = FALSE))
        }
        return(character(0))
    }

    if (length(trackstr) > 1) {
        if (!dataframe) {
            return(unlist(lapply(trackstr, gtrack.dbs, dataframe = FALSE), use.names = TRUE))
        }
        res <- lapply(trackstr, function(t) gtrack.dbs(t, dataframe = TRUE))
        return(do.call(rbind, res))
    }

    if (!(trackstr %in% get("GTRACKS", envir = .misha))) {
        .gstop_track_not_found(trackstr)
    }

    track_dbs <- get("GTRACK_DBS", envir = .misha)
    if (is.null(track_dbs) || !(trackstr %in% names(track_dbs))) {
        dbs <- gtrack.db(trackstr)
    } else {
        dbs <- track_dbs[[trackstr]]
    }

    if (!dataframe) {
        names(dbs) <- rep(trackstr, length(dbs))
        return(dbs)
    }

    data.frame(
        track = rep(trackstr, length(dbs)),
        db = dbs,
        stringsAsFactors = FALSE
    )
}

#' Returns information about a track
#'
#' Returns information about a track.
#'
#' Returns information about the track (type, dimensions, size in bytes, etc.).
#' The fields in the returned value vary depending on the type of the track.
#'
#' @param track track name
#' @param validate if TRUE, validates the track index file integrity (for indexed tracks). Default: FALSE
#' @return A list that contains track properties
#' @seealso \code{\link{gtrack.exists}}, \code{\link{gtrack.ls}}
#' @keywords ~track ~info ~property
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gtrack.info("dense_track")
#' gtrack.info("rects_track")
#'
#' @export gtrack.info
gtrack.info <- function(track = NULL, validate = FALSE) {
    if (is.null(substitute(track))) {
        stop("Usage: gtrack.info(track)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    # Use context switching to ensure GWD points to the correct database
    .with_track_context(trackstr, function() {
        .gcall("gtrackinfo", trackstr, validate, .misha_env())
    })
}

#' Returns a list of track names
#'
#' Returns a list of track names in Genomic Database.
#'
#' This function returns a list of tracks whose name or track attribute value
#' match a pattern (see 'grep'). If called without any arguments all tracks are
#' returned.
#'
#' If pattern is specified without a track attribute (i.e. in the form of
#' 'pattern') then filtering is applied to the track names. If pattern is
#' supplied with a track attribute (i.e. in the form of 'name = pattern') then
#' track attribute is matched against the pattern.
#'
#' Multiple patterns are applied one after another. The resulted list of tracks
#' should match all the patterns.
#'
#' When multiple databases are connected, the 'db' parameter can be used to
#' filter tracks to only those from a specific database. For prefixed databases
#' (those with a `.misha` config file containing a `prefix` field), track names
#' are returned with their prefix (e.g., "at@my.track").
#'
#' @param ... these arguments are of either form 'pattern' or 'attribute =
#' pattern'
#' @param db optional database path or prefix string to filter tracks. If specified,
#' only tracks from that database are returned. Can be a path (e.g., "/path/to/db")
#' or a prefix (e.g., "at" for tracks from the database with prefix "at").
#' @param ignore.case,perl,fixed,useBytes see 'grep'
#' @return An array that contains the names of tracks that match the supplied
#' patterns.
#' @seealso \code{\link{grep}}, \code{\link{gtrack.exists}},
#' \code{\link{gtrack.create}}, \code{\link{gtrack.rm}}, \code{\link{gtrack.db}}
#' @keywords ~intervals ~ls
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' # get all track names
#' gtrack.ls()
#'
#' # get track names that match the pattern "den*"
#' gtrack.ls("den*")
#'
#' # get track names whose "created.by" attribute match the pattern
#' # "create_sparse"
#' gtrack.ls(created.by = "create_sparse")
#'
#' # get track names whose names match the pattern "den*" and whose
#' # "created.by" attribute match the pattern "track"
#' gtrack.ls("den*", created.by = "track")
#'
#' @export gtrack.ls
gtrack.ls <- function(..., db = NULL, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
    .gcheckroot()

    args <- as.list(substitute(list(...)))[-1L]
    args <- list(...)

    tracks <- get("GTRACKS", envir = .misha)

    if (is.null(tracks) || !length(tracks)) {
        return(NULL)
    }

    # Filter by database if specified
    if (!is.null(db)) {
        # Try to resolve db argument as prefix first, then as path
        db_path <- .gresolve_db_arg(db)
        if (is.null(db_path)) {
            # Fallback: try normalizing as path
            db_path <- tryCatch(
                normalizePath(db, mustWork = FALSE),
                error = function(e) db
            )
        }

        track_db <- get("GTRACK_DB", envir = .misha)
        if (is.null(track_db)) {
            return(NULL)
        }
        track_db_vec <- unlist(track_db, use.names = TRUE)
        db_by_track <- track_db_vec[tracks]
        tracks <- tracks[!is.na(db_by_track) & db_by_track == db_path]
        if (length(tracks) == 0) {
            return(NULL)
        }
    }

    if (length(args) >= 1) {
        attrs <- c()
        patterns <- c()

        # first filter out file names (this filtering is faster than filtering by track variable)
        for (i in 1:length(args)) {
            arg <- as.character(args[[i]])
            if (is.null(names(args)) || names(args)[i] == "") {
                tracks <- grep(arg, tracks, value = TRUE, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
            } else {
                attrs <- c(attrs, names(args)[i])
                patterns <- c(patterns, as.character(args[[i]]))
            }
        }

        # filter out by attributes
        if (length(attrs)) {
            attrs_table <- .gcall("gget_tracks_attrs", tracks, attrs, .misha_env())
            if (is.null(attrs_table)) {
                return(NULL)
            }

            cols <- colnames(attrs_table)
            for (i in 1:length(attrs)) {
                idx <- which(cols == attrs[i])[1]
                if (!is.na(idx)) {
                    attrs_table <- subset(attrs_table, grepl(patterns[i], attrs_table[, idx], ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes))
                    if (!nrow(attrs_table)) {
                        return(NULL)
                    }
                }
            }
            tracks <- rownames(attrs_table)
        }
    }

    tracks
}

#' Renames or moves a track
#'
#' Renames a track or moves it to a different namespace within the same database,
#' or moves it between databases when using prefixed names.
#'
#' This function renames a track or moves it to a different namespace (directory).
#' When both source and destination have the same prefix (or both are unprefixed),
#' a simple file rename is performed. When different prefixes are used
#' (e.g., "src@track" to "dst@track"), the track is copied to the destination
#' database and then deleted from the source.
#'
#' @param src source track name, optionally with prefix (e.g., "at@my.track")
#' @param dest destination track name, optionally with prefix (e.g., "al@my.track")
#' @return None.
#' @seealso \code{\link{gtrack.copy}}, \code{\link{gtrack.rm}},
#' \code{\link{gtrack.exists}}, \code{\link{gtrack.ls}}
#' @keywords ~track
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gtrack.create_sparse("test_track", "Test", gintervals(1, 0, 100), 1)
#' gtrack.mv("test_track", "renamed_track")
#' gtrack.exists("renamed_track")
#' gtrack.rm("renamed_track", force = TRUE)
#'
#' @export gtrack.mv
gtrack.mv <- function(src = NULL, dest = NULL) {
    if (is.null(substitute(src)) || is.null(substitute(dest))) {
        stop("Usage: gtrack.mv(src, dest)", call. = FALSE)
    }
    .gcheckroot()

    srcname <- do.call(.gexpr2str, list(substitute(src)), envir = parent.frame())
    destname <- do.call(.gexpr2str, list(substitute(dest)), envir = parent.frame())

    if (srcname == destname) {
        stop("Source and destination track names are the same", call. = FALSE)
    }

    # Check source exists
    if (!(srcname %in% get("GTRACKS", envir = .misha))) {
        .gstop_track_not_found(srcname)
    }

    # Check destination doesn't exist
    if (destname %in% get("GTRACKS", envir = .misha)) {
        stop(sprintf("Track %s already exists", destname), call. = FALSE)
    }

    # Parse source and destination prefixes
    src_parsed <- .gparse_prefixed_name(srcname)
    dest_parsed <- .gparse_prefixed_name(destname)

    src_dir <- .track_dir(srcname)
    src_db <- .gtrack_db_path(srcname)

    # Determine destination database
    if (!is.null(dest_parsed$prefix)) {
        dest_db <- .gprefix_to_db(dest_parsed$prefix)
        if (is.null(dest_db)) {
            .gstop_unknown_prefix(dest_parsed$prefix)
        }
    } else {
        # Unprefixed dest - use source database for same-db rename, or GWD
        dest_db <- if (!is.null(src_db)) src_db else .gdb.resolve_db_for_path(get("GWD", envir = .misha))
    }

    dest_base <- file.path(dest_db, "tracks")
    dest_dir <- file.path(dest_base, paste0(gsub("\\.", "/", dest_parsed$name), ".track"))

    # Check if cross-database move
    same_db <- identical(src_db, dest_db)

    # Check write permission
    .gcheck_write_permission(src_dir, "move track from")

    if (same_db) {
        # Same database - use file rename
        dest_parent <- dirname(dest_dir)
        if (!dir.exists(dest_parent)) {
            dir.create(dest_parent, recursive = TRUE, showWarnings = FALSE)
        }

        success <- file.rename(src_dir, dest_dir)
        if (!success) {
            stop(sprintf("Failed to rename track %s to %s", srcname, destname), call. = FALSE)
        }

        # Clean up empty parent directories from source
        .cleanup_empty_dirs(dirname(src_dir))
    } else {
        # Cross-database move - copy then delete
        .gcheck_write_permission(dest_base, "create track in")

        dest_parent <- dirname(dest_dir)
        if (!dir.exists(dest_parent)) {
            dir.create(dest_parent, recursive = TRUE, showWarnings = FALSE)
        }

        # Copy directory
        success <- .copy_dir_recursive(src_dir, dest_dir)
        if (!success) {
            stop(sprintf("Failed to copy track %s to %s", srcname, destname), call. = FALSE)
        }

        # Delete source
        unlink(src_dir, recursive = TRUE)

        # Clean up empty parent directories from source
        .cleanup_empty_dirs(dirname(src_dir))
    }

    # Qualify the destination name
    qualified_destname <- .gqualify_name(dest_parsed$name, dest_db)

    # Update cache
    .gdb.rm_track(srcname)
    .gdb.add_track(qualified_destname, dest_db)

    invisible()
}

#' Copies a track
#'
#' Creates a copy of an existing track.
#'
#' This function creates a copy of a track. The destination database is determined
#' by the prefix in the destination name. If the destination has a prefix
#' (e.g., "al@track"), the track is created in that database. If no prefix is
#' specified, the track is created in the current working directory (.misha$GWD).
#'
#' This allows copying tracks between databases:
#' \code{gtrack.copy("src@track", "dst@track")} copies from the "src" database
#' to the "dst" database.
#'
#' @param src source track name, optionally with prefix (e.g., "at@my.track")
#' @param dest destination track name, optionally with prefix (e.g., "al@my.track")
#' @return None.
#' @seealso \code{\link{gtrack.mv}}, \code{\link{gtrack.rm}},
#' \code{\link{gtrack.exists}}, \code{\link{gtrack.ls}}
#' @keywords ~track
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gtrack.copy("dense_track", "dense_track_copy")
#' gtrack.exists("dense_track_copy")
#' gtrack.rm("dense_track_copy", force = TRUE)
#'
#' @export gtrack.copy
gtrack.copy <- function(src = NULL, dest = NULL) {
    if (is.null(substitute(src)) || is.null(substitute(dest))) {
        stop("Usage: gtrack.copy(src, dest)", call. = FALSE)
    }
    .gcheckroot()

    srcname <- do.call(.gexpr2str, list(substitute(src)), envir = parent.frame())
    destname <- do.call(.gexpr2str, list(substitute(dest)), envir = parent.frame())

    if (srcname == destname) {
        stop("Source and destination track names are the same", call. = FALSE)
    }

    # Check source exists
    if (!(srcname %in% get("GTRACKS", envir = .misha))) {
        .gstop_track_not_found(srcname)
    }

    # Check destination doesn't exist
    if (destname %in% get("GTRACKS", envir = .misha)) {
        stop(sprintf("Track %s already exists", destname), call. = FALSE)
    }

    src_dir <- .track_dir(srcname)

    # Parse destination prefix
    dest_parsed <- .gparse_prefixed_name(destname)

    # Determine destination database
    if (!is.null(dest_parsed$prefix)) {
        dest_db <- .gprefix_to_db(dest_parsed$prefix)
        if (is.null(dest_db)) {
            .gstop_unknown_prefix(dest_parsed$prefix)
        }
        dest_gwd <- file.path(dest_db, "tracks")
    } else {
        dest_gwd <- get("GWD", envir = .misha)
        dest_db <- .gdb.resolve_db_for_path(dest_gwd)
    }

    dest_dir <- file.path(dest_gwd, paste0(gsub("\\.", "/", dest_parsed$name), ".track"))

    # Check write permission for destination
    .gcheck_write_permission(dest_gwd, "copy track to")

    # Create destination parent directory if needed
    dest_parent <- dirname(dest_dir)
    if (!dir.exists(dest_parent)) {
        dir.create(dest_parent, recursive = TRUE, showWarnings = FALSE)
    }

    # Copy track contents using helper
    if (!.copy_dir_recursive(src_dir, dest_dir)) {
        stop(sprintf("Failed to copy track %s to %s", srcname, destname), call. = FALSE)
    }

    # Qualify the destination name
    qualified_destname <- .gqualify_name(dest_parsed$name, dest_db)

    # Update cache
    .gdb.add_track(qualified_destname, dest_db)

    invisible()
}

#' Deletes a track
#'
#' Deletes a track.
#'
#' This function deletes a track from the Genomic Database. By default
#' 'gtrack.rm' requires the user to interactively confirm the deletion. Set
#' 'force' to 'TRUE' to suppress the user prompt.
#'
#' For prefixed databases, you can specify the track with its prefix
#' (e.g., "at@my.track") to delete from a specific database.
#'
#' @param track track name, optionally with a prefix (e.g., "at@my.track")
#' @param force if 'TRUE', suppresses user confirmation of a named track removal
#' @param db optional database path or prefix string to delete the track from when
#' multiple databases are connected and the track name is not prefixed
#' @return None.
#' @seealso \code{\link{gtrack.exists}}, \code{\link{gtrack.ls}},
#' \code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.smooth}}
#' @keywords ~track
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gtrack.create("new_track", "Test track", "2 * dense_track")
#' gtrack.exists("new_track")
#' gtrack.rm("new_track", force = TRUE)
#' gtrack.exists("new_track")
#'
#' @export gtrack.rm
gtrack.rm <- function(track = NULL, force = FALSE, db = NULL) {
    if (is.null(substitute(track))) {
        stop("Usage: gtrack.rm(track, force = FALSE, db = NULL)", call. = FALSE)
    }
    .gcheckroot()

    trackname <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())

    # Parse prefix from track name if present
    parsed <- .gparse_prefixed_name(trackname)
    base_name <- parsed$name

    if (!is.null(db)) {
        # db parameter can be a path or a prefix string
        db_path <- .gresolve_db_arg(db)
        if (is.null(db_path)) {
            db_path <- tryCatch(normalizePath(db, mustWork = FALSE), error = function(e) db)
        }
        groots <- get("GROOTS", envir = .misha)
        if (is.null(groots)) {
            groots <- get("GROOT", envir = .misha)
        }
        if (!is.null(groots) && !(db_path %in% groots)) {
            stop(sprintf("Database %s is not connected", db), call. = FALSE)
        }
        dirname <- file.path(db_path, "tracks", paste0(gsub("\\.", "/", base_name), ".track"))
    } else {
        dirname <- .track_dir(trackname)
    }

    # check whether track appears among GTRACKS
    if (!(trackname %in% get("GTRACKS", envir = .misha))) {
        if (force) {
            unlink(dirname, recursive = TRUE)
            .gdb.rm_track(trackname, trackdir = dirname, db = if (!is.null(db)) db_path else NULL)
            return(invisible())
        }
        .gstop_track_not_found(trackname)
    }

    if (!is.null(db) && !file.exists(dirname)) {
        if (force) {
            return(invisible())
        }
        .gstop_track_not_found(trackname, sprintf("Track %s does not exist in database %s", trackname, db))
    }

    # Check write permission before attempting to delete
    .gcheck_write_permission(dirname, "delete track from")

    answer <- "N"
    if (force) {
        answer <- "Y"
    } else {
        str <- sprintf("Are you sure you want to delete track %s (Y/N)? ", trackname)
        message(str)
        answer <- toupper(readLines(n = 1))
    }

    if (answer == "Y" || answer == "YES") {
        # remove the track
        unlink(dirname, recursive = TRUE)

        if (dir.exists(dirname)) {
            message(sprintf("Failed to delete track %s", trackname))
        } else {
            # refresh the list of GTRACKS, etc.
            .gdb.rm_track(trackname, trackdir = dirname, db = if (!is.null(db)) db_path else NULL)
        }
    }
}
