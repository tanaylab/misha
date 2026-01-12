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
#' filter tracks to only those from a specific database.
#'
#' @param ... these arguments are of either form 'pattern' or 'attribute =
#' pattern'
#' @param db optional database path to filter tracks. If specified, only tracks
#' from that database are returned.
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
        db <- normalizePath(db, mustWork = FALSE)
        track_db <- get("GTRACK_DB", envir = .misha)
        if (is.null(track_db)) {
            return(NULL)
        }
        track_db_vec <- unlist(track_db, use.names = TRUE)
        db_by_track <- track_db_vec[tracks]
        tracks <- tracks[!is.na(db_by_track) & db_by_track == db]
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
#' Renames a track or moves it to a different namespace within the same database.
#'
#' This function renames a track or moves it to a different namespace (directory)
#' within the same database. The track cannot be moved to a different database.
#' Use \code{\link{gtrack.copy}} followed by \code{\link{gtrack.rm}} if you need
#' to move a track between databases.
#'
#' @param src source track name
#' @param dest destination track name
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
        stop(sprintf("Track %s does not exist", srcname), call. = FALSE)
    }

    # Check destination doesn't exist
    if (destname %in% get("GTRACKS", envir = .misha)) {
        stop(sprintf("Track %s already exists", destname), call. = FALSE)
    }

    src_dir <- .track_dir(srcname)
    src_db <- .gtrack_db_path(srcname)

    # Destination must be in the same database
    dest_base <- if (!is.null(src_db)) file.path(src_db, "tracks") else get("GWD", envir = .misha)
    dest_dir <- file.path(dest_base, paste0(gsub("\\.", "/", destname), ".track"))

    # Check write permission
    .gcheck_write_permission(src_dir, "rename track from")

    # Create destination parent directory if needed
    dest_parent <- dirname(dest_dir)
    if (!dir.exists(dest_parent)) {
        dir.create(dest_parent, recursive = TRUE, showWarnings = FALSE)
    }

    # Move the track directory
    success <- file.rename(src_dir, dest_dir)
    if (!success) {
        stop(sprintf("Failed to rename track %s to %s", srcname, destname), call. = FALSE)
    }

    # Clean up empty parent directories from source
    .cleanup_empty_dirs(dirname(src_dir))

    # Update cache
    .gdb.rm_track(srcname)
    .gdb.add_track(destname, src_db)

    invisible()
}

#' Copies a track
#'
#' Creates a copy of an existing track.
#'
#' This function creates a copy of a track. The new track is created in the
#' current working directory (.misha$GWD), which may be in a different database than
#' the source track when multiple databases are connected.
#'
#' @param src source track name
#' @param dest destination track name
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
        stop(sprintf("Track %s does not exist", srcname), call. = FALSE)
    }

    # Check destination doesn't exist
    if (destname %in% get("GTRACKS", envir = .misha)) {
        stop(sprintf("Track %s already exists", destname), call. = FALSE)
    }

    src_dir <- .track_dir(srcname)

    # Destination is in current GWD (may be different database)
    gwd <- get("GWD", envir = .misha)
    dest_dir <- file.path(gwd, paste0(gsub("\\.", "/", destname), ".track"))

    # Check write permission for destination
    .gcheck_write_permission(gwd, "copy track to")

    # Create destination directory
    if (!dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)) {
        if (!dir.exists(dest_dir)) {
            stop(sprintf("Failed to create destination directory for track %s", destname), call. = FALSE)
        }
    }

    # Copy track contents
    src_files <- list.files(src_dir, full.names = TRUE, recursive = FALSE, all.files = TRUE)
    for (src_file in src_files) {
        # Skip . and .. entries
        if (basename(src_file) %in% c(".", "..")) next

        dest_file <- file.path(dest_dir, basename(src_file))
        if (dir.exists(src_file)) {
            # Copy subdirectory
            dir.create(dest_file, recursive = TRUE, showWarnings = FALSE)
            sub_files <- list.files(src_file, full.names = TRUE, recursive = TRUE, all.files = TRUE)
            for (sf in sub_files) {
                rel_path <- sub(paste0("^", src_file, "/"), "", sf)
                df <- file.path(dest_file, rel_path)
                dir.create(dirname(df), recursive = TRUE, showWarnings = FALSE)
                file.copy(sf, df, copy.mode = TRUE)
            }
        } else {
            file.copy(src_file, dest_file, copy.mode = TRUE)
        }
    }

    # Update cache - determine which database the dest is in
    dest_db <- NULL
    groots <- get("GROOTS", envir = .misha)
    for (groot in groots) {
        if (startsWith(gwd, file.path(groot, "tracks"))) {
            dest_db <- groot
            break
        }
    }

    .gdb.add_track(destname, dest_db)

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
#' @param track track name
#' @param force if 'TRUE', suppresses user confirmation of a named track removal
#' @param db optional database path to delete the track from when multiple
#' databases are connected
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
    if (!is.null(db)) {
        db <- normalizePath(db, mustWork = FALSE)
        groots <- get("GROOTS", envir = .misha)
        if (is.null(groots)) {
            groots <- get("GROOT", envir = .misha)
        }
        if (!is.null(groots) && !(db %in% groots)) {
            stop(sprintf("Database %s is not connected", db), call. = FALSE)
        }
        dirname <- file.path(db, "tracks", paste0(gsub("\\.", "/", trackname), ".track"))
    } else {
        dirname <- .track_dir(trackname)
    }

    # check whether track appears among GTRACKS
    if (!(trackname %in% get("GTRACKS", envir = .misha))) {
        if (force) {
            unlink(dirname, recursive = TRUE)
            .gdb.rm_track(trackname, trackdir = dirname, db = db)
            return(invisible())
        }
        stop(sprintf("Track %s does not exist", trackname), call. = FALSE)
    }

    if (!is.null(db) && !file.exists(dirname)) {
        if (force) {
            return(invisible())
        }
        stop(sprintf("Track %s does not exist in database %s", trackname, db), call. = FALSE)
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
            .gdb.rm_track(trackname, trackdir = dirname, db = db)
        }
    }
}
