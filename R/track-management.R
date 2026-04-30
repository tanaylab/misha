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

#' Returns the database/dataset path for a track
#'
#' Returns the path of the database or dataset containing a track.
#'
#' When datasets are loaded, tracks can come from either the working database
#' or from loaded datasets. This function returns the source path for each track.
#'
#' @param track track name or a vector of track names
#' @return Character vector of database/dataset paths. Returns NA for non-existent tracks.
#' @seealso \code{\link{gtrack.dbs}}, \code{\link{gtrack.exists}},
#' \code{\link{gtrack.ls}}, \code{\link{gdataset.ls}}
#' @keywords ~track ~database ~path
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gtrack.dataset("dense_track")
#'
#' @export gtrack.dataset
gtrack.dataset <- function(track = NULL) {
    if (is.null(substitute(track))) {
        stop("Usage: gtrack.dataset(track)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    if (length(trackstr) == 0) {
        return(character(0))
    }

    track_db <- get("GTRACK_DATASET", envir = .misha)
    if (is.null(track_db) || length(track_db) == 0) {
        tracks <- get("GTRACKS", envir = .misha)
        groot <- get("GROOT", envir = .misha)
        return(ifelse(trackstr %in% tracks, groot, NA_character_))
    }

    unname(track_db[trackstr])
}

#' Returns the database paths that contain track(s)
#'
#' Returns all database paths that contain a version of a track.
#'
#' When datasets are loaded, a track may exist in multiple locations (working
#' database and/or datasets). This function computes on-demand and returns all
#' such paths, which is useful for debugging when using \code{force=TRUE} with
#' \code{gdataset.load()}.
#'
#' @param track track name or a vector of track names
#' @param dataframe return a data frame with columns \code{track} and \code{db}
#' instead of a named character vector.
#' @return A named character vector of database paths for each track. If
#' \code{dataframe} is TRUE, returns a data frame with columns \code{track} and
#' \code{db}, with multiple rows per track when it appears in multiple databases.
#' @seealso \code{\link{gtrack.dataset}}, \code{\link{gtrack.exists}},
#' \code{\link{gtrack.ls}}, \code{\link{gdataset.ls}}
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
    .gdb.resource_dbs_impl(trackstr, ".track", "track", dataframe, gtrack.dbs)
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
#' \code{\link{gtrack.create}}, \code{\link{gtrack.rm}}, \code{\link{gtrack.dataset}}
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
        track_db <- get("GTRACK_DATASET", envir = .misha)
        if (is.null(track_db)) {
            if (!identical(db, get("GROOT", envir = .misha))) {
                return(NULL)
            }
        }
        db_by_track <- track_db[tracks]
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

#' Copies one or more tracks
#'
#' Creates a copy of an existing track, optionally to a different database.
#' Transparently handles format mismatches (per-chromosome vs indexed) and
#' chromosome-order differences between source and destination databases.
#'
#' Chromosomes that exist in the source database but not in the destination
#' are dropped with a warning.
#'
#' @param src source track name(s). Either a single name or a character
#'   vector of names.
#' @param dest destination name. If \code{src} is a single name, this is the
#'   destination track name (defaults to \code{src}). If \code{src} is a
#'   vector, \code{dest} is treated as a namespace prefix (e.g. \code{"ns"}
#'   produces \code{"ns.track1"}, \code{"ns.track2"}, ...). NULL keeps each
#'   track's name.
#' @param db destination database root. Must be the current \code{GROOT} or a
#'   member of \code{GDATASETS}. NULL means the current working directory db.
#' @param overwrite if TRUE, replace an existing destination track.
#'
#' @return invisibly, the character vector of created track names.
#' @seealso \code{\link{gtrack.mv}}, \code{\link{gtrack.rm}}
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
gtrack.copy <- function(src = NULL, dest = NULL, db = NULL, overwrite = FALSE) {
    if (is.null(substitute(src))) {
        stop("Usage: gtrack.copy(src, dest = NULL, db = NULL, overwrite = FALSE)", call. = FALSE)
    }
    .gcheckroot()

    # Resolve src names: support both unquoted single name (back-compat) and
    # character vectors.
    if (is.character(src) && length(src) > 1) {
        srcnames <- src
    } else if (is.character(src)) {
        srcnames <- src
    } else {
        srcnames <- do.call(.gexpr2str, list(substitute(src)), envir = parent.frame())
    }

    if (is.null(dest)) {
        destnames <- srcnames
    } else if (length(srcnames) == 1) {
        destnames <- if (is.character(dest)) {
            dest
        } else {
            do.call(.gexpr2str, list(substitute(dest)), envir = parent.frame())
        }
    } else {
        if (!is.character(dest) || length(dest) != 1) {
            stop("When copying multiple tracks, 'dest' must be a single namespace prefix or NULL.",
                call. = FALSE
            )
        }
        destnames <- paste(dest, srcnames, sep = ".")
    }

    dest_db <- .gtrack.copy.resolve_dest_db(db)
    .gcheck_write_permission(file.path(dest_db, "tracks"), "copy track to")

    created <- character(0)
    for (i in seq_along(srcnames)) {
        created <- c(created, .gtrack.copy.one(srcnames[i], destnames[i], dest_db, overwrite))
    }

    invisible(created)
}

# Resolve the destination db path. NULL -> the db that owns the current GWD.
# Otherwise validate that the explicit path is loaded.
.gtrack.copy.resolve_dest_db <- function(db) {
    if (is.null(db)) {
        gwd <- get("GWD", envir = .misha)
        groot <- get("GROOT", envir = .misha)
        gdatasets <- get("GDATASETS", envir = .misha)
        if (is.null(gdatasets)) gdatasets <- character(0)
        for (g in c(groot, gdatasets)) {
            if (.gpath_is_within(gwd, file.path(g, "tracks"))) {
                return(g)
            }
        }
        return(groot)
    }
    db <- normalizePath(db, mustWork = TRUE)
    groot <- get("GROOT", envir = .misha)
    gdatasets <- get("GDATASETS", envir = .misha)
    if (is.null(gdatasets)) gdatasets <- character(0)
    if (!(db %in% c(groot, gdatasets))) {
        stop(sprintf(
            "Destination db %s is not the current GROOT and not a loaded dataset; load it with gdataset.load() first.",
            db
        ), call. = FALSE)
    }
    db
}

# Copy a single track. Returns the destination track name on success.
.gtrack.copy.one <- function(srcname, destname, dest_db, overwrite) {
    if (!(srcname %in% get("GTRACKS", envir = .misha))) {
        stop(sprintf("Track %s does not exist", srcname), call. = FALSE)
    }

    src_db <- .gtrack_db_path(srcname)
    if (is.null(src_db)) src_db <- get("GROOT", envir = .misha)

    if (srcname == destname && identical(src_db, dest_db)) {
        stop(sprintf("Source and destination are the same track: %s", srcname), call. = FALSE)
    }

    # Check destination existence
    existing_db <- .gtrack_db_path(destname)
    if (!is.null(existing_db) && identical(existing_db, dest_db)) {
        if (!overwrite) {
            stop(sprintf(
                "Track %s already exists in %s; use overwrite=TRUE to replace.",
                destname, dest_db
            ), call. = FALSE)
        }
        gtrack.rm(destname, force = TRUE, db = dest_db)
    }

    src_dir <- .track_dir(srcname)
    dest_dir <- file.path(dest_db, "tracks", paste0(gsub("\\.", "/", destname), ".track"))
    dest_parent <- dirname(dest_dir)
    if (!dir.exists(dest_parent)) {
        dir.create(dest_parent, recursive = TRUE, showWarnings = FALSE)
    }

    src_indexed <- file.exists(file.path(src_dir, "track.idx"))
    dest_indexed <- .gdb.is_indexed_at(dest_db)
    src_chroms <- .gdb.chrom_names_at(src_db)
    dest_chroms <- .gdb.chrom_names_at(dest_db)

    info <- gtrack.info(srcname)

    # 2D track guard
    if (info$type %in% c("rectangles", "points") &&
        !identical(src_chroms, dest_chroms)) {
        stop(sprintf(
            "Cross-db copy of 2D track %s requires identical chromosome order in source and destination.",
            srcname
        ), call. = FALSE)
    }

    # Legacy single-file (1D) tracks: refuse if dest is indexed
    if (info$type %in% c("dense", "sparse", "array") && !dir.exists(src_dir)) {
        if (dest_indexed) {
            stop(sprintf(
                "Track %s is in legacy single-file format; convert with gtrack.convert(\"%s\") first.",
                srcname, srcname
            ), call. = FALSE)
        }
        if (!file.copy(src_dir, dest_dir, copy.mode = TRUE)) {
            stop(sprintf("Failed to copy %s to %s", srcname, destname), call. = FALSE)
        }
        .gdb.add_track(destname, dest_db)
        return(destname)
    }

    .gtrack.copy.pipeline(
        src_dir, dest_dir, src_chroms, dest_chroms,
        src_indexed, dest_indexed, info$type, destname, dest_db
    )
    .gdb.add_track(destname, dest_db)
    destname
}

# The full per-track copy pipeline:
#   copy dir -> [decode if src indexed] -> [drop unmapped chroms] -> [encode if dest indexed]
.gtrack.copy.pipeline <- function(src_dir, dest_dir,
                                  src_chroms, dest_chroms,
                                  src_indexed, dest_indexed,
                                  track_type, destname, dest_db) {
    same_order <- identical(src_chroms, dest_chroms)

    # Fast path: same format + same chrom order
    if (same_order && (src_indexed == dest_indexed)) {
        .gtrack.copy.raw_dir(src_dir, dest_dir)
        return(invisible())
    }

    # 2D track: at this point we know either format or chrom order differs.
    # If chrom orders differ, we already errored in .gtrack.copy.one.
    if (track_type %in% c("rectangles", "points")) {
        if (src_indexed && !dest_indexed) {
            stop(sprintf(
                "Cross-db copy of indexed 2D track %s into a per-chromosome database is not yet supported.",
                destname
            ), call. = FALSE)
        }
        .gtrack.copy.raw_dir(src_dir, dest_dir)
        if (dest_indexed && !src_indexed) {
            .with_db_context(dest_db, function() {
                gtrack.2d.convert_to_indexed(destname, remove.old = TRUE)
            })
        }
        return(invisible())
    }

    # 1D pipeline
    .gtrack.copy.raw_dir(src_dir, dest_dir)

    if (src_indexed) {
        .gtrack.split_indexed_to_per_chrom(dest_dir, src_chroms, remove_indexed = TRUE)
    }

    # Drop per-chrom files for chroms not in dest_chroms.
    # Tolerate the chr-prefix variant the indexed-conversion code already tolerates
    # (see src/GenomeTrackIndexedFormat.cpp:218-234).
    files_in_dir <- list.files(dest_dir, full.names = FALSE)
    dest_with_variants <- unique(c(
        dest_chroms,
        ifelse(startsWith(dest_chroms, "chr"), substr(dest_chroms, 4, nchar(dest_chroms)),
            paste0("chr", dest_chroms)
        )
    ))
    internal <- c("track.idx", "track.dat", ".attrs", ".vars", ".meta")
    candidates_for_drop <- setdiff(files_in_dir, internal)
    dropped <- candidates_for_drop[!(candidates_for_drop %in% dest_with_variants)]
    if (length(dropped) > 0) {
        warning(sprintf(
            "gtrack.copy(%s): dropped chromosomes not present in destination: %s",
            destname, paste(dropped, collapse = ", ")
        ), call. = FALSE)
        for (f in dropped) unlink(file.path(dest_dir, f))
    }

    if (dest_indexed) {
        .with_db_context(dest_db, function() gtrack.convert_to_indexed(destname))
    }
    invisible()
}

# Raw file copy of an entire .track directory.
.gtrack.copy.raw_dir <- function(src_dir, dest_dir) {
    if (!dir.create(dest_dir, showWarnings = FALSE) && !dir.exists(dest_dir)) {
        stop(sprintf("Failed to create %s", dest_dir), call. = FALSE)
    }
    contents <- list.files(src_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    for (item in contents) {
        if (!file.copy(item, dest_dir, recursive = TRUE, copy.mode = TRUE)) {
            unlink(dest_dir, recursive = TRUE)
            stop(sprintf("Failed to copy %s into %s", item, dest_dir), call. = FALSE)
        }
    }
}

# Temporarily switch GWD to the given db's tracks/ for the duration of fn().
# Mirrors .with_track_context (R/db-cache.R:297) but for an explicit dest db path.
.with_db_context <- function(dest_db, fn) {
    correct_gwd <- file.path(dest_db, "tracks")
    current_gwd <- get("GWD", envir = .misha)
    if (correct_gwd == current_gwd) {
        return(fn())
    }
    assign("GWD", correct_gwd, envir = .misha)
    on.exit(assign("GWD", current_gwd, envir = .misha), add = TRUE)
    fn()
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
        groot <- get("GROOT", envir = .misha)
        gdatasets <- get("GDATASETS", envir = .misha)
        if (is.null(gdatasets)) gdatasets <- character(0)
        groots <- c(groot, gdatasets)
        if (!(db %in% groots)) {
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
