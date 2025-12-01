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

    dir <- get("GWD", envir = .misha)

    res <- ""

    if (get("GWD", envir = .misha) != paste(get("GROOT", envir = .misha), "tracks", sep = "/")) {
        rescan <- TRUE
    }

    db.filename <- paste(get("GROOT", envir = .misha), ".db.cache", sep = "/")

    suppressWarnings({ # disable warnings since dir() on non dir or non existing dir produces warnings
        if (!rescan) {
            retv <- try(
                {
                    f <- file(db.filename, "rb")
                    res <- unserialize(f)
                    close(f)
                },
                silent = TRUE
            )

            if (inherits(retv, "try-error")) {
                rescan <- TRUE
            }
        }

        if (rescan) {
            res <- .gcall("gfind_tracks_n_intervals", dir, .misha_env())
            if (get("GWD", envir = .misha) == paste(get("GROOT", envir = .misha), "tracks", sep = "/")) {
                try(
                    {
                        f <- file(db.filename, "wb")
                        serialize(res, f)
                        close(f)
                        .gdb.cache_clear_dirty(get("GROOT", envir = .misha))
                    },
                    silent = TRUE
                )
            } else {
                unlink(db.filename, recursive = TRUE)
                .gdb.cache_clear_dirty(get("GROOT", envir = .misha))
            }
        }
    })

    tracks <- res[[1]]
    intervals <- res[[2]]

    tracks <- sort(tracks)
    intervals <- sort(intervals)

    res <- intersect(tracks, intervals)
    if (length(res) > 0) {
        stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))
    }

    assign("GTRACKS", tracks, envir = .misha)
    assign("GINTERVS", intervals, envir = .misha)
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
    if (!is.na(match(track, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Track %s already exists", track), call. = FALSE)
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

.track_dir <- function(trackname) {
    sprintf("%s.track", paste(get("GWD", envir = .misha), gsub("\\.", "/", trackname), sep = "/"))
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

.gdb.cache_update_lists <- function() {
    if (!exists("GROOT", envir = .misha, inherits = FALSE) ||
        !exists("GTRACKS", envir = .misha, inherits = FALSE) ||
        !exists("GINTERVS", envir = .misha, inherits = FALSE)) {
        return(invisible(FALSE))
    }

    groot <- get("GROOT", envir = .misha)
    if (is.null(groot) || groot == "") {
        return(invisible(FALSE))
    }

    tracks <- get("GTRACKS", envir = .misha)
    intervals <- get("GINTERVS", envir = .misha)
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
