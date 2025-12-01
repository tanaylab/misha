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
    .gcall("gtrackinfo", trackstr, validate, .misha_env())
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
#' @param ... these arguments are of either form 'pattern' or 'attribute =
#' pattern'
#' @param ignore.case,perl,fixed,useBytes see 'grep'
#' @return An array that contains the names of tracks that match the supplied
#' patterns.
#' @seealso \code{\link{grep}}, \code{\link{gtrack.exists}},
#' \code{\link{gtrack.create}}, \code{\link{gtrack.rm}}
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
gtrack.ls <- function(..., ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
    .gcheckroot()

    args <- as.list(substitute(list(...)))[-1L]
    args <- list(...)

    tracks <- get("GTRACKS", envir = .misha)

    if (is.null(tracks) || !length(tracks)) {
        return(NULL)
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
gtrack.rm <- function(track = NULL, force = FALSE) {
    if (is.null(substitute(track))) {
        stop("Usage: gtrack.rm(track, force = FALSE)", call. = FALSE)
    }
    .gcheckroot()

    trackname <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    dirname <- .track_dir(trackname)

    # check whether track appears among GTRACKS
    if (!(trackname %in% get("GTRACKS", envir = .misha))) {
        if (force) {
            .rm_track_dir(trackname)
            return(invisible())
        }
        stop(sprintf("Track %s does not exist", trackname), call. = FALSE)
    }

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
        .rm_track_dir(trackname)

        if (dir.exists(dirname)) {
            message(sprintf("Failed to delete track %s", trackname))
        } else {
            # refresh the list of GTRACKS, etc.
            .gdb.rm_track(trackname)
        }
    }
}
