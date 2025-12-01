# Track conversion functions

#' Converts a track to the most current format
#'
#' Converts a track (if needed) to the most current format.
#'
#' This function converts a track to the most current format. It should be used
#' if a track created by an old version of the library cannot be read anymore
#' by the newer version. The old track is given by 'src.track'. After
#' conversion a new track 'tgt.track' is created. If 'tgt.track' is 'NULL' the
#' source track is overwritten.
#'
#' @param src.track source track name
#' @param tgt.track target track name. If 'NULL' the source track is
#' overwritten.
#' @return None
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.create_sparse}}
#' @keywords ~track ~convert
#' @export gtrack.convert
gtrack.convert <- function(src.track = NULL, tgt.track = NULL) {
    if (is.null(substitute(src.track))) {
        stop("Usage: gtrack.convert(src.track, tgt.track = NULL)", call. = FALSE)
    }
    .gcheckroot()

    src.trackstr <- do.call(.gexpr2str, list(substitute(src.track)), envir = parent.frame())
    if (is.na(match(src.trackstr, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Track %s does not exist", src.trackstr), call. = FALSE)
    }

    tgt.trackstr <- ""
    if (is.null(substitute(tgt.track))) {
        tgt.trackstr <- paste(src.trackstr, "_converted", sep = "")
        counter <- 2
        while (!is.na(match(tgt.trackstr, get("GTRACKS", envir = .misha)))) {
            tgt.trackstr <- paste(src.trackstr, "_converted", counter, sep = "")
            counter <- counter + 1
        }
    } else {
        tgt.trackstr <- do.call(.gexpr2str, list(substitute(tgt.track)), envir = parent.frame())
        .gconfirmtrackcreate(tgt.trackstr)
    }

    src.dirname <- sprintf("%s.track", paste(get("GWD", envir = .misha), gsub("\\.", "/", src.trackstr), sep = "/"))
    tgt.dirname <- sprintf("%s.track", paste(get("GWD", envir = .misha), gsub("\\.", "/", tgt.trackstr), sep = "/"))
    success <- FALSE

    tryCatch(
        {
            .gcall("gtrackconvert", src.trackstr, tgt.trackstr, .misha_env())

            # copy all supplimentary data of a track (vars, etc.)
            if (!system(sprintf("cp -r -u %s/. %s", src.dirname, tgt.dirname))) {
                # if tgt track is null move it to the source track
                if (is.null(substitute(tgt.track))) {
                    unlink(src.dirname, recursive = TRUE)
                    success <- TRUE
                    file.rename(tgt.dirname, src.dirname)
                }
            } else {
                msg <- sprintf("Failed to copy some or all track supplementary data from %s to %s", src.dirname, tgt.dirname)
                if (is.null(substitute(tgt.track))) {
                    msg <- paste(msg,
                        sprintf(
                            "Track %s will remain unchanged.\nA new converted track named %s was created without supplementary data.",
                            src.trackstr, tgt.trackstr
                        ),
                        sep = "\n"
                    )
                }
                warning(msg, call. = FALSE)
            }
            success <- TRUE

            # If database is indexed, automatically convert the track to indexed format
            if (.gdb.is_indexed()) {
                gtrack.convert_to_indexed(tgt.trackstr)
            }
        },
        finally = {
            if (!success) {
                unlink(tgt.dirname, recursive = TRUE)
            }
            .gdb.rm_track(tgt.trackstr)
        }
    )
    retv <- 0 # suppress return value
}
