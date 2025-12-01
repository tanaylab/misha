# Track creation functions (create, create_pwm_energy, create_sparse)

#' Creates a track from a track expression
#'
#' Creates a track from a track expression.
#'
#' This function creates a new track named track. The values of the track are
#' determined by evaluation of 'expr' - a numeric track expression. The type of
#' the new track is determined by the type of the iterator. 'Fixed bin',
#' 'Sparse' or 'Rectangles' track can be created accordingly. 'description' is
#' added as a track attribute.
#'
#' @param track track name
#' @param description a character string description
#' @param expr track expression
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @return None.
#' @seealso \code{\link{gtrack.2d.create}}, \code{\link{gtrack.create_sparse}},
#' \code{\link{gtrack.smooth}}, \code{\link{gtrack.modify}},
#' \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}
#' @keywords ~create ~track
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' ## Creates a new track that is a sum of values from 'dense' and
#' ## 2 * non-nan values of 'sparse' track. The new track type is
#' ## Dense with a bin size that equals to '70'.
#' gtrack.create("mixed_track", "Test track",
#'     "dense_track +
#'               replace(sparse_track, is.nan(sparse_track), 0) * 2",
#'     iterator = 70
#' )
#' gtrack.info("mixed_track")
#' gtrack.rm("mixed_track", force = TRUE)
#'
#' @export gtrack.create
gtrack.create <- function(track = NULL, description = NULL, expr = NULL, iterator = NULL, band = NULL) {
    if (is.null(substitute(track)) || is.null(description) || is.null(substitute(expr))) {
        stop("Usage: gtrack.create(track, description, expr, iterator = NULL, band = NULL)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
    trackdir <- .track_dir(trackstr)

    direxisted <- file.exists(trackdir)

    if (!is.na(match(trackstr, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Track %s already exists", trackstr), call. = FALSE)
    }

    .gconfirmtrackcreate(trackstr)
    success <- FALSE
    tryCatch(
        {
            if (.ggetOption("gmultitasking")) {
                .gcall("gtrackcreate_multitask", trackstr, exprstr, .iterator, band, .misha_env())
            } else {
                .gcall("gtrackcreate", trackstr, exprstr, .iterator, band, .misha_env())
            }
            .gdb.add_track(trackstr)
            .gtrack.attr.set(
                trackstr, "created.by",
                sprintf("gtrack.create(%s, description, %s, iterator=%s)", trackstr, exprstr, deparse(substitute(iterator), width.cutoff = 500)[1]), TRUE
            )
            .gtrack.attr.set(trackstr, "created.date", date(), TRUE)
            .gtrack.attr.set(trackstr, "created.user", Sys.getenv("USER"), TRUE)
            .gtrack.attr.set(trackstr, "description", description, TRUE)
            success <- TRUE

            # If database is indexed, automatically convert the track to indexed format
            # Only convert 1D tracks (dense, sparse, array) - 2D tracks cannot be converted
            if (.gdb.is_indexed()) {
                track_info <- gtrack.info(trackstr)
                if (track_info$type %in% c("dense", "sparse", "array")) {
                    gtrack.convert_to_indexed(trackstr)
                }
            }
        },
        finally = {
            if (!success && !direxisted) {
                unlink(trackdir, recursive = TRUE)
                .gdb.rm_track(trackstr)
            }
        }
    )
    retv <- 0 # suppress return value
}


#' Creates a new track from PSSM energy function
#'
#' Creates a new track from PSSM energy function.
#'
#' This function creates a new track with values of a PSSM energy function.
#' PSSM parameters (nucleotide probability per position and pluralization) are
#' determined by 'pssmset' key and data files ('pssmset.key' and
#' 'pssmset.data'). These two files must be located in 'GROOT/pssms' directory.
#' The type of the created track is determined by the type of the iterator.
#' 'description' is added as a track attribute.
#'
#' @param track track name
#' @param description a character string description
#' @param pssmset name of PSSM set: 'pssmset.key' and 'pssmset.data' must be
#' presented in 'GROOT/pssms' directory
#' @param pssmid PSSM id
#' @param prior prior
#' @param iterator track expression iterator for the newly created track
#' @return None.
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.smooth}},
#' \code{\link{gtrack.modify}}, \code{\link{gtrack.rm}},
#' \code{\link{gtrack.info}}, \code{\link{gdir.create}}
#' @keywords ~energy ~pssm ~pwm ~track
#' @examples
#' \donttest{
#' gdb.init_examples()
#' gtrack.create_pwm_energy("pwm_energy_track", "Test track", "pssm",
#'     3, 0.01,
#'     iterator = 100
#' )
#' gextract("pwm_energy_track", gintervals(1, 0, 1000))
#' }
#'
#' @export gtrack.create_pwm_energy
gtrack.create_pwm_energy <- function(track = NULL, description = NULL, pssmset = NULL, pssmid = NULL, prior = NULL, iterator = NULL) {
    if (is.null(substitute(track)) || is.null(description) || is.null(pssmset) || is.null(pssmid) || is.null(prior) || is.null(iterator)) {
        stop("Usage: gtrack.create_pwm_energy(track, description, pssmset, pssmid, prior, iterator)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    trackdir <- .track_dir(trackstr)
    direxisted <- file.exists(trackdir)
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    if (!is.na(match(trackstr, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Track %s already exists", trackstr), call. = FALSE)
    }

    .gconfirmtrackcreate(trackstr)
    success <- FALSE
    tryCatch(
        {
            if (.ggetOption("gmultitasking")) {
                .gcall("gcreate_pwm_energy_multitask", trackstr, pssmset, pssmid, prior, .iterator, .misha_env())
            } else {
                .gcall("gcreate_pwm_energy", trackstr, pssmset, pssmid, prior, .iterator, .misha_env())
            }
            .gdb.add_track(trackstr)
            .gtrack.attr.set(
                trackstr, "created.by",
                sprintf(
                    "gtrack.create_pwm_energy(%s, description, \"%s\", %g, %g, iterator=%s)",
                    trackstr, pssmset, as.numeric(pssmid), as.numeric(prior), deparse(substitute(iterator), width.cutoff = 500)[1]
                ), TRUE
            )
            .gtrack.attr.set(trackstr, "created.date", date(), TRUE)
            .gtrack.attr.set(trackstr, "created.user", Sys.getenv("USER"), TRUE)
            .gtrack.attr.set(trackstr, "description", description, TRUE)
            success <- TRUE

            # If database is indexed, automatically convert the track to indexed format
            # Only convert 1D tracks (dense, sparse, array) - 2D tracks cannot be converted
            if (.gdb.is_indexed()) {
                track_info <- gtrack.info(trackstr)
                if (track_info$type %in% c("dense", "sparse", "array")) {
                    gtrack.convert_to_indexed(trackstr)
                }
            }
        },
        finally = {
            if (!success && !direxisted) {
                unlink(trackdir, recursive = TRUE)
                .gdb.rm_track(trackstr)
            }
        }
    )
    retv <- 0 # suppress return value
}


#' Creates a 'Sparse' track from intervals and values
#'
#' Creates a 'Sparse' track from intervals and values.
#'
#' This function creates a new 'Sparse' track with values at given intervals.
#' 'description' is added as a track attribute.
#'
#' @param track track name
#' @param description a character string description
#' @param intervals a set of one-dimensional intervals
#' @param values an array of numeric values - one for each interval
#' @return None.
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.smooth}}, \code{\link{gtrack.modify}},
#' \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}
#' @keywords ~create ~sparse ~track
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' intervs <- gintervals.load("annotations")
#' gtrack.create_sparse(
#'     "test_sparse", "Test track", intervs,
#'     1:dim(intervs)[1]
#' )
#' gextract("test_sparse", .misha$ALLGENOME)
#' gtrack.rm("test_sparse", force = TRUE)
#'
#' @export gtrack.create_sparse
gtrack.create_sparse <- function(track = NULL, description = NULL, intervals = NULL, values = NULL) {
    if (is.null(substitute(track)) || is.null(description) || is.null(intervals) || is.null(values)) {
        stop("Usage: gtrack.create_sparse(track, description, intervals, values)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    trackdir <- .track_dir(trackstr)

    direxisted <- file.exists(trackdir)

    if (!is.na(match(trackstr, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Track %s already exists", trackstr), call. = FALSE)
    }

    .gconfirmtrackcreate(trackstr)
    success <- FALSE
    tryCatch(
        {
            .gcall("gtrack_create_sparse", trackstr, intervals, values, .misha_env())
            .gdb.add_track(trackstr)
            .gtrack.attr.set(trackstr, "created.by", sprintf("gtrack.create_sparse(%s, description, intervals, values)", trackstr), TRUE)
            .gtrack.attr.set(trackstr, "created.date", date(), TRUE)
            .gtrack.attr.set(trackstr, "created.user", Sys.getenv("USER"), TRUE)
            .gtrack.attr.set(trackstr, "description", description, TRUE)
            success <- TRUE

            # If database is indexed, automatically convert the track to indexed format
            # Only convert 1D tracks (dense, sparse, array) - 2D tracks cannot be converted
            if (.gdb.is_indexed()) {
                track_info <- gtrack.info(trackstr)
                if (track_info$type %in% c("dense", "sparse", "array")) {
                    gtrack.convert_to_indexed(trackstr)
                }
            }
        },
        finally = {
            if (!success && !direxisted) {
                unlink(trackdir, recursive = TRUE)
                .gdb.rm_track(trackstr)
            }
        }
    )
    retv <- 0 # suppress return value
}
# Dense track creation function

#' Creates a 'Dense' track from intervals and values
#'
#' Creates a 'Dense' track from intervals and values.
#'
#' This function creates a new 'Dense' track with values at given intervals.
#' 'description' is added as a track attribute.
#'
#' @param track track name
#' @param description a character string description
#' @param intervals a set of one-dimensional intervals
#' @param values an array of numeric values - one for each interval
#' @param binsize bin size of the newly created 'Dense' track
#' @param defval default track value for genomic regions not covered by the intervals
#' @return None.
#' @seealso \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.import}},
#' \code{\link{gtrack.modify}}, \code{\link{gtrack.rm}},
#' \code{\link{gtrack.info}}
#' @keywords ~create ~dense ~track
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' intervs <- gintervals.load("annotations")
#' gtrack.create_dense(
#'     "test_dense", "Test dense track", intervs,
#'     1:dim(intervs)[1], 50, 0
#' )
#' gextract("test_dense", .misha$ALLGENOME)
#' gtrack.rm("test_dense", force = TRUE)
#'
#' @export gtrack.create_dense
gtrack.create_dense <- function(track = NULL, description = NULL, intervals = NULL, values = NULL,
                                binsize = NULL, defval = NaN) {
    if (is.null(substitute(track)) || is.null(description) || is.null(intervals) || is.null(values) || is.null(binsize)) {
        stop("Usage: gtrack.create_dense(track, description, intervals, values, binsize, defval = NaN)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    trackdir <- .track_dir(trackstr)

    direxisted <- file.exists(trackdir)

    if (!is.na(match(trackstr, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Track %s already exists", trackstr), call. = FALSE)
    }

    # Create a data frame from intervals and values
    if (length(values) != nrow(intervals)) {
        stop("Length of values must match the number of intervals", call. = FALSE)
    }

    intervalData <- data.frame(
        chrom = as.character(intervals$chrom),
        start = as.numeric(intervals$start),
        end = as.numeric(intervals$end),
        value = as.numeric(values)
    )

    .gconfirmtrackcreate(trackstr)
    success <- FALSE

    tryCatch(
        {
            # Call the C++ function with the data frame
            .gcall("gtrack_create_dense", trackstr, intervalData, binsize, defval, .misha_env())

            # Add the track to the database
            .gdb.add_track(trackstr)

            # Set track attributes
            .gtrack.attr.set(
                trackstr, "created.by",
                sprintf("gtrack.create_dense(%s, description, intervals, values, %d, %g)", trackstr, binsize, defval), TRUE
            )
            .gtrack.attr.set(trackstr, "created.date", date(), TRUE)
            .gtrack.attr.set(trackstr, "created.user", Sys.getenv("USER"), TRUE)
            .gtrack.attr.set(trackstr, "description", description, TRUE)

            # Set the track type to "dense"
            .gtrack.attr.set(trackstr, "type", "dense", TRUE)

            # Set the binsize attribute which is required for the tests
            .gtrack.attr.set(trackstr, "binsize", binsize, TRUE)

            success <- TRUE

            # If database is indexed, automatically convert the track to indexed format
            # Only convert 1D tracks (dense, sparse, array) - 2D tracks cannot be converted
            if (.gdb.is_indexed()) {
                track_info <- gtrack.info(trackstr)
                if (track_info$type %in% c("dense", "sparse", "array")) {
                    gtrack.convert_to_indexed(trackstr)
                }
            }
        },
        finally = {
            if (!success && !direxisted) {
                unlink(trackdir, recursive = TRUE)
                .gdb.rm_track(trackstr)
            }
        }
    )

    retv <- 0 # suppress return value
}
