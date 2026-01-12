# Track modification functions (lookup, modify, smooth)

#' Creates a new track from a lookup table based on track expression
#'
#' Evaluates track expression and translates the values into bin indices that
#' are used in turn to retrieve values from a lookup table and create a track.
#'
#' This function evaluates the track expression for all iterator intervals and
#' translates this value into an index based on the breaks. This index is then
#' used to address the lookup table and create with its values a new track.
#' More than one 'expr'-'breaks' pair can be used. In that case 'lookup_table'
#' is addressed in a multidimensional manner, i.e. 'lookup_table[i1, i2, ...]'.
#'
#' The range of bins is determined by 'breaks' argument. For example: 'breaks =
#' c(x1, x2, x3, x4)' represents three different intervals (bins): (x1, x2],
#' (x2, x3], (x3, x4].
#'
#' If 'include.lowest' is 'TRUE' the the lowest value is included in the first
#' interval, i.e. in [x1, x2].
#'
#' 'force.binning' parameter controls what should be done when the value of
#' 'expr' exceeds the range determined by 'breaks'. If 'force.binning' is
#' 'TRUE' then values smaller than the minimal break will be translated to
#' index 1, and the values exceeding the maximal break will be translated to
#' index 'M-1' where 'M' is the number of breaks. If 'force.binning' is 'FALSE'
#' the out-of-range values will produce 'NaN' values.
#'
#' Regardless of 'force.binning' value if the value of 'expr' is 'NaN' then the
#' value in the track would be 'NaN' too.
#'
#' 'description' is added as a track attribute.
#'
#' @param track track name
#' @param description a character string description
#' @param lookup_table a multi-dimensional array containing the values that are
#' returned by the function
#' @param ... pairs of track expressions and breaks
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param force.binning if 'TRUE', the values smaller than the minimal break
#' will be translated to index 1, and the values that exceed the maximal break
#' will be translated to index N-1 where N is the number of breaks. If 'FALSE'
#' the out-of-range values will produce NaN values.
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @return None.
#' @seealso \code{\link{glookup}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.smooth}},
#' \code{\link{gtrack.modify}}, \code{\link{gtrack.rm}},
#' \code{\link{gtrack.info}}, \code{\link{gdir.create}}
#' @keywords ~lookup ~track
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' ## one-dimensional example
#' breaks1 <- seq(0.1, 0.2, length.out = 6)
#' gtrack.lookup(
#'     "lookup_track", "Test track", 1:5, "dense_track",
#'     breaks1
#' )
#' gtrack.rm("lookup_track", force = TRUE)
#'
#' ## two-dimensional example
#' t <- array(1:15, dim = c(5, 3))
#' breaks2 <- seq(0.31, 0.37, length.out = 4)
#' gtrack.lookup(
#'     "lookup_track", "Test track", t, "dense_track",
#'     breaks1, "2 * dense_track", breaks2
#' )
#' gtrack.rm("lookup_track", force = TRUE)
#'
#' @export gtrack.lookup
gtrack.lookup <- function(track = NULL, description = NULL, lookup_table = NULL, ..., include.lowest = FALSE, force.binning = TRUE, iterator = NULL, band = NULL) {
    args <- as.list(substitute(list(...)))[-1L]
    if (is.null(substitute(track)) || is.null(description) || is.null(lookup_table) || length(args) < 2 || length(args) %% 2 != 0) {
        stop("Usage: gtrack.lookup(track, description, lookup_table, [expr, breaks]+, include.lowest = FALSE, force.binning = TRUE, iterator = NULL, band = NULL)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())

    exprs <- c()
    breaks <- list()

    for (i in (0:(length(args) / 2 - 1))) {
        exprs <- append(exprs, do.call(.gexpr2str, list(args[[i * 2 + 1]]), envir = parent.frame()))
        breaks[[length(breaks) + 1]] <- eval.parent(args[[i * 2 + 2]])
    }

    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    .gconfirmtrackcreate(trackstr)
    trackdir <- .track_dir(trackstr)
    direxisted <- file.exists(trackdir)
    success <- FALSE
    tryCatch(
        {
            .gcall("gtrack_bintransform", trackstr, exprs, breaks, include.lowest, force.binning, lookup_table, .iterator, band, .misha_env())
            .gdb.add_track(trackstr)
            created.by <- sprintf("gtrack.lookup(%s, description, lookup_table", trackstr)
            for (i in (1:length(exprs))) {
                created.by <- sprintf("%s, %s, breaks%d", created.by, exprs[i], i)
            }
            created.by <- sprintf("%s, include.lowest = %s, force.binning = %s)", created.by, include.lowest, force.binning)
            .gtrack.attr.set(trackstr, "created.by", created.by, TRUE)
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

#' Modifies track contents
#'
#' Modifies 'Dense' track contents.
#'
#' This function modifies the contents of a 'Dense' track by the values of
#' 'expr'. 'intervals' argument controls which portion of the track is
#' modified. The iterator policy is set internally to the bin size of the
#' track.
#'
#' @param track track name
#' @param expr track expression
#' @param intervals genomic scope for which track is modified
#' @return None.
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.rm}}
#' @keywords ~modify ~track
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' intervs <- gintervals(1, 300, 800)
#' gextract("dense_track", intervs)
#' gtrack.modify("dense_track", "dense_track * 2", intervs)
#' gextract("dense_track", intervs)
#' gtrack.modify("dense_track", "dense_track / 2", intervs)
#'
#' @export gtrack.modify
gtrack.modify <- function(track = NULL, expr = NULL, intervals = NULL) {
    if (is.null(substitute(track)) || is.null(substitute(expr))) {
        stop("Usage: gtrack.modify(track, expr, intervals = .misha$ALLGENOME)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (is.null(intervals)) {
        intervals <- get("ALLGENOME", envir = .misha)
    }

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())

    .gcall("gtrack_modify", trackstr, exprstr, intervals, iterator = trackstr, .misha_env())

    str <- sprintf("gtrack.modify(%s, %s, intervs)", trackstr, exprstr)
    created.by.str <- gtrack.attr.export(trackstr, "created.by")[1, 1]
    if (is.null(created.by.str)) {
        created.by.str <- str
    } else {
        created.by.str <- paste(created.by.str, str, sep = "\n")
    }

    .gtrack.attr.set(trackstr, "created.by", created.by.str, TRUE)

    retv <- 0 # suppress return value
}

#' Creates a new track from smoothed values of track expression
#'
#' Creates a new track from smoothed values of track expression.
#'
#' This function creates a new 'Dense' track named 'track'. The values of the
#' track are results of smoothing the values of 'expr'.
#'
#' Each track value at coordinate 'C' is determined by smoothing non 'NaN'
#' values of 'expr' over the window around 'C'. The window size is controlled
#' by 'winsize' and is given in coordinate units (not in number of bins),
#' defining the total regions to be considered when smoothing (on both sides of
#' the central point). Two different algorithms can be used for smoothing:
#'
#' "MEAN" - an arithmetic average.
#'
#' "LINEAR_RAMP" - a weighted arithmetic average, where the weights linearly
#' decrease as the distance from the center of the window increases.
#'
#' 'weight_thr' determines the function behavior when some of the values in the
#' window are missing or 'NaN' (missing values may occur at the edges of each
#' chromosome when the window covers an area beyond chromosome boundaries).
#' 'weight_thr' sets the weight sum threshold below which smoothing algorithm
#' returns 'NaN' rather than a smoothing value based on non 'NaN' values in the
#' window.
#'
#' 'smooth_nans' controls what would be the smoothed value if the central value
#' in the window is 'NaN'. If 'smooth_nans' is 'FALSE' then the smoothed value
#' is set to 'NaN' regardless of 'weight_thr' parameter. Otherwise it is
#' calculated normally.
#'
#' 'description' is added as a track attribute.
#'
#' Iterator policy must be of "fixed bin" type.
#'
#' @param track track name
#' @param description a character string description
#' @param expr track expression
#' @param winsize size of smoothing window
#' @param weight_thr smoothing weight threshold
#' @param smooth_nans if 'FALSE' track value is always set to 'NaN' if central
#' window value is 'NaN', otherwise it is calculated from the rest of non 'NaN'
#' values
#' @param alg smoothing algorithm - "MEAN" or "LINEAR_RAMP"
#' @param iterator track expression iterator of 'Fixed bin' type
#' @return None.
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.modify}},
#' \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}
#' @keywords ~smooth ~track
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gtrack.smooth("smoothed_track", "Test track", "dense_track", 500)
#' gextract("dense_track", "smoothed_track", gintervals(1, 0, 1000))
#' gtrack.rm("smoothed_track", force = TRUE)
#'
#' @export gtrack.smooth
gtrack.smooth <- function(track = NULL, description = NULL, expr = NULL, winsize = NULL, weight_thr = 0, smooth_nans = FALSE, alg = "LINEAR_RAMP", iterator = NULL) {
    if (is.null(substitute(track)) || is.null(description) || is.null(substitute(expr)) || is.null(winsize)) {
        stop("Usage: gtrack.smooth(track, description, expr, winsize, weight_thr = 0, smooth_nans = FALSE, alg = \"LINEAR_RAMP\" (\"LINEAR_RAMP\" | \"MEAN\"), iterator = NULL)",
            call. = FALSE
        )
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    .gconfirmtrackcreate(trackstr)
    trackdir <- .track_dir(trackstr)
    direxisted <- file.exists(trackdir)
    success <- FALSE
    tryCatch(
        {
            .gcall("gsmooth", trackstr, exprstr, winsize, weight_thr, smooth_nans, alg, .iterator, .misha_env())
            .gdb.add_track(trackstr)
            .gtrack.attr.set(
                trackstr, "created.by",
                sprintf("gtrack.smooth(%s, description, %s, %s, %s, %s, %s)", trackstr, exprstr, as.character(winsize), as.character(weight_thr), as.character(smooth_nans), alg), TRUE
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
