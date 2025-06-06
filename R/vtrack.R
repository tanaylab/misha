.gvtrack <- function(vtrack) {
    vtrackstr <- do.call(.gexpr2str, list(substitute(vtrack)), envir = parent.frame())
    if (!is.character(vtrackstr) || length(vtrackstr) != 1) {
        stop(sprintf("Virtual track must be specified as a character string"), call. = FALSE)
    }

    if (is.na(match(vtrackstr, gvtrack.ls()))) {
        stop(sprintf("Virtual track %s does not exist", vtrackstr), call. = FALSE)
    }
    vtrackstr
}

.gvtrack.get <- function(vtrackstr) {
    if (is.na(match(vtrackstr, gvtrack.ls()))) {
        stop(sprintf("Virtual track %s does not exist", vtrackstr), call. = FALSE)
    }

    gwd <- get("GWD", envir = .misha)
    get("GVTRACKS", envir = .misha)[[gwd]][[vtrackstr]]
}

.gvtrack.set <- function(vtrackstr, var) {
    if (exists("GVTRACKS", envir = .misha)) {
        gvtracks <- get("GVTRACKS", envir = .misha)
    } else {
        gvtracks <- list()
    }

    gwds <- names(gvtracks)
    if (!is.list(gvtracks) || (length(gvtracks) && !is.character(gwds)) || length(gvtracks) != length(gwds)) {
        stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = FALSE)
    }

    gwd <- get("GWD", envir = .misha)
    idx1 <- match(gwd, gwds)
    if (is.na(idx1)) {
        gwds <- c(gwds, gwd)
        idx1 <- length(gwds)
        gvtracks[[idx1]] <- list()
    }

    vtracks <- gvtracks[[idx1]]
    names(gvtracks) <- gwds

    vtracknames <- names(vtracks)
    if (!is.list(vtracks) || (length(vtracks) && !is.character(vtracknames)) || length(vtracks) != length(vtracknames)) {
        stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = FALSE)
    }

    idx2 <- match(vtrackstr, vtracknames)
    if (is.na(idx2)) {
        if (!is.na(match(vtrackstr, get("GTRACKS", envir = .misha)))) {
            stop(sprintf("Track %s already exists", vtrackstr), call. = FALSE)
        }

        if (!is.na(match(vtrackstr, get("GINTERVS", envir = .misha)))) {
            stop(sprintf("Interval %s already exists", vtrackstr), call. = FALSE)
        }

        if (.ggetOption(".gautocompletion", FALSE) && exists(vtrackstr, envir = .misha)) {
            stop(sprintf("Variable \"%s\" shadows the name of identically named virtual track.\nPlease remove this variable from the environment or switch off autocompletion mode.", vtrackstr), call. = FALSE)
        }

        vtracknames <- c(vtracknames, vtrackstr)
        idx2 <- length(vtracknames)
    }

    gvtracks[[idx1]][[idx2]] <- var

    if (!is.null(var)) {
        names(gvtracks[[idx1]]) <- vtracknames

        envir <- .misha_env()
        assign("GVTRACKS", gvtracks, envir$.misha)
        .gcall("gcheck_vtrack", vtrackstr, envir)
    }

    success <- FALSE
    old.gvtracks <- NULL
    if (exists("GVTRACKS", envir = .misha)) {
        old.gvtracks <- get("GVTRACKS", envir = .misha)
    }

    success <- FALSE
    tryCatch({
        assign("GVTRACKS", gvtracks, envir = .misha)
        success <- TRUE
    })
}



#' Creates a new virtual track
#'
#' Creates a new virtual track.
#'
#' This function creates a new virtual track named 'vtrack' with the given
#' source, function and parameters. 'src' can be either a track or intervals
#' (1D or 2D). Use the following table for a reference of all valid source,
#' function and parameters combinations:
#'
#' \emph{src = [Track], func = "avg", params = NULL} \cr Average track value in
#' iterator interval.
#'
#' \emph{src = [Track], func = "max", params = NULL} \cr Maximal track value in
#' iterator interval.
#'
#' \emph{src = [Track], func = "min", params = NULL} \cr Minimal track value in
#' iterator interval.
#'
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "nearest", params =
#' NULL} \cr Mean track value in iterator interval. If there are no track
#' values covered by an iterator interator (can occur only in 'Sparse' track),
#' the nearest track value is returned.
#'
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "stddev", params =
#' NULL} \cr Unbiased standard deviation of track values in iterator interval.
#'
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "sum", params =
#' NULL} \cr Sum of track values in iterator interval.
#'
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "quantile", params
#' = [Percentile in the range of [0, 1]]} \cr Quantile of track values in
#' iterator interval.
#'
#' \emph{src = ['Dense' track], func = "global.percentile", params = NULL} \cr
#' Percentile of an average track value in iterator interval relatively to all
#' values of the track.
#'
#' \emph{src = ['Dense' track], func = "global.percentile.max", params = NULL}
#' \cr Percentile of a maximal track value in iterator interval relatively to
#' all values of the track.
#'
#' \emph{src = ['Dense' track], func = "global.percentile.min", params = NULL}
#' \cr Percentile of a minimal track value in iterator interval relatively to
#' all values of the track.
#'
#' \emph{src = [2D track], func = "area", params = NULL} \cr Area covered by
#' iterator interval.
#'
#' \emph{src = [2D track], func = "weighted.sum", params = NULL} \cr Weighted
#' sum of values where each weight equals to the intersection area between the
#' iterator interval and the rectangle containing the value.
#'
#' \emph{src = [1D intervals], func = "distance", params = [Minimal distance
#' from center (default: 0)]} \cr Given the center 'C' of the current iterator
#' interval returns 'DC * X/2', where 'DC' is the normalized distance to the
#' center of the interval that contains 'C', and 'X' is the value of the
#' parameter. If no interval contains 'C' the resulted value is 'D + XXX/2'
#' where 'D' is the distance between 'C' and the edge of the closest interval.
#' Distance can be positive or negative depending on the position of the
#' coordinate relative to the interval and the strand (-1 or 1) of the
#' interval. Distance is always positive if 'strand' is '0' or if 'strand'
#' column is missing. Distance is 'NA' if no intervals exist for the current
#' chromosome.
#'
#' \emph{src = [1D intervals], func = "distance.center", params = NULL} \cr
#' Given the center 'C' of the current iterator interval returns 'NaN' if 'C'
#' is outside of the intervals, otherwise returns the distance between 'C' and
#' the center of the closest interval. Distance can be positive or negative
#' depending on the position of the coordinate relative to the interval and the
#' strand (-1 or 1) of the interval. Distance is always positive if 'strand' is
#' '0' or if 'strand' column is missing.
#'
#' \emph{src = [1D intervals], func = "coverage", params = NULL} \cr
#' For each iterator interval, calculates the fraction of its length that is covered by the
#' source intervals. Returns a value between 0 and 1. For example, if an iterator interval is [100,200]
#' and the source intervals cover positions 120-140 and 160-170, the coverage would be 0.3
#' ((20 + 10) / 100 = 0.3). Overlapping source intervals are first unified.
#'
#' \emph{func = "pwm", params = list(pssm = matrix, bidirect = TRUE,
#' prior = 0.01, extend = TRUE)} \cr
#' Calculates total log-likelihood score of DNA sequence against PSSM.
#' Uses log-sum-exp over all positions. For bidirect=TRUE, scans both
#' strands. Prior adds pseudocounts, extend=TRUE allows scoring at boundaries.
#'
#' \emph{func = "pwm.max", params = list(pssm = matrix, bidirect = TRUE,
#' prior = 0.01, extend = TRUE)} \cr
#' Returns maximum log-likelihood score of best PSSM match. bidirect=TRUE
#' checks both strands. Prior adds pseudocounts, extend=TRUE allows boundary
#' scoring.
#'
#' \emph{func = "pwm.max.pos", params = list(pssm = matrix, bidirect = TRUE,
#' prior = 0.01, extend = TRUE)} \cr
#' Returns 1-based position of best PSSM match.
#' If bidirect=TRUE, the position would be positive if the best hit was at the
#' forward strand, and negative if it was at the reverse strand. When strand is
#' -1 the position is still according to the forward strand, but the hit is at
#' the end of the match.
#' Prior adds pseudocounts, extend=TRUE allows boundary scoring.
#'
#' For all PWM functions:
#' \itemize{
#'   \item pssm: Position-specific scoring matrix (A,C,G,T frequencies)
#'   \item bidirect: If TRUE, scans both strands; if FALSE, forward only
#'   \item prior: Pseudocount for frequencies (default: 0.01)
#'   \item extend: If TRUE, computes boundary scores
#'   \item strand: If 1, scans forward strand; if -1, scans reverse strand.
#' For strand == 1, the energy (and position of the best match) would be at
#' the beginning of the match, for strand == -1, the energy (and position of
#' the best match) would be at the end of the match.
#' }
#'
#' PWM parameters are accepted as list or individual parameters (see examples).
#'
#' \emph{func = "kmer.count", params = list(kmer = "ACGT", extend = TRUE, strand = 0)} \cr
#' Counts occurrences of the specified kmer in each interval.
#'
#' \emph{func = "kmer.frac", params = list(kmer = "ACGT", extend = TRUE, strand = 0)} \cr
#' Calculates the fraction of possible positions in each interval that contain the specified kmer.
#'
#' For kmer functions:
#' \itemize{
#'   \item kmer: The DNA sequence to count (case-insensitive)
#'   \item extend: If TRUE (default), considers k-mers starting at any position within the interval, which may require looking at sequence data beyond the interval's end. If FALSE, the sequence is not extended, so k-mers starting near the end of the interval that would cross the boundary are not counted.
#'   \item strand: If 1, counts kmers on forward strand; if -1, counts kmers on reverse strand. If
#'  0, counts kmers on both strands. Default is 0.
#' }
#'
#' Kmer parameters are accepted as list or individual parameters (see examples).
#' Note that for palindromic kmers, setting strand to 1 or -1 is recommended to avoid double counting.
#'
#' \emph{func = "kmer.fft", params = list(kmer = "CG", freq = 10.2, extend = TRUE, window = "hann")} \cr
#' Performs FFT analysis on kmer occurrence signal and returns the power at the specified frequency.
#' The frequency is specified in cycles per base. Window functions reduce spectral leakage.
#'
#' \emph{func = "kmer.fft.peak", params = list(kmer = "CG", extend = TRUE, window = "hann")} \cr
#' Returns the dominant frequency (in cycles per base) from FFT analysis of kmer occurrences.
#' Useful for detecting periodic patterns in kmer distribution.
#'
#' \emph{func = "kmer.fft.peak.power", params = list(kmer = "CG", extend = TRUE, window = "hann")} \cr
#' Returns the power at the dominant frequency from FFT analysis of kmer occurrences.
#' Indicates the strength of the most prominent periodic pattern.
#'
#' For kmer FFT functions:
#' \itemize{
#'   \item kmer: The DNA sequence to analyze (case-insensitive)
#'   \item freq: The frequency to evaluate in cycles per base (only for kmer.fft)
#'   \item extend: If TRUE (default), considers k-mers starting at any position within the interval, which may require looking at sequence data beyond the interval's end. The FFT is computed on a signal of the interval's length. If FALSE, the sequence is not extended, so k-mers starting near the end of the interval that would cross the boundary are not counted.
#'   \item window: Window function to apply ("none", "hann", "hamming", "blackman", default: "hann"). These functions are used to reduce spectral leakage from the FFT.
#'     \itemize{
#'       \item "none": No window.
#'       \item "hann": Hann window (default). Good for general purpose.
#'       \item "hamming": Hamming window. Minimizes the nearest side lobe.
#'       \item "blackman": Blackman window. Better stop-band attenuation but wider main lobe.
#'     }
#'   For more information, see: \url{https://en.wikipedia.org/wiki/Window_function}
#' }
#'
#' Modify iterator behavior with 'gvtrack.iterator' or 'gvtrack.iterator.2d'.
#'
#' @param vtrack virtual track name
#' @param src source (track/intervals). NULL for PWM functions
#' @param func function name (see above)
#' @param params function parameters (see above)
#' @param ... additional PWM parameters
#' @return None.
#' @seealso \code{\link{gvtrack.info}}, \code{\link{gvtrack.iterator}},
#' \code{\link{gvtrack.iterator.2d}}, \code{\link{gvtrack.array.slice}},
#' \code{\link{gvtrack.ls}}, \code{\link{gvtrack.rm}}
#' @keywords ~virtual
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
#' gextract("dense_track", "vtrack1", "vtrack2",
#'     gintervals(1, 0, 10000),
#'     iterator = 1000
#' )
#'
#' gvtrack.create("vtrack3", "dense_track", "global.percentile")
#' gvtrack.create("vtrack4", "annotations", "distance")
#' gdist(
#'     "vtrack3", seq(0, 1, l = 10), "vtrack4",
#'     seq(-500, 500, 200)
#' )
#'
#' gvtrack.create("cov", "annotations", "coverage")
#' gextract("cov", gintervals(1, 0, 1000), iterator = 100)
#'
#' pssm <- matrix(
#'     c(
#'         0.7, 0.1, 0.1, 0.1, # Example PSSM
#'         0.1, 0.7, 0.1, 0.1,
#'         0.1, 0.1, 0.7, 0.1,
#'         0.1, 0.1, 0.7, 0.1,
#'         0.1, 0.1, 0.7, 0.1,
#'         0.1, 0.1, 0.7, 0.1
#'     ),
#'     ncol = 4, byrow = TRUE
#' )
#' colnames(pssm) <- c("A", "C", "G", "T")
#' gvtrack.create(
#'     "motif_score", NULL, "pwm",
#'     list(pssm = pssm, bidirect = TRUE, prior = 0.01)
#' )
#' gvtrack.create("max_motif_score", NULL, "pwm.max",
#'     pssm = pssm, bidirect = TRUE, prior = 0.01
#' )
#' gvtrack.create("max_motif_pos", NULL, "pwm.max.pos",
#'     pssm = pssm
#' )
#' gextract(
#'     c(
#'         "dense_track", "motif_score", "max_motif_score",
#'         "max_motif_pos"
#'     ),
#'     gintervals(1, 0, 10000),
#'     iterator = 500
#' )
#'
#' # Kmer counting examples
#' gvtrack.create("cg_count", NULL, "kmer.count", kmer = "CG", strand = 1)
#' gvtrack.create("cg_frac", NULL, "kmer.frac", kmer = "CG", strand = 1)
#' gextract(c("cg_count", "cg_frac"), gintervals(1, 0, 10000), iterator = 1000)
#'
#' gvtrack.create("at_pos", NULL, "kmer.count", kmer = "AT", strand = 1)
#' gvtrack.create("at_neg", NULL, "kmer.count", kmer = "AT", strand = -1)
#' gvtrack.create("at_both", NULL, "kmer.count", kmer = "AT", strand = 0)
#' gextract(c("at_pos", "at_neg", "at_both"), gintervals(1, 0, 10000), iterator = 1000)
#'
#' # GC content
#' gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
#' gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
#' gextract("g_frac + c_frac", gintervals(1, 0, 10000),
#'     iterator = 1000,
#'     colnames = "gc_content"
#' )
#'
#' # FFT analysis examples
#' gvtrack.create("cg_fft", NULL, "kmer.fft",
#'     params = list(kmer = "CG", freq = 0.1, extend = TRUE)
#' )
#' gvtrack.create("cg_peak_freq", NULL, "kmer.fft.peak",
#'     params = list(kmer = "CG", extend = TRUE)
#' )
#' gvtrack.create("cg_peak_power", NULL, "kmer.fft.peak.power",
#'     params = list(kmer = "CG", extend = TRUE)
#' )
#' gextract(c("cg_fft", "cg_peak_freq", "cg_peak_power"),
#'     gintervals(1, 0, 10000),
#'     iterator = 1000
#' )
#' @export gvtrack.create
gvtrack.create <- function(vtrack = NULL, src = NULL, func = NULL, params = NULL, ...) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.create(vtrack, src, func = NULL, params = NULL, ...)", call. = FALSE)
    }
    if (is.null(substitute(src)) && !(func %in% c("pwm", "pwm.max", "pwm.max.pos", "kmer.count", "kmer.frac", "kmer.fft", "kmer.fft.peak", "kmer.fft.peak.power"))) {
        stop("Usage: gvtrack.create(vtrack, src, func = NULL, params = NULL, ...)", call. = FALSE)
    }

    .gcheckroot()

    if (!is.null(func) && func %in% c("pwm", "pwm.max", "pwm.max.pos")) {
        dots <- list(...)

        if (!is.null(params)) {
            # params list
            if (!is.list(params) || !("pssm" %in% names(params))) {
                stop("pwm function requires a list with at least 'pssm' matrix parameter")
            }
            dots <- params
        }

        if (!("pssm" %in% names(dots))) {
            stop("pwm function requires a 'pssm' matrix parameter")
        }
        pssm <- dots$pssm
        bidirect <- if (!is.null(dots$bidirect)) dots$bidirect else TRUE
        prior <- if (!is.null(dots$prior)) dots$prior else 0.01
        extend <- if (!is.null(dots$extend)) dots$extend else TRUE
        strand <- if (!is.null(dots$strand)) dots$strand else 1


        if (!all(c("A", "C", "G", "T") %in% colnames(pssm))) {
            stop("PSSM must be a nx4 matrix with colnames A, C, G, T")
        }

        pssm <- pssm[, c("A", "C", "G", "T")]

        if (is.data.frame(pssm)) {
            pssm <- as.matrix(pssm)
        }

        if (!is.numeric(prior) || prior < 0 || prior > 1) {
            stop("prior must be a number between 0 and 1")
        }

        if (!is.logical(bidirect)) {
            stop("bidirect must be TRUE or FALSE")
        }

        if (!is.logical(extend)) {
            stop("extend must be TRUE or FALSE")
        }

        if (strand != 1 && strand != -1) {
            stop("strand must be 1 or -1")
        }

        # Normalize PSSM and add prior
        pssm <- sweep(pssm, 1, rowSums(pssm), "/") # Normalize rows
        if (prior > 0) {
            pssm <- pssm + prior
            pssm <- sweep(pssm, 1, rowSums(pssm), "/") # Renormalize after adding prior
        }

        # Set params with processed values
        params <- list(
            pssm = pssm,
            bidirect = bidirect,
            prior = prior,
            extend = extend,
            strand = strand
        )
    } else if (!is.null(func) && func %in% c("kmer.count", "kmer.frac")) {
        # Check for kmer parameter
        dots <- list(...)

        if (!is.null(params)) {
            # Handle as list or string parameter
            if (is.list(params)) {
                # params list
                if (!is.list(params) || !("kmer" %in% names(params))) {
                    stop("kmer function requires a list with at least 'kmer' parameter")
                }
                kmer_params <- params
            }
        } else if ("kmer" %in% names(dots)) {
            # Use named parameters
            kmer_params <- dots
        } else {
            stop("kmer functions require a 'kmer' parameter")
        }

        kmer_params$extend <- if (!is.null(kmer_params$extend)) kmer_params$extend else TRUE
        kmer_params$strand <- if (!is.null(kmer_params$strand)) kmer_params$strand else 0

        # Validate required kmer parameter
        if (!("kmer" %in% names(kmer_params)) || !is.character(kmer_params$kmer) || length(kmer_params$kmer) != 1) {
            stop("kmer parameter must be a single string")
        }

        if (nchar(kmer_params$kmer) == 0) {
            stop("kmer sequence cannot be empty")
        }

        if (kmer_params$kmer == grevcomp(kmer_params$kmer) && kmer_params$strand == 0) {
            warning(paste0("kmer sequence '", kmer_params$kmer, "' is palindromic, please set strand to 1 or -1 to avoid double counting"))
        }

        params <- kmer_params
    } else if (!is.null(func) && func %in% c("kmer.fft", "kmer.fft.peak", "kmer.fft.peak.power")) {
        dots <- list(...)

        if (!is.null(params)) {
            if (!is.list(params) || !("kmer" %in% names(params))) {
                stop("FFT functions require a list with at least 'kmer' parameter")
            }
            kmer_params <- params
        } else if ("kmer" %in% names(dots)) {
            kmer_params <- dots
        } else {
            stop("FFT functions require a 'kmer' parameter")
        }

        # Validate kmer
        if (!is.character(kmer_params$kmer) || length(kmer_params$kmer) != 1) {
            stop("kmer parameter must be a single string")
        }

        # Set defaults
        kmer_params$extend <- if (!is.null(kmer_params$extend)) kmer_params$extend else TRUE
        kmer_params$window <- if (!is.null(kmer_params$window)) kmer_params$window else "hann"

        # Validate freq for kmer.fft
        if (func == "kmer.fft") {
            if (is.null(kmer_params$freq) || !is.numeric(kmer_params$freq) || length(kmer_params$freq) != 1) {
                stop("kmer.fft requires a numeric 'freq' parameter")
            }
            if (kmer_params$freq < 0 || kmer_params$freq > 0.5) {
                stop("freq parameter must be between 0 and 0.5")
            }
        }

        # Validate window parameter
        if (!is.null(kmer_params$window)) {
            if (!is.character(kmer_params$window) || length(kmer_params$window) != 1) {
                stop("window parameter must be a single string")
            }
            if (!(kmer_params$window %in% c("none", "hann", "hamming", "blackman"))) {
                stop("window parameter must be one of: 'none', 'hann', 'hamming', 'blackman'")
            }
        }

        params <- kmer_params
    }

    vtrackstr <- do.call(.gexpr2str, list(substitute(vtrack)), envir = parent.frame())
    srcstr <- do.call(.gexpr2str, list(substitute(src)), envir = parent.frame())

    if (!is.na(match(vtrackstr, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Cannot create virtual track: regular track named %s already exists", vtrackstr), call. = FALSE)
    }

    if (!is.na(match(vtrackstr, get("GINTERVS", envir = .misha)))) {
        stop(sprintf("Cannot create virtual track: intervals named %s already exists", vtrackstr), call. = FALSE)
    }

    if (vtrackstr != make.names(vtrackstr)) {
        stop(sprintf("\"%s\" is not a syntactically valid name for a variable", vtrackstr), call. = FALSE)
    }

    var <- list()
    if (is.character(srcstr) && !is.na(match(srcstr, get("GTRACKS", envir = .misha)))) {
        var$src <- srcstr
    } else {
        var$src <- src
    }
    var$func <- func
    var$params <- params

    .gvtrack.set(vtrackstr, var)

    retv <- NULL
}



#' Returns the definition of a virtual track
#'
#' Returns the definition of a virtual track.
#'
#' This function returns the internal representation of a virtual track.
#'
#' @param vtrack virtual track name
#' @return Internal representation of a virtual track.
#' @seealso \code{\link{gvtrack.create}}
#' @keywords ~virtual
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.info("vtrack1")
#'
#' @export gvtrack.info
gvtrack.info <- function(vtrack = NULL) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.info(vtrack)", call. = FALSE)
    }
    .gcheckroot()

    vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())
    .gvtrack.get(vtrackstr)
}



#' Defines modification rules for a one-dimensional iterator in a virtual track
#'
#' Defines modification rules for a one-dimensional iterator in a virtual
#' track.
#'
#' This function defines modification rules for one-dimensional iterator
#' intervals in a virtual track.
#'
#' 'dim' converts a 2D iterator interval (chrom1, start1, end1, chrom2, start2,
#' end2) to a 1D interval. If 'dim' is '1' the interval is converted to
#' (chrom1, start1, end1). If 'dim' is '2' the interval is converted to
#' (chrom2, start2, end2). If 1D iterator is used 'dim' must be set to 'NULL'
#' or '0' (meaning: no conversion is made).
#'
#' Iterator interval's 'start' coordinate is modified by adding 'sshift'.
#' Similarly 'end' coordinate is altered by adding 'eshift'.
#'
#' @param vtrack virtual track name
#' @param dim use 'NULL' or '0' for 1D iterators. '1' converts 2D iterator to
#' (chrom1, start1, end1) , '2' converts 2D iterator to (chrom2, start2, end2)
#' @param sshift shift of 'start' coordinate
#' @param eshift shift of 'end' coordinate
#' @return None.
#' @seealso \code{\link{gvtrack.create}}, \code{\link{gvtrack.iterator.2d}}
#' @keywords ~virtual
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' gvtrack.create("vtrack1", "dense_track")
#' gvtrack.iterator("vtrack1", sshift = 200, eshift = 200)
#' gextract("dense_track", "vtrack1", gintervals(1, 0, 500))
#'
#' gvtrack.create("vtrack2", "dense_track")
#' gvtrack.iterator("vtrack2", dim = 1)
#' gextract("vtrack2", gintervals.2d(1, 0, 1000, 1, 0, -1),
#'     iterator = "rects_track"
#' )
#'
#' @export gvtrack.iterator
gvtrack.iterator <- function(vtrack = NULL, dim = NULL, sshift = 0, eshift = 0) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.iterator(vtrack, dim = NULL, sshift = 0, eshift = 0)", call. = FALSE)
    }
    .gcheckroot()

    vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())

    var <- .gvtrack.get(vtrackstr)
    itr <- list()
    itr$type <- "1d"
    itr$dim <- dim
    itr$sshift <- sshift
    itr$eshift <- eshift
    var$itr <- itr
    .gvtrack.set(vtrackstr, var)
    retv <- NULL
}



#' Defines modification rules for a two-dimensional iterator in a virtual track
#'
#' Defines modification rules for a two-dimensional iterator in a virtual
#' track.
#'
#' This function defines modification rules for one-dimensional iterator
#' intervals in a virtual track.
#'
#' Iterator interval's 'start1' coordinate is modified by adding 'sshift1'.
#' Similarly 'end1', 'start2', 'end2' coordinates are altered by adding
#' 'eshift1', 'sshift2' and 'eshift2' accordingly.
#'
#' @param vtrack virtual track name
#' @param sshift1 shift of 'start1' coordinate
#' @param eshift1 shift of 'end1' coordinate
#' @param sshift2 shift of 'start2' coordinate
#' @param eshift2 shift of 'end2' coordinate
#' @return None.
#' @seealso \code{\link{gvtrack.create}}, \code{\link{gvtrack.iterator}}
#' @keywords ~virtual
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "rects_track")
#' gvtrack.iterator.2d("vtrack1", sshift1 = 1000, eshift1 = 2000)
#' gextract(
#'     "rects_track", "vtrack1",
#'     gintervals.2d(1, 0, 5000, 2, 0, 5000)
#' )
#'
#' @export gvtrack.iterator.2d
gvtrack.iterator.2d <- function(vtrack = NULL, sshift1 = 0, eshift1 = 0, sshift2 = 0, eshift2 = 0) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.iterator.2d(vtrack, sshift1 = 0, eshift1 = 0, sshift2 = 0, eshift2 = 0)", call. = FALSE)
    }
    .gcheckroot()

    vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())

    var <- .gvtrack.get(vtrackstr)
    itr <- list()
    itr$type <- "2d"
    itr$sshift1 <- sshift1
    itr$eshift1 <- eshift1
    itr$sshift2 <- sshift2
    itr$eshift2 <- eshift2
    var$itr <- itr
    .gvtrack.set(vtrackstr, var)
    retv <- NULL
}



#' Returns a list of virtual track names
#'
#' Returns a list of virtual track names.
#'
#' This function returns a list of virtual tracks that exist in current R
#' environment that match the pattern (see 'grep'). If called without any
#' arguments all virtual tracks are returned.
#'
#' @param pattern,ignore.case,perl,fixed,useBytes see 'grep'
#' @return An array that contains the names of virtual tracks.
#' @seealso \code{\link{grep}}, \code{\link{gvtrack.create}},
#' \code{\link{gvtrack.rm}}
#' @keywords ~virtual ~ls
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
#' gvtrack.ls()
#' gvtrack.ls(pattern = "*2")
#'
#' @export gvtrack.ls
gvtrack.ls <- function(pattern = "", ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
    .gcheckroot()

    if (!exists("GVTRACKS", envir = .misha)) {
        return(NULL)
    }

    gvtracks <- get("GVTRACKS", envir = .misha)
    gwds <- names(gvtracks)
    if (!is.list(gvtracks) || (length(gvtracks) && !is.character(gwds)) || length(gvtracks) != length(gwds)) {
        stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = FALSE)
    }

    gwd <- get("GWD", envir = .misha)
    idx <- match(gwd, gwds)
    if (is.na(idx)) {
        return(NULL)
    }

    vtracks <- gvtracks[[idx]]
    vtracknames <- names(vtracks)
    if (!is.list(vtracks) || (length(vtracks) && !is.character(vtracknames)) || length(vtracks) != length(vtracknames)) {
        stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = FALSE)
    }

    if (!length(vtracks)) {
        return(NULL)
    }

    if (pattern != "") {
        grep(pattern, vtracknames, value = TRUE, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
    } else {
        vtracknames
    }
}



#' Deletes a virtual track
#'
#' Deletes a virtual track.
#'
#' This function deletes a virtual track from current R environment.
#'
#' @param vtrack virtual track name
#' @return None.
#' @seealso \code{\link{gvtrack.create}}, \code{\link{gvtrack.ls}}
#' @keywords ~virtual
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
#' gvtrack.ls()
#' gvtrack.rm("vtrack1")
#' gvtrack.ls()
#'
#' @export gvtrack.rm
gvtrack.rm <- function(vtrack = NULL) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.rm(vtrack)", call. = FALSE)
    }
    .gcheckroot()

    vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())
    .gvtrack.set(vtrackstr, NULL)
    retv <- NULL
}



#' Defines rules for a single value calculation of a virtual 'Array' track
#'
#' Defines how a single value within an interval is achieved for a virtual
#' track based on 'Array' track.
#'
#' A track (regular or virtual) used in a track expression is expected to
#' return one value for each track interval. 'Array' tracks store multiple
#' values per interval (one for each 'column') and hence if used in a track
#' expression one must define the way of how a single value should be deduced
#' from several ones.
#'
#' By default if an 'Array' track is used in a track expressions, its interval
#' value would be the average of all column values that are not NaN.
#' 'gvtrack.array.slice' allows to select specific columns and to specify the
#' function applied to their values.
#'
#' 'slice' parameter allows to choose the columns. Columns can be indicated by
#' their names or their indices. If 'slice' is 'NULL' the non-NaN values of all
#' track columns are used.
#'
#' 'func' parameter determines the function applied to the columns' values. Use
#' the following table for a reference of all valid functions and parameters
#' combinations:
#'
#' \emph{func = "avg", params = NULL} \cr Average of columns' values.
#'
#' \emph{func = "max", params = NULL} \cr Maximum of columns' values.
#'
#' \emph{func = "min", params = NULL} \cr Minimum of columns' values.
#'
#' \emph{func = "stdev", params = NULL} \cr Unbiased standard deviation of
#' columns' values.
#'
#' \emph{func = "sum", params = NULL} \cr Sum of columns' values.
#'
#' \emph{func = "quantile", params = [Percentile in the range of [0, 1]]} \cr
#' Quantile of columns' values.
#'
#' @param vtrack virtual track name
#' @param slice a vector of column names or column indices or 'NULL'
#' @param func,params see below
#' @return None.
#' @seealso \code{\link{gvtrack.create}},
#' \code{\link{gtrack.array.get_colnames}}, \code{\link{gtrack.array.extract}}
#' @keywords ~virtual ~array
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "array_track")
#' gvtrack.array.slice("vtrack1", c("col2", "col4"), "max")
#' gextract("vtrack1", gintervals(1, 0, 1000))
#'
#' @export gvtrack.array.slice
gvtrack.array.slice <- function(vtrack = NULL, slice = NULL, func = "avg", params = NULL) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.array.slice(vtrack, slice = NULL, func = \"avg\", params = NULL)", call. = FALSE)
    }
    .gcheckroot()

    vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())

    var <- .gvtrack.get(vtrackstr)
    .slice <- list()
    .slice$slice <- .gslice(var$src, slice)$slice
    .slice$func <- func
    .slice$params <- params
    var$slice <- .slice
    .gvtrack.set(vtrackstr, var)

    retv <- NULL
}
