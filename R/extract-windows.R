#' One-hot encoding of DNA sequences
#'
#' Returns the one-hot encoding of genomic intervals of equal width, read
#' directly from the genome sequence files, or of character strings of equal
#' length.
#'
#' The result is an array of dimensions \code{c(N, 4, L)}: \code{x[n, c, p]} is
#' 1 when position \code{p} of sequence \code{n} is base \code{c}, and 0
#' otherwise. The channels are A, C, G and T, in that order (the names of the
#' second dimension). Lowercase (soft-masked) bases are encoded like uppercase
#' ones. Any other character (N, IUPAC ambiguity codes, gaps) and, with
#' \code{pad = TRUE}, every position outside the chromosome has all four
#' channels 0.
#'
#' The axes mean the same as in pymisha's \code{gseq_onehot}, whose numpy array
#' has shape \code{(N, 4, L)} as well: element \code{[n, c, p]} is the same
#' value in both. (The memory order differs: R arrays are column-major, numpy's
#' are row-major. Converting to a torch tensor in either language keeps the
#' indexing, i.e. \code{(batch, channels, length)}.)
#'
#' For intervals with a 'strand' column, a row with strand -1 (or "-") is
#' reverse complemented, as in \code{\link{gseq.extract}}: its encoding is that
#' of the reverse complement of the forward sequence. The encoding is computed
#' in C++ straight from the sequence store; no character strings are made.
#'
#' The array takes \code{N * 4 * L} bytes for \code{dtype = "raw"}, 4 times
#' that for "integer" and 8 times for "double".
#'
#' @param x intervals of equal width (a data frame or the name of an intervals
#' set), or a character vector of sequences of equal length
#' @param pad for intervals: if TRUE, an interval may run past either end of its
#' chromosome and the positions outside it are all 0, as the "N" of
#' \code{gseq.extract(pad = TRUE)}. If FALSE (the default), such an interval is
#' an error.
#' @param dtype storage type of the result: "raw" (one byte per value, the
#' default), "integer" or "double"
#' @return An array of dimensions \code{c(N, 4, L)}.
#' @seealso \code{\link{gseq.extract}}, \code{\link{gextract.windows}}
#' @keywords ~extract ~DNA ~sequence
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' x <- gseq.onehot(gintervals(1, c(10000, 20000), c(10010, 20010)))
#' dim(x)
#' x[1, , ]
#'
#' gseq.onehot(c("ACGTN", "acgtn"), dtype = "integer")[2, , ]
#'
#' @export gseq.onehot
gseq.onehot <- function(x = NULL, pad = FALSE, dtype = c("raw", "integer", "double")) {
    if (is.null(x)) {
        stop("Usage: gseq.onehot(x, pad = FALSE, dtype = \"raw\")", call. = FALSE)
    }
    if (!is.logical(pad) || length(pad) != 1 || is.na(pad)) {
        stop("pad must be TRUE or FALSE", call. = FALSE)
    }
    dtype <- .gonehot_dtype(dtype)

    if (is.character(x) && !(length(x) == 1 && .gonehot_is_intervals_set(x))) {
        res <- .gcall("C_gseq_onehot_strings", x, FALSE, dtype)
        return(.gonehot_dimnames(res))
    }

    .gcheckroot()
    if (is.character(x)) {
        x <- gintervals.load(x)
    }
    if (!is.data.frame(x) || ncol(x) < 3 || !identical(names(x)[1:3], c("chrom", "start", "end"))) {
        stop("x must be 1D intervals or a character vector of sequences", call. = FALSE)
    }
    if (!nrow(x)) {
        stop("x is empty", call. = FALSE)
    }
    w <- .gonehot_width(x$start, x$end)
    chroms <- as.character(.gchroms(as.character(x$chrom)))
    res <- .gcall(
        "C_gseq_onehot", chroms, as.numeric(x$start), as.numeric(w), .gonehot_minus(x$strand, nrow(x)),
        pad, dtype, .misha_env()
    )
    .gonehot_dimnames(res)
}

.gonehot_dtype <- function(dtype) {
    dtype <- match.arg(dtype, c("raw", "integer", "double"))
    match(dtype, c("raw", "integer", "double")) - 1L
}

.gonehot_is_intervals_set <- function(x) {
    exists("GROOT", envir = .misha, inherits = FALSE) &&
        x %in% get("GINTERVS", envir = .misha)
}

.gonehot_dimnames <- function(res) {
    dimnames(res) <- list(NULL, c("A", "C", "G", "T"), NULL)
    res
}

.gonehot_width <- function(starts, ends) {
    if (!is.numeric(starts) || !is.numeric(ends) || anyNA(starts) || anyNA(ends) ||
        any(!is.finite(starts)) || any(!is.finite(ends))) {
        stop("Interval coordinates must be finite numbers", call. = FALSE)
    }
    if (any(starts != round(starts)) || any(ends != round(ends))) {
        stop("Interval coordinates must be whole numbers", call. = FALSE)
    }
    w <- ends[1] - starts[1]
    if (w <= 0 || any(ends - starts != w)) {
        stop("All intervals must have the same positive width", call. = FALSE)
    }
    w
}

# TRUE for the rows with strand -1 or "-"
.gonehot_minus <- function(strand, n = length(strand)) {
    if (is.null(strand)) {
        return(logical(n))
    }
    neg <- if (is.numeric(strand)) strand == -1 else as.character(strand) == "-"
    neg %in% TRUE
}

#' Fixed-width windows as arrays for machine learning
#'
#' Extracts fixed-width windows around a set of loci and returns their DNA as a
#' one-hot array and the values of track expressions over them as a numeric
#' array, ready to feed a sequence model. One call replaces
#' \code{\link{gseq.extract}}, a one-hot encoding of the strings and
#' \code{\link{gextract.matrix}}; the one-hot array is encoded in C++ straight
#' from the sequence files.
#'
#' Each locus (a row of 'intervals') is centred at \code{mid = start + (end -
#' start) \%/\% 2}. Its sequence window is \code{[mid - width \%/\% 2 - jitter,
#' mid - width \%/\% 2 + width + jitter)} and its signal window is the same
#' with \code{signal_width}: an odd width reaches one base further right than
#' left. A positive \code{jitter} widens both windows by \code{jitter} on each
#' side, so that a data loader can crop randomly shifted windows of the
#' requested widths later without reading the genome again (as in tangermeme's
#' \code{extract_loci}).
#'
#' Rows are dropped, and listed in \code{dropped} with the first reason that
#' applies, when
#' \itemize{
#' \item{"off_chrom": the sequence window or (with tracks) the signal window
#' runs past either end of the chromosome;}
#' \item{"excluded": either window overlaps an interval of 'exclude' (by at
#' least one base; intervals of 'exclude' on chromosomes the genome does not have
#' are ignored);}
#' \item{"min_counts" / "max_counts": the sum of the values of track expression
#' 'count_track' over the signal window (as returned, i.e. after binning, NaN
#' counted as 0) is below 'min_counts' or above 'max_counts'.}
#' }
#'
#' With \code{use_strand = TRUE} and a 'strand' column in 'intervals', the rows
#' with strand -1 (or "-") are reverse complemented: their one-hot encoding is
#' that of the reverse complement, and their signals are reversed along the
#' length axis.
#'
#' The arrays have the same axes as pymisha's \code{gextract_windows}:
#' \code{x[n, c, p]} and \code{y[n, t, b]}, i.e. \code{(batch, channels,
#' length)}. (R arrays are column-major, so the memory order differs from
#' numpy's; indexing, and conversion to a torch tensor, are the same.)
#'
#' Memory: 'x' takes \code{N * 4 * (width + 2 * jitter)} bytes with
#' \code{dtype = "raw"} (4 and 8 times that for "integer" and "double"), and
#' 'y' takes \code{N * T * (signal_width + 2 * jitter) / binsize * 8} bytes for
#' T track expressions. The signals are extracted with
#' \code{\link{gextract.matrix}} in chunks of about 1e7 values, so the
#' extraction itself adds little to that. For 100,000 windows of 2 kb, 'x'
#' takes 800 MB; 'y' of three tracks at 1 bp resolution takes 4.8 GB, or 150 MB
#' in 32 bp bins.
#'
#' @param intervals the loci: 1D intervals (a data frame or the name of an
#' intervals set), or a list of such, e.g. \code{list(pos = peaks, neg =
#' background)}, whose rows are concatenated and labelled by a 'set' column
#' (the names of the list, or their indices)
#' @param tracks track expressions (a character vector) for 'y', or NULL for
#' sequence only. Its names, if any, name the second dimension of 'y'.
#' @param width width of the sequence window
#' @param signal_width width of the signal window
#' @param jitter extra margin added to each side of both windows
#' @param binsize bin size of the signals: a positive integer that divides
#' \code{signal_width + 2 * jitter}. Each value is the one
#' \code{\link{gextract.matrix}} returns for the bin (e.g. the mean of a dense
#' track; use a virtual track with \code{func = "sum"} for sums).
#' @param use_strand if TRUE, rows with strand -1 are reverse complemented (see
#' Details)
#' @param exclude intervals whose overlapping windows are dropped, or NULL
#' @param min_counts,max_counts bounds on the sum of 'count_track' over the
#' signal window (inclusive), or NULL
#' @param count_track the track expression the count filters use: an index into
#' 'tracks', or one of its names
#' @param sequence if FALSE, 'x' is not computed (NULL)
#' @param dtype storage type of 'x': "raw" (one byte per value, the default),
#' "integer" or "double"
#' @return A list with
#' \item{x}{the one-hot sequences, an array of dimensions \code{c(N, 4, width
#' + 2 * jitter)} with channels A, C, G, T (see \code{\link{gseq.onehot}}), or
#' NULL}
#' \item{y}{the signals, a numeric array of dimensions \code{c(N, T,
#' (signal_width + 2 * jitter) / binsize)}, or NULL without 'tracks'}
#' \item{intervals}{the kept rows, one per window, in their original order:
#' 'chrom', 'start' and 'end' of the sequence window (jitter included), the
#' 'strand' column if any, 'signal_start' and 'signal_end' with tracks, the
#' other columns of 'intervals', 'index' (the row in 'intervals', 1-based) and,
#' for a list, 'set'}
#' \item{dropped}{the dropped rows: 'index' and 'reason'}
#' @seealso \code{\link{gseq.onehot}}, \code{\link{gextract.matrix}},
#' \code{\link{gseq.extract}}
#' @keywords ~extract ~DNA ~sequence
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' loci <- data.frame(chrom = "chr1", start = c(5000, 20000, 100), end = c(5001, 20001, 101))
#' w <- gextract.windows(loci, tracks = c(dense = "dense_track"), width = 500, binsize = 50)
#' dim(w$x)
#' dim(w$y)
#' w$dropped # the third locus is too close to the chromosome start
#'
#' @export gextract.windows
gextract.windows <- function(intervals = NULL,
                             tracks = NULL,
                             width = NULL,
                             signal_width = width,
                             jitter = 0,
                             binsize = 1,
                             use_strand = TRUE,
                             exclude = NULL,
                             min_counts = NULL,
                             max_counts = NULL,
                             count_track = 1,
                             sequence = TRUE,
                             dtype = c("raw", "integer", "double")) {
    if (is.null(intervals) || is.null(width)) {
        stop("Usage: gextract.windows(intervals, tracks = NULL, width, signal_width = width, jitter = 0, binsize = 1, use_strand = TRUE, exclude = NULL, min_counts = NULL, max_counts = NULL, count_track = 1, sequence = TRUE, dtype = \"raw\")", call. = FALSE)
    }
    .gcheckroot()
    envir <- parent.frame()
    dtype <- .gonehot_dtype(dtype)

    is_count <- function(v, name, min = 0) {
        if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v != round(v) || v < min) {
            stop(sprintf("%s must be a whole number of at least %d", name, min), call. = FALSE)
        }
        as.numeric(v)
    }
    width <- is_count(width, "width", 1)
    jitter <- is_count(jitter, "jitter", 0)
    binsize <- is_count(binsize, "binsize", 1)
    for (flag in c("use_strand", "sequence")) {
        v <- get(flag)
        if (!is.logical(v) || length(v) != 1 || is.na(v)) {
            stop(sprintf("%s must be TRUE or FALSE", flag), call. = FALSE)
        }
    }

    if (!is.null(tracks)) {
        if (!is.character(tracks) || !length(tracks) || anyNA(tracks)) {
            stop("tracks must be a character vector of track expressions", call. = FALSE)
        }
        signal_width <- is_count(signal_width, "signal_width", 1)
        signal_len <- signal_width + 2 * jitter
        if (signal_len %% binsize != 0) {
            stop(sprintf("binsize (%g) does not divide signal_width + 2 * jitter (%g)", binsize, signal_len), call. = FALSE)
        }
        track_names <- if (is.null(names(tracks))) unname(tracks) else ifelse(names(tracks) == "", tracks, names(tracks))
    }
    if (!is.null(min_counts) || !is.null(max_counts)) {
        if (is.null(tracks)) {
            stop("min_counts and max_counts are measured on tracks, so tracks must be given", call. = FALSE)
        }
        for (v in list(min_counts, max_counts)) {
            if (!is.null(v) && (!is.numeric(v) || length(v) != 1 || is.na(v))) {
                stop("min_counts and max_counts must be numbers or NULL", call. = FALSE)
            }
        }
        ct <- if (is.character(count_track)) match(count_track, track_names) else count_track
        if (length(ct) != 1 || is.na(ct) || ct != round(ct) || ct < 1 || ct > length(tracks)) {
            stop("count_track must be an index into tracks or one of their names", call. = FALSE)
        }
    }

    # the loci
    load1 <- function(iv) {
        if (is.character(iv) && length(iv) == 1) {
            iv <- gintervals.load(iv)
        }
        if (!is.data.frame(iv) || ncol(iv) < 3 || !identical(names(iv)[1:3], c("chrom", "start", "end"))) {
            stop("intervals must be 1D intervals or a list of 1D intervals", call. = FALSE)
        }
        iv
    }
    set <- NULL
    if (is.list(intervals) && !is.data.frame(intervals)) {
        sets <- lapply(intervals, load1)
        labels <- if (is.null(names(sets))) seq_along(sets) else names(sets)
        set <- rep(labels, vapply(sets, nrow, integer(1)))
        cols <- Reduce(intersect, lapply(sets, names))
        intervals <- do.call(rbind, lapply(sets, function(iv) as.data.frame(iv)[, cols, drop = FALSE]))
        rownames(intervals) <- NULL
    } else {
        intervals <- load1(intervals)
    }
    n0 <- nrow(intervals)
    starts <- intervals$start
    ends <- intervals$end
    if (!is.numeric(starts) || !is.numeric(ends) || any(!is.finite(starts)) || any(!is.finite(ends)) ||
        any(starts != round(starts)) || any(ends != round(ends))) {
        stop("Interval coordinates must be finite whole numbers", call. = FALSE)
    }
    if (any(ends <= starts)) {
        stop("Interval start must be lesser than end", call. = FALSE)
    }
    chroms <- as.character(.gchroms(as.character(intervals$chrom)))
    allgenome <- get("ALLGENOME", envir = .misha)[[1]]
    chrom_lens <- allgenome$end[match(chroms, as.character(allgenome$chrom))]

    mids <- starts + (ends - starts) %/% 2
    win_len <- width + 2 * jitter
    seq_start <- mids - width %/% 2 - jitter
    lo <- seq_start
    hi <- seq_start + win_len
    if (!is.null(tracks)) {
        sig_start <- mids - signal_width %/% 2 - jitter
        lo <- pmin(lo, sig_start)
        hi <- pmax(hi, sig_start + signal_len)
    }

    reason <- rep(NA_character_, n0)
    reason[lo < 0 | hi > chrom_lens] <- "off_chrom"
    if (!is.null(exclude)) {
        hit <- .gwindows_overlap(chroms, lo, hi, load1(exclude))
        reason[is.na(reason) & hit] <- "excluded"
    }

    # signals of the rows left so far
    y <- NULL
    keep <- which(is.na(reason))
    if (!is.null(tracks)) {
        num_bins <- signal_len / binsize
        y <- array(NA_real_, dim = c(length(keep), length(tracks), num_bins))
        if (length(keep)) {
            sig_iv <- data.frame(chrom = chroms[keep], start = sig_start[keep], end = sig_start[keep] + signal_len)
            chunk <- max(1, floor(1e7 / (num_bins * length(tracks))))
            for (first in seq(1, length(keep), by = chunk)) {
                rows <- first:min(first + chunk - 1, length(keep))
                m <- do.call(gextract.matrix, c(as.list(unname(tracks)), list(intervals = sig_iv[rows, ], iterator = binsize)), envir = envir)
                m <- array(m, dim = c(length(rows), num_bins, length(tracks)))
                for (t in seq_along(tracks)) {
                    y[rows, t, ] <- m[, , t]
                }
                rm(m)
            }
            if (use_strand) {
                neg <- which(.gonehot_minus(intervals$strand, n0)[keep])
                if (length(neg)) {
                    y[neg, , ] <- y[neg, , num_bins:1, drop = FALSE]
                }
            }
        }
        dimnames(y) <- list(NULL, track_names, NULL)

        if (!is.null(min_counts) || !is.null(max_counts)) {
            counts <- rowSums(matrix(y[, ct, ], nrow = length(keep)), na.rm = TRUE)
            failed <- rep(NA_character_, length(keep))
            if (!is.null(max_counts)) {
                failed[counts > max_counts] <- "max_counts"
            }
            if (!is.null(min_counts)) {
                failed[counts < min_counts] <- "min_counts"
            }
            reason[keep] <- failed
            ok <- is.na(failed)
            if (!all(ok)) {
                y <- y[ok, , , drop = FALSE]
            }
            keep <- keep[ok]
        }
    }

    x <- NULL
    if (sequence) {
        x <- .gcall(
            "C_gseq_onehot", chroms[keep], as.numeric(seq_start[keep]), win_len,
            .gonehot_minus(if (use_strand) intervals$strand else NULL, n0)[keep], FALSE, dtype, .misha_env()
        )
        x <- .gonehot_dimnames(x)
    }

    out <- data.frame(
        chrom = intervals$chrom[keep], start = seq_start[keep], end = seq_start[keep] + win_len,
        stringsAsFactors = FALSE
    )
    if (!is.null(intervals$strand)) {
        out$strand <- intervals$strand[keep]
    }
    if (!is.null(tracks)) {
        out$signal_start <- sig_start[keep]
        out$signal_end <- sig_start[keep] + signal_len
    }
    extra <- setdiff(names(intervals), c("chrom", "start", "end", "strand", names(out), "index", "set"))
    for (col in extra) {
        out[[col]] <- intervals[[col]][keep]
    }
    out$index <- keep
    if (!is.null(set)) {
        out$set <- set[keep]
    }
    dropped_idx <- which(!is.na(reason))
    dropped <- data.frame(index = dropped_idx, reason = reason[dropped_idx], stringsAsFactors = FALSE)

    list(x = x, y = y, intervals = out, dropped = dropped)
}

# TRUE for the windows [lo, hi) on chroms that overlap an interval of exclude by at least one base
.gwindows_overlap <- function(chroms, lo, hi, exclude) {
    hit <- logical(length(chroms))
    if (!nrow(exclude)) {
        return(hit)
    }
    ex <- data.frame(chrom = as.character(exclude$chrom), start = exclude$start, end = exclude$end, stringsAsFactors = FALSE)
    ex <- ex[ex$end > ex$start, , drop = FALSE]
    # chromosomes the genome does not have are ignored
    uchroms <- unique(ex$chrom)
    norm <- vapply(uchroms, function(ch) {
        tryCatch(as.character(.gchroms(ch)), error = function(e) NA_character_)
    }, character(1))
    ex$chrom <- unname(norm[ex$chrom])
    ex <- ex[!is.na(ex$chrom), , drop = FALSE]

    for (ch in intersect(unique(ex$chrom), unique(chroms))) {
        e <- ex[ex$chrom == ch, , drop = FALSE]
        e <- e[order(e$start), , drop = FALSE]
        max_end <- cummax(e$end)
        w <- which(chroms == ch)
        # the exclusion intervals that start before each window's end; one of them overlaps the
        # window iff the furthest end among them is past the window's start
        j <- findInterval(hi[w], e$start, left.open = TRUE)
        hit[w] <- j > 0 & max_end[pmax(j, 1)] > lo[w]
    }
    hit
}
