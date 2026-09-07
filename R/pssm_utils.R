# Internal helpers for PSSM handling shared across sequence utilities and virtual tracks

.coerce_pssm_matrix <- function(pssm,
                                numeric_msg = "pssm must be a numeric matrix",
                                ncol_msg = "pssm must have columns named A, C, G, T",
                                colnames_msg = "pssm columns must be named A, C, G, T") {
    # Handle data frames by extracting required columns first
    if (is.data.frame(pssm)) {
        cols <- colnames(pssm)
        if (is.null(cols) || !all(c("A", "C", "G", "T") %in% cols)) {
            stop(colnames_msg)
        }
        # Extract only the required columns before converting to matrix
        # This avoids issues with non-numeric columns
        pssm <- as.matrix(pssm[, c("A", "C", "G", "T"), drop = FALSE])
    }

    if (!is.matrix(pssm) || !is.numeric(pssm)) {
        stop(numeric_msg)
    }

    cols <- colnames(pssm)
    if (is.null(cols) || !all(c("A", "C", "G", "T") %in% cols)) {
        stop(colnames_msg)
    }

    pssm[, c("A", "C", "G", "T"), drop = FALSE]
}

# One validator for both Potts entry points - gseq.potts() and
# .vtrack_params_potts(). It returns the model in the exact shape C++ takes:
# `e` a W x 4 double matrix, `J` an npair x 16 double matrix, `pairs` an
# npair x 2 integer matrix (1-based), `intercept` a single double.
#
# It is deliberately tolerant of the extra elements a motifmodel Potts carries
# (`width`, `pair_strength`, `attr`, `link`) so that a fitted model can be
# passed verbatim, and it CROSS-CHECKS `width` against nrow(e) when present,
# which catches a hand-edited model.
.coerce_potts_model <- function(model, what) {
    if (!is.list(model)) {
        stop(sprintf("%s: `model` must be a list with 'e', 'J', 'pairs' and 'intercept'; it is a %s", what, class(model)[1L]), call. = FALSE)
    }
    if (!("e" %in% names(model))) {
        stop(sprintf("%s: `model` requires an 'e' matrix - one row per position, columns A, C, G, T", what), call. = FALSE)
    }

    e <- .coerce_pssm_matrix(
        model$e,
        numeric_msg = sprintf("%s: 'e' must be a numeric matrix or data frame with numeric columns", what),
        ncol_msg = sprintf("%s: 'e' must have columns named A, C, G, T", what),
        colnames_msg = sprintf("%s: 'e' must have columns named A, C, G, T", what)
    )
    if (nrow(e) < 1L) {
        stop(sprintf("%s: 'e' has no rows, so the model covers no positions", what), call. = FALSE)
    }
    if (any(!is.finite(e))) {
        stop(sprintf("%s: 'e' has %d non-finite entries", what, sum(!is.finite(e))), call. = FALSE)
    }
    W <- nrow(e)

    if (!is.null(model$width)) {
        if (!is.numeric(model$width) || length(model$width) != 1L || model$width != W) {
            stop(sprintf("%s: 'width' says %s but 'e' has %d rows", what, paste(model$width, collapse = ", "), W), call. = FALSE)
        }
    }

    pairs <- model$pairs
    if (is.null(pairs)) {
        pairs <- matrix(integer(0), 0L, 2L)
    }
    if (is.data.frame(pairs)) {
        pairs <- as.matrix(pairs)
    }
    if (!is.matrix(pairs) || !is.numeric(pairs) || ncol(pairs) != 2L) {
        stop(sprintf("%s: 'pairs' must be a numeric matrix with 2 columns, one row per coupling", what), call. = FALSE)
    }
    storage.mode(pairs) <- "integer"
    dimnames(pairs) <- NULL
    npair <- nrow(pairs)
    if (npair) {
        if (anyNA(pairs) || any(pairs < 1L) || any(pairs > W)) {
            stop(sprintf("%s: every 'pairs' index must be in 1..%d", what, W), call. = FALSE)
        }
        if (any(pairs[, 1L] >= pairs[, 2L])) {
            stop(sprintf("%s: every row of 'pairs' must have the lower position first (i < j); row %d does not", what, which(pairs[, 1L] >= pairs[, 2L])[1L]), call. = FALSE)
        }
        if (anyDuplicated(paste(pairs[, 1L], pairs[, 2L]))) {
            stop(sprintf("%s: 'pairs' names the same pair of positions twice; a pair carries one coupling table", what), call. = FALSE)
        }
    }

    J <- model$J
    if (is.null(J)) {
        J <- list()
    }
    if (is.list(J)) {
        if (length(J) != npair) {
            stop(sprintf("%s: 'J' holds %d tables and 'pairs' names %d", what, length(J), npair), call. = FALSE)
        }
        for (k in seq_along(J)) {
            Jk <- J[[k]]
            if (!is.matrix(Jk) || !is.numeric(Jk) || !identical(dim(Jk), c(4L, 4L))) {
                stop(sprintf("%s: 'J' table %d must be a 4 x 4 numeric matrix", what, k), call. = FALSE)
            }
            if (any(!is.finite(Jk))) {
                stop(sprintf("%s: 'J' table %d has non-finite entries", what, k), call. = FALSE)
            }
        }
        # as.numeric() of a 4x4 is column-major, so element [a, b] lands at
        # index (b - 1) * 4 + a - which is what the C++ kernel indexes.
        Jflat <- if (npair) {
            matrix(unlist(lapply(J, as.numeric), use.names = FALSE), nrow = npair, byrow = TRUE)
        } else {
            matrix(numeric(0), 0L, 16L)
        }
    } else {
        if (!is.matrix(J) || !is.numeric(J) || nrow(J) != npair || (npair && ncol(J) != 16L)) {
            stop(sprintf("%s: a flat 'J' must be an %d x 16 numeric matrix", what, npair), call. = FALSE)
        }
        if (any(!is.finite(J))) {
            stop(sprintf("%s: 'J' has non-finite entries", what), call. = FALSE)
        }
        Jflat <- J
    }
    storage.mode(Jflat) <- "double"
    dimnames(Jflat) <- NULL

    intercept <- if (is.null(model$intercept)) 0 else model$intercept
    if (!is.numeric(intercept) || length(intercept) != 1L || !is.finite(intercept)) {
        stop(sprintf("%s: 'intercept' must be a single finite number", what), call. = FALSE)
    }

    storage.mode(e) <- "double"
    list(e = e, J = Jflat, pairs = pairs, intercept = as.numeric(intercept))
}

# The C++-ready parameter list, built one way for both entry points. Its shape
# is what PottsParams::parse() reads, so a field renamed here must be renamed
# there; that is the price of having one parser instead of two.
.potts_params <- function(model, bidirect = TRUE, extend = TRUE, strand = 1L,
                          score.thresh = 0, what) {
    p <- .coerce_potts_model(model, what)

    if (!is.logical(bidirect) || length(bidirect) != 1L || is.na(bidirect)) {
        stop(sprintf("%s: bidirect must be TRUE or FALSE", what), call. = FALSE)
    }
    if (!is.logical(extend) || length(extend) != 1L || is.na(extend)) {
        stop(sprintf("%s: extend must be TRUE or FALSE", what), call. = FALSE)
    }
    if (!is.numeric(strand) || length(strand) != 1L || !(strand %in% c(1, -1))) {
        stop(sprintf("%s: strand must be 1 or -1", what), call. = FALSE)
    }
    if (!is.numeric(score.thresh) || length(score.thresh) != 1L || is.na(score.thresh)) {
        stop(sprintf("%s: score.thresh must be a single number", what), call. = FALSE)
    }

    # Under bidirect both strands are scored at every anchor, so `strand` cannot
    # change which anchors are scored or what any of them scores. All it changes
    # is the orientation the scan target is fetched in, and hence the ORDER the
    # anchors are visited in - and the *.pos functions break an exact tie by
    # keeping the first anchor they see. So an unclamped strand = -1 reports a
    # different, equally maximal, position for the same interval wherever two
    # anchors tie, which a repeat region supplies readily. Clamped here rather
    # than at either entry point so both get it: the pwm family clamps the same
    # combination, for the same reason.
    if (bidirect) {
        strand <- 1L
    }

    list(
        e = p$e,
        J = p$J,
        pairs = p$pairs,
        intercept = p$intercept,
        bidirect = bidirect,
        extend = extend,
        strand = as.integer(strand),
        score.thresh = as.numeric(score.thresh)
    )
}

# Normalize a `score.thresh` to a single double, for every function that takes
# one: pwm.count, gseq.pwm(mode = "count"), the pwm.edit_distance family,
# pwm.n_mutations and gseq.pwm_edits.
#
# They all compare the threshold against a PWM log-likelihood - the LSE
# variants against a log-sum-exp of those scores, which can be positive where a
# single score cannot, but nothing here or in C++ constrains the range - so one
# contract covers all of them.
#
# Numbers, character spellings of numbers and factors of those are accepted:
# thresholds routinely arrive from a config file or a read.csv column, and
# gseq.pwm's as.numeric() in the .Call always accepted them.
#
# Everything else is rejected, and the rejections matter as much as the
# coercions. The two gseq entry points read the threshold with an
# unconditional Rf_asReal() - GseqString.cpp:631 and GseqPwmEdits.cpp:1097 - so
# a vector silently collapses to its first element there. The virtual-track
# path errors instead of truncating (TrackExpressionParams.h:233 checks
# Rf_isReal and Rf_length before reading REAL(rthresh)[0]), but out of C++ and
# about the wrong thing. A logical would coerce to 1 or 0 - a threshold no PSSM
# reaches with a non-zero prior, so it would count nothing, or need an
# impossible number of edits, forever. Both are the silent-wrong-answer class
# that making this parameter mandatory exists to close, so neither may pass
# through the coercion that closes it.
.coerce_score_thresh <- function(score.thresh) {
    if (length(score.thresh) == 0) {
        stop("score.thresh must be a single value, and this one is empty: NULL, or a zero-length vector.", call. = FALSE)
    }

    if (length(score.thresh) != 1) {
        stop(sprintf(
            "score.thresh must be a single value, got a vector of length %d. Only one threshold is used, so pass the single value you mean.",
            length(score.thresh)
        ), call. = FALSE)
    }

    # A factor from read.csv() is a character value as far as the caller is
    # concerned; as.numeric() on it would silently return the level index.
    if (is.factor(score.thresh)) {
        score.thresh <- as.character(score.thresh)
    }

    if (!is.numeric(score.thresh) && !is.character(score.thresh)) {
        stop(sprintf(
            "score.thresh must be a single number, or a character spelling of one, not a %s.",
            class(score.thresh)[1]
        ), call. = FALSE)
    }

    out <- suppressWarnings(as.numeric(score.thresh))
    if (is.na(out)) {
        stop(sprintf(
            "score.thresh must be a single number, or something coercible to one; '%s' is not.",
            paste(as.character(score.thresh), collapse = "")
        ), call. = FALSE)
    }

    out
}
