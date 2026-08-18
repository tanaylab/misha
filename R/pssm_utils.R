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

# Normalize a `score.thresh` to a single double, for every function that takes
# one: pwm.count, gseq.pwm(mode = "count"), the pwm.edit_distance family and
# pwm.n_mutations.
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
# coercions. gseq.pwm(mode = "count") reads the threshold with an
# unconditional Rf_asReal() (GseqString.cpp:631), so a vector silently
# collapses to its first element there. The virtual-track
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
