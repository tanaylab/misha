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

# Normalize a `score.thresh` (pwm.count / gseq.pwm(mode = "count")) to a single
# double.
#
# A length-1 character is coerced: thresholds routinely arrive from a config
# file or a read.csv column, and gseq.pwm's as.numeric() in the .Call always
# accepted them. A vector is rejected instead of silently collapsing to its
# first element, which is the same silent-wrong-answer class the mandatory
# threshold closes.
.coerce_score_thresh <- function(score.thresh) {
    if (length(score.thresh) != 1) {
        stop(sprintf(
            "score.thresh must be a single value, got a vector of length %d. Counting uses one threshold, so the remaining elements would be dropped silently.",
            length(score.thresh)
        ), call. = FALSE)
    }

    # A factor from read.csv() is a character value as far as the caller is
    # concerned; as.numeric() on it would silently return the level index.
    if (is.factor(score.thresh)) {
        score.thresh <- as.character(score.thresh)
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
