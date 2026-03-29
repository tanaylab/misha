remove_all_vtracks <- function() {
    vtracks <- gvtrack.ls()
    for (vtrack in vtracks) {
        do.call(gvtrack.rm, list(vtrack = vtrack))
    }
}

# Helper function to compute exact log-sum-exp as done in C++
log_sum_exp <- function(x) {
    x <- x[is.finite(x)] # Remove any NA/NaN values
    if (length(x) == 0) {
        return(-Inf)
    }
    if (length(x) == 1) {
        return(x)
    }

    # Sort in descending order to match C++ algorithm
    x <- sort(x, decreasing = TRUE)
    sum <- x[1]

    # Iteratively add using same approach as C++
    for (i in 2:length(x)) {
        if (sum > x[i]) {
            sum <- sum + log1p(exp(x[i] - sum))
        } else {
            sum <- x[i] + log1p(exp(sum - x[i]))
        }
    }
    return(sum)
}

# Helper function for manual PWM calculation matching C++ exactly
manual_pwm_scores_single_strand <- function(seq, pssm, prior = 0.01) {
    motif_length <- nrow(pssm)
    scores <- numeric()

    # Add prior and normalize exactly as in C++
    if (prior > 0) {
        normalized_pssm <- matrix(0, nrow = nrow(pssm), ncol = ncol(pssm))
        for (i in 1:nrow(pssm)) {
            row <- pssm[i, ] + prior
            sum <- sum(row)
            normalized_pssm[i, ] <- row / sum
        }
        pssm <- normalized_pssm
    }

    # Calculate scores for each possible position
    for (i in 1:(nchar(seq) - motif_length + 1)) {
        subseq <- substr(seq, i, i + motif_length - 1)
        score <- 0
        valid_score <- TRUE

        for (j in 1:motif_length) {
            base <- substr(subseq, j, j)
            base_idx <- switch(base,
                "A" = 1,
                "C" = 2,
                "G" = 3,
                "T" = 4,
                NA
            )

            if (is.na(base_idx)) {
                valid_score <- FALSE
                break
            }

            prob <- pssm[j, base_idx]
            if (prob == 0) {
                score <- -Inf
                valid_score <- FALSE
                break
            } else {
                score <- score + log(prob)
            }
        }

        if (!valid_score) score <- -Inf
        scores <- c(scores, score)
    }
    return(scores)
}


# Helper function to create a standard PSSM for testing
create_test_pssm <- function() {
    pssm <- matrix(c(
        1.0, 0.0, 0.0, 0.0, # Only A
        0.0, 1.0, 0.0, 0.0 # Only C
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")
    return(pssm)
}

# Helper function to compute PWM edit distance manually (for testing)
# When scan_all = TRUE, returns the minimum edits across every start position in seq.
# When scan_all = FALSE, seq is treated as a single motif-length window.
manual_pwm_edit_distance <- function(seq, pssm, threshold, max_edits = NULL, scan_all = TRUE) {
    motif_len <- nrow(pssm)

    if (nchar(seq) < motif_len) {
        return(NA_real_)
    }

    log_pssm <- log(pssm)
    col_max <- apply(log_pssm, 1, max)

    score_window <- function(window_seq) {
        adjusted_score <- 0
        mandatory_edits <- 0
        gains <- numeric(0)

        for (i in seq_len(motif_len)) {
            base <- substr(window_seq, i, i)
            base_idx <- switch(base,
                "A" = 1,
                "C" = 2,
                "G" = 3,
                "T" = 4,
                NA
            )

            if (is.na(base_idx)) {
                base_score <- min(log_pssm[i, ])
            } else {
                base_score <- log_pssm[i, base_idx]
            }

            if (!is.finite(base_score)) {
                mandatory_edits <- mandatory_edits + 1
                adjusted_score <- adjusted_score + col_max[i]
            } else {
                adjusted_score <- adjusted_score + base_score
                gains <- c(gains, col_max[i] - base_score)
            }
        }

        deficit <- threshold - adjusted_score
        if (deficit <= 0) {
            if (!is.null(max_edits) && mandatory_edits > max_edits) {
                return(NA_real_)
            }
            return(mandatory_edits)
        }

        total_max_gain <- sum(col_max) - adjusted_score
        if (total_max_gain < deficit - 1e-12) {
            return(NA_real_)
        }

        gains_sorted <- sort(gains, decreasing = TRUE)

        if (!is.null(max_edits)) {
            remaining_budget <- max_edits - mandatory_edits
            if (remaining_budget < 0) {
                return(NA_real_)
            }
            if (remaining_budget < length(gains_sorted)) {
                gains_sorted <- gains_sorted[seq_len(remaining_budget)]
            }
        }

        acc <- 0
        edits <- mandatory_edits
        for (gain in gains_sorted) {
            edits <- edits + 1
            acc <- acc + gain
            if (acc >= deficit) {
                return(edits)
            }
        }

        return(NA_real_)
    }

    if (!scan_all) {
        return(score_window(seq))
    }

    seq_len_total <- nchar(seq)
    best <- NA_real_
    for (start_idx in seq_len(seq_len_total - motif_len + 1)) {
        window_seq <- substr(seq, start_idx, start_idx + motif_len - 1)
        cand <- score_window(window_seq)
        if (is.na(best) || (!is.na(cand) && cand < best)) {
            best <- cand
        }
    }
    best
}

# Helper function to compute PWM edit distance in "below" direction manually.
# Returns the minimum number of substitutions to bring the best window's score
# BELOW the threshold (i.e., score <= threshold).
# When scan_all = TRUE, returns the minimum edits across every start position.
# When scan_all = FALSE, seq is treated as a single motif-length window.
manual_pwm_edit_distance_below <- function(seq, pssm, threshold, max_edits = NULL, scan_all = TRUE) {
    motif_len <- nrow(pssm)

    if (nchar(seq) < motif_len) {
        return(NA_real_)
    }

    log_pssm <- log(pssm)
    col_min <- apply(log_pssm, 1, min)

    score_window <- function(window_seq) {
        current_score <- 0
        has_neg_inf <- FALSE
        losses <- numeric(0)

        for (i in seq_len(motif_len)) {
            base <- substr(window_seq, i, i)
            base_idx <- switch(base,
                "A" = 1,
                "C" = 2,
                "G" = 3,
                "T" = 4,
                NA
            )

            if (is.na(base_idx)) {
                base_score <- min(log_pssm[i, ])
            } else {
                base_score <- log_pssm[i, base_idx]
            }

            if (!is.finite(base_score)) {
                # Score is -Inf: total score is -Inf, already below any threshold
                has_neg_inf <- TRUE
                break
            }

            current_score <- current_score + base_score
            loss <- base_score - col_min[i]
            losses <- c(losses, loss)
        }

        if (has_neg_inf) {
            # Score is -Inf, which is <= any finite threshold
            return(0)
        }

        # surplus = how much the current score exceeds the threshold
        surplus <- current_score - threshold
        if (surplus <= 0) {
            # Already at or below threshold
            return(0)
        }

        # Check if even switching all positions to worst can cover the surplus
        total_possible_loss <- sum(losses)
        if (total_possible_loss < surplus - 1e-12) {
            return(NA_real_)
        }

        # Sort losses descending and greedily accumulate
        losses_sorted <- sort(losses, decreasing = TRUE)

        if (!is.null(max_edits)) {
            if (max_edits < length(losses_sorted)) {
                losses_sorted <- losses_sorted[seq_len(max_edits)]
            }
        }

        acc <- 0
        edits <- 0
        for (loss in losses_sorted) {
            edits <- edits + 1
            acc <- acc + loss
            if (acc >= surplus) {
                return(edits)
            }
        }

        return(NA_real_)
    }

    if (!scan_all) {
        return(score_window(seq))
    }

    seq_len_total <- nchar(seq)
    best <- NA_real_
    for (start_idx in seq_len(seq_len_total - motif_len + 1)) {
        window_seq <- substr(seq, start_idx, start_idx + motif_len - 1)
        cand <- score_window(window_seq)
        if (is.na(best) || (!is.na(cand) && cand < best)) {
            best <- cand
        }
    }
    best
}
