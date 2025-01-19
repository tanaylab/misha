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
