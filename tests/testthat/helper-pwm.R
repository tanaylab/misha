remove_all_vtracks <- function() {
    vtracks <- gvtrack.ls()
    for (vtrack in vtracks) {
        do.call(gvtrack.rm, list(vtrack = vtrack))
    }
}

# Helper function for computing log sum exp
log_sum_exp <- function(x, na.rm = TRUE) {
    max_x <- max(x, na.rm = TRUE)
    max_x + log(sum(exp(x - max_x), na.rm = TRUE))
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

# Helper function to calculate PWM scores manually
manual_pwm_score <- function(seq, pssm) {
    motif_length <- nrow(pssm)
    score <- 0
    for (i in 1:motif_length) {
        base <- substr(seq, i, i)
        base_idx <- switch(base,
            "A" = 1,
            "C" = 2,
            "G" = 3,
            "T" = 4
        )
        score <- score + log(pssm[i, base_idx])
    }
    return(score)
}

# Helper function to calculate scores for all possible positions
manual_pwm_scores <- function(seq, pssm) {
    motif_length <- nrow(pssm)
    scores <- numeric()
    for (i in 1:(nchar(seq) - motif_length + 1)) {
        subseq <- substr(seq, i, i + motif_length - 1)
        scores <- c(scores, manual_pwm_score(subseq, pssm))
    }
    return(scores)
}
