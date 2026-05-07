# Oracle helpers for PWM gradient vtrack tests.
#
# These are pure R reference implementations for the linearized gradient (B)
# and ISM gradient (A) defined in the PWM gradient vtrack design.
# They are intentionally simple (correctness over speed) and reuse the same
# prior + per-row normalization as `manual_pwm_scores_single_strand` in
# helper-pwm.R, so engine outputs can be validated against these in later
# tasks.

# Map base character -> 1-based index (A=1, C=2, G=3, T=4). NA for non-ACGT.
.base_to_idx <- function(base) {
    switch(base,
        "A" = 1L,
        "C" = 2L,
        "G" = 3L,
        "T" = 4L,
        NA_integer_
    )
}

# Complement of an index (1<->4, 2<->3).
.complement_idx <- function(idx) c(4L, 3L, 2L, 1L)[idx]

# Reverse-complement of a sequence string.
.rc_seq <- function(seq) {
    comp <- chartr("ACGT", "TGCA", seq)
    paste(rev(strsplit(comp, "")[[1]]), collapse = "")
}

# Apply prior + per-row normalization, then take log. Matches the path used in
# manual_pwm_scores_single_strand (same prior treatment).
.normalize_pssm <- function(pssm, prior = 0.01) {
    if (prior > 0) {
        normalized <- matrix(0, nrow = nrow(pssm), ncol = ncol(pssm))
        for (i in seq_len(nrow(pssm))) {
            row <- pssm[i, ] + prior
            normalized[i, ] <- row / sum(row)
        }
        pssm <- normalized
    }
    pssm
}

# Per-anchor scores in a window.
#   $fwd[a] : log-prob of fwd PSSM applied to seq[a..a+L-1]
#   $rc[a]  : log-prob of fwd PSSM applied to RC of seq[a..a+L-1]
#             (== rc PSSM applied to seq[a..a+L-1])
#   $comb[a]: log(exp(fwd[a]) + exp(rc[a]))
.oracle_anchor_scores <- function(seq, pssm, prior = 0.01) {
    L <- nrow(pssm)
    n <- nchar(seq)
    stopifnot(n >= L)

    # Fwd uses the same code path as manual_pwm_scores_single_strand.
    fwd <- manual_pwm_scores_single_strand(seq, pssm, prior = prior)

    # Rc: apply fwd PSSM to RC of each window. Easiest: build rc PSSM
    # = pssm[L:1, c(4,3,2,1)] and run on the original seq. This is
    # equivalent (since prior + normalization is row-wise) to running the
    # fwd PSSM on the RC of the seq.
    rc_pssm <- pssm[L:1, c(4L, 3L, 2L, 1L), drop = FALSE]
    if (!is.null(colnames(pssm))) colnames(rc_pssm) <- colnames(pssm)
    rc <- manual_pwm_scores_single_strand(seq, rc_pssm, prior = prior)

    comb <- vapply(seq_along(fwd), function(i) {
        log_sum_exp(c(fwd[i], rc[i]))
    }, numeric(1))

    list(fwd = fwd, rc = rc, comb = comb)
}

# Aggregate a vector of per-anchor scores via "lse" or "max".
.aggregate_scores <- function(s, aggregate) {
    aggregate <- match.arg(aggregate, c("lse", "max"))
    if (aggregate == "max") {
        return(max(s))
    }
    log_sum_exp(s)
}

# Pick the per-anchor score vector to use given (bidirect, strand).
#   bidirect = TRUE  -> comb
#   bidirect = FALSE, strand = +1 -> fwd
#   bidirect = FALSE, strand = -1 -> rc
.pick_scores <- function(scores, bidirect, strand) {
    if (bidirect) {
        return(scores$comb)
    }
    if (strand == 1L || strand == 1) {
        return(scores$fwd)
    }
    if (strand == -1L || strand == -1) {
        return(scores$rc)
    }
    stop("strand must be +1 or -1 when bidirect = FALSE")
}

# Apply a per-anchor spatial log-factor to a vector of per-anchor scores.
# `spat_factor` is a vector of multiplicative factors (in linear scale); we
# convert to log and add. `spat_bin_size` controls how many anchors share each
# factor (matches the engine's get_spatial_log_factor convention). When
# `spat_factor` is NULL or empty, scores are returned unchanged.
.apply_spat <- function(scores, spat_factor, spat_bin_size = 1L) {
    if (is.null(spat_factor) || length(spat_factor) == 0) {
        return(scores)
    }
    n <- length(scores)
    bin_sz <- max(1L, as.integer(spat_bin_size))
    bins <- ((seq_len(n) - 1L) %/% bin_sz) + 1L
    bins <- pmin(bins, length(spat_factor))
    log_factors <- log(pmax(1e-30, spat_factor))[bins]
    scores + log_factors
}

# Linearized gradient at p = interval start, head anchor, column 0 (= base b_p).
#
# Returns a single non-negative scalar. For an all-N seq or otherwise invalid
# head base, returns 0.
manual_pwm_grad <- function(pssm, seq, aggregate = c("lse", "max"),
                            bidirect = TRUE, strand = 1L, prior = 0.01,
                            spat_factor = NULL, spat_bin_size = 1L) {
    aggregate <- match.arg(aggregate)
    L <- nrow(pssm)

    # Normalized log PSSM (used for diff lookups).
    M <- log(.normalize_pssm(pssm, prior = prior))

    # Per-anchor scores (with spatial weighting applied).
    scores <- .oracle_anchor_scores(seq, pssm, prior = prior)
    scores$fwd <- .apply_spat(scores$fwd, spat_factor, spat_bin_size)
    scores$rc <- .apply_spat(scores$rc, spat_factor, spat_bin_size)
    # Recompute comb after spatial weighting (spatial is the same for fwd and
    # rc at a given anchor, so log(exp(fwd_spat) + exp(rc_spat)) = log(exp(fwd)
    # + exp(rc)) + spat, but recomputing keeps the code symmetric and handles
    # the all-Inf cases cleanly).
    scores$comb <- vapply(seq_along(scores$fwd), function(i) {
        log_sum_exp(c(scores$fwd[i], scores$rc[i]))
    }, numeric(1))

    s <- .pick_scores(scores, bidirect, strand)
    f <- .aggregate_scores(s, aggregate)

    # Head base.
    b_p_char <- substr(seq, 1L, 1L)
    b_p <- .base_to_idx(b_p_char)
    if (is.na(b_p)) {
        return(0)
    }

    # Diffs at column 0 of head anchor for each strand.
    diff_fwd <- M[1L, b_p] - min(M[1L, ])
    diff_rc <- M[L, .complement_idx(b_p)] - min(M[L, ])

    if (!bidirect) {
        diff <- if (strand == 1L || strand == 1) diff_fwd else diff_rc
        if (aggregate == "max") {
            # 1-based index of head anchor is 1.
            return(if (which.max(s) == 1L) diff else 0)
        }
        # LSE: weight by softmax of head anchor.
        w_p <- exp(s[1L] - f)
        return(w_p * diff)
    }

    # Bidirect: per-strand softmax at the head anchor, combined diff.
    fwd_p <- scores$fwd[1L]
    rc_p <- scores$rc[1L]
    comb_p <- scores$comb[1L]
    sm_fwd <- exp(fwd_p - comb_p)
    sm_rc <- 1 - sm_fwd
    combined_diff <- sm_fwd * diff_fwd + sm_rc * diff_rc

    if (aggregate == "max") {
        return(if (which.max(s) == 1L) combined_diff else 0)
    }
    w_p_comb <- exp(comb_p - f)
    w_p_comb * combined_diff
}

# In-silico mutagenesis at p = interval start.
# g = f_actual - min over b' != b_p of f(seq with seq[1] := b').
# If b_p is not ACGT, returns 0.
manual_pwm_grad_ism <- function(pssm, seq, aggregate = c("lse", "max"),
                                bidirect = TRUE, strand = 1L, prior = 0.01,
                                spat_factor = NULL, spat_bin_size = 1L) {
    aggregate <- match.arg(aggregate)

    aggregate_for_seq <- function(s_seq) {
        sc <- .oracle_anchor_scores(s_seq, pssm, prior = prior)
        sc$fwd <- .apply_spat(sc$fwd, spat_factor, spat_bin_size)
        sc$rc <- .apply_spat(sc$rc, spat_factor, spat_bin_size)
        sc$comb <- vapply(seq_along(sc$fwd), function(i) {
            log_sum_exp(c(sc$fwd[i], sc$rc[i]))
        }, numeric(1))
        s <- .pick_scores(sc, bidirect, strand)
        .aggregate_scores(s, aggregate)
    }

    f_actual <- aggregate_for_seq(seq)

    b_p_char <- substr(seq, 1L, 1L)
    b_p <- .base_to_idx(b_p_char)
    if (is.na(b_p)) {
        return(0)
    }

    bases <- c("A", "C", "G", "T")
    others <- bases[-b_p]
    f_alt <- vapply(others, function(b) {
        modified <- paste0(b, substr(seq, 2L, nchar(seq)))
        aggregate_for_seq(modified)
    }, numeric(1))

    f_actual - min(f_alt)
}
