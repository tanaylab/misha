# The pure-R Potts oracle. Deliberately slow, deliberately non-reusing: it is
# the ground truth for src/PottsModel.h, so it must not share a line of logic
# with it. Do not "optimise" it and do not fix it to match the kernel.
#
# Base order is A, C, G, T everywhere. A model here is a plain list, which is
# also the shape gvtrack.create() takes as `params` and the shape a fitted
# Potts model already has.

POTTS_BASES <- c("A", "C", "G", "T")

# A reproducible test model. `npair_mode` is "full" (every pair), "sparse"
# (every third pair) or "none" (order 1).
potts_ref_model <- function(W, npair_mode = c("full", "sparse", "none"), seed = 1L) {
    npair_mode <- match.arg(npair_mode)
    set.seed(seed)
    e <- matrix(round(rnorm(W * 4), 4),
        nrow = W, ncol = 4,
        dimnames = list(NULL, POTTS_BASES)
    )
    full <- if (W >= 2L) t(utils::combn(W, 2)) else matrix(integer(0), 0L, 2L)
    pairs <- switch(npair_mode,
        full = full,
        sparse = full[seq(1L, nrow(full), by = 3L), , drop = FALSE],
        none = matrix(integer(0), 0L, 2L)
    )
    storage.mode(pairs) <- "integer"
    dimnames(pairs) <- NULL
    J <- lapply(seq_len(nrow(pairs)), function(k) matrix(round(rnorm(16), 4), 4L, 4L))
    list(e = e, J = J, pairs = pairs, intercept = round(rnorm(1), 4))
}

# Score ONE window of exactly nrow(e) bases. NA for any non-ACGT base, for NA,
# and for the wrong width - a Potts has no prior, so an N has no energy.
potts_ref_window <- function(win, model) {
    W <- nrow(model$e)
    if (length(win) != 1L || is.na(win) || nchar(win) != W) {
        return(NA_real_)
    }
    codes <- match(strsplit(toupper(win), "", fixed = TRUE)[[1L]], POTTS_BASES)
    if (anyNA(codes)) {
        return(NA_real_)
    }
    s <- model$intercept
    for (i in seq_len(W)) {
        s <- s + model$e[i, codes[i]]
    }
    for (k in seq_len(nrow(model$pairs))) {
        i <- model$pairs[k, 1L]
        j <- model$pairs[k, 2L]
        s <- s + model$J[[k]][codes[i], codes[j]]
    }
    s
}

# The reverse-complement twin, from the transform in the design doc:
#   e'[i][b]           = e[W+1-i][comp(b)]
#   pairs'             = (W+1-p2, W+1-p1), re-sorted
#   J'_{(u,v)}[a][b]   = J_k[comp(b)][comp(a)]
# so that potts_ref_window(w, potts_ref_rc(m)) == potts_ref_window(revcomp(w), m).
potts_ref_rc <- function(model) {
    W <- nrow(model$e)
    comp <- c(4L, 3L, 2L, 1L) # A<->T, C<->G in A,C,G,T order
    e2 <- matrix(0, nrow = W, ncol = 4L, dimnames = list(NULL, POTTS_BASES))
    for (i in seq_len(W)) {
        for (b in 1:4) {
            e2[i, b] <- model$e[W + 1L - i, comp[b]]
        }
    }
    np <- nrow(model$pairs)
    if (np) {
        u <- W + 1L - model$pairs[, 2L]
        v <- W + 1L - model$pairs[, 1L]
        ord <- order(u, v)
        pairs2 <- cbind(u, v)[ord, , drop = FALSE]
        storage.mode(pairs2) <- "integer"
        dimnames(pairs2) <- NULL
        J2 <- lapply(ord, function(k) {
            Jk <- model$J[[k]]
            M <- matrix(0, 4L, 4L)
            for (a in 1:4) {
                for (b in 1:4) {
                    M[a, b] <- Jk[comp[b], comp[a]]
                }
            }
            M
        })
    } else {
        pairs2 <- matrix(integer(0), 0L, 2L)
        J2 <- list()
    }
    list(e = e2, J = J2, pairs = pairs2, intercept = model$intercept)
}

# Per-anchor values over `seq`, in the original genome's forward orientation.
# `union` says how the two strands combine at one anchor: "lse" for potts,
# potts.max and potts.count; "max" for potts.max.pos. See the design doc's
# strand table - this asymmetry is inherited from the pwm family on purpose.
potts_ref_anchors <- function(seq, model, bidirect = TRUE, strand = 1L,
                              union = c("lse", "max")) {
    union <- match.arg(union)
    W <- nrow(model$e)
    n <- nchar(seq) - W + 1L
    if (n < 1L) {
        return(numeric(0))
    }
    rc <- potts_ref_rc(model)
    vapply(seq_len(n), function(i) {
        win <- substr(seq, i, i + W - 1L)
        fwd <- potts_ref_window(win, model)
        rev <- potts_ref_window(win, rc)
        if (!bidirect) {
            return(if (strand == -1L) rev else fwd)
        }
        if (is.na(fwd) || is.na(rev)) {
            return(NA_real_)
        }
        if (union == "max") max(fwd, rev) else log_sum_exp(c(fwd, rev))
    }, numeric(1))
}
