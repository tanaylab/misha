# Helpers for test-gextract-variants.R.

VARIANT_BASES <- c("A", "C", "G", "T")

variant_test_pssm <- function() {
    matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.4, 0.1, 0.4, 0.1
    ), ncol = 4, byrow = TRUE, dimnames = list(NULL, VARIANT_BASES))
}

# One vtrack of each sequence family, with shifts. Vtracks belong to the
# working database, so call this again after every gsetroot().
define_variant_vtracks <- function() {
    pssm <- variant_test_pssm()
    gvtrack.create("v_pwm", NULL, "pwm.max", params = list(pssm = pssm, bidirect = TRUE, extend = TRUE, prior = 0.01))
    gvtrack.iterator("v_pwm", sshift = -5, eshift = 3)
    m <- potts_ref_model(W = 6L, npair_mode = "full", seed = 3L)
    gvtrack.create("v_potts", NULL, "potts.max", params = c(m, list(bidirect = FALSE, strand = -1)))
    gvtrack.iterator("v_potts", sshift = 4, eshift = 4)
    gvtrack.create("v_kmer", NULL, "kmer.count", params = list(kmer = "ACG", extend = TRUE, strand = 0))
    gvtrack.iterator("v_kmer", sshift = -2, eshift = 2)
    gvtrack.create("v_masked", NULL, "masked.count")
    gvtrack.iterator("v_masked", sshift = -3, eshift = 3)
    gvtrack.create("v_ed", NULL, "pwm.edit_distance", params = list(pssm = pssm, score.thresh = -2, max_indels = 1))
}

VARIANT_VTRACKS <- c("v_pwm", "v_potts", "v_kmer", "v_masked", "v_ed")

# A random chromosome with some soft-masked (lowercase) bases.
variant_test_seq <- function(n = 400L, seed = 7L) {
    set.seed(seed)
    paste(sample(c(VARIANT_BASES, tolower(VARIANT_BASES)), n,
        replace = TRUE, prob = c(rep(0.2, 4), rep(0.05, 4))
    ), collapse = "")
}

# A one-chromosome database ("chr1") built from `seq`; returns its root.
make_seq_db <- function(seq, env = parent.frame()) {
    dir <- tempfile("gvariants_")
    dir.create(dir)
    withr::defer(unlink(dir, recursive = TRUE), envir = env)
    fa <- file.path(dir, "g.fa")
    cat(">chr1\n", seq, "\n", sep = "", file = fa)
    db <- file.path(dir, "db")
    suppressMessages(gdb.create(groot = db, fasta = fa, verbose = FALSE))
    db
}

# `seq` with the 0-based positions `pos` replaced by `base`, each keeping the
# case of the base it replaces.
substitute_bases <- function(seq, pos, base) {
    chars <- strsplit(seq, "")[[1]]
    for (i in seq_along(pos)) {
        old <- chars[pos[i] + 1]
        chars[pos[i] + 1] <- if (old == tolower(old)) tolower(base[i]) else toupper(base[i])
    }
    paste(chars, collapse = "")
}

# The value of `vtrack` on 1 bp bins `bins` of chr1 in database `db`.
values_on_db <- function(db, vtrack, bins) {
    suppressMessages(gsetroot(db))
    define_variant_vtracks()
    gextract(vtrack, gintervals("chr1", min(bins), max(bins) + 1), iterator = 1)[[vtrack]][bins - min(bins) + 1]
}
