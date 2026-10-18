create_isolated_test_db()

test_that("reach: a bin [s, e) of each sequence vtrack reads [s + lo, e + hi)", {
    define_variant_vtracks()
    pssm <- variant_test_pssm()
    gvtrack.create("r_noext", NULL, "pwm.max", params = list(pssm = pssm, bidirect = TRUE, extend = FALSE, prior = 0.01))
    gvtrack.create("r_kmer_fwd", NULL, "kmer.count", params = list(kmer = "ACG", extend = TRUE, strand = 1))
    gvtrack.create("r_kmer_noext", NULL, "kmer.count", params = list(kmer = "ACG", extend = FALSE, strand = 0))
    reach <- function(expr) misha:::.gvariants_reach(expr)
    expect_equal(reach("v_pwm"), cbind(lo = -5, hi = 3 + 3))
    expect_equal(reach("r_noext"), cbind(lo = 0, hi = 0))
    expect_equal(reach("v_potts"), cbind(lo = 4, hi = 4 + 5))
    expect_equal(reach("v_kmer"), cbind(lo = -2 - 2, hi = 2 + 2))
    expect_equal(reach("r_kmer_fwd"), cbind(lo = 0, hi = 2))
    expect_equal(reach("r_kmer_noext"), cbind(lo = 0, hi = 0))
    expect_equal(reach("v_masked"), cbind(lo = -3, hi = 3))
    expect_equal(reach("v_ed"), cbind(lo = 0, hi = 4 + 1 - 1))
    expect_equal(nrow(reach("v_pwm + 2 * v_potts + test.fixedbin")), 2L)
})

test_that("reach: an expression with no sequence vtrack is an error", {
    expect_error(misha:::.gvariants_reach("test.fixedbin * 2"), "no sequence virtual track")
})

test_that("engine: a substituted row equals the same bin on a really edited genome", {
    local_db_state()
    seq <- variant_test_seq()
    db <- make_seq_db(seq)
    subs <- data.frame(pos = c(100, 101, 150, 151, 152, 200), base = c("A", "C", "G", "T", "A", "C"))
    subs <- subs[toupper(substr(seq, subs$pos + 1, subs$pos + 1)) != subs$base, ]
    # A for loop, not lapply: make_seq_db() deletes its database when the
    # calling frame exits, and an lapply FUN frame exits at once.
    edited <- list()
    for (s in seq_len(nrow(subs))) {
        edited[[s]] <- make_seq_db(substitute_bases(seq, subs$pos[s], subs$base[s]))
    }
    bins <- 90:210
    # Bin-major: consecutive rows share a bin and differ in the substitution,
    # which is where a cache that ignores the substitution would leak.
    rows <- expand.grid(sub = seq_len(nrow(subs)), bin = bins)
    for (vt in VARIANT_VTRACKS) {
        suppressMessages(gsetroot(db))
        define_variant_vtracks()
        got <- misha:::.gvariants_eval(
            vt, rep("chr1", nrow(rows)), rows$bin, rows$bin + 1,
            subs$pos[rows$sub], match(subs$base[rows$sub], VARIANT_BASES) - 1L
        )
        for (s in seq_len(nrow(subs))) {
            expect_equal(got[rows$sub == s], values_on_db(edited[[s]], vt, bins),
                info = sprintf("%s, substitution %d", vt, s)
            )
        }
    }
})

test_that("engine: a row with no substitution after substituted rows reads the reference", {
    local_db_state()
    seq <- variant_test_seq(seed = 8L)
    db <- make_seq_db(seq)
    b <- 120
    alt <- setdiff(VARIANT_BASES, toupper(substr(seq, b + 1, b + 1)))[1]
    for (vt in VARIANT_VTRACKS) {
        want <- values_on_db(db, vt, b)
        got <- misha:::.gvariants_eval(
            vt, rep("chr1", 3), rep(b, 3), rep(b + 1, 3), rep(b, 3),
            c(match(alt, VARIANT_BASES) - 1L, -1L, -1L)
        )
        expect_equal(got[2:3], rep(want, 2), info = vt)
    }
})

test_that("engine: unsorted rows are an error", {
    local_db_state()
    db <- make_seq_db(variant_test_seq())
    suppressMessages(gsetroot(db))
    define_variant_vtracks()
    expect_error(misha:::.gvariants_eval("v_pwm", c("chr1", "chr1"), c(50, 40), c(51, 41), c(50, 40), c(0L, 0L)), "sorted")
})
