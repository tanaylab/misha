create_isolated_test_db()

test_that("the pure-R oracle agrees with motifmodel::motif_score()", {
    skip_if_not(requireNamespace("motifmodel", quietly = TRUE) &&
        is.function(getExportedValue("motifmodel", "motif_score")))

    for (mode in c("full", "sparse", "none")) {
        m <- potts_ref_model(W = 8L, npair_mode = mode, seed = 7L)
        mm <- motifmodel::Potts(m$e, m$J, if (nrow(m$pairs)) m$pairs else NULL,
            intercept = m$intercept
        )
        set.seed(11L)
        seqs <- vapply(1:20, function(i) {
            paste(sample(POTTS_BASES, 8L, replace = TRUE), collapse = "")
        }, character(1))

        ref <- vapply(seqs, potts_ref_window, numeric(1), model = m)
        got <- motifmodel::motif_score(mm, seqs)

        # motifmodel's Potts() projects into the zero-sum gauge, which changes
        # the parameters but not the function, so the SCORES must match even
        # though the tables do not. That is the invariant worth pinning.
        expect_equal(as.numeric(got), as.numeric(ref),
            tolerance = 1e-8,
            info = paste("npair_mode =", mode)
        )
    }
})

test_that("the oracle's RC transform is self-consistent", {
    for (mode in c("full", "sparse", "none")) {
        m <- potts_ref_model(W = 6L, npair_mode = mode, seed = 3L)
        rc <- potts_ref_rc(m)
        set.seed(5L)
        seqs <- vapply(1:50, function(i) {
            paste(sample(POTTS_BASES, 6L, replace = TRUE), collapse = "")
        }, character(1))
        for (s in seqs) {
            expect_equal(potts_ref_window(s, rc),
                potts_ref_window(grevcomp(s), m),
                tolerance = 1e-12,
                ignore_attr = TRUE,
                info = paste(mode, s)
            )
        }
        # and it is an involution
        expect_equal(potts_ref_rc(rc)$e, m$e, tolerance = 1e-12)
        expect_equal(potts_ref_rc(rc)$intercept, m$intercept)
    }
})
