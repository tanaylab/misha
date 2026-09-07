create_isolated_test_db()

test_that("the potts params handler validates and normalises", {
    m <- potts_ref_model(W = 5L, npair_mode = "sparse", seed = 89L)

    p <- misha:::.vtrack_params_potts("potts", m, list())
    expect_equal(dim(p$e), c(5L, 4L))
    expect_equal(colnames(p$e), c("A", "C", "G", "T"))
    expect_equal(dim(p$J), c(nrow(m$pairs), 16L))
    expect_equal(dim(p$pairs), dim(m$pairs))
    expect_true(is.integer(p$pairs))
    expect_equal(p$intercept, m$intercept)
    expect_true(p$bidirect)
    expect_true(p$extend)

    # the flat J is column-major per block, so J[[k]][a, b] lands at (b-1)*4+a
    k <- 1L
    expect_equal(p$J[k, (3L - 1L) * 4L + 2L], m$J[[k]][2L, 3L])
})

test_that("the potts params handler tolerates a whole motifmodel Potts", {
    skip_if_not(requireNamespace("motifmodel", quietly = TRUE))
    m <- potts_ref_model(W = 5L, npair_mode = "full", seed = 97L)
    mm <- motifmodel::Potts(m$e, m$J, m$pairs, intercept = m$intercept)

    # width, pair_strength, attr and link ride along and are ignored
    expect_silent(p <- misha:::.vtrack_params_potts("potts", mm, list()))
    expect_equal(nrow(p$e), 5L)
})

test_that("the potts params handler rejects what it should", {
    m <- potts_ref_model(W = 5L, npair_mode = "full", seed = 101L)

    expect_error(misha:::.vtrack_params_potts("potts", list(J = m$J), list()), "'e'")

    bad <- m
    bad$width <- 7L
    expect_error(misha:::.vtrack_params_potts("potts", bad, list()), "width")

    bad <- m
    bad$pairs[1L, ] <- c(3L, 2L)
    expect_error(misha:::.vtrack_params_potts("potts", bad, list()), "lower position first")

    bad <- m
    bad$pairs[2L, ] <- bad$pairs[1L, ]
    expect_error(misha:::.vtrack_params_potts("potts", bad, list()), "twice")

    bad <- m
    bad$J <- bad$J[-1L]
    expect_error(misha:::.vtrack_params_potts("potts", bad, list()), "tables")

    # spatial parameters are not part of the potts family and must be named
    expect_error(
        misha:::.vtrack_params_potts("potts", NULL, c(m, list(spat_factor = c(1, 2)))),
        "spat_factor"
    )
    # and potts.count has no default threshold
    expect_error(misha:::.vtrack_params_potts("potts.count", m, list()), "score.thresh")
})
