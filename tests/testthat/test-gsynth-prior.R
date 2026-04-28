# Tests for gsynth.train(prior = ...) — informative Dirichlet prior added in misha 5.7.0.

test_that("gsynth.train accepts prior = 'marginal' (default) and stores resolved prior", {
    gdb.init_examples()
    if ("g_frac" %in% gvtrack.ls()) gvtrack.rm("g_frac")
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gsynth.train(
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )

    expect_equal(model$prior_mode, "marginal")
    expect_true(is.matrix(model$prior))
    expect_equal(dim(model$prior), c(model$total_bins, 4L))
    row_sums <- rowSums(model$prior)
    # Either sums to 1 (observed) or exactly 0.25 each (uniform fallback)
    expect_true(all(abs(row_sums - 1) < 1e-9))

    gvtrack.rm("g_frac")
})

test_that("gsynth.train with prior = NULL produces uniform 1/4 rows", {
    gdb.init_examples()
    if ("g_frac" %in% gvtrack.ls()) gvtrack.rm("g_frac")
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gsynth.train(
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200,
        prior = NULL
    )

    expect_equal(model$prior_mode, "uniform")
    expect_true(all(abs(model$prior - 0.25) < 1e-12))

    gvtrack.rm("g_frac")
})
