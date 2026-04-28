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

test_that("gsynth.save/load round-trips the prior matrix", {
    gdb.init_examples()
    if ("g_frac" %in% gvtrack.ls()) gvtrack.rm("g_frac")
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gsynth.train(
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200,
        prior = c(A = 0.4, C = 0.1, G = 0.1, T = 0.4)
    )

    out_dir <- file.path(tempdir(), "test_prior_roundtrip")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
    gsynth.save(model, out_dir)
    loaded <- gsynth.load(out_dir)

    expect_equal(loaded$prior_mode, model$prior_mode)
    expect_equal(loaded$prior, model$prior, tolerance = 1e-9)

    gvtrack.rm("g_frac")
})

test_that("loading a .gsm without a prior block defaults to uniform prior", {
    gdb.init_examples()
    if ("g_frac" %in% gvtrack.ls()) gvtrack.rm("g_frac")
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gsynth.train(
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )

    out_dir <- file.path(tempdir(), "test_legacy_no_prior")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
    gsynth.save(model, out_dir)

    # Strip prior block from metadata to simulate an old .gsm
    meta <- yaml::read_yaml(file.path(out_dir, "metadata.yaml"))
    meta$prior <- NULL
    meta$prior_mode <- NULL
    yaml::write_yaml(meta, file.path(out_dir, "metadata.yaml"), precision = 15)

    loaded <- gsynth.load(out_dir)
    expect_equal(loaded$prior_mode, "uniform")
    expect_true(all(abs(loaded$prior - 0.25) < 1e-12))

    gvtrack.rm("g_frac")
})
