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

test_that("prior = 'marginal' makes unobserved contexts return per-bin base composition", {
    gdb.init_examples()
    if ("g_frac" %in% gvtrack.ls()) gvtrack.rm("g_frac")
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")

    test_intervals <- gintervals(1, 0, 100000)

    model <- gsynth.train(
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200,
        prior = "marginal",
        pseudocount = 1
    )

    counts_b1 <- model$model_data$counts[[1]]
    cdf_b1 <- model$model_data$cdf[[1]]
    row_sums <- rowSums(counts_b1)
    unobs_idx <- which(row_sums == 0)

    skip_if(length(unobs_idx) == 0L, "no unobserved context in bin 1 for this fixture")

    ctx <- unobs_idx[1]
    p_obs <- diff(c(0, cdf_b1[ctx, ]))

    pi <- model$prior[1, ]
    expect_equal(p_obs, pi, tolerance = 1e-5)

    sums <- colSums(counts_b1)
    expected_pi <- sums / sum(sums)
    expect_equal(pi, expected_pi, tolerance = 1e-9)

    gvtrack.rm("g_frac")
})

test_that("explicit length-4 prior overrides unobserved context probabilities", {
    gdb.init_examples()
    if ("g_frac" %in% gvtrack.ls()) gvtrack.rm("g_frac")
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")

    test_intervals <- gintervals(1, 0, 50000)
    target <- c(0.10, 0.40, 0.40, 0.10)

    model <- gsynth.train(
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200,
        prior = target,
        pseudocount = 1
    )

    expect_equal(as.numeric(model$prior[1, ]), target, tolerance = 1e-9)

    found <- FALSE
    for (b in seq_len(model$total_bins)) {
        counts_b <- model$model_data$counts[[b]]
        cdf_b <- model$model_data$cdf[[b]]
        unobs <- which(rowSums(counts_b) == 0)
        if (length(unobs) == 0L) next
        ctx <- unobs[1]
        p_obs <- diff(c(0, cdf_b[ctx, ]))
        expect_equal(p_obs, target, tolerance = 1e-5)
        found <- TRUE
        break
    }
    skip_if(!found, "no unobserved context in any bin for this fixture")

    gvtrack.rm("g_frac")
})

test_that("prior = NULL preserves uniform 1/4 at unobserved contexts (regression)", {
    gdb.init_examples()
    if ("g_frac" %in% gvtrack.ls()) gvtrack.rm("g_frac")
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gsynth.train(
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200,
        prior = NULL,
        pseudocount = 1
    )

    counts_b1 <- model$model_data$counts[[1]]
    cdf_b1 <- model$model_data$cdf[[1]]
    unobs_idx <- which(rowSums(counts_b1) == 0)
    skip_if(length(unobs_idx) == 0L, "no unobserved context in bin 1 for this fixture")

    ctx <- unobs_idx[1]
    p_obs <- diff(c(0, cdf_b1[ctx, ]))
    expect_equal(p_obs, rep(0.25, 4L), tolerance = 1e-5)

    gvtrack.rm("g_frac")
})

test_that("prior validation rejects malformed input", {
    expect_error(.gsynth_resolve_prior_arg("bogus", 4L), "Unknown prior")
    expect_error(
        .gsynth_resolve_prior_arg(c(-0.1, 0.5, 0.3, 0.3), 4L),
        "non-negative"
    )
    expect_error(.gsynth_resolve_prior_arg(c(0, 0, 0, 0), 4L), "sum to > 0")
    # Mismatched matrix dims
    expect_error(
        .gsynth_resolve_prior_arg(matrix(0.25, nrow = 3, ncol = 4), 4L),
        "must be 4 x 4"
    )
    # Length-4 named vector with wrong names
    expect_error(
        .gsynth_resolve_prior_arg(c(X = 0.25, Y = 0.25, Z = 0.25, W = 0.25), 4L),
        "names A, C, G, T"
    )
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
