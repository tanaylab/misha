test_that(".gsynth_build_log_p recovers log probabilities from cdf", {
    gdb.init_examples()

    if ("g_frac" %in% gvtrack.ls()) gvtrack.rm("g_frac")
    if ("c_frac" %in% gvtrack.ls()) gvtrack.rm("c_frac")
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    # Tiny k=2 unstratified model on a small chr1 region.
    model <- gsynth.train(
        intervals = gintervals(1, 0, 50000),
        iterator = 200,
        k = 2L,
        pseudocount = 1
    )

    log_p <- .gsynth_build_log_p(model)

    # One bin (unstratified)
    expect_equal(length(log_p), 1L)
    # 16 contexts (4^2) x 4 bases
    expect_equal(dim(log_p[[1]]), c(16L, 4L))

    # Recover p from cdf, take log, compare against helper.
    cdf <- model$model_data$cdf[[1]]
    p_from_cdf <- cbind(
        cdf[, 1],
        cdf[, 2] - cdf[, 1],
        cdf[, 3] - cdf[, 2],
        cdf[, 4] - cdf[, 3]
    )
    expect_equal(log_p[[1]], log(p_from_cdf), tolerance = 1e-10)

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gsynth.score validates required arguments", {
    gdb.init_examples()

    # Build a tiny model to use as input.
    model <- gsynth.train(
        intervals = gintervals(1, 0, 50000),
        iterator = 200,
        k = 2L
    )

    expect_error(
        gsynth.score(model = "not a model", track = "x"),
        "model must be a gsynth.model"
    )
    expect_error(
        gsynth.score(model = model),
        "track is required"
    )
    expect_error(
        gsynth.score(model = model, track = "x", resolution = -1),
        "resolution must be"
    )
    expect_error(
        gsynth.score(model = model, track = "x", resolution = 1.5),
        "resolution must be"
    )
})

test_that("gsynth.score writes a queryable track end-to-end", {
    gdb.init_examples()
    on.exit(
        {
            if (gtrack.exists("test_score_e2e")) {
                gtrack.rm("test_score_e2e", force = TRUE)
            }
        },
        add = TRUE
    )

    model <- gsynth.train(
        intervals = gintervals(1, 0, 50000),
        iterator = 200,
        k = 2L
    )

    gsynth.score(
        model = model,
        track = "test_score_e2e",
        intervals = gintervals(1, 0, 50000),
        resolution = 200,
        overwrite = TRUE
    )

    expect_true(gtrack.exists("test_score_e2e"))

    vals <- gextract("test_score_e2e",
        intervals = gintervals(1, 0, 50000),
        iterator = 200
    )
    # First 200bp bin overlaps the chromosome boundary (positions
    # 0..199 contain the first k=2 bp which are NA), so it must be NA.
    expect_true(is.na(vals$test_score_e2e[1]))
    # All other bins must have a finite negative value (negative log-p
    # over 200 bp).
    rest <- vals$test_score_e2e[-1]
    expect_true(all(is.finite(rest)))
    expect_true(all(rest < 0))
    # Value should be roughly bounded: -log(4) * 200 = -277 is the
    # near-uniform ceiling; the model is much better than uniform but
    # also much worse than perfect, so allow a wide range.
    expect_true(all(rest > -1000 & rest < -50))
})
