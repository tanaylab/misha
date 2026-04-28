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
