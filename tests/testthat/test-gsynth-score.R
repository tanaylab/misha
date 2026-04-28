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

test_that("gsynth.score: resolution invariance (sum of bp track == windowed track)", {
    gdb.init_examples()
    on.exit(
        {
            for (t in c("ts_bp", "ts_w200")) {
                if (gtrack.exists(t)) gtrack.rm(t, force = TRUE)
            }
        },
        add = TRUE
    )

    model <- gsynth.train(
        intervals = gintervals(1, 0, 50000),
        iterator = 200,
        k = 2L
    )

    test_iv <- gintervals(1, 1000, 5000) # avoid chrom-start NA region

    gsynth.score(
        model = model, track = "ts_bp",
        intervals = test_iv, resolution = 1,
        overwrite = TRUE
    )
    gsynth.score(
        model = model, track = "ts_w200",
        intervals = test_iv, resolution = 200,
        overwrite = TRUE
    )

    bp <- gextract("ts_bp", intervals = test_iv, iterator = 1)
    w200 <- gextract("ts_w200", intervals = test_iv, iterator = 200)

    # Sum the bp track in non-overlapping width-200 windows aligned to
    # multiples of 200 from chrom start. The interval starts at 1000
    # (multiple of 200) so windows align perfectly.
    expect_equal(nrow(bp), test_iv$end - test_iv$start)

    expected <- tapply(
        bp$ts_bp,
        (bp$start %/% 200) * 200,
        function(v) {
            if (any(is.na(v))) NA_real_ else sum(v)
        }
    )

    cmp <- merge(
        data.frame(
            start = as.integer(names(expected)),
            expected = as.numeric(expected)
        ),
        data.frame(start = w200$start, w200 = w200$ts_w200),
        by = "start"
    )
    # Sums of 200 floats can accumulate ~1e-5 relative error.
    expect_equal(cmp$w200, cmp$expected, tolerance = 1e-5)
})

test_that("gsynth.score: first k bp of chromosome are NA", {
    gdb.init_examples()
    on.exit(
        {
            if (gtrack.exists("ts_chrom_start")) {
                gtrack.rm("ts_chrom_start", force = TRUE)
            }
        },
        add = TRUE
    )

    model <- gsynth.train(
        intervals = gintervals(1, 0, 50000),
        iterator = 200,
        k = 5L # bigger k makes the boundary effect easier to see
    )

    gsynth.score(
        model = model, track = "ts_chrom_start",
        intervals = gintervals(1, 0, 1000),
        resolution = 1,
        overwrite = TRUE
    )

    vals <- gextract("ts_chrom_start",
        intervals = gintervals(1, 0, 1000),
        iterator = 1
    )

    # First k=5 bp must be NA; the rest must be finite.
    expect_true(all(is.na(vals$ts_chrom_start[1:5])))
    expect_true(all(is.finite(vals$ts_chrom_start[6:nrow(vals)])))
})

test_that("gsynth.score: sparse_policy='uniform' replaces sparse-bin NaN with log(1/4)", {
    gdb.init_examples()
    on.exit(
        {
            for (t in c("ts_sparse_na", "ts_sparse_uniform")) {
                if (gtrack.exists(t)) gtrack.rm(t, force = TRUE)
            }
            for (vt in c("g_frac_sp", "c_frac_sp")) {
                if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
            }
        },
        add = TRUE
    )

    gvtrack.create("g_frac_sp", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac_sp", NULL, "kmer.frac", kmer = "C")

    # Train with high min_obs and a tiny region so most strata fall
    # below the threshold and are marked NA.
    model <- suppressWarnings(gsynth.train(
        list(expr = "g_frac_sp + c_frac_sp", breaks = seq(0, 1, 0.05)),
        intervals = gintervals(1, 0, 30000),
        iterator = 200,
        k = 2L,
        min_obs = 1e6 # forces nearly all bins to be sparse
    ))
    expect_true(length(model$sparse_bins) > 0L)

    iv <- gintervals(1, 1000, 5000) # avoid chrom-start NA

    gsynth.score(
        model = model, track = "ts_sparse_na",
        intervals = iv, resolution = 1,
        sparse_policy = "NA", overwrite = TRUE
    )
    gsynth.score(
        model = model, track = "ts_sparse_uniform",
        intervals = iv, resolution = 1,
        sparse_policy = "uniform", overwrite = TRUE
    )

    na_vals <- gextract("ts_sparse_na", intervals = iv, iterator = 1)
    unif_vals <- gextract("ts_sparse_uniform", intervals = iv, iterator = 1)

    # In NA mode, every position whose stratum bin is sparse is NA.
    # Most positions should be NA given the aggressive min_obs.
    expect_true(mean(is.na(na_vals$ts_sparse_na)) > 0.5)

    # In uniform mode, those same positions must be log(1/4).
    sparse_pos <- which(is.na(na_vals$ts_sparse_na) &
        is.finite(unif_vals$ts_sparse_uniform))
    expect_true(length(sparse_pos) > 0L)
    take <- head(sparse_pos, 10)
    expect_equal(unif_vals$ts_sparse_uniform[take],
        rep(log(0.25), length(take)),
        tolerance = 1e-6
    )
})
