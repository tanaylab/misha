# Tests for .gsm format (gsynth.save / gsynth.load / gsynth.convert)

# Re-roots into a database of its own; hand this process's overlay back to the next
# file in the parallel worker (helper-test_db.R).
restore_groot_on_exit()


test_that("gsynth.save creates directory .gsm format with correct files", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gsynth.train(
        list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.1)),
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )

    out_dir <- file.path(tempdir(), "test_gsm_dir")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    gsynth.save(model, out_dir)

    expect_true(dir.exists(out_dir))
    expect_true(file.exists(file.path(out_dir, "metadata.yaml")))
    expect_true(file.exists(file.path(out_dir, "counts.bin")))
    expect_true(file.exists(file.path(out_dir, "cdf.bin")))

    # Check binary file sizes: total_bins * 1024 * 4 * 8 bytes
    expected_bytes <- model$total_bins * 1024L * 4L * 8L
    expect_equal(file.info(file.path(out_dir, "counts.bin"))$size, expected_bytes)
    expect_equal(file.info(file.path(out_dir, "cdf.bin"))$size, expected_bytes)

    # Check metadata
    meta <- yaml::read_yaml(file.path(out_dir, "metadata.yaml"))
    expect_equal(meta$format, "gsynth_model")
    expect_equal(meta$version, 1L)
    expect_equal(meta$markov_order, 5L)
    expect_equal(meta$n_dims, 2L)
    expect_equal(meta$total_bins, as.integer(model$total_bins))

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("directory .gsm round-trip preserves all model fields", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gsynth.train(
        list(
            expr = "g_frac + c_frac",
            breaks = seq(0, 1, 0.1),
            bin_merge = list(list(from = 0.8, to = c(0.7, 0.8)))
        ),
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )

    out_dir <- file.path(tempdir(), "test_gsm_roundtrip")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    gsynth.save(model, out_dir)
    loaded <- gsynth.load(out_dir)

    expect_s3_class(loaded, "gsynth.model")

    # Scalar fields
    expect_equal(loaded$n_dims, model$n_dims)
    expect_equal(loaded$total_bins, model$total_bins)
    expect_equal(loaded$total_kmers, model$total_kmers)
    expect_equal(loaded$total_masked, model$total_masked)
    expect_equal(loaded$total_n, model$total_n)
    expect_equal(loaded$min_obs, model$min_obs)

    # dim_sizes
    expect_equal(loaded$dim_sizes, model$dim_sizes)

    # per_bin_kmers
    expect_equal(loaded$per_bin_kmers, model$per_bin_kmers)

    # dim_specs
    for (d in seq_len(model$n_dims)) {
        expect_equal(loaded$dim_specs[[d]]$expr, model$dim_specs[[d]]$expr)
        expect_equal(loaded$dim_specs[[d]]$breaks, model$dim_specs[[d]]$breaks)
        expect_equal(loaded$dim_specs[[d]]$num_bins, model$dim_specs[[d]]$num_bins)
        expect_equal(loaded$dim_specs[[d]]$bin_map, model$dim_specs[[d]]$bin_map)
    }

    # Model data — counts and CDF matrices
    expect_equal(length(loaded$model_data$counts), length(model$model_data$counts))
    expect_equal(length(loaded$model_data$cdf), length(model$model_data$cdf))

    for (i in seq_len(model$total_bins)) {
        expect_equal(loaded$model_data$counts[[i]], model$model_data$counts[[i]],
            tolerance = 1e-15,
            info = sprintf("counts matrix mismatch at bin %d", i)
        )
        expect_equal(loaded$model_data$cdf[[i]], model$model_data$cdf[[i]],
            tolerance = 1e-15,
            info = sprintf("cdf matrix mismatch at bin %d", i)
        )
    }

    # Sparse bins should match
    expect_equal(loaded$sparse_bins, model$sparse_bins)

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("ZIP .gsm round-trip preserves all model fields", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gsynth.train(
        list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.1)),
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )

    zip_path <- file.path(tempdir(), "test_model.gsm.zip")
    on.exit(unlink(zip_path), add = TRUE)

    gsynth.save(model, zip_path, compress = TRUE)
    expect_true(file.exists(zip_path))

    loaded <- gsynth.load(zip_path)
    expect_s3_class(loaded, "gsynth.model")

    # Scalar fields
    expect_equal(loaded$n_dims, model$n_dims)
    expect_equal(loaded$total_bins, model$total_bins)
    expect_equal(loaded$total_kmers, model$total_kmers)
    expect_equal(loaded$per_bin_kmers, model$per_bin_kmers)

    # Matrices preserved exactly
    for (i in seq_len(model$total_bins)) {
        expect_equal(loaded$model_data$counts[[i]], model$model_data$counts[[i]],
            tolerance = 1e-15
        )
        expect_equal(loaded$model_data$cdf[[i]], model$model_data$cdf[[i]],
            tolerance = 1e-15
        )
    }

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("legacy RDS backward compatibility", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) gvtrack.rm("test_vt")
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gsynth.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Save as legacy RDS
    rds_file <- tempfile(fileext = ".rds")
    on.exit(unlink(rds_file), add = TRUE)
    saveRDS(model, rds_file)

    # Load should detect and handle legacy RDS
    loaded <- gsynth.load(rds_file)
    expect_s3_class(loaded, "gsynth.model")
    expect_equal(loaded$n_dims, model$n_dims)
    expect_equal(loaded$total_bins, model$total_bins)
    expect_equal(loaded$total_kmers, model$total_kmers)

    gvtrack.rm("test_vt")
})

test_that("gsynth.convert converts RDS to directory .gsm", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) gvtrack.rm("test_vt")
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gsynth.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    rds_file <- tempfile(fileext = ".rds")
    saveRDS(model, rds_file)
    on.exit(unlink(rds_file), add = TRUE)

    gsm_dir <- file.path(tempdir(), "test_convert_dir")
    on.exit(unlink(gsm_dir, recursive = TRUE), add = TRUE)

    gsynth.convert(rds_file, gsm_dir)
    expect_true(dir.exists(gsm_dir))
    expect_true(file.exists(file.path(gsm_dir, "metadata.yaml")))

    loaded <- gsynth.load(gsm_dir)
    expect_s3_class(loaded, "gsynth.model")
    expect_equal(loaded$total_bins, model$total_bins)
    expect_equal(loaded$total_kmers, model$total_kmers)

    gvtrack.rm("test_vt")
})

test_that("gsynth.convert works with compress=TRUE", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) gvtrack.rm("test_vt")
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gsynth.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    rds_file <- tempfile(fileext = ".rds")
    saveRDS(model, rds_file)
    on.exit(unlink(rds_file), add = TRUE)

    zip_file <- file.path(tempdir(), "test_convert.gsm.zip")
    on.exit(unlink(zip_file), add = TRUE)

    gsynth.convert(rds_file, zip_file, compress = TRUE)
    expect_true(file.exists(zip_file))

    loaded <- gsynth.load(zip_file)
    expect_s3_class(loaded, "gsynth.model")
    expect_equal(loaded$total_bins, model$total_bins)

    gvtrack.rm("test_vt")
})

test_that("0D (unstratified) model .gsm round-trip", {
    gdb.init_examples()

    test_intervals <- gintervals(1, 0, 50000)

    model <- gsynth.train(
        intervals = test_intervals,
        iterator = 200
    )

    # Directory round-trip
    out_dir <- file.path(tempdir(), "test_gsm_0d")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    gsynth.save(model, out_dir)
    loaded <- gsynth.load(out_dir)

    expect_s3_class(loaded, "gsynth.model")
    expect_equal(loaded$n_dims, 0L)
    expect_equal(loaded$total_bins, 1L)
    expect_equal(length(loaded$dim_specs), 0)
    expect_equal(length(loaded$dim_sizes), 0)
    expect_equal(loaded$total_kmers, model$total_kmers)
    expect_equal(loaded$total_masked, model$total_masked)
    expect_equal(loaded$total_n, model$total_n)

    # Counts and CDF preserved
    expect_equal(loaded$model_data$counts[[1]], model$model_data$counts[[1]],
        tolerance = 1e-15
    )
    expect_equal(loaded$model_data$cdf[[1]], model$model_data$cdf[[1]],
        tolerance = 1e-15
    )

    # ZIP round-trip for 0D
    zip_path <- file.path(tempdir(), "test_0d.gsm.zip")
    on.exit(unlink(zip_path), add = TRUE)
    gsynth.save(model, zip_path, compress = TRUE)
    loaded_zip <- gsynth.load(zip_path)
    expect_equal(loaded_zip$n_dims, 0L)
    expect_equal(loaded_zip$total_bins, 1L)
    expect_equal(loaded_zip$model_data$counts[[1]], model$model_data$counts[[1]],
        tolerance = 1e-15
    )
})

test_that("dim_specs with non-identity bin_map are preserved", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gsynth.train(
        list(
            expr = "g_frac + c_frac",
            breaks = seq(0, 1, 0.1),
            bin_merge = list(list(from = 0.8, to = c(0.7, 0.8)))
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # The first dim should have a non-identity bin_map
    bm <- model$dim_specs[[1]]$bin_map
    expect_false(all(bm == seq_along(bm)),
        info = "Expected non-identity bin_map from bin_merge"
    )

    out_dir <- file.path(tempdir(), "test_gsm_binmap")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    gsynth.save(model, out_dir)
    loaded <- gsynth.load(out_dir)

    expect_equal(loaded$dim_specs[[1]]$bin_map, model$dim_specs[[1]]$bin_map)

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gsynth.save errors on non-model input", {
    expect_error(gsynth.save(list(a = 1), tempdir()), "gsynth.model")
})

test_that("gsynth.load errors on missing file", {
    expect_error(gsynth.load("/nonexistent/path/model.gsm"), "not found")
})

test_that("gsynth.convert errors on non-model RDS", {
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp))
    saveRDS(list(a = 1), tmp)
    expect_error(gsynth.convert(tmp, tempfile()), "valid gsynth.model")
})

test_that("metadata.yaml fields match spec", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gsynth.train(
        list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.1)),
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )

    out_dir <- file.path(tempdir(), "test_gsm_meta")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    gsynth.save(model, out_dir)
    meta <- yaml::read_yaml(file.path(out_dir, "metadata.yaml"))

    # Check required top-level fields
    expect_true("format" %in% names(meta))
    expect_true("version" %in% names(meta))
    expect_true("markov_order" %in% names(meta))
    expect_true("n_dims" %in% names(meta))
    expect_true("dim_sizes" %in% names(meta))
    expect_true("total_bins" %in% names(meta))
    expect_true("pseudocount" %in% names(meta))
    expect_true("min_obs" %in% names(meta))
    expect_true("total_kmers" %in% names(meta))
    expect_true("total_masked" %in% names(meta))
    expect_true("total_n" %in% names(meta))
    expect_true("per_bin_kmers" %in% names(meta))
    expect_true("dim_specs" %in% names(meta))
    expect_true("data" %in% names(meta))

    # Check data section
    expect_equal(meta$data$counts$dtype, "float64")
    expect_equal(meta$data$counts$order, "C")
    expect_equal(meta$data$counts$file, "counts.bin")
    expect_equal(meta$data$cdf$dtype, "float64")
    expect_equal(meta$data$cdf$order, "C")
    expect_equal(meta$data$cdf$file, "cdf.bin")

    # Check dim_specs entries
    expect_equal(length(meta$dim_specs), 2)
    for (d in seq_len(2)) {
        s <- meta$dim_specs[[d]]
        expect_true("expr" %in% names(s))
        expect_true("breaks" %in% names(s))
        expect_true("num_bins" %in% names(s))
    }

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("sparse bins are preserved in .gsm round-trip", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    # Use min_obs to create sparse bins
    model <- suppressWarnings(gsynth.train(
        list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.05)),
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.05)),
        intervals = test_intervals,
        iterator = 200,
        min_obs = 1000
    ))

    expect_true(length(model$sparse_bins) > 0, info = "Need sparse bins for this test")

    out_dir <- file.path(tempdir(), "test_gsm_sparse")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    gsynth.save(model, out_dir)
    loaded <- gsynth.load(out_dir)

    expect_equal(loaded$sparse_bins, model$sparse_bins)
    expect_equal(loaded$min_obs, model$min_obs)

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("a loaded model samples and scores like the model it was saved from", {
    gdb.init_examples()
    if ("g_frac_it" %in% gvtrack.ls()) gvtrack.rm("g_frac_it")
    gvtrack.create("g_frac_it", NULL, "kmer.frac", kmer = "G")
    withr::defer({
        gvtrack.rm("g_frac_it")
        for (tr in c("zz_gsm_score_a", "zz_gsm_score_b")) {
            if (gtrack.exists(tr)) gtrack.rm(tr, force = TRUE)
        }
    })

    model <- gsynth.train(
        list(expr = "g_frac_it", breaks = c(0, 0.2, 0.25, 0.3, 1)),
        intervals = gintervals(1, 0, 50000),
        iterator = 500,
        k = 3
    )

    paths <- file.path(tempdir(), c("test_gsm_iter", "test_gsm_iter.zip"))
    withr::defer(unlink(paths, recursive = TRUE))
    for (compress in c(FALSE, TRUE)) {
        f <- paths[compress + 1]
        gsynth.save(model, f, compress = compress)
        loaded <- gsynth.load(f)
        expect_equal(loaded$iterator, 500)

        iv <- gintervals(1, 1234, 9876)
        expect_identical(
            gsynth.sample(loaded, output_format = "vector", intervals = iv, seed = 17),
            gsynth.sample(model, output_format = "vector", intervals = iv, seed = 17)
        )
    }

    iv <- gintervals(1, 1000, 6000)
    gsynth.score(model, "zz_gsm_score_a", intervals = iv, resolution = 100)
    gsynth.score(loaded, "zz_gsm_score_b", intervals = iv, resolution = 100)
    a <- gextract("zz_gsm_score_a", intervals = iv, iterator = 100)$zz_gsm_score_a
    b <- gextract("zz_gsm_score_b", intervals = iv, iterator = 100)$zz_gsm_score_b
    expect_true(any(is.finite(a)))
    expect_identical(a, b)

    # pymisha writes the iterator under the same key
    meta <- yaml::read_yaml(file.path(paths[1], "metadata.yaml"))
    expect_equal(meta$iterator, 500L)
})

test_that("a model saved without its iterator asks for it, and an unstratified one samples", {
    gdb.init_examples()
    if ("g_frac_noit" %in% gvtrack.ls()) gvtrack.rm("g_frac_noit")
    gvtrack.create("g_frac_noit", NULL, "kmer.frac", kmer = "G")
    withr::defer(gvtrack.rm("g_frac_noit"))

    # a .gsm written by a gsynth.save that did not store the iterator
    old_save <- function(model, f) {
        gsynth.save(model, f)
        meta <- yaml::read_yaml(file.path(f, "metadata.yaml"))
        meta$iterator <- NULL
        yaml::write_yaml(meta, file.path(f, "metadata.yaml"), precision = 15)
    }

    model <- gsynth.train(
        list(expr = "g_frac_noit", breaks = c(0, 0.25, 1)),
        intervals = gintervals(1, 0, 20000),
        iterator = 200,
        k = 2
    )
    f <- file.path(tempdir(), "test_gsm_noiter")
    withr::defer(unlink(f, recursive = TRUE))
    old_save(model, f)
    loaded <- gsynth.load(f)
    expect_null(loaded$iterator)
    iv <- gintervals(1, 0, 1000)
    expect_error(
        gsynth.sample(loaded, output_format = "vector", intervals = iv, seed = 1),
        "model\\$iterator <- "
    )
    expect_error(gsynth.score(loaded, "zz_gsm_noiter", intervals = iv), "model\\$iterator <- ")
    expect_false(gtrack.exists("zz_gsm_noiter"))
    loaded$iterator <- 200
    expect_identical(
        gsynth.sample(loaded, output_format = "vector", intervals = iv, seed = 1),
        gsynth.sample(model, output_format = "vector", intervals = iv, seed = 1)
    )

    # every position of an unstratified model is in its single bin, whatever the iterator
    model0 <- gsynth.train(intervals = gintervals(1, 0, 20000), iterator = 200, k = 2)
    f0 <- file.path(tempdir(), "test_gsm_noiter0")
    withr::defer(unlink(f0, recursive = TRUE))
    old_save(model0, f0)
    loaded0 <- gsynth.load(f0)
    expect_null(loaded0$iterator)
    expect_identical(
        gsynth.sample(loaded0, output_format = "vector", intervals = iv, seed = 1),
        gsynth.sample(model0, output_format = "vector", intervals = iv, seed = 1)
    )
    expect_error(gsynth.score(loaded0, "zz_gsm_noiter0", intervals = iv), "or pass resolution")
})

test_that("a model trained with the implicit iterator of a dense track samples, scores and saves", {
    gdb.init_examples()
    if ("dense_avg_it" %in% gvtrack.ls()) gvtrack.rm("dense_avg_it")
    gvtrack.create("dense_avg_it", "dense_track", "avg")
    withr::defer({
        gvtrack.rm("dense_avg_it")
        for (tr in c("zz_impl_a", "zz_impl_b")) if (gtrack.exists(tr)) gtrack.rm(tr, force = TRUE)
    })
    rng <- gsummary("dense_track", intervals = gintervals(1, 0, 50000))
    breaks <- seq(rng["Min"], rng["Max"], length.out = 4)
    binsize <- gtrack.info("dense_track")$bin.size

    # no iterator: the bins are those of dense_track
    model <- gsynth.train(list(expr = "dense_avg_it", breaks = breaks), intervals = gintervals(1, 0, 50000), k = 2)
    expect_equal(model$iterator, binsize)
    explicit <- gsynth.train(list(expr = "dense_avg_it", breaks = breaks), intervals = gintervals(1, 0, 50000), iterator = binsize, k = 2)
    expect_equal(model$model_data$counts, explicit$model_data$counts)

    iv <- gintervals(1, 1000, 5000)
    s <- gsynth.sample(model, output_format = "vector", intervals = iv, seed = 2)
    expect_identical(s, gsynth.sample(explicit, output_format = "vector", intervals = iv, seed = 2))

    f <- file.path(tempdir(), "test_gsm_implicit_iter")
    withr::defer(unlink(f, recursive = TRUE))
    gsynth.save(model, f)
    loaded <- gsynth.load(f)
    expect_equal(loaded$iterator, binsize)
    expect_identical(gsynth.sample(loaded, output_format = "vector", intervals = iv, seed = 2), s)

    gsynth.score(model, "zz_impl_a", intervals = iv)
    gsynth.score(loaded, "zz_impl_b", intervals = iv)
    a <- gextract("zz_impl_a", intervals = iv, iterator = binsize)$zz_impl_a
    expect_true(any(is.finite(a)))
    expect_identical(gextract("zz_impl_b", intervals = iv, iterator = binsize)$zz_impl_b, a)
})

test_that("a .gsm without a prior key loads (metadata$prior partially matched prior_mode)", {
    gdb.init_examples()
    m <- gsynth.train(intervals = gintervals(1, 0, 20000), iterator = 200, k = 2)
    f <- file.path(withr::local_tempdir(), "m_noprior")
    gsynth.save(m, f)
    meta <- yaml::read_yaml(file.path(f, "metadata.yaml"))
    meta$prior <- NULL
    yaml::write_yaml(meta, file.path(f, "metadata.yaml"))
    loaded <- gsynth.load(f)
    expect_s3_class(loaded, class(m)[1])
    expect_equal(loaded$k, m$k)
})

test_that("an RDS model saved before priors existed loads with a uniform prior and saves again", {
    gdb.init_examples()
    m <- gsynth.train(intervals = gintervals(1, 0, 20000), iterator = 200, k = 2)
    m$prior <- NULL
    m$prior_mode <- NULL
    old <- file.path(withr::local_tempdir(), "old_model.rds")
    saveRDS(m, old)
    loaded <- gsynth.load(old)
    expect_identical(loaded$prior_mode, "uniform")
    expect_true(is.matrix(loaded[["prior"]]))
    expect_equal(ncol(loaded[["prior"]]), 4L)
    expect_no_error(gsynth.save(loaded, file.path(withr::local_tempdir(), "resaved")))
})

test_that("printing a model with a non-uniform prior_mode but no prior does not fail", {
    gdb.init_examples()
    m <- gsynth.train(intervals = gintervals(1, 0, 20000), iterator = 200, k = 2)
    m$prior <- NULL
    m$prior_mode <- "custom"
    expect_no_error(capture.output(print(m)))
})
