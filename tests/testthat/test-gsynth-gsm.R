# Tests for .gsm format (gsynth.save / gsynth.load / gsynth.convert)

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
