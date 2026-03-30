# Tests for variable Markov order k parameter in gsynth functions

test_that("gsynth.train errors on invalid k values", {
    gdb.init("/home/aviezerl/hg38")

    test_intervals <- gintervals(1, 0, 100000)

    # k = 0 should error (below minimum of 1)
    expect_error(
        gsynth.train(intervals = test_intervals, iterator = 200, k = 0),
        "k must be a single integer between 1 and 8"
    )

    # k = 9 should error (above maximum of 8)
    expect_error(
        gsynth.train(intervals = test_intervals, iterator = 200, k = 9),
        "k must be a single integer between 1 and 8"
    )

    # k = 3.5 should error (non-integer numeric)
    expect_error(
        gsynth.train(intervals = test_intervals, iterator = 200, k = 3.5),
        "k must be a single integer between 1 and 8"
    )

    # k = "abc" should error (as.integer("abc") produces NA)
    expect_error(
        suppressWarnings(
            gsynth.train(intervals = test_intervals, iterator = 200, k = "abc")
        ),
        "k must be a single integer between 1 and 8"
    )
})

test_that("gsynth.train with k=3 basic training and sampling", {
    gdb.init("/home/aviezerl/hg38")

    test_intervals <- gintervals(1, 0, 100000)

    # Train a 0D model with k=3
    model <- gsynth.train(
        intervals = test_intervals,
        iterator = 200,
        k = 3
    )

    expect_s3_class(model, "gsynth.model")
    expect_equal(model$k, 3L)
    expect_equal(model$num_kmers, 64L) # 4^3

    # CDF matrices should be 64 x 4
    cdf_mat <- model$model_data$cdf[[1]]
    expect_equal(dim(cdf_mat), c(64, 4))

    # Sample and verify valid DNA
    output_fasta <- tempfile(fileext = ".fa")
    on.exit(unlink(output_fasta), add = TRUE)

    gsynth.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = test_intervals,
        seed = 60427
    )

    expect_true(file.exists(output_fasta))
    lines <- readLines(output_fasta)
    seq_lines <- lines[!grepl("^>", lines)]
    full_seq <- paste(seq_lines, collapse = "")

    # Should only contain A, C, G, T
    expect_true(grepl("^[ACGT]+$", full_seq))

    # Should have correct length (100kb)
    expect_equal(nchar(full_seq), 100000)
})

test_that("gsynth.train with k=1 minimal model", {
    gdb.init("/home/aviezerl/hg38")

    test_intervals <- gintervals(1, 0, 100000)

    model <- gsynth.train(
        intervals = test_intervals,
        iterator = 200,
        k = 1
    )

    expect_s3_class(model, "gsynth.model")
    expect_equal(model$k, 1L)
    expect_equal(model$num_kmers, 4L) # 4^1

    # CDF matrices should be 4 x 4
    cdf_mat <- model$model_data$cdf[[1]]
    expect_equal(dim(cdf_mat), c(4, 4))

    # CDF should be valid
    expect_true(all(cdf_mat >= 0))
    expect_true(all(cdf_mat <= 1))
    expect_true(all(abs(cdf_mat[, 4] - 1) < 1e-5))

    # Sample and verify valid DNA
    output_fasta <- tempfile(fileext = ".fa")
    on.exit(unlink(output_fasta), add = TRUE)

    gsynth.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = gintervals(1, 0, 10000),
        seed = 42
    )

    lines <- readLines(output_fasta)
    seq_lines <- lines[!grepl("^>", lines)]
    full_seq <- paste(seq_lines, collapse = "")
    expect_true(grepl("^[ACGT]+$", full_seq))
    expect_equal(nchar(full_seq), 10000)
})

test_that("gsynth.train with k=7 larger model", {
    gdb.init("/home/aviezerl/hg38")

    test_intervals <- gintervals(1, 0, 100000)

    model <- gsynth.train(
        intervals = test_intervals,
        iterator = 200,
        k = 7
    )

    expect_s3_class(model, "gsynth.model")
    expect_equal(model$k, 7L)
    expect_equal(model$num_kmers, 16384L) # 4^7

    # CDF matrices should be 16384 x 4
    cdf_mat <- model$model_data$cdf[[1]]
    expect_equal(dim(cdf_mat), c(16384, 4))

    # Sample and verify
    output_fasta <- tempfile(fileext = ".fa")
    on.exit(unlink(output_fasta), add = TRUE)

    gsynth.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = gintervals(1, 0, 10000),
        seed = 99
    )

    lines <- readLines(output_fasta)
    seq_lines <- lines[!grepl("^>", lines)]
    full_seq <- paste(seq_lines, collapse = "")
    expect_true(grepl("^[ACGT]+$", full_seq))
    expect_equal(nchar(full_seq), 10000)
})

test_that("save/load round-trip with k != 5", {
    gdb.init("/home/aviezerl/hg38")

    test_intervals <- gintervals(1, 0, 100000)

    model <- gsynth.train(
        intervals = test_intervals,
        iterator = 200,
        k = 3
    )

    # Save as .gsm directory
    out_dir <- file.path(tempdir(), "test_gsm_k3")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    gsynth.save(model, out_dir)

    # Verify .gsm metadata has version 2 for k != 5
    meta <- yaml::read_yaml(file.path(out_dir, "metadata.yaml"))
    expect_equal(meta$version, 2L)
    expect_equal(meta$markov_order, 3L)

    loaded <- gsynth.load(out_dir)

    expect_s3_class(loaded, "gsynth.model")
    expect_equal(loaded$k, 3L)
    expect_equal(loaded$num_kmers, 64L)

    # Verify CDF dimensions match
    expect_equal(dim(loaded$model_data$cdf[[1]]), c(64, 4))

    # Verify model data matches exactly
    expect_equal(loaded$model_data$counts[[1]], model$model_data$counts[[1]],
        tolerance = 1e-15
    )
    expect_equal(loaded$model_data$cdf[[1]], model$model_data$cdf[[1]],
        tolerance = 1e-15
    )

    # Verify sampling is identical with same seed.
    # Loaded model has iterator=NULL (not stored in .gsm format),
    # so restore it from the original model before sampling.
    loaded$iterator <- model$iterator

    out1 <- tempfile(fileext = ".fa")
    out2 <- tempfile(fileext = ".fa")
    on.exit(unlink(c(out1, out2)), add = TRUE)

    sample_intervals <- gintervals(1, 0, 10000)

    gsynth.sample(model, out1,
        output_format = "fasta",
        intervals = sample_intervals, seed = 60427
    )
    gsynth.sample(loaded, out2,
        output_format = "fasta",
        intervals = sample_intervals, seed = 60427
    )

    expect_identical(readLines(out1), readLines(out2))
})

test_that("k=5 produces identical results to default", {
    gdb.init("/home/aviezerl/hg38")

    test_intervals <- gintervals(1, 0, 100000)

    # Train with explicit k=5
    model_k5 <- gsynth.train(
        intervals = test_intervals,
        iterator = 200,
        k = 5
    )

    # Train with default (no k specified)
    model_default <- gsynth.train(
        intervals = test_intervals,
        iterator = 200
    )

    expect_equal(model_k5$k, 5L)
    expect_equal(model_default$k, 5L)
    expect_equal(model_k5$num_kmers, 1024L)
    expect_equal(model_default$num_kmers, 1024L)

    # CDF matrices should be identical
    expect_equal(
        model_k5$model_data$cdf[[1]],
        model_default$model_data$cdf[[1]],
        tolerance = 1e-15
    )

    # Counts should be identical
    expect_equal(
        model_k5$model_data$counts[[1]],
        model_default$model_data$counts[[1]],
        tolerance = 1e-15
    )
})

test_that("k with 1D stratification", {
    gdb.init("/home/aviezerl/hg38")

    for (vt in c("g_frac_k", "c_frac_k")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac_k", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac_k", NULL, "kmer.frac", kmer = "C")
    on.exit({
        if ("g_frac_k" %in% gvtrack.ls()) gvtrack.rm("g_frac_k")
        if ("c_frac_k" %in% gvtrack.ls()) gvtrack.rm("c_frac_k")
    }, add = TRUE)

    test_intervals <- gintervals(1, 0, 100000)

    # Train 1D model with k=3
    model <- gsynth.train(
        list(expr = "g_frac_k + c_frac_k", breaks = seq(0, 1, 0.2)), # 5 bins
        intervals = test_intervals,
        iterator = 200,
        k = 3
    )

    expect_s3_class(model, "gsynth.model")
    expect_equal(model$k, 3L)
    expect_equal(model$num_kmers, 64L)
    expect_equal(model$n_dims, 1)
    expect_equal(model$dim_sizes, 5)
    expect_equal(model$total_bins, 5)

    # All CDF matrices should be 64 x 4
    for (bin in seq_len(model$total_bins)) {
        expect_equal(dim(model$model_data$cdf[[bin]]), c(64, 4))
    }

    # Sample and verify
    output_fasta <- tempfile(fileext = ".fa")
    on.exit(unlink(output_fasta), add = TRUE)

    gsynth.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = gintervals(1, 0, 10000),
        seed = 60427
    )

    lines <- readLines(output_fasta)
    seq_lines <- lines[!grepl("^>", lines)]
    full_seq <- paste(seq_lines, collapse = "")
    expect_true(grepl("^[ACGT]+$", full_seq))
    expect_equal(nchar(full_seq), 10000)
})
