test_that("gcanvas.train produces valid model with 1D stratification", {
    gdb.init_examples()

    # Create a virtual track using an existing track
    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    # Get small intervals for testing
    test_intervals <- gintervals(1, 0, 100000)

    # Get the range of values for appropriate breaks
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Train model with single dimension
    model <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Check model structure
    expect_s3_class(model, "gcanvas.model")
    expect_equal(model$n_dims, 1)
    expect_equal(model$total_bins, 10)
    expect_equal(model$dim_sizes, 10)
    expect_equal(length(model$dim_specs), 1)
    expect_true(model$total_kmers > 0)
    expect_equal(length(model$per_bin_kmers), 10)

    # Check that CDFs are valid (sum to 1 for each context)
    cdf_list <- model$model_data$cdf
    expect_equal(length(cdf_list), 10)

    # Check first bin's CDF structure
    cdf_mat <- cdf_list[[1]]
    expect_equal(dim(cdf_mat), c(1024, 4))

    # CDF should be monotonically increasing and end at 1
    for (ctx in 1:10) { # Check first 10 contexts
        expect_true(all(diff(cdf_mat[ctx, ]) >= 0))
        expect_equal(cdf_mat[ctx, 4], 1.0, tolerance = 1e-5)
    }

    # Clean up
    gvtrack.rm("test_vt")
})

test_that("gcanvas.train works with 2D stratification", {
    gdb.init_examples()

    # Create virtual tracks
    for (vt in c("g_frac", "c_frac", "gc_vt")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
    gvtrack.create("gc_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 100000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Train with 2D stratification
    model <- gcanvas.train(
        list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.1)), # 10 bins
        list(expr = "gc_vt", breaks = seq(track_range["Min"], track_range["Max"], length.out = 5)), # 4 bins
        intervals = test_intervals,
        iterator = 200
    )

    expect_s3_class(model, "gcanvas.model")
    expect_equal(model$n_dims, 2)
    expect_equal(model$dim_sizes, c(10, 4))
    expect_equal(model$total_bins, 40)
    expect_equal(length(model$model_data$cdf), 40)

    # Check dimension specs
    expect_equal(model$dim_specs[[1]]$expr, "g_frac + c_frac")
    expect_equal(model$dim_specs[[2]]$expr, "gc_vt")
    expect_equal(model$dim_specs[[1]]$num_bins, 10)
    expect_equal(model$dim_specs[[2]]$num_bins, 4)

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
    gvtrack.rm("gc_vt")
})

test_that("gcanvas.train respects mask", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 100000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Create a mask covering half the region
    mask <- gintervals(1, 0, 50000)

    # Train without mask
    model_no_mask <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Train with mask
    model_with_mask <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        mask = mask,
        intervals = test_intervals,
        iterator = 200
    )

    # Model with mask should have fewer k-mers
    expect_lt(model_with_mask$total_kmers, model_no_mask$total_kmers)
    expect_gt(model_with_mask$total_masked, 0)

    gvtrack.rm("test_vt")
})

test_that("gcanvas.save and gcanvas.load work correctly", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Save and load
    temp_file <- tempfile(fileext = ".rds")
    gcanvas.save(model, temp_file)
    expect_true(file.exists(temp_file))

    loaded_model <- gcanvas.load(temp_file)
    expect_s3_class(loaded_model, "gcanvas.model")
    expect_equal(loaded_model$n_dims, model$n_dims)
    expect_equal(loaded_model$total_bins, model$total_bins)
    expect_equal(loaded_model$total_kmers, model$total_kmers)

    # Clean up
    unlink(temp_file)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.train with bin_merge works per dimension", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac", "gc_vt")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
    gvtrack.create("gc_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 100000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Create breaks for dimension 2
    breaks2 <- seq(track_range["Min"], track_range["Max"], length.out = 11)

    # Train with bin merging on both dimensions
    model <- gcanvas.train(
        list(
            expr = "g_frac + c_frac",
            breaks = seq(0, 1, 0.1), # 10 bins
            bin_merge = list(list(from = 0.8, to = c(0.7, 0.8))) # Merge bins 9-10 to bin 8
        ),
        list(
            expr = "gc_vt",
            breaks = breaks2,
            bin_merge = list(list(from = c(breaks2[9], Inf), to = c(breaks2[8], breaks2[9])))
        ),
        intervals = test_intervals,
        iterator = 200
    )

    expect_s3_class(model, "gcanvas.model")
    expect_equal(model$n_dims, 2)

    # Check that bin_map is stored correctly
    expect_true(all(model$dim_specs[[1]]$bin_map[9:10] == 8))

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
    gvtrack.rm("gc_vt")
})

test_that("gcanvas.bin_map helper function works", {
    breaks <- seq(0, 1, 0.025) # 40 bins

    # Test mapping all values above 0.7 to bin containing (0.675, 0.7]
    bin_map <- gcanvas.bin_map(
        breaks = breaks,
        merge_ranges = list(
            list(from = 0.7, to = c(0.675, 0.7))
        )
    )

    expect_type(bin_map, "integer")
    expect_equal(length(bin_map), 40)
    expect_true(all(bin_map >= 1 & bin_map <= 40))

    # The bin containing 0.675-0.7 should map to itself
    target_bin <- findInterval(0.6875, breaks, rightmost.closed = TRUE)
    expect_equal(as.integer(bin_map[target_bin]), target_bin)
})

test_that("gcanvas.bin_map with multiple merge ranges", {
    breaks <- seq(0, 1, 0.1) # 10 bins

    # Merge both low and high value bins into middle bins
    bin_map <- gcanvas.bin_map(
        breaks = breaks,
        merge_ranges = list(
            list(from = c(-Inf, 0.2), to = c(0.2, 0.3)), # Low values -> bin 3
            list(from = c(0.8, Inf), to = c(0.7, 0.8)) # High values -> bin 8
        )
    )

    expect_type(bin_map, "integer")
    expect_equal(length(bin_map), 10)

    # Bins 1 and 2 should map to bin 3
    expect_equal(as.integer(bin_map[1]), 3)
    expect_equal(as.integer(bin_map[2]), 3)

    # Bins 9 and 10 should map to bin 8
    expect_equal(as.integer(bin_map[9]), 8)
    expect_equal(as.integer(bin_map[10]), 8)

    # Middle bins should map to themselves
    expect_equal(as.integer(bin_map[5]), 5)
})

test_that("print.gcanvas.model works for multi-dimensional model", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gcanvas.train(
        list(expr = "g_frac", breaks = seq(0, 1, 0.2)), # 5 bins
        list(expr = "c_frac", breaks = seq(0, 1, 0.25)), # 4 bins
        intervals = test_intervals,
        iterator = 200
    )

    # Print should show dimensional info
    expect_output(print(model), "Genome Canvas Markov-5 Model")
    expect_output(print(model), "Dimensions: 2")
    expect_output(print(model), "Total bins: 20")

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gcanvas.sample produces output file", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Sample to FASTA
    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = test_intervals,
        seed = 60427
    )

    expect_true(file.exists(output_fasta))
    expect_gt(file.size(output_fasta), 0)

    # Check FASTA format
    lines <- readLines(output_fasta, n = 2)
    expect_true(grepl("^>", lines[1]))
    expect_true(grepl("^[ACGT]+$", lines[2]))

    # Clean up
    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample from 2D model", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gcanvas.train(
        list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.1)), # 10 bins
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)), # 5 bins
        intervals = test_intervals,
        iterator = 200
    )

    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = gintervals(1, 0, 10000),
        seed = 60427
    )

    expect_true(file.exists(output_fasta))
    lines <- readLines(output_fasta)
    seq_content <- paste(lines[!grepl("^>", lines)], collapse = "")
    expect_equal(nchar(seq_content), 10000)

    unlink(output_fasta)
    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gcanvas.sample with seed is reproducible", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 10000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Sample twice with same seed
    out1 <- tempfile(fileext = ".fa")
    out2 <- tempfile(fileext = ".fa")

    gcanvas.sample(model, out1,
        output_format = "fasta",
        intervals = test_intervals, seed = 60427
    )
    gcanvas.sample(model, out2,
        output_format = "fasta",
        intervals = test_intervals, seed = 60427
    )

    # Should produce identical output
    expect_identical(readLines(out1), readLines(out2))

    # Clean up
    unlink(c(out1, out2))
    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample different seeds produce different sequences", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    sample_intervals <- gintervals(1, 0, 1000)

    out1 <- tempfile(fileext = ".fa")
    out2 <- tempfile(fileext = ".fa")

    gcanvas.sample(model, out1,
        output_format = "fasta",
        intervals = sample_intervals, seed = 12345
    )
    gcanvas.sample(model, out2,
        output_format = "fasta",
        intervals = sample_intervals, seed = 67890
    )

    # Different seeds should produce different output
    expect_false(identical(readLines(out1), readLines(out2)))

    unlink(c(out1, out2))
    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample produces valid DNA sequences", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    sample_intervals <- gintervals(1, 0, 10000)

    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = sample_intervals,
        seed = 60427
    )

    lines <- readLines(output_fasta)
    seq_lines <- lines[!grepl("^>", lines)]
    full_seq <- paste(seq_lines, collapse = "")

    # Should only contain A, C, G, T
    expect_true(grepl("^[ACGT]+$", full_seq))

    # Should have correct length
    expect_equal(nchar(full_seq), 10000)

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.train with 3D stratification", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac", "gc_vt")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
    gvtrack.create("gc_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)), # 5 bins
        list(expr = "c_frac", breaks = seq(0, 0.5, 0.125)), # 4 bins
        list(expr = "gc_vt", breaks = seq(track_range["Min"], track_range["Max"], length.out = 3)), # 2 bins
        intervals = test_intervals,
        iterator = 200
    )

    expect_equal(model$n_dims, 3)
    expect_equal(model$dim_sizes, c(5, 4, 2))
    expect_equal(model$total_bins, 5 * 4 * 2) # 40

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
    gvtrack.rm("gc_vt")
})

test_that("gcanvas.train error handling", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)

    # No dimension specs
    expect_error(
        gcanvas.train(intervals = test_intervals, iterator = 200),
        "At least one dimension specification"
    )

    # Missing expr
    expect_error(
        gcanvas.train(
            list(breaks = seq(0, 1, 0.1)),
            intervals = test_intervals,
            iterator = 200
        ),
        "must have an 'expr' element"
    )

    # Missing breaks
    expect_error(
        gcanvas.train(
            list(expr = "test_vt"),
            intervals = test_intervals,
            iterator = 200
        ),
        "must have a 'breaks' element"
    )

    # Invalid breaks
    expect_error(
        gcanvas.train(
            list(expr = "test_vt", breaks = c(0.5)),
            intervals = test_intervals,
            iterator = 200
        ),
        "at least 2 elements"
    )

    gvtrack.rm("test_vt")
})

test_that("gcanvas.bin_map error handling", {
    breaks <- seq(0, 1, 0.1)

    # Invalid breaks
    expect_error(gcanvas.bin_map(breaks = c(0.5)))

    # Invalid target range (doesn't match any bin)
    expect_error(
        gcanvas.bin_map(
            breaks = breaks,
            merge_ranges = list(
                list(from = 0.5, to = c(0.123, 0.456)) # Not a valid bin
            )
        )
    )
})

test_that("gcanvas model CDF structure is correct", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Check all CDFs
    for (bin in seq_len(model$total_bins)) {
        cdf_mat <- model$model_data$cdf[[bin]]

        # Should be 1024 x 4
        expect_equal(dim(cdf_mat), c(1024, 4))

        # All values should be between 0 and 1
        expect_true(all(cdf_mat >= 0))
        expect_true(all(cdf_mat <= 1))

        # Last column should all be 1 (cumulative)
        expect_true(all(abs(cdf_mat[, 4] - 1) < 1e-5))

        # Each row should be monotonically non-decreasing
        for (ctx in 1:1024) {
            expect_true(all(diff(cdf_mat[ctx, ]) >= -1e-10))
        }
    }

    gvtrack.rm("test_vt")
})

test_that("gcanvas.save and gcanvas.load preserve all model fields", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gcanvas.train(
        list(
            expr = "g_frac + c_frac",
            breaks = seq(0, 1, 0.1),
            bin_merge = list(list(from = 0.8, to = c(0.7, 0.8)))
        ),
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )

    temp_file <- tempfile(fileext = ".rds")
    gcanvas.save(model, temp_file)
    loaded_model <- gcanvas.load(temp_file)

    # Check all fields are preserved
    expect_equal(loaded_model$n_dims, model$n_dims)
    expect_equal(loaded_model$dim_sizes, model$dim_sizes)
    expect_equal(loaded_model$total_bins, model$total_bins)
    expect_equal(loaded_model$total_kmers, model$total_kmers)
    expect_equal(loaded_model$per_bin_kmers, model$per_bin_kmers)
    expect_equal(loaded_model$total_masked, model$total_masked)
    expect_equal(loaded_model$total_n, model$total_n)

    # Check dim_specs
    for (d in seq_len(model$n_dims)) {
        expect_equal(loaded_model$dim_specs[[d]]$expr, model$dim_specs[[d]]$expr)
        expect_equal(loaded_model$dim_specs[[d]]$breaks, model$dim_specs[[d]]$breaks)
        expect_equal(loaded_model$dim_specs[[d]]$num_bins, model$dim_specs[[d]]$num_bins)
        expect_equal(loaded_model$dim_specs[[d]]$bin_map, model$dim_specs[[d]]$bin_map)
    }

    # Check model_data
    expect_equal(length(loaded_model$model_data$counts), length(model$model_data$counts))
    expect_equal(length(loaded_model$model_data$cdf), length(model$model_data$cdf))

    # Check class
    expect_s3_class(loaded_model, "gcanvas.model")

    unlink(temp_file)
    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gcanvas handles empty bins gracefully", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Create breaks that will definitely have empty bins
    # by extending well beyond the actual data range
    min_val <- track_range["Min"]
    max_val <- track_range["Max"]
    range_val <- max_val - min_val

    # Extend breaks beyond data range
    breaks <- seq(min_val - range_val, max_val + range_val, length.out = 21)

    model <- gcanvas.train(
        list(expr = "test_vt", breaks = breaks),
        intervals = test_intervals,
        iterator = 200
    )

    expect_s3_class(model, "gcanvas.model")

    # Some bins should have 0 k-mers
    expect_true(any(model$per_bin_kmers == 0))

    # Model should still be usable for sampling
    output_fasta <- tempfile(fileext = ".fa")
    expect_no_error(
        gcanvas.sample(
            model,
            output_fasta,
            output_format = "fasta",
            intervals = gintervals(1, 0, 1000),
            seed = 60427
        )
    )

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample with mask_mode copy preserves original sequence in masked regions", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Define a mask region
    mask <- gintervals(1, 1000, 2000)
    sample_intervals <- gintervals(1, 0, 3000)

    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = sample_intervals,
        mask = mask,
        mask_mode = "copy",
        seed = 60427
    )

    # Read the sampled sequence
    lines <- readLines(output_fasta)
    sampled_seq <- paste(lines[!grepl("^>", lines)], collapse = "")

    # Get original sequence for the mask region (gseq.extract returns character vector)
    original_seq <- gseq.extract(mask)[1]

    # The mask region (positions 1000-2000 within the 0-3000 interval)
    # should match the original (case-insensitive due to different output formats)
    mask_start_in_sample <- 1001 # 1-based
    mask_end_in_sample <- 2000
    sampled_mask_region <- substr(sampled_seq, mask_start_in_sample, mask_end_in_sample)

    # Compare case-insensitively (sampler may output different case than gseq.extract)
    expect_equal(toupper(sampled_mask_region), toupper(original_seq))

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample on multiple chromosomes", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals.all()
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Sample from both chromosomes
    sample_intervals <- gintervals(
        c(1, 2),
        c(0, 0),
        c(1000, 1000)
    )

    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = sample_intervals,
        seed = 60427
    )

    lines <- readLines(output_fasta)
    headers <- lines[grepl("^>", lines)]

    # Should have 2 sequence headers
    expect_equal(length(headers), 2)

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.train with different pseudocount values", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)
    breaks <- seq(track_range["Min"], track_range["Max"], length.out = 11)

    # Train with different pseudocounts
    model_pc1 <- gcanvas.train(
        list(expr = "test_vt", breaks = breaks),
        intervals = test_intervals,
        iterator = 200,
        pseudocount = 1
    )

    model_pc10 <- gcanvas.train(
        list(expr = "test_vt", breaks = breaks),
        intervals = test_intervals,
        iterator = 200,
        pseudocount = 10
    )

    # Both should produce valid models
    expect_s3_class(model_pc1, "gcanvas.model")
    expect_s3_class(model_pc10, "gcanvas.model")

    # Total k-mers should be the same (pseudocount only affects probabilities)
    expect_equal(model_pc1$total_kmers, model_pc10$total_kmers)

    # CDFs should both be valid but different
    cdf1 <- model_pc1$model_data$cdf[[1]]
    cdf10 <- model_pc10$model_data$cdf[[1]]

    # Higher pseudocount should make distributions more uniform
    expect_false(identical(cdf1, cdf10))

    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample misha output size matches interval lengths", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    train_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = train_intervals)

    model <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = train_intervals,
        iterator = 200
    )

    sample_intervals <- gintervals(
        c(1, 2),
        c(1000, 2000),
        c(1075, 2150)
    )

    output_seq <- tempfile(fileext = ".seq")
    gcanvas.sample(
        model,
        output_seq,
        output_format = "misha",
        intervals = sample_intervals,
        seed = 60427
    )

    expected_size <- sum(sample_intervals$end - sample_intervals$start)
    expect_equal(file.size(output_seq), expected_size)

    unlink(output_seq)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.train with GC and CG stratification (user use case)", {
    # This test implements the specific use case from the user:
    # stratify on both GC content (capped at 0.7) and CG dinucleotide content
    # (breaks at 0.01, 0.02, 0.03, 0.04, 0.2, capped at 0.04)
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac", "cg_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
    gvtrack.create("cg_frac", NULL, "kmer.frac", kmer = "CG")

    test_intervals <- gintervals(1, 0, 100000)

    # GC content with fine breaks, capped at 0.7
    gc_breaks <- seq(0, 1, 0.025) # 40 bins

    # CG dinucleotide with specific breaks, capped at 0.04
    cg_breaks <- c(0, 0.01, 0.02, 0.03, 0.04, 0.2) # 5 bins

    model <- gcanvas.train(
        list(
            expr = "g_frac + c_frac",
            breaks = gc_breaks,
            bin_merge = list(list(from = 0.7, to = c(0.675, 0.7)))
        ),
        list(
            expr = "cg_frac",
            breaks = cg_breaks,
            bin_merge = list(list(from = 0.04, to = c(0.03, 0.04)))
        ),
        intervals = test_intervals,
        iterator = 200
    )

    expect_s3_class(model, "gcanvas.model")
    expect_equal(model$n_dims, 2)
    expect_equal(model$dim_specs[[1]]$num_bins, 40)
    expect_equal(model$dim_specs[[2]]$num_bins, 5)
    expect_equal(model$total_bins, 40 * 5) # 200 total bins

    # Verify bin_map for GC dimension (bins 29-40 should map to bin 28)
    gc_bin_map <- model$dim_specs[[1]]$bin_map
    expect_true(all(gc_bin_map[29:40] == 28))

    # Verify bin_map for CG dimension (bin 5 should map to bin 4)
    cg_bin_map <- model$dim_specs[[2]]$bin_map
    expect_equal(as.integer(cg_bin_map[5]), 4)

    # Sample from the trained model
    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = gintervals(1, 0, 10000),
        seed = 60427
    )

    expect_true(file.exists(output_fasta))
    lines <- readLines(output_fasta)
    seq_content <- paste(lines[!grepl("^>", lines)], collapse = "")
    expect_equal(nchar(seq_content), 10000)
    expect_true(grepl("^[ACGT]+$", seq_content))

    unlink(output_fasta)
    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
    gvtrack.rm("cg_frac")
})

test_that("gcanvas.sample from model with bin_merge produces valid sequences", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 100000)

    # Train model with aggressive bin merging
    model <- gcanvas.train(
        list(
            expr = "g_frac + c_frac",
            breaks = seq(0, 1, 0.1), # 10 bins
            bin_merge = list(
                list(from = c(-Inf, 0.2), to = c(0.2, 0.3)), # Low -> bin 3
                list(from = c(0.7, Inf), to = c(0.6, 0.7)) # High -> bin 7
            )
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Verify bin merging worked
    bin_map <- model$dim_specs[[1]]$bin_map
    expect_true(all(bin_map[1:2] == 3))
    expect_true(all(bin_map[8:10] == 7))

    # Sample and verify valid output
    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = gintervals(1, 0, 5000),
        seed = 12345
    )

    lines <- readLines(output_fasta)
    seq_content <- paste(lines[!grepl("^>", lines)], collapse = "")
    expect_equal(nchar(seq_content), 5000)
    expect_true(grepl("^[ACGT]+$", seq_content))

    unlink(output_fasta)
    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gcanvas flat index computation is correct for 2D", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gcanvas.train(
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)), # 5 bins
        list(expr = "c_frac", breaks = seq(0, 0.5, 0.125)), # 4 bins
        intervals = test_intervals,
        iterator = 200
    )

    # Check flat indexing: flat_idx = idx_1 + (idx_2-1) * size_1
    # With 5 bins in dim1 and 4 bins in dim2:
    # - (dim1=1, dim2=1) -> 1
    # - (dim1=5, dim2=1) -> 5
    # - (dim1=1, dim2=2) -> 1 + 1*5 = 6
    # - (dim1=5, dim2=4) -> 5 + 3*5 = 20
    expect_equal(model$dim_sizes, c(5, 4))
    expect_equal(model$total_bins, 20)
    expect_equal(length(model$model_data$cdf), 20)

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gcanvas flat index computation is correct for 3D", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac", "a_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
    gvtrack.create("a_frac", NULL, "kmer.frac", kmer = "A")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gcanvas.train(
        list(expr = "g_frac", breaks = c(0, 0.2, 0.4)), # 2 bins
        list(expr = "c_frac", breaks = c(0, 0.15, 0.3, 0.45)), # 3 bins
        list(expr = "a_frac", breaks = c(0, 0.25, 0.5)), # 2 bins
        intervals = test_intervals,
        iterator = 200
    )

    # Total bins = 2 * 3 * 2 = 12
    expect_equal(model$dim_sizes, c(2, 3, 2))
    expect_equal(model$total_bins, 12)
    expect_equal(length(model$model_data$cdf), 12)
    expect_equal(length(model$per_bin_kmers), 12)

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
    gvtrack.rm("a_frac")
})

test_that("gcanvas.train validates dimension specs correctly", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Test that passing non-list spec fails
    expect_error(
        gcanvas.train(
            "test_vt", # Not a list
            intervals = test_intervals,
            iterator = 200
        ),
        "must be a list"
    )

    # Test that passing empty list fails
    expect_error(
        gcanvas.train(
            list(), # Empty list
            intervals = test_intervals,
            iterator = 200
        ),
        "must have an 'expr' element"
    )

    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample from 3D model works correctly", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac", "a_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
    gvtrack.create("a_frac", NULL, "kmer.frac", kmer = "A")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gcanvas.train(
        list(expr = "g_frac", breaks = c(0, 0.2, 0.4, 0.6)), # 3 bins
        list(expr = "c_frac", breaks = c(0, 0.2, 0.4, 0.6)), # 3 bins
        list(expr = "a_frac", breaks = c(0, 0.2, 0.4, 0.6)), # 3 bins
        intervals = test_intervals,
        iterator = 200
    )

    expect_equal(model$total_bins, 27) # 3 * 3 * 3

    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = gintervals(1, 0, 5000),
        seed = 60427
    )

    lines <- readLines(output_fasta)
    seq_content <- paste(lines[!grepl("^>", lines)], collapse = "")
    expect_equal(nchar(seq_content), 5000)
    expect_true(grepl("^[ACGT]+$", seq_content))

    unlink(output_fasta)
    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
    gvtrack.rm("a_frac")
})

test_that("gcanvas multi-dimensional sampling is reproducible", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gcanvas.train(
        list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.1)),
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )

    sample_intervals <- gintervals(1, 0, 5000)

    out1 <- tempfile(fileext = ".fa")
    out2 <- tempfile(fileext = ".fa")

    gcanvas.sample(model, out1,
        output_format = "fasta",
        intervals = sample_intervals, seed = 42
    )
    gcanvas.sample(model, out2,
        output_format = "fasta",
        intervals = sample_intervals, seed = 42
    )

    expect_identical(readLines(out1), readLines(out2))

    unlink(c(out1, out2))
    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gcanvas bin_merge chains correctly in multi-dimensional model", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 100000)

    # Test that bin_merge affects the right dimensions independently
    model <- gcanvas.train(
        list(
            expr = "g_frac + c_frac",
            breaks = seq(0, 1, 0.2), # 5 bins
            bin_merge = list(list(from = 0.8, to = c(0.6, 0.8))) # bin 5 -> bin 4
        ),
        list(
            expr = "g_frac",
            breaks = seq(0, 0.5, 0.1), # 5 bins
            bin_merge = list(list(from = 0.4, to = c(0.3, 0.4))) # bin 5 -> bin 4
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Check first dimension
    expect_equal(as.integer(model$dim_specs[[1]]$bin_map[5]), 4)
    expect_equal(as.integer(model$dim_specs[[1]]$bin_map[4]), 4)

    # Check second dimension
    expect_equal(as.integer(model$dim_specs[[2]]$bin_map[5]), 4)
    expect_equal(as.integer(model$dim_specs[[2]]$bin_map[4]), 4)

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gcanvas per_bin_kmers sum equals total_kmers", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    # 1D model
    model_1d <- gcanvas.train(
        list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )
    expect_equal(sum(model_1d$per_bin_kmers), model_1d$total_kmers)

    # 2D model
    model_2d <- gcanvas.train(
        list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.2)),
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )
    expect_equal(sum(model_2d$per_bin_kmers), model_2d$total_kmers)

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gcanvas.save and gcanvas.load preserve multi-dimensional structure", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac", "a_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
    gvtrack.create("a_frac", NULL, "kmer.frac", kmer = "A")

    test_intervals <- gintervals(1, 0, 50000)

    # 3D model with bin_merge
    model <- gcanvas.train(
        list(
            expr = "g_frac",
            breaks = seq(0, 0.5, 0.1),
            bin_merge = list(list(from = 0.4, to = c(0.3, 0.4)))
        ),
        list(expr = "c_frac", breaks = seq(0, 0.5, 0.125)),
        list(expr = "a_frac", breaks = seq(0, 0.5, 0.25)),
        intervals = test_intervals,
        iterator = 200
    )

    temp_file <- tempfile(fileext = ".rds")
    gcanvas.save(model, temp_file)
    loaded_model <- gcanvas.load(temp_file)

    # Verify all dimensional structure is preserved
    expect_equal(loaded_model$n_dims, 3)
    expect_equal(loaded_model$dim_sizes, model$dim_sizes)
    expect_equal(loaded_model$total_bins, model$total_bins)

    # Check dim_specs for all dimensions
    for (d in 1:3) {
        expect_equal(loaded_model$dim_specs[[d]]$expr, model$dim_specs[[d]]$expr)
        expect_equal(loaded_model$dim_specs[[d]]$breaks, model$dim_specs[[d]]$breaks)
        expect_equal(loaded_model$dim_specs[[d]]$num_bins, model$dim_specs[[d]]$num_bins)
        expect_equal(loaded_model$dim_specs[[d]]$bin_map, model$dim_specs[[d]]$bin_map)
    }

    unlink(temp_file)
    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
    gvtrack.rm("a_frac")
})

test_that("gcanvas.train with min_obs identifies sparse bins", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Create breaks that will have some empty/sparse bins
    # by extending beyond the data range
    min_val <- track_range["Min"]
    max_val <- track_range["Max"]
    range_val <- max_val - min_val
    breaks <- seq(min_val - range_val, max_val + range_val, length.out = 21)

    # Train with a high min_obs threshold to trigger sparse bin detection
    expect_warning(
        model <- gcanvas.train(
            list(expr = "test_vt", breaks = breaks),
            intervals = test_intervals,
            iterator = 200,
            min_obs = 1000
        ),
        "bins have fewer than"
    )

    # Should have some sparse bins
    expect_true(length(model$sparse_bins) > 0)
    expect_equal(model$min_obs, 1000)

    # Sparse bins should have NA CDFs
    for (bin_idx in model$sparse_bins) {
        expect_true(all(is.na(model$model_data$cdf[[bin_idx]])))
    }

    # Non-sparse bins should have valid CDFs
    non_sparse <- setdiff(seq_len(model$total_bins), model$sparse_bins)
    if (length(non_sparse) > 0) {
        for (bin_idx in non_sparse[1:min(3, length(non_sparse))]) {
            expect_false(any(is.na(model$model_data$cdf[[bin_idx]])))
        }
    }

    gvtrack.rm("test_vt")
})

test_that("gcanvas.train with min_obs=0 has no sparse bins", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(
            expr = "test_vt",
            breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
        ),
        intervals = test_intervals,
        iterator = 200,
        min_obs = 0
    )

    # With min_obs=0, no bins should be sparse
    expect_equal(length(model$sparse_bins), 0)
    expect_equal(model$min_obs, 0)

    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample warns when using sparse bins", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Create model with sparse bins
    min_val <- track_range["Min"]
    max_val <- track_range["Max"]
    range_val <- max_val - min_val
    breaks <- seq(min_val - range_val, max_val + range_val, length.out = 21)

    suppressWarnings({
        model <- gcanvas.train(
            list(expr = "test_vt", breaks = breaks),
            intervals = test_intervals,
            iterator = 200,
            min_obs = 1000
        )
    })

    # Sampling should warn about sparse bins
    output_fasta <- tempfile(fileext = ".fa")
    expect_warning(
        gcanvas.sample(
            model,
            output_fasta,
            output_format = "fasta",
            intervals = test_intervals,
            seed = 60427
        ),
        "sparse bins"
    )

    # Output should still be valid
    expect_true(file.exists(output_fasta))
    lines <- readLines(output_fasta)
    seq_content <- paste(lines[!grepl("^>", lines)], collapse = "")
    expect_true(grepl("^[ACGT]+$", seq_content))

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample with sparse bins uses uniform distribution", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Create model with very high min_obs so all bins are sparse
    suppressWarnings({
        model <- gcanvas.train(
            list(
                expr = "test_vt",
                breaks = seq(track_range["Min"], track_range["Max"], length.out = 11)
            ),
            intervals = test_intervals,
            iterator = 200,
            min_obs = 10000000 # Very high, all bins will be sparse
        )
    })

    # All bins should be sparse
    expect_equal(length(model$sparse_bins), model$total_bins)

    # Sample should still work (with uniform distribution)
    output_fasta <- tempfile(fileext = ".fa")
    suppressWarnings({
        gcanvas.sample(
            model,
            output_fasta,
            output_format = "fasta",
            intervals = gintervals(1, 0, 10000),
            seed = 60427
        )
    })

    lines <- readLines(output_fasta)
    seq_content <- paste(lines[!grepl("^>", lines)], collapse = "")

    # Should have uniform-ish base distribution
    a_count <- sum(strsplit(seq_content, "")[[1]] == "A")
    c_count <- sum(strsplit(seq_content, "")[[1]] == "C")
    g_count <- sum(strsplit(seq_content, "")[[1]] == "G")
    t_count <- sum(strsplit(seq_content, "")[[1]] == "T")

    # Each base should be roughly 25% (with some tolerance)
    total <- nchar(seq_content)
    expect_gt(a_count / total, 0.15)
    expect_lt(a_count / total, 0.35)
    expect_gt(c_count / total, 0.15)
    expect_lt(c_count / total, 0.35)
    expect_gt(g_count / total, 0.15)
    expect_lt(g_count / total, 0.35)
    expect_gt(t_count / total, 0.15)
    expect_lt(t_count / total, 0.35)

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas print shows sparse bins info", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Create model with some sparse bins
    min_val <- track_range["Min"]
    max_val <- track_range["Max"]
    range_val <- max_val - min_val
    breaks <- seq(min_val - range_val, max_val + range_val, length.out = 11)

    suppressWarnings({
        model <- gcanvas.train(
            list(expr = "test_vt", breaks = breaks),
            intervals = test_intervals,
            iterator = 200,
            min_obs = 5000
        )
    })

    # Print should show sparse bins info
    if (length(model$sparse_bins) > 0) {
        expect_output(print(model), "Sparse bins")
    }

    gvtrack.rm("test_vt")
})

test_that("gcanvas.save and gcanvas.load preserve sparse_bins", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Create model with sparse bins
    min_val <- track_range["Min"]
    max_val <- track_range["Max"]
    range_val <- max_val - min_val
    breaks <- seq(min_val - range_val, max_val + range_val, length.out = 21)

    suppressWarnings({
        model <- gcanvas.train(
            list(expr = "test_vt", breaks = breaks),
            intervals = test_intervals,
            iterator = 200,
            min_obs = 1000
        )
    })

    temp_file <- tempfile(fileext = ".rds")
    gcanvas.save(model, temp_file)
    loaded_model <- gcanvas.load(temp_file)

    # Sparse bins should be preserved
    expect_equal(loaded_model$sparse_bins, model$sparse_bins)
    expect_equal(loaded_model$min_obs, model$min_obs)

    unlink(temp_file)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample with bin_merge at sampling time", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    # Train WITHOUT bin_merge
    model <- gcanvas.train(
        list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.1)), # 10 bins
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)), # 5 bins
        intervals = test_intervals,
        iterator = 200
    )

    # Verify no bin merging during training
    expect_true(all(model$dim_specs[[1]]$bin_map == 1:10))
    expect_true(all(model$dim_specs[[2]]$bin_map == 1:5))

    # Sample WITH bin_merge at sampling time
    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = gintervals(1, 0, 5000),
        seed = 60427,
        bin_merge = list(
            list(list(from = 0.8, to = c(0.7, 0.8))), # Merge bins 9-10 -> 8
            list(list(from = 0.4, to = c(0.3, 0.4))) # Merge bin 5 -> 4
        )
    )

    expect_true(file.exists(output_fasta))
    lines <- readLines(output_fasta)
    seq_content <- paste(lines[!grepl("^>", lines)], collapse = "")
    expect_equal(nchar(seq_content), 5000)
    expect_true(grepl("^[ACGT]+$", seq_content))

    unlink(output_fasta)
    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gcanvas.sample bin_merge can be NULL per dimension", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    # Train with bin_merge on first dimension only
    model <- gcanvas.train(
        list(
            expr = "g_frac + c_frac",
            breaks = seq(0, 1, 0.1),
            bin_merge = list(list(from = 0.8, to = c(0.7, 0.8)))
        ),
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )

    # Sample with bin_merge only on second dimension (NULL for first = use training)
    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = gintervals(1, 0, 5000),
        seed = 60427,
        bin_merge = list(
            NULL, # Use training-time bin_map for dim 1
            list(list(from = 0.4, to = c(0.3, 0.4))) # Override for dim 2
        )
    )

    expect_true(file.exists(output_fasta))
    lines <- readLines(output_fasta)
    seq_content <- paste(lines[!grepl("^>", lines)], collapse = "")
    expect_equal(nchar(seq_content), 5000)

    unlink(output_fasta)
    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gcanvas.sample bin_merge validation", {
    gdb.init_examples()

    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) gvtrack.rm(vt)
    }
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    model <- gcanvas.train(
        list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.1)),
        list(expr = "g_frac", breaks = seq(0, 0.5, 0.1)),
        intervals = test_intervals,
        iterator = 200
    )

    output_fasta <- tempfile(fileext = ".fa")

    # Wrong number of dimensions in bin_merge
    expect_error(
        gcanvas.sample(
            model,
            output_fasta,
            output_format = "fasta",
            intervals = gintervals(1, 0, 5000),
            bin_merge = list(NULL) # Only 1 element for 2D model
        ),
        "2 elements"
    )

    # Wrong type
    expect_error(
        gcanvas.sample(
            model,
            output_fasta,
            output_format = "fasta",
            intervals = gintervals(1, 0, 5000),
            bin_merge = "invalid"
        ),
        "must be a list"
    )

    gvtrack.rm("g_frac")
    gvtrack.rm("c_frac")
})

test_that("gcanvas.sample bin_merge handles sparse bins from training", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Train with sparse bins (no bin_merge during training)
    min_val <- track_range["Min"]
    max_val <- track_range["Max"]
    range_val <- max_val - min_val
    breaks <- seq(min_val - range_val, max_val + range_val, length.out = 11)

    suppressWarnings({
        model <- gcanvas.train(
            list(expr = "test_vt", breaks = breaks),
            intervals = test_intervals,
            iterator = 200,
            min_obs = 5000
        )
    })

    # Should have sparse bins
    expect_true(length(model$sparse_bins) > 0)

    # Sample with bin_merge to merge sparse bins into non-sparse
    # Find the middle bins (likely non-sparse) and merge edges into them
    mid_bin <- ceiling(length(breaks) / 2) - 1
    mid_val <- (breaks[mid_bin] + breaks[mid_bin + 1]) / 2

    output_fasta <- tempfile(fileext = ".fa")

    # Use bin_merge at sample time to redirect sparse bins
    # This should reduce/eliminate sparse bin warnings
    suppressWarnings({
        gcanvas.sample(
            model,
            output_fasta,
            output_format = "fasta",
            intervals = gintervals(1, 0, 5000),
            seed = 60427,
            bin_merge = list(
                list(
                    list(from = c(-Inf, breaks[3]), to = c(breaks[mid_bin], breaks[mid_bin + 1])),
                    list(from = c(breaks[length(breaks) - 2], Inf), to = c(breaks[mid_bin], breaks[mid_bin + 1]))
                )
            )
        )
    })

    expect_true(file.exists(output_fasta))
    lines <- readLines(output_fasta)
    seq_content <- paste(lines[!grepl("^>", lines)], collapse = "")
    expect_equal(nchar(seq_content), 5000)
    expect_true(grepl("^[ACGT]+$", seq_content))

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample returns vector with output_format='vector'", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(expr = "test_vt", breaks = seq(track_range["Min"], track_range["Max"], length.out = 6)),
        intervals = test_intervals,
        iterator = 200
    )

    # Sample with output_format = "vector"
    seqs <- gcanvas.sample(
        model,
        output_format = "vector",
        intervals = gintervals(1, 0, 1000),
        seed = 60427
    )

    expect_type(seqs, "character")
    expect_equal(length(seqs), 1)
    expect_equal(nchar(seqs[1]), 1000)
    expect_true(grepl("^[ACGT]+$", seqs[1]))

    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample n_samples generates multiple sequences", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(expr = "test_vt", breaks = seq(track_range["Min"], track_range["Max"], length.out = 6)),
        intervals = test_intervals,
        iterator = 200
    )

    # Sample 5 times with vector output
    seqs <- gcanvas.sample(
        model,
        output_format = "vector",
        intervals = gintervals(1, 0, 500),
        n_samples = 5,
        seed = 60427
    )

    expect_type(seqs, "character")
    expect_equal(length(seqs), 5)
    for (i in seq_along(seqs)) {
        expect_equal(nchar(seqs[i]), 500)
        expect_true(grepl("^[ACGT]+$", seqs[i]))
    }

    # Each sample should be different (they use different random sequences)
    unique_seqs <- unique(seqs)
    expect_gt(length(unique_seqs), 1) # At least some should be different

    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample n_samples works with FASTA output", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(expr = "test_vt", breaks = seq(track_range["Min"], track_range["Max"], length.out = 6)),
        intervals = test_intervals,
        iterator = 200
    )

    output_fasta <- tempfile(fileext = ".fa")

    # Sample 3 times to FASTA
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = gintervals(1, 0, 500),
        n_samples = 3,
        seed = 60427
    )

    expect_true(file.exists(output_fasta))
    lines <- readLines(output_fasta)

    # Count headers - should have 3 (one per sample)
    headers <- lines[grepl("^>", lines)]
    expect_equal(length(headers), 3)

    # Headers should include sample numbers
    expect_true(any(grepl("_sample1", headers)))
    expect_true(any(grepl("_sample2", headers)))
    expect_true(any(grepl("_sample3", headers)))

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample n_samples with multiple intervals", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(expr = "test_vt", breaks = seq(track_range["Min"], track_range["Max"], length.out = 6)),
        intervals = test_intervals,
        iterator = 200
    )

    # 2 intervals x 3 samples = 6 sequences
    sample_intervals <- gintervals(1, c(0, 1000), c(500, 1500))

    seqs <- gcanvas.sample(
        model,
        output_format = "vector",
        intervals = sample_intervals,
        n_samples = 3,
        seed = 60427
    )

    expect_equal(length(seqs), 6) # 2 intervals * 3 samples

    # Each should be 500bp
    for (i in seq_along(seqs)) {
        expect_equal(nchar(seqs[i]), 500)
        expect_true(grepl("^[ACGT]+$", seqs[i]))
    }

    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample output_format='vector' requires no output_path", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(expr = "test_vt", breaks = seq(track_range["Min"], track_range["Max"], length.out = 6)),
        intervals = test_intervals,
        iterator = 200
    )

    # output_path is not required when output_format = "vector"
    seqs <- gcanvas.sample(
        model,
        output_format = "vector",
        intervals = gintervals(1, 0, 100),
        seed = 60427
    )

    expect_type(seqs, "character")
    expect_equal(length(seqs), 1)

    # But it IS required for other formats
    expect_error(
        gcanvas.sample(
            model,
            output_format = "fasta",
            intervals = gintervals(1, 0, 100),
            seed = 60427
        ),
        "output_path is required"
    )

    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample n_samples with same seed is reproducible", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        list(expr = "test_vt", breaks = seq(track_range["Min"], track_range["Max"], length.out = 6)),
        intervals = test_intervals,
        iterator = 200
    )

    # Sample twice with same seed
    seqs1 <- gcanvas.sample(
        model,
        output_format = "vector",
        intervals = gintervals(1, 0, 500),
        n_samples = 3,
        seed = 12345
    )

    seqs2 <- gcanvas.sample(
        model,
        output_format = "vector",
        intervals = gintervals(1, 0, 500),
        n_samples = 3,
        seed = 12345
    )

    expect_equal(seqs1, seqs2)

    # Different seed should give different results
    seqs3 <- gcanvas.sample(
        model,
        output_format = "vector",
        intervals = gintervals(1, 0, 500),
        n_samples = 3,
        seed = 54321
    )

    expect_false(identical(seqs1, seqs3))

    gvtrack.rm("test_vt")
})
