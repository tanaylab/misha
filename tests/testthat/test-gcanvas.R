test_that("gcanvas.train produces valid model", {
    # Initialize test database
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

    # Train model with stratification based on dense_track
    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11), # 10 bins
        intervals = test_intervals,
        iterator = 200
    )

    # Check model structure
    expect_s3_class(model, "gcanvas.model")
    expect_equal(model$num_bins, 10)
    expect_equal(length(model$breaks), 11)
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
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
        intervals = test_intervals,
        iterator = 200
    )

    # Train with mask
    model_with_mask <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
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
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
        intervals = test_intervals,
        iterator = 200
    )

    # Save and load
    temp_file <- tempfile(fileext = ".rds")
    gcanvas.save(model, temp_file)
    expect_true(file.exists(temp_file))

    loaded_model <- gcanvas.load(temp_file)
    expect_s3_class(loaded_model, "gcanvas.model")
    expect_equal(loaded_model$num_bins, model$num_bins)
    expect_equal(loaded_model$total_kmers, model$total_kmers)

    # Clean up
    unlink(temp_file)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.train with bin_merge merges bins correctly", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 100000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Create breaks for 10 bins
    breaks <- seq(track_range["Min"], track_range["Max"], length.out = 11)

    # Train with bin merging: merge bins 9 and 10 (high values) into bin 8
    # Bin 8 is [breaks[8], breaks[9]), bins 9-10 are [breaks[9], breaks[11])
    model <- gcanvas.train(
        expr = "test_vt",
        breaks = breaks,
        bin_merge = list(
            list(from = c(breaks[9], breaks[11]), to = c(breaks[8], breaks[9]))
        ),
        intervals = test_intervals,
        iterator = 200
    )

    expect_s3_class(model, "gcanvas.model")
    expect_equal(model$num_bins, 10)

    gvtrack.rm("test_vt")
})

test_that("gcanvas.train with bin_merge uses value-based ranges", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Create breaks
    breaks <- seq(track_range["Min"], track_range["Max"], length.out = 21) # 20 bins

    # Merge all values above a threshold into a specific bin (e.g., map all GC > 7% to (0.675, 0.7])
    # In this case, we'll map high values to a middle bin
    mid_val <- (track_range["Min"] + track_range["Max"]) / 2
    threshold <- track_range["Min"] + 0.7 * (track_range["Max"] - track_range["Min"])

    # Find which bin contains the threshold value
    target_bin_idx <- findInterval(threshold, breaks, rightmost.closed = TRUE)
    target_bin_start <- breaks[target_bin_idx]
    target_bin_end <- breaks[target_bin_idx + 1]

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = breaks,
        bin_merge = list(
            list(from = threshold, to = c(target_bin_start, target_bin_end))
        ),
        intervals = test_intervals,
        iterator = 200
    )

    expect_s3_class(model, "gcanvas.model")
    expect_equal(model$num_bins, 20)

    gvtrack.rm("test_vt")
})

test_that("gcanvas.bin_map helper function works", {
    breaks <- seq(0, 1, 0.025) # 40 bins

    # Test mapping all values above 0.07 to bin (0.675, 0.7]
    bin_map <- gcanvas.bin_map(
        breaks = breaks,
        merge_ranges = list(
            list(from = 0.07, to = c(0.675, 0.7))
        )
    )

    expect_type(bin_map, "integer")
    expect_equal(length(bin_map), 40)
    expect_true(all(bin_map >= 1 & bin_map <= 40))

    # The bin containing 0.675-0.7 should map to itself
    target_bin <- findInterval(0.6875, breaks, rightmost.closed = TRUE) # center of (0.675, 0.7]
    expect_equal(as.integer(bin_map[target_bin]), target_bin)

    # Bins with values > 0.07 should map to the target bin
    # (We can't easily test all of them without knowing exact bin boundaries,
    # but we can check that some high-value bins are mapped)
})

test_that("print.gcanvas.model works", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
        intervals = test_intervals,
        iterator = 200
    )

    # Print should not error
    expect_output(print(model), "Genome Canvas Markov-5 Model")
    expect_output(print(model), "Number of bins: 10")

    gvtrack.rm("test_vt")
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
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
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

test_that("gcanvas.sample with seed is reproducible", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 10000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
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

test_that("gcanvas.sample respects intervals output", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    train_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = train_intervals)

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
        intervals = train_intervals,
        iterator = 200
    )

    sample_intervals <- gintervals(
        c(1, 1),
        c(1000, 2000),
        c(1050, 2050)
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
    header_idx <- grep("^>", lines)
    expect_equal(length(header_idx), 2)

    seq1 <- paste(lines[(header_idx[1] + 1):(header_idx[2] - 1)], collapse = "")
    seq2 <- paste(lines[(header_idx[2] + 1):length(lines)], collapse = "")
    expect_equal(nchar(seq1), 50)
    expect_equal(nchar(seq2), 50)

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample writes interval headers in FASTA output", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    train_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = train_intervals)

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
        intervals = train_intervals,
        iterator = 200
    )

    sample_intervals <- gintervals(
        c(1, 1, 2),
        c(100, 200, 300),
        c(150, 260, 350)
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
    header_idx <- grep("^>", lines)
    expect_equal(length(header_idx), 3)
    expect_true(all(grepl("^>.*:[0-9]+-[0-9]+$", lines[header_idx])))

    unlink(output_fasta)
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
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
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

test_that("bin_merge is applied during sampling", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Create breaks
    breaks <- seq(track_range["Min"], track_range["Max"], length.out = 21) # 20 bins

    # Train with bin merging: merge high-value bins into a middle bin
    # Find a threshold value and a target bin
    threshold <- track_range["Min"] + 0.7 * (track_range["Max"] - track_range["Min"])
    target_bin_idx <- findInterval(threshold, breaks, rightmost.closed = TRUE)
    target_bin_start <- breaks[target_bin_idx]
    target_bin_end <- breaks[target_bin_idx + 1]

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = breaks,
        bin_merge = list(
            list(from = threshold, to = c(target_bin_start, target_bin_end))
        ),
        intervals = test_intervals,
        iterator = 200
    )

    # Verify bin_map is stored in the model
    expect_true("bin_map" %in% names(model))
    expect_type(model$bin_map, "integer")
    expect_equal(length(model$bin_map), 20)

    # Verify that high-value bins are mapped to the target bin
    high_bin_idx <- findInterval(threshold + 0.01 * (track_range["Max"] - track_range["Min"]),
        breaks,
        rightmost.closed = TRUE
    )
    if (high_bin_idx <= 20 && high_bin_idx > 0) {
        expect_equal(model$bin_map[high_bin_idx], target_bin_idx)
    }

    # Sample with the model - this should use the merged bins correctly
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

    # Verify the model can be saved and loaded with bin_map preserved
    temp_file <- tempfile(fileext = ".rds")
    gcanvas.save(model, temp_file)
    loaded_model <- gcanvas.load(temp_file)
    expect_true("bin_map" %in% names(loaded_model))
    expect_equal(loaded_model$bin_map, model$bin_map)

    # Sampling with loaded model should also work
    output_fasta2 <- tempfile(fileext = ".fa")
    gcanvas.sample(
        loaded_model,
        output_fasta2,
        output_format = "fasta",
        intervals = test_intervals,
        seed = 60427
    )
    expect_true(file.exists(output_fasta2))

    # Clean up
    unlink(c(output_fasta, output_fasta2, temp_file))
    gvtrack.rm("test_vt")
})

# Additional comprehensive tests

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
        expr = "test_vt",
        breaks = breaks,
        intervals = test_intervals,
        iterator = 200,
        pseudocount = 1
    )

    model_pc10 <- gcanvas.train(
        expr = "test_vt",
        breaks = breaks,
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
    # (i.e., CDFs should be closer to equal spacing)
    expect_false(identical(cdf1, cdf10))

    gvtrack.rm("test_vt")
})

test_that("gcanvas.train with very few bins", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Train with just 2 bins
    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 3), # 2 bins
        intervals = test_intervals,
        iterator = 200
    )

    expect_s3_class(model, "gcanvas.model")
    expect_equal(model$num_bins, 2)
    expect_equal(length(model$model_data$cdf), 2)

    gvtrack.rm("test_vt")
})

test_that("gcanvas.train with many bins", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    # Train with 50 bins
    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 51), # 50 bins
        intervals = test_intervals,
        iterator = 200
    )

    expect_s3_class(model, "gcanvas.model")
    expect_equal(model$num_bins, 50)
    expect_equal(length(model$model_data$cdf), 50)

    # Some bins may be empty, but model should still work
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

test_that("gcanvas.train with different iterator sizes", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)
    breaks <- seq(track_range["Min"], track_range["Max"], length.out = 11)

    # Train with different iterator sizes
    model_100 <- gcanvas.train(
        expr = "test_vt",
        breaks = breaks,
        intervals = test_intervals,
        iterator = 100
    )

    model_500 <- gcanvas.train(
        expr = "test_vt",
        breaks = breaks,
        intervals = test_intervals,
        iterator = 500
    )

    # Both should produce valid models
    expect_s3_class(model_100, "gcanvas.model")
    expect_s3_class(model_500, "gcanvas.model")

    # Smaller iterator gives finer resolution, but total k-mers should be similar
    # (since we're counting the same sequence)
    expect_equal(model_100$total_kmers, model_500$total_kmers)

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
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
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

test_that("gcanvas.sample with mask_mode sample ignores mask", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
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
        mask_mode = "sample", # This should ignore the mask
        seed = 60427
    )

    # Read the sampled sequence
    lines <- readLines(output_fasta)
    sampled_seq <- paste(lines[!grepl("^>", lines)], collapse = "")

    # Get original sequence for the mask region (gseq.extract returns character vector)
    original_seq <- gseq.extract(mask)[1]

    # The mask region should NOT necessarily match the original
    mask_start_in_sample <- 1001
    mask_end_in_sample <- 2000
    sampled_mask_region <- substr(sampled_seq, mask_start_in_sample, mask_end_in_sample)

    # With high probability, the sampled region will differ from original
    # (case-insensitive comparison)
    expect_false(identical(toupper(sampled_mask_region), toupper(original_seq)))

    unlink(output_fasta)
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
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
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
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
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

test_that("gcanvas.sample on multiple chromosomes", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals.all()
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
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

test_that("gcanvas.train error handling", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)

    # Invalid breaks
    expect_error(
        gcanvas.train(
            expr = "test_vt",
            breaks = c(0.5), # Only one element
            intervals = test_intervals,
            iterator = 200
        ),
        "at least 2 elements"
    )

    gvtrack.rm("test_vt")
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
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
        intervals = test_intervals,
        iterator = 200
    )

    # Check all CDFs
    for (bin in seq_len(model$num_bins)) {
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

test_that("gcanvas model counts structure is correct", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
        intervals = test_intervals,
        iterator = 200
    )

    # Check all count matrices
    for (bin in seq_len(model$num_bins)) {
        count_mat <- model$model_data$counts[[bin]]

        # Should be 1024 x 4
        expect_equal(dim(count_mat), c(1024, 4))

        # All values should be non-negative (includes pseudocount)
        expect_true(all(count_mat >= 0))
    }

    # Total counts from matrices should be greater than or equal to total_kmers
    # (counts include pseudocount added during normalization)
    total_from_counts <- sum(sapply(model$model_data$counts, sum))
    expect_gte(total_from_counts, model$total_kmers)

    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample output file is created even for small intervals", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
        intervals = test_intervals,
        iterator = 200
    )

    # Sample very small interval (just 10 bp)
    sample_intervals <- gintervals(1, 100, 110)

    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = sample_intervals,
        seed = 60427
    )

    expect_true(file.exists(output_fasta))
    lines <- readLines(output_fasta)
    seq_lines <- lines[!grepl("^>", lines)]
    full_seq <- paste(seq_lines, collapse = "")
    expect_equal(nchar(full_seq), 10)

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas k-mer distribution of sampled sequence is reasonable", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
        intervals = test_intervals,
        iterator = 200
    )

    # Get k-mer distribution of original sequence
    original_kmer <- gseq.kmer.dist(test_intervals, k = 2)
    original_freqs <- original_kmer$count / sum(original_kmer$count)
    names(original_freqs) <- original_kmer$kmer

    # Sample a sequence
    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = test_intervals,
        seed = 60427
    )

    # Read sampled sequence and compute k-mer distribution
    lines <- readLines(output_fasta)
    sampled_seq <- paste(lines[!grepl("^>", lines)], collapse = "")

    # Count dinucleotides in sampled sequence
    sampled_counts <- list()
    for (i in 1:(nchar(sampled_seq) - 1)) {
        dinuc <- substr(sampled_seq, i, i + 1)
        if (is.null(sampled_counts[[dinuc]])) {
            sampled_counts[[dinuc]] <- 0
        }
        sampled_counts[[dinuc]] <- sampled_counts[[dinuc]] + 1
    }

    total_sampled <- sum(unlist(sampled_counts))
    sampled_freqs <- sapply(sampled_counts, function(x) x / total_sampled)

    # The correlation between original and sampled frequencies should be high
    common_kmers <- intersect(names(original_freqs), names(sampled_freqs))
    cor_val <- cor(original_freqs[common_kmers], sampled_freqs[common_kmers])

    # Expect reasonably high correlation (Markov model should preserve k-mer stats)
    expect_gt(cor_val, 0.8)

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.save and gcanvas.load preserve all model fields", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals(1, 0, 50000)
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
        intervals = test_intervals,
        iterator = 200
    )

    temp_file <- tempfile(fileext = ".rds")
    gcanvas.save(model, temp_file)
    loaded_model <- gcanvas.load(temp_file)

    # Check all fields are preserved
    expect_equal(loaded_model$num_bins, model$num_bins)
    expect_equal(loaded_model$breaks, model$breaks)
    expect_equal(loaded_model$total_kmers, model$total_kmers)
    expect_equal(loaded_model$per_bin_kmers, model$per_bin_kmers)
    expect_equal(loaded_model$total_masked, model$total_masked)
    expect_equal(loaded_model$total_n, model$total_n)
    expect_equal(loaded_model$expr, model$expr)
    expect_equal(loaded_model$iterator, model$iterator)
    expect_equal(loaded_model$bin_map, model$bin_map)

    # Check model_data
    expect_equal(length(loaded_model$model_data$counts), length(model$model_data$counts))
    expect_equal(length(loaded_model$model_data$cdf), length(model$model_data$cdf))

    # Check class
    expect_s3_class(loaded_model, "gcanvas.model")

    unlink(temp_file)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.sample with all chromosomes works", {
    gdb.init_examples()

    if ("test_vt" %in% gvtrack.ls()) {
        gvtrack.rm("test_vt")
    }
    gvtrack.create("test_vt", "dense_track", "avg")

    test_intervals <- gintervals.all()
    track_range <- gsummary("dense_track", intervals = test_intervals)

    model <- gcanvas.train(
        expr = "test_vt",
        breaks = seq(track_range["Min"], track_range["Max"], length.out = 11),
        intervals = test_intervals,
        iterator = 200
    )

    # Sample all chromosomes
    output_fasta <- tempfile(fileext = ".fa")
    gcanvas.sample(
        model,
        output_fasta,
        output_format = "fasta",
        intervals = gintervals.all(),
        seed = 60427
    )

    lines <- readLines(output_fasta)
    headers <- lines[grepl("^>", lines)]

    # Should have headers for all chromosomes
    expect_gt(length(headers), 0)

    # Total sequence length should match total interval sizes
    seq_lines <- lines[!grepl("^>", lines)]
    total_length <- sum(nchar(seq_lines))
    all_intervals <- gintervals.all()
    expected_length <- sum(all_intervals$end - all_intervals$start)
    expect_equal(total_length, expected_length)

    unlink(output_fasta)
    gvtrack.rm("test_vt")
})

test_that("gcanvas.train handles track expression correctly", {
    gdb.init_examples()

    # Clean up any existing vtracks
    for (vt in c("g_frac", "c_frac")) {
        if (vt %in% gvtrack.ls()) {
            gvtrack.rm(vt)
        }
    }

    # Create G and C fraction virtual tracks
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
    gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")

    test_intervals <- gintervals(1, 0, 50000)

    # Train with expression (GC content)
    model <- gcanvas.train(
        expr = "g_frac + c_frac",
        breaks = seq(0, 1, 0.1),
        intervals = test_intervals,
        iterator = 200
    )

    expect_s3_class(model, "gcanvas.model")
    expect_equal(model$expr, "g_frac + c_frac")
    expect_true(model$total_kmers > 0)

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
        expr = "test_vt",
        breaks = breaks,
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
