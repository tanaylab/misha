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
