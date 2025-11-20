load_test_db()
# Helper function to manually count kmer occurrences in a sequence
# This matches the C++ implementation's approach
count_kmers_manually <- function(sequence, kmer, original_interval, extended_interval = NULL) {
    kmer <- toupper(kmer)
    sequence <- toupper(sequence)

    # Calculate where in the extended sequence the original interval starts
    original_offset <- 0
    if (!is.null(extended_interval) && extended_interval$start < original_interval$start) {
        original_offset <- original_interval$start - extended_interval$start
    }

    # Calculate the end position in the extended sequence that corresponds
    # to the end of the original interval
    original_end_pos <- nchar(sequence)
    if (!is.null(extended_interval) && extended_interval$end > original_interval$end) {
        original_end_pos <- original_end_pos - (extended_interval$end - original_interval$end)
    }

    # Count kmers whose starting position falls within the original interval
    count <- 0
    kmer_length <- nchar(kmer)

    # Only search positions within the original interval bounds
    for (i in (original_offset + 1):(original_offset + original_interval$end - original_interval$start)) {
        # Make sure we don't go out of bounds
        if (i + kmer_length - 1 <= nchar(sequence)) {
            if (substr(sequence, i, i + kmer_length - 1) == kmer) {
                count <- count + 1
            }
        }
    }

    # Calculate possible positions for fraction
    possible_positions <- original_interval$end - original_interval$start - kmer_length + 1
    possible_positions <- max(0, possible_positions)

    return(list(
        count = count,
        possible_positions = possible_positions,
        fraction = if (possible_positions > 0) count / possible_positions else 0
    ))
}

test_that("kmer.count and kmer.frac vtracks work with basic inputs", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals)) # CCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCC

    # Create virtual tracks with different k-mers
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA", strand = 1)
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA", strand = 1)
    gvtrack.create("count_ccc", NULL, "kmer.count", kmer = "CCC", strand = 1)
    gvtrack.create("frac_ccc", NULL, "kmer.frac", kmer = "CCC", strand = 1)

    # Extract scores for the entire interval
    scores <- gextract(c("count_ta", "frac_ta", "count_ccc", "frac_ccc"), test_intervals, iterator = test_intervals)

    # Retrieve extended sequence to account for extension
    extended_interval <- test_intervals
    extended_interval$start <- extended_interval$start - 1 # For 2-letter kmer, extend by 1 on each side
    extended_interval$end <- extended_interval$end + 1
    extended_seq <- toupper(gseq.extract(extended_interval))

    # Use our helper function for "TA"
    ta_results <- count_kmers_manually(extended_seq, "TA", test_intervals, extended_interval)

    # Use our helper function for "CCC"
    ccc_results <- count_kmers_manually(extended_seq, "CCC", test_intervals, extended_interval)

    # Test that counts match expected values
    expect_equal(scores$count_ta, ta_results$count)
    expect_equal(scores$frac_ta, ta_results$fraction, tolerance = 1e-5)
    expect_equal(scores$count_ccc, ccc_results$count)
    expect_equal(scores$frac_ccc, ccc_results$fraction, tolerance = 1e-5)

    # Manually verify the "CCC" count by looking at positions in the sequence
    # This ensures that "CCCC" is correctly counted as two occurrences of "CCC"
    manual_ccc_count <- 0
    for (i in 1:(nchar(seq) - 2)) {
        if (substr(seq, i, i + 2) == "CCC") {
            manual_ccc_count <- manual_ccc_count + 1
        }
    }
    expect_equal(scores$count_ccc, manual_ccc_count)
})

test_that("kmer honors gvtrack.iterator shifts for count and frac", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 2000, 2040)
    k <- "TA" # 2-mer

    # Create vtracks (forward strand) and apply different iterator shifts
    gvtrack.create("count_small", NULL, "kmer.count", kmer = k, strand = 1)
    gvtrack.create("count_large", NULL, "kmer.count", kmer = k, strand = 1)
    gvtrack.create("frac_small", NULL, "kmer.frac", kmer = k, strand = 1)
    gvtrack.create("frac_large", NULL, "kmer.frac", kmer = k, strand = 1)

    gvtrack.iterator("count_small", sshift = -10, eshift = 10)
    gvtrack.iterator("frac_small", sshift = -10, eshift = 10)
    gvtrack.iterator("count_large", sshift = -1000, eshift = 1000)
    gvtrack.iterator("frac_large", sshift = -1000, eshift = 1000)

    scores <- gextract(c("count_small", "count_large", "frac_small", "frac_large"), test_intervals, iterator = test_intervals)

    # Build intervals matching how the engine applies per-vtrack iterator modifiers
    # Small window: iterator shift of 10 on each side => original (per-vtrack) interval is shifted
    small_orig <- test_intervals
    small_orig$start <- pmax(0, small_orig$start - 10)
    small_orig$end <- small_orig$end + 10
    # Fetch interval adds extension by (k-1) at the end for forward strand
    small_fetch <- small_orig
    small_fetch$end <- small_fetch$end + (nchar(k) - 1)
    seq_small <- toupper(gseq.extract(small_fetch))
    res_small <- count_kmers_manually(seq_small, k, small_orig, small_fetch)

    expect_equal(scores$count_small, res_small$count)
    expect_equal(scores$frac_small, res_small$fraction, tolerance = 1e-5)

    # Large window: iterator shift 1000 each side
    large_orig <- test_intervals
    large_orig$start <- pmax(0, large_orig$start - 1000)
    large_orig$end <- large_orig$end + 1000
    large_fetch <- large_orig
    large_fetch$end <- large_fetch$end + (nchar(k) - 1)
    seq_large <- toupper(gseq.extract(large_fetch))
    res_large <- count_kmers_manually(seq_large, k, large_orig, large_fetch)

    expect_equal(scores$count_large, res_large$count)
    expect_equal(scores$frac_large, res_large$fraction, tolerance = 1e-5)

    # The small and large windows should generally give different totals
    expect_true(scores$count_small != scores$count_large || abs(scores$frac_small - scores$frac_large) > 1e-12)

    # Asymmetric shift case: start-only and end-only shifts
    gvtrack.create("count_asym", NULL, "kmer.count", kmer = k, strand = 1)
    gvtrack.iterator("count_asym", sshift = -25, eshift = 5)
    asym_scores <- gextract("count_asym", test_intervals, iterator = test_intervals)

    asym_orig <- test_intervals
    asym_orig$start <- pmax(0, asym_orig$start - 25)
    asym_orig$end <- asym_orig$end + 5
    asym_fetch <- asym_orig
    asym_fetch$end <- asym_fetch$end + (nchar(k) - 1)
    seq_asym <- toupper(gseq.extract(asym_fetch))
    res_asym <- count_kmers_manually(seq_asym, k, asym_orig, asym_fetch)
    expect_equal(asym_scores$count_asym, res_asym$count)
})

test_that("kmer.count and kmer.frac work with smaller iterator intervals", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals)) # CCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCC

    # Create virtual tracks
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA", strand = 1)
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA", strand = 1)

    # Extract scores with 10bp iterator
    scores_10bp <- gextract(c("count_ta", "frac_ta"), test_intervals, iterator = 10)

    # Test each 10bp segment separately using our helper function
    expected_counts <- numeric(4)
    expected_fracs <- numeric(4)

    for (i in 1:4) {
        bin_interval <- gintervals(1, 200 + (i - 1) * 10, 200 + i * 10)

        # Create extended interval for correct comparison
        ext_bin_interval <- bin_interval
        ext_bin_interval$start <- ext_bin_interval$start - 1
        ext_bin_interval$end <- ext_bin_interval$end + 1
        ext_seq <- toupper(gseq.extract(ext_bin_interval))

        # Count kmers in this segment
        results <- count_kmers_manually(ext_seq, "TA", bin_interval, ext_bin_interval)
        expected_counts[i] <- results$count
        expected_fracs[i] <- results$fraction
    }

    # Test count values for each bin
    expect_equal(scores_10bp$count_ta, expected_counts)
    expect_equal(scores_10bp$frac_ta, expected_fracs, tolerance = 1e-5)
})

test_that("kmer.count and kmer.frac handle case sensitivity correctly", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals))

    # Create virtual tracks with different case patterns
    gvtrack.create("count_ta_upper", NULL, "kmer.count", kmer = "TA", strand = 1)
    gvtrack.create("count_ta_lower", NULL, "kmer.count", kmer = "ta", strand = 1)
    gvtrack.create("count_ta_mixed", NULL, "kmer.count", kmer = "Ta", strand = 1)

    # Extract scores
    scores <- gextract(c("count_ta_upper", "count_ta_lower", "count_ta_mixed"),
        test_intervals,
        iterator = test_intervals
    )

    # Extended interval for comparison
    extended_interval <- test_intervals
    extended_interval$start <- extended_interval$start - 1
    extended_interval$end <- extended_interval$end + 1
    extended_seq <- toupper(gseq.extract(extended_interval))

    # Use helper function to calculate expected count
    expected <- count_kmers_manually(extended_seq, "TA", test_intervals, extended_interval)

    # All should give the same results (case-insensitive matching)
    expect_equal(scores$count_ta_upper, expected$count)
    expect_equal(scores$count_ta_lower, expected$count)
    expect_equal(scores$count_ta_mixed, expected$count)
})

test_that("kmer.count and kmer.frac handle longer k-mers correctly", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)

    # Create virtual tracks with longer k-mers
    gvtrack.create("count_taac", NULL, "kmer.count", kmer = "TAAC", strand = 1)
    gvtrack.create("frac_taac", NULL, "kmer.frac", kmer = "TAAC", strand = 1)
    gvtrack.create("count_ccctaa", NULL, "kmer.count", kmer = "CCCTAA", strand = 1)
    gvtrack.create("frac_ccctaa", NULL, "kmer.frac", kmer = "CCCTAA", strand = 1)

    # Extract scores
    scores <- gextract(c("count_taac", "frac_taac", "count_ccctaa", "frac_ccctaa"),
        test_intervals,
        iterator = test_intervals
    )

    # Create extended intervals with appropriate extension for each kmer length
    extended_interval_taac <- test_intervals
    extended_interval_taac$start <- extended_interval_taac$start - 3 # extension for 4-letter kmer
    extended_interval_taac$end <- extended_interval_taac$end + 3
    extended_seq_taac <- toupper(gseq.extract(extended_interval_taac))

    extended_interval_ccctaa <- test_intervals
    extended_interval_ccctaa$start <- extended_interval_ccctaa$start - 5 # extension for 6-letter kmer
    extended_interval_ccctaa$end <- extended_interval_ccctaa$end + 5
    extended_seq_ccctaa <- toupper(gseq.extract(extended_interval_ccctaa))

    # Calculate expected results using our helper function
    taac_results <- count_kmers_manually(extended_seq_taac, "TAAC", test_intervals, extended_interval_taac)
    ccctaa_results <- count_kmers_manually(extended_seq_ccctaa, "CCCTAA", test_intervals, extended_interval_ccctaa)

    # Test that counts match expected values
    expect_equal(scores$count_taac, taac_results$count)
    expect_equal(scores$frac_taac, taac_results$fraction, tolerance = 1e-5)
    expect_equal(scores$count_ccctaa, ccctaa_results$count)
    expect_equal(scores$frac_ccctaa, ccctaa_results$fraction, tolerance = 1e-5)
})

test_that("kmer.count and kmer.frac handle edge cases correctly", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals))

    # Test k-mer length equal to sequence length
    gvtrack.create("count_full", NULL, "kmer.count", kmer = seq, strand = 1)
    gvtrack.create("frac_full", NULL, "kmer.frac", kmer = seq, strand = 1)

    # Test 1bp k-mer
    gvtrack.create("count_a", NULL, "kmer.count", kmer = "A", strand = 1)
    gvtrack.create("frac_a", NULL, "kmer.frac", kmer = "A", strand = 1)

    # Extract scores
    scores <- gextract(c("count_full", "frac_full", "count_a", "frac_a"),
        test_intervals,
        iterator = test_intervals
    )

    # For full-length k-mer - there should be only one occurrence and one possible position
    # For single base k-mer - use our helper function with appropriate extension
    extended_interval_a <- test_intervals
    extended_interval_a$start <- extended_interval_a$start - 0 # No need to extend for 1-letter kmer
    extended_interval_a$end <- extended_interval_a$end + 0
    extended_seq_a <- toupper(gseq.extract(extended_interval_a))

    a_results <- count_kmers_manually(extended_seq_a, "A", test_intervals, extended_interval_a)

    # Test results - full sequence should have exactly one occurrence
    expect_equal(scores$count_full, 1)
    expect_equal(scores$frac_full, 1)

    # Single letter kmer
    expect_equal(scores$count_a, a_results$count)
    expect_equal(scores$frac_a, a_results$fraction, tolerance = 1e-5)
})

test_that("kmer.count and kmer.frac handle no matches correctly", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)

    # K-mer that doesn't exist in the sequence
    gvtrack.create("count_xyz", NULL, "kmer.count", kmer = "XYZ", strand = 1)
    gvtrack.create("frac_xyz", NULL, "kmer.frac", kmer = "XYZ", strand = 1)

    # Extract scores
    scores <- gextract(c("count_xyz", "frac_xyz"), test_intervals, iterator = test_intervals)

    # Should return 0 for count and 0 for fraction
    expect_equal(scores$count_xyz, 0)
    expect_equal(scores$frac_xyz, 0)
})

test_that("kmer.count and kmer.frac work with different iterator sizes", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)

    # Create virtual tracks
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA", strand = 1)
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA", strand = 1)

    # Extract scores with different iterator sizes
    scores_5bp <- gextract(c("count_ta", "frac_ta"), test_intervals, iterator = 5)
    scores_10bp <- gextract(c("count_ta", "frac_ta"), test_intervals, iterator = 10)
    scores_20bp <- gextract(c("count_ta", "frac_ta"), test_intervals, iterator = 20)

    # We can't easily calculate exact expected values for each bin size with extensions
    # So instead check that the sum of counts is consistent
    expect_equal(sum(scores_5bp$count_ta), sum(scores_10bp$count_ta))
    expect_equal(sum(scores_10bp$count_ta), sum(scores_20bp$count_ta))

    # The average of fractions should be similar (but not exactly equal due to different bin sizes)
    # but with the extension, we just check they're all reasonable values
    expect_true(all(scores_5bp$frac_ta >= 0 & scores_5bp$frac_ta <= 1))
    expect_true(all(scores_10bp$frac_ta >= 0 & scores_10bp$frac_ta <= 1))
    expect_true(all(scores_20bp$frac_ta >= 0 & scores_20bp$frac_ta <= 1))
})

test_that("kmer.count and kmer.frac work with overlapping k-mers", {
    remove_all_vtracks()

    # Create a test interval with repeating sequence
    test_intervals <- gintervals(1, 200, 220)

    # Create virtual tracks for overlapping k-mer
    gvtrack.create("count_taa", NULL, "kmer.count", kmer = "TAA", strand = 1)
    gvtrack.create("frac_taa", NULL, "kmer.frac", kmer = "TAA", strand = 1)

    # Extract scores
    scores <- gextract(c("count_taa", "frac_taa"), test_intervals, iterator = test_intervals)

    # Get extended sequence to account for kmer extension
    extended_interval <- test_intervals
    extended_interval$start <- extended_interval$start - 2 # For 3-letter kmer
    extended_interval$end <- extended_interval$end + 2
    extended_seq <- toupper(gseq.extract(extended_interval))

    # Original sequence
    seq <- toupper(gseq.extract(test_intervals))

    # Manual count with extension
    taa_count <- 0
    original_offset <- 2
    for (i in original_offset:(original_offset + nchar(seq) - nchar("TAA"))) {
        if (substr(extended_seq, i, i + nchar("TAA") - 1) == "TAA") {
            taa_count <- taa_count + 1
        }
    }
    taa_possible_positions <- nchar(seq) - nchar("TAA") + 1

    # Test counts
    expect_equal(scores$count_taa, taa_count)
    expect_equal(scores$frac_taa, taa_count / taa_possible_positions, tolerance = 1e-5)
})

test_that("kmer.count and kmer.frac handle interval boundaries correctly", {
    remove_all_vtracks()

    # Create a longer interval and sub-intervals
    full_interval <- gintervals(1, 200, 240)
    sub_interval1 <- gintervals(1, 200, 220)
    sub_interval2 <- gintervals(1, 220, 240)

    # Create virtual tracks
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA", strand = 1)
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA", strand = 1)

    # Extract scores for all intervals
    full_scores <- gextract(c("count_ta", "frac_ta"), full_interval, iterator = full_interval)
    sub1_scores <- gextract(c("count_ta", "frac_ta"), sub_interval1, iterator = sub_interval1)
    sub2_scores <- gextract(c("count_ta", "frac_ta"), sub_interval2, iterator = sub_interval2)

    # The sum of sub-interval counts doesn't necessarily equal the full interval count
    # because of the extension at boundaries, but they should be close
    # We'll check the total count is within a reasonable range
    expect_true(abs(full_scores$count_ta - (sub1_scores$count_ta + sub2_scores$count_ta)) <= 2)
})

test_that("kmer.count and kmer.frac handle multiple chromosomes", {
    remove_all_vtracks()

    # Create intervals on different chromosomes
    interval_chr1 <- gintervals(1, 200, 220)
    interval_chr2 <- gintervals(2, 200, 220)

    # Create virtual tracks
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA", strand = 1)
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA", strand = 1)

    # Extract scores for both chromosomes
    scores_chr1 <- gextract(c("count_ta", "frac_ta"), interval_chr1, iterator = interval_chr1)
    scores_chr2 <- gextract(c("count_ta", "frac_ta"), interval_chr2, iterator = interval_chr2)

    # Both should return valid results (not checking specific values,
    # just ensuring proper functionality across chromosomes)
    expect_true(is.numeric(scores_chr1$count_ta))
    expect_true(is.numeric(scores_chr1$frac_ta))
    expect_true(is.numeric(scores_chr2$count_ta))
    expect_true(is.numeric(scores_chr2$frac_ta))
})

test_that("kmer.count and kmer.frac return correct results with gdist", {
    remove_all_vtracks()

    # Create test interval
    test_interval <- gintervals(1, 200, 300)

    # Create virtual tracks
    gvtrack.create("count_a", NULL, "kmer.count", kmer = "A", strand = 1)
    gvtrack.create("frac_a", NULL, "kmer.frac", kmer = "A", strand = 1)

    # Use gdist to bin k-mer counts
    count_dist <- gdist("count_a",
        breaks = 0:10,
        intervals = test_interval, iterator = 10
    )

    # Use gdist to bin k-mer fractions
    frac_dist <- gdist("frac_a",
        breaks = seq(0, 0.5, by = 0.1),
        intervals = test_interval, iterator = 10
    )

    # Verify distribution shape (specific values not tested)
    expect_true(all(count_dist >= 0))
    expect_true(all(frac_dist >= 0))
    expect_equal(length(count_dist), 10)
    expect_equal(length(frac_dist), 5)
})

test_that("kmer.count and kmer.frac preserve position information", {
    remove_all_vtracks()

    # Create test interval
    test_interval <- gintervals(1, 200, 300)

    # Create virtual tracks
    gvtrack.create("count_ccc", NULL, "kmer.count", kmer = "CCC", strand = 1)

    # Extract with position information
    scores <- gextract("count_ccc",
        test_interval,
        iterator = 10
    ) %>%
        arrange(intervalID)


    # Verify that position columns are preserved
    expect_true("chrom" %in% colnames(scores))
    expect_true("start" %in% colnames(scores))
    expect_true("end" %in% colnames(scores))
    expect_true("count_ccc" %in% colnames(scores))

    # Check that positions make sense
    expect_equal(scores$start[1], 200)
    expect_equal(scores$end[10], 300)
    expect_equal(nrow(scores), 10)
})

test_that("kmer.count and kmer.frac work with other virtual track functions", {
    remove_all_vtracks()

    # Create test interval
    test_interval <- gintervals(1, 200, 300)

    # Create different virtual tracks on the same sequence
    gvtrack.create("count_a", NULL, "kmer.count", kmer = "A", strand = 1)
    gvtrack.create("frac_a", NULL, "kmer.frac", kmer = "A", strand = 1)

    # Create a distance virtual track
    annotations <- gintervals(1, 220, 240)
    gvtrack.create("dist_track", annotations, "distance", strand = 1)

    # Extract all tracks together
    scores <- gextract(c("count_a", "frac_a", "dist_track"),
        test_interval,
        iterator = 20
    )

    # Verify that all tracks return results
    expect_true(all(c("count_a", "frac_a", "dist_track") %in% colnames(scores)))
    expect_equal(nrow(scores), 5)

    # The count and fraction should be correlated
    correlation <- cor(scores$count_a, scores$frac_a, use = "complete.obs")
    expect_true(correlation > 0.8)
})

test_that("kmer.count and kmer.frac parameter validation works", {
    remove_all_vtracks()

    # Empty k-mer should error
    expect_error(gvtrack.create("empty_kmer", NULL, "kmer.count", kmer = "", strand = 1))

    # Non-string parameter should error
    expect_error(gvtrack.create("numeric_kmer", NULL, "kmer.count", kmer = 123, strand = 1))

    # No parameter should error
    expect_error(gvtrack.create("no_param", NULL, "kmer.count", strand = 1))

    # Multiple k-mers should error
    expect_error(gvtrack.create("multi_kmer", NULL, "kmer.count", kmer = c("AA", "TT"), strand = 1))
})

test_that("kmer.count and kmer.frac work with explicit extension parameters", {
    remove_all_vtracks()

    # Create a test interval
    test_interval <- gintervals(1, 200, 220)

    # Create virtual tracks with different extension settings
    gvtrack.create("count_default", NULL, "kmer.count", kmer = "TA", strand = 1) # Default extension=TRUE
    gvtrack.create("count_extend", NULL, "kmer.count", kmer = "TA", extend = TRUE, strand = 1)
    gvtrack.create("count_no_extend", NULL, "kmer.count", kmer = "TA", extend = FALSE, strand = 1)

    # Create fraction tracks with the same settings
    gvtrack.create("frac_default", NULL, "kmer.frac", kmer = "TA", strand = 1)
    gvtrack.create("frac_extend", NULL, "kmer.frac", kmer = "TA", extend = TRUE, strand = 1)
    gvtrack.create("frac_no_extend", NULL, "kmer.frac", kmer = "TA", extend = FALSE, strand = 1)

    # Extract scores
    results <- gextract(
        c(
            "count_default", "count_extend", "count_no_extend",
            "frac_default", "frac_extend", "frac_no_extend"
        ),
        test_interval,
        iterator = test_interval
    )

    # Default and explicit extend should be the same
    expect_equal(results$count_default, results$count_extend)
    expect_equal(results$frac_default, results$frac_extend)

    # No extension should potentially return different results
    # Depending on the test interval, they might be different or the same
    # Here we just test that they're all valid numeric values
    expect_true(is.numeric(results$count_no_extend))
    expect_true(is.numeric(results$frac_no_extend))
})

test_that("kmer.count and kmer.frac extension handles interval boundaries properly", {
    remove_all_vtracks()

    # Create test intervals that span sequence boundaries
    start_interval <- gintervals(1, 0, 10) # At the beginning of chromosome
    end_interval <- gintervals(1, 9990, 10000) # Near the end of chromosome

    # Create virtual tracks with different extension modes
    gvtrack.create("count_with_ext", NULL, "kmer.count", kmer = "ACGT", extend = TRUE, strand = 1)
    gvtrack.create("count_no_ext", NULL, "kmer.count", kmer = "ACGT", extend = FALSE, strand = 1)

    # Extract scores
    start_results <- gextract(c("count_with_ext", "count_no_ext"), start_interval, iterator = start_interval)
    end_results <- gextract(c("count_with_ext", "count_no_ext"), end_interval, iterator = end_interval)

    # We expect valid results even at chromosome boundaries
    expect_true(is.numeric(start_results$count_with_ext))
    expect_true(is.numeric(start_results$count_no_ext))
    expect_true(is.numeric(end_results$count_with_ext))
    expect_true(is.numeric(end_results$count_no_ext))

    # Extension should be properly limited by chromosome boundaries
    # This test just ensures we don't crash or get NA values
    expect_true(!is.na(start_results$count_with_ext))
    expect_true(!is.na(end_results$count_with_ext))
})

test_that("kmer.count and kmer.frac performance is acceptable", {
    skip_on_cran()

    # Create large test interval
    test_interval <- gintervals(1, 1, 10000)

    # Create virtual tracks with different k-mer lengths
    gvtrack.create("count_2mer", NULL, "kmer.count", kmer = "AT", strand = 1)
    gvtrack.create("count_5mer", NULL, "kmer.count", kmer = "ATGCG", strand = 1)
    gvtrack.create("count_10mer", NULL, "kmer.count", kmer = "ATGCGCGTAA", strand = 1)

    # Measure timing with different k-mer lengths
    time_2mer <- system.time(gextract("count_2mer", test_interval, iterator = 1000))
    time_5mer <- system.time(gextract("count_5mer", test_interval, iterator = 1000))
    time_10mer <- system.time(gextract("count_10mer", test_interval, iterator = 1000))

    # Performance expectations (adjust threshold as needed)
    expect_true(time_2mer["elapsed"] < 10)
    expect_true(time_5mer["elapsed"] < 10)
    expect_true(time_10mer["elapsed"] < 10)

    # Longer k-mers should not be dramatically slower
    expect_true(time_10mer["elapsed"] < time_2mer["elapsed"] * 3)
})

test_that("kmer.count and kmer.frac handle overlapping k-mers correctly", {
    remove_all_vtracks()

    # Create a test interval with known content
    test_interval <- gintervals(1, 200, 210)
    seq <- toupper(gseq.extract(test_interval)) # Should be "CCCTAACCCT"

    # Check for overlapping CCC occurrences (there should be 2 in "CCCTAACCCT")
    gvtrack.create("count_ccc", NULL, "kmer.count", kmer = "CCC", strand = 1)

    # Extract score
    scores <- gextract("count_ccc", test_interval, iterator = test_interval)

    # Create extended interval
    extended_interval <- test_interval
    extended_interval$start <- extended_interval$start - 2
    extended_interval$end <- extended_interval$end + 2
    extended_seq <- toupper(gseq.extract(extended_interval))

    # Manually count CCC occurrences
    manual_count <- 0
    for (i in 3:(3 + nchar(seq) - 3)) {
        if (substr(extended_seq, i, i + 2) == "CCC") {
            manual_count <- manual_count + 1
        }
    }

    # This should be 2
    expect_equal(scores$count_ccc, manual_count)
})

test_that("kmer.count and kmer.frac handle strand parameter correctly", {
    remove_all_vtracks()

    # Create a test interval with a known palindromic sequence
    # Let's use "ACGT" which is its own reverse complement
    test_intervals <- gintervals(1, 200, 220)

    # Create virtual tracks for different strand settings
    expect_warning(gvtrack.create("count_both", NULL, "kmer.count", list(kmer = "ACGT", strand = 0)))
    gvtrack.create("count_fwd", NULL, "kmer.count", list(kmer = "ACGT", strand = 1))
    gvtrack.create("count_rev", NULL, "kmer.count", list(kmer = "ACGT", strand = -1))

    # For a non-palindromic sequence like "AAGT", the reverse complement is "ACTT"
    gvtrack.create("count_aagt_both", NULL, "kmer.count", list(kmer = "AAGT", strand = 0))
    gvtrack.create("count_aagt_fwd", NULL, "kmer.count", list(kmer = "AAGT", strand = 1))
    gvtrack.create("count_aagt_rev", NULL, "kmer.count", list(kmer = "AAGT", strand = -1))

    # Extract scores
    scores <- gextract(
        c(
            "count_both", "count_fwd", "count_rev",
            "count_aagt_both", "count_aagt_fwd", "count_aagt_rev"
        ),
        test_intervals,
        iterator = test_intervals
    )

    # For palindromic sequence, forward-only and reverse-only should be the same
    # But both should be approximately double either one
    expect_true(scores$count_fwd == scores$count_rev)
    expect_true(scores$count_both == 2 * scores$count_fwd)

    # For non-palindromic, forward and reverse could differ based on the sequence content
    # But the total should equal the sum of forward and reverse
    expect_equal(scores$count_aagt_both, scores$count_aagt_fwd + scores$count_aagt_rev)

    # Test with fraction too
    expect_warning(gvtrack.create("frac_both", NULL, "kmer.frac", list(kmer = "ACGT", strand = 0)))
    gvtrack.create("frac_fwd", NULL, "kmer.frac", list(kmer = "ACGT", strand = 1))
    gvtrack.create("frac_rev", NULL, "kmer.frac", list(kmer = "ACGT", strand = -1))

    frac_scores <- gextract(c("frac_both", "frac_fwd", "frac_rev"),
        test_intervals,
        iterator = test_intervals
    )

    # For palindromic sequence, forward-only and reverse-only fractions should be the same
    expect_equal(frac_scores$frac_fwd, frac_scores$frac_rev, tolerance = 1e-5)

    # Both strands fraction might be up to 2x either strand fraction, but could be less
    # if there's overlap in positions (can't exceed 1.0)
    expect_true(frac_scores$frac_both >= frac_scores$frac_fwd)
    expect_true(frac_scores$frac_both <= min(1.0, 2 * frac_scores$frac_fwd))
})

test_that("kmer functions properly handle strand parameter for forward strand", {
    remove_all_vtracks()

    # Create test interval with a known sequence
    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval)) # Example: CCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCC

    # Create virtual tracks with explicit forward strand
    gvtrack.create("count_ta_fwd", NULL, "kmer.count", list(kmer = "TA", strand = 1))
    gvtrack.create("frac_ta_fwd", NULL, "kmer.frac", list(kmer = "TA", strand = 1))


    # Extract scores
    scores <- gextract(c("count_ta_fwd", "frac_ta_fwd"), test_interval, iterator = test_interval)

    # Manually count occurrences of "TA" in forward orientation
    ta_count <- 0
    for (i in 1:(nchar(seq) - 1)) {
        if (substr(seq, i, i + 1) == "TA") {
            ta_count <- ta_count + 1
        }
    }

    # Test that only forward strand matches are counted
    expect_true(scores$count_ta_fwd > 0) # Should find some occurrences

    # Sum of counts in small chunks should equal total count
    expect_warning(gvtrack.create("count_ta_both", NULL, "kmer.count", list(kmer = "TA", strand = 0)))
    scores_total <- gextract("count_ta_both", test_interval, iterator = test_interval)
    scores_chunks <- gextract("count_ta_both", test_interval, iterator = 10)
    expect_equal(scores_total$count_ta_both, sum(scores_chunks$count_ta_both))
})
test_that("kmer functions properly handle strand parameter for reverse strand", {
    remove_all_vtracks()

    # Create test interval with a known sequence
    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval))
    revcomp_seq <- grevcomp(seq)

    # We need a non-palindromic kmer for this test
    # "TA" is palindromic (its revcomp is also "TA")
    # Let's use "AG" instead, whose revcomp is "CT"

    # Create virtual tracks with explicit reverse strand
    gvtrack.create("count_ag_rev", NULL, "kmer.count", list(kmer = "AG", strand = -1, extend = FALSE))
    gvtrack.create("frac_ag_rev", NULL, "kmer.frac", list(kmer = "AG", strand = -1))

    # Extract scores
    scores <- gextract(c("count_ag_rev", "frac_ag_rev"), test_interval, iterator = test_interval)

    # The reverse strand search for "AG" is equivalent to finding "CT" in the forward strand
    gvtrack.create("count_ct_fwd", NULL, "kmer.count", list(kmer = "CT", strand = 1, extend = FALSE))
    scores_ct_fwd <- gextract("count_ct_fwd", test_interval, iterator = test_interval)

    # Count of "AG" on reverse strand should equal count of its reverse complement "CT" on forward strand
    expect_equal(scores$count_ag_rev, scores_ct_fwd$count_ct_fwd)

    # Manually verify: counts of "AG" on reverse strand should match counts of "CT" in forward sequence
    ct_count <- 0
    for (i in 1:(nchar(seq) - 1)) {
        if (substr(seq, i, i + 1) == "CT") {
            ct_count <- ct_count + 1
        }
    }

    ag_rev_count <- 0
    for (i in 1:(nchar(revcomp_seq) - 1)) {
        if (substr(revcomp_seq, i, i + 1) == "AG") {
            ag_rev_count <- ag_rev_count + 1
        }
    }

    expect_equal(ct_count, ag_rev_count)
    expect_equal(scores$count_ag_rev, ct_count)
})

test_that("kmer functions on both strands equal the sum of forward and reverse strands", {
    remove_all_vtracks()

    # Create test interval
    test_interval <- gintervals(1, 200, 240)

    # Create non-palindromic kmer tracks for all strand options
    kmer <- "AGTC" # Non-palindromic kmer
    gvtrack.create("count_both", NULL, "kmer.count", list(kmer = kmer, strand = 0))
    gvtrack.create("count_fwd", NULL, "kmer.count", list(kmer = kmer, strand = 1))
    gvtrack.create("count_rev", NULL, "kmer.count", list(kmer = kmer, strand = -1))

    gvtrack.create("frac_both", NULL, "kmer.frac", list(kmer = kmer, strand = 0))
    gvtrack.create("frac_fwd", NULL, "kmer.frac", list(kmer = kmer, strand = 1))
    gvtrack.create("frac_rev", NULL, "kmer.frac", list(kmer = kmer, strand = -1))

    # Extract scores
    scores <- gextract(
        c("count_both", "count_fwd", "count_rev", "frac_both", "frac_fwd", "frac_rev"),
        test_interval,
        iterator = test_interval
    )

    # Count for both strands should equal sum of forward and reverse counts
    expect_equal(scores$count_both, scores$count_fwd + scores$count_rev)

    # Fraction for both strands should reflect the combined possible positions
    # For a kmer of length k in a sequence of length n:
    # - Forward strand has (n-k+1) possible positions
    # - Reverse strand has (n-k+1) possible positions
    # - Total positions considering both strands is 2*(n-k+1)

    # The fraction when considering both strands can be calculated as:
    # (forward_count + reverse_count) / (2 * possible_positions)
    # This should equal the reported fraction for both strands

    seq_length <- test_interval$end - test_interval$start
    kmer_length <- nchar(kmer)
    possible_positions <- seq_length - kmer_length + 1

    expected_both_fraction <- (scores$count_fwd + scores$count_rev) / (2 * possible_positions)
    # Due to possible implementation differences, use a tolerance
    expect_equal(scores$frac_both, expected_both_fraction, tolerance = 0.01)
})

test_that("kmer functions handle palindromic sequences correctly with strand parameter", {
    remove_all_vtracks()

    # Create test interval
    test_interval <- gintervals(1, 200, 240)

    # Palindromic kmers (reverse complement equals original)
    palindromic_kmers <- c("AT", "TA", "GC", "CG", "ATAT", "CGCG", "ACGT")

    # Test with "ACGT" which is its own reverse complement
    palindrome <- "ACGT"

    # Test with warning for palindromic kmer with strand=0
    expect_warning(
        gvtrack.create("count_pal_both", NULL, "kmer.count", list(kmer = palindrome, strand = 0))
    )

    # Create tracks for forward and reverse separately (should be identical for palindromes)
    gvtrack.create("count_pal_fwd", NULL, "kmer.count", list(kmer = palindrome, strand = 1))
    gvtrack.create("count_pal_rev", NULL, "kmer.count", list(kmer = palindrome, strand = -1))

    # Extract scores
    scores <- gextract(c("count_pal_fwd", "count_pal_rev"), test_interval, iterator = test_interval)

    # For palindromic sequences, forward and reverse counts should be identical
    expect_equal(scores$count_pal_fwd, scores$count_pal_rev)

    # Verify with multiple palindromic kmers of different lengths
    for (pal in palindromic_kmers[3:5]) { # Test a few different ones
        gvtrack.rm("count_pal_fwd")
        gvtrack.rm("count_pal_rev")

        gvtrack.create("count_pal_fwd", NULL, "kmer.count", list(kmer = pal, strand = 1))
        gvtrack.create("count_pal_rev", NULL, "kmer.count", list(kmer = pal, strand = -1))

        pal_scores <- gextract(c("count_pal_fwd", "count_pal_rev"), test_interval, iterator = test_interval)

        # Should be equal for palindromes
        expect_equal(pal_scores$count_pal_fwd, pal_scores$count_pal_rev)
    }
})

test_that("kmer functions correctly handle strand with overlapping intervals", {
    remove_all_vtracks()

    # Create test intervals that overlap
    interval1 <- gintervals(1, 200, 220)
    interval2 <- gintervals(1, 210, 230)
    combined <- gintervals(1, 200, 230)

    # Non-palindromic kmer to test with
    kmer <- "AGTC"

    # Create tracks for each strand parameter
    gvtrack.create("count_fwd", NULL, "kmer.count", list(kmer = kmer, strand = 1))
    gvtrack.create("count_rev", NULL, "kmer.count", list(kmer = kmer, strand = -1))
    gvtrack.create("count_both", NULL, "kmer.count", list(kmer = kmer, strand = 0))

    # Extract scores for each interval
    scores1 <- gextract("count_fwd", interval1, iterator = interval1)
    scores2 <- gextract("count_fwd", interval2, iterator = interval2)
    scores_combined <- gextract("count_fwd", combined, iterator = combined)

    # Test reverse strand
    rev_scores1 <- gextract("count_rev", interval1, iterator = interval1)
    rev_scores2 <- gextract("count_rev", interval2, iterator = interval2)
    rev_scores_combined <- gextract("count_rev", combined, iterator = combined)

    # Test both strands
    both_scores1 <- gextract("count_both", interval1, iterator = interval1)
    both_scores2 <- gextract("count_both", interval2, iterator = interval2)
    both_scores_combined <- gextract("count_both", combined, iterator = combined)

    # Due to extension, exact equality checks might not work, but we can check relationships
    # The combined interval should have at least as many kmers as each individual interval
    expect_true(scores_combined$count_fwd >= scores1$count_fwd)
    expect_true(scores_combined$count_fwd >= scores2$count_fwd)

    expect_true(rev_scores_combined$count_rev >= rev_scores1$count_rev)
    expect_true(rev_scores_combined$count_rev >= rev_scores2$count_rev)

    expect_true(both_scores_combined$count_both >= both_scores1$count_both)
    expect_true(both_scores_combined$count_both >= both_scores2$count_both)

    # The combined counts on both strands should equal sum of forward and reverse
    all_scores_combined <- gextract(c("count_fwd", "count_rev", "count_both"), combined, iterator = combined)
    expect_equal(all_scores_combined$count_both, all_scores_combined$count_fwd + all_scores_combined$count_rev)
})

test_that("kmer functions handle strand parameter with different kmer lengths", {
    remove_all_vtracks()

    # Create test interval
    test_interval <- gintervals(1, 200, 240)

    # Test with different kmer lengths
    kmer_lengths <- c(2, 4, 6, 8)

    for (len in kmer_lengths) {
        # Create a non-palindromic kmer of specified length
        # Instead of random generation, use specific kmers that we know are not palindromic
        if (len == 2) {
            kmer <- "AG"
        } # AG -> CT
        else if (len == 4) {
            kmer <- "AGCT"
        } # AGCT -> AGCT (palindromic, skip test)
        else if (len == 6) {
            kmer <- "AGCTAG"
        } # AGCTAG -> CTAGCT
        else if (len == 8) {
            kmer <- "AGCTAGCT"
        } # AGCTAGCT -> AGCTAGCT (palindromic, skip test)
        else {
            kmer <- paste0(rep("AG", len / 2), collapse = "")
        } # Non-palindromic for even lengths

        # Skip if palindromic
        if (kmer == grevcomp(kmer)) {
            next
        }

        gvtrack.create("count_fwd", NULL, "kmer.count", list(kmer = kmer, strand = 1, extend = FALSE))
        gvtrack.create("count_rev", NULL, "kmer.count", list(kmer = kmer, strand = -1, extend = FALSE))
        gvtrack.create("count_both", NULL, "kmer.count", list(kmer = kmer, strand = 0, extend = FALSE))

        # Test the relationship between forward, reverse, and both strands
        scores <- gextract(c("count_fwd", "count_rev", "count_both"), test_interval, iterator = test_interval)

        # Both strands should be the sum of forward and reverse
        expect_equal(scores$count_both, scores$count_fwd + scores$count_rev)

        # Create a track for the reverse complement kmer on forward strand
        rc_kmer <- grevcomp(kmer)
        gvtrack.create("count_rc_fwd", NULL, "kmer.count", list(kmer = rc_kmer, strand = 1, extend = FALSE))

        # Reverse complement on forward strand should equal original kmer on reverse strand
        rc_scores <- gextract("count_rc_fwd", test_interval, iterator = test_interval)
        expect_equal(scores$count_rev, rc_scores$count_rc_fwd)

        # Clean up for next iteration
        gvtrack.rm("count_fwd")
        gvtrack.rm("count_rev")
        gvtrack.rm("count_both")
        gvtrack.rm("count_rc_fwd")
    }
})

test_that("kmer functions correctly handle negative strand with extension at boundaries", {
    remove_all_vtracks()

    # Create test intervals at chromosome boundaries
    start_interval <- gintervals(1, 0, 20)
    end_interval <- gintervals(1, 9980, 10000)

    # Use a kmer that's likely to appear at chromosome boundaries
    kmer <- "ACGT"

    # Create tracks with different extension and strand settings
    gvtrack.create("count_fwd_ext", NULL, "kmer.count", list(kmer = kmer, strand = 1, extend = TRUE))
    gvtrack.create("count_fwd_noext", NULL, "kmer.count", list(kmer = kmer, strand = 1, extend = FALSE))
    gvtrack.create("count_rev_ext", NULL, "kmer.count", list(kmer = kmer, strand = -1, extend = TRUE))
    gvtrack.create("count_rev_noext", NULL, "kmer.count", list(kmer = kmer, strand = -1, extend = FALSE))

    # Extract scores at boundaries
    start_scores <- gextract(c("count_fwd_ext", "count_fwd_noext", "count_rev_ext", "count_rev_noext"),
        start_interval,
        iterator = start_interval
    )
    end_scores <- gextract(c("count_fwd_ext", "count_fwd_noext", "count_rev_ext", "count_rev_noext"),
        end_interval,
        iterator = end_interval
    )

    # Extension should have valid results (not NA) even at chromosome boundaries
    expect_true(!is.na(start_scores$count_fwd_ext))
    expect_true(!is.na(start_scores$count_rev_ext))
    expect_true(!is.na(end_scores$count_fwd_ext))
    expect_true(!is.na(end_scores$count_rev_ext))

    # Non-extension should have valid results
    expect_true(!is.na(start_scores$count_fwd_noext))
    expect_true(!is.na(start_scores$count_rev_noext))
    expect_true(!is.na(end_scores$count_fwd_noext))
    expect_true(!is.na(end_scores$count_rev_noext))

    # At the start boundary, reverse strand with extension should still work
    # At the end boundary, forward strand with extension should still work
    # The extension should be properly limited by chromosome boundaries
})

test_that("kmer functions handle strand with mixed-case sequences", {
    remove_all_vtracks()

    # Create test interval
    test_interval <- gintervals(1, 200, 240)

    # Use non-palindromic kmers for clearer testing
    # "AG" -> "CT"

    # Create virtual tracks with different case patterns
    gvtrack.create("count_upper_fwd", NULL, "kmer.count", list(kmer = "AG", strand = 1, extend = FALSE))
    gvtrack.create("count_lower_fwd", NULL, "kmer.count", list(kmer = "ag", strand = 1, extend = FALSE))
    gvtrack.create("count_mixed_fwd", NULL, "kmer.count", list(kmer = "Ag", strand = 1, extend = FALSE))

    gvtrack.create("count_upper_rev", NULL, "kmer.count", list(kmer = "AG", strand = -1, extend = FALSE))
    gvtrack.create("count_lower_rev", NULL, "kmer.count", list(kmer = "ag", strand = -1, extend = FALSE))
    gvtrack.create("count_mixed_rev", NULL, "kmer.count", list(kmer = "Ag", strand = -1, extend = FALSE))

    # Extract scores
    scores_fwd <- gextract(c("count_upper_fwd", "count_lower_fwd", "count_mixed_fwd"),
        test_interval,
        iterator = test_interval
    )
    scores_rev <- gextract(c("count_upper_rev", "count_lower_rev", "count_mixed_rev"),
        test_interval,
        iterator = test_interval
    )

    # Case should not matter for forward strand
    expect_equal(scores_fwd$count_upper_fwd, scores_fwd$count_lower_fwd)
    expect_equal(scores_fwd$count_lower_fwd, scores_fwd$count_mixed_fwd)

    # Case should not matter for reverse strand either
    expect_equal(scores_rev$count_upper_rev, scores_rev$count_lower_rev)
    expect_equal(scores_rev$count_lower_rev, scores_rev$count_mixed_rev)

    # Verify that reverse complement works correctly regardless of case
    # "AG" reverse complement is "CT"
    gvtrack.create("count_ct_upper_fwd", NULL, "kmer.count", list(kmer = "CT", strand = 1, extend = FALSE))
    gvtrack.create("count_ct_lower_fwd", NULL, "kmer.count", list(kmer = "ct", strand = 1, extend = FALSE))
    scores_ct <- gextract(c("count_ct_upper_fwd", "count_ct_lower_fwd"),
        test_interval,
        iterator = test_interval
    )

    # CT forward should match AG reverse
    expect_equal(scores_ct$count_ct_upper_fwd, scores_rev$count_upper_rev)
    expect_equal(scores_ct$count_ct_lower_fwd, scores_rev$count_lower_rev)
})

test_that("sum of fraction of T, C, G, A equals 1", {
    remove_all_vtracks()

    test_interval <- gintervals(1, 200, 240)
    gvtrack.create("frac_t", NULL, "kmer.frac", kmer = "T", strand = 1)
    gvtrack.create("frac_c", NULL, "kmer.frac", kmer = "C", strand = 1)
    gvtrack.create("frac_g", NULL, "kmer.frac", kmer = "G", strand = 1)
    gvtrack.create("frac_a", NULL, "kmer.frac", kmer = "A", strand = 1)

    scores <- gextract("frac_a + frac_c + frac_g + frac_t", test_interval, iterator = test_interval, colnames = "s")
    expect_equal(sum(scores$s), 1, tolerance = 1e-5)

    gvtrack.create("frac_t", NULL, "kmer.frac", kmer = "T", strand = -1)
    gvtrack.create("frac_c", NULL, "kmer.frac", kmer = "C", strand = -1)
    gvtrack.create("frac_g", NULL, "kmer.frac", kmer = "G", strand = -1)
    gvtrack.create("frac_a", NULL, "kmer.frac", kmer = "A", strand = -1)

    scores <- gextract("frac_a + frac_c + frac_g + frac_t", test_interval, iterator = test_interval, colnames = "s")
    expect_equal(sum(scores$s), 1, tolerance = 1e-5)

    gvtrack.create("frac_t", NULL, "kmer.frac", kmer = "T", strand = 0)
    gvtrack.create("frac_c", NULL, "kmer.frac", kmer = "C", strand = 0)
    gvtrack.create("frac_g", NULL, "kmer.frac", kmer = "G", strand = 0)
    gvtrack.create("frac_a", NULL, "kmer.frac", kmer = "A", strand = 0)

    scores <- gextract("frac_a + frac_c + frac_g + frac_t", test_interval, iterator = test_interval, colnames = "s")
    expect_equal(sum(scores$s), 1, tolerance = 1e-5)
})

test_that("kmer.frac correctly doubles possible positions when strand=0", {
    remove_all_vtracks()

    # Create test interval with known content
    test_interval <- gintervals(1, 200, 210)
    seq <- toupper(gseq.extract(test_interval)) # Should be "CCCTAACCCT"

    # Create non-palindromic kmer
    kmer <- "AG"
    rc_kmer <- "CT" # Reverse complement of AG

    # Define the same kmer with different strand settings
    gvtrack.create("frac_fwd", NULL, "kmer.frac", list(kmer = kmer, strand = 1, extend = FALSE))
    gvtrack.create("frac_rev", NULL, "kmer.frac", list(kmer = kmer, strand = -1, extend = FALSE))
    gvtrack.create("frac_both", NULL, "kmer.frac", list(kmer = kmer, strand = 0, extend = FALSE))

    # Also create complementary forward track to verify rev behavior
    gvtrack.create("frac_rc_fwd", NULL, "kmer.frac", list(kmer = rc_kmer, strand = 1, extend = FALSE))

    # Extract scores
    scores <- gextract(
        c("frac_fwd", "frac_rev", "frac_both", "frac_rc_fwd"),
        test_interval,
        iterator = test_interval
    )

    # Count AG and CT occurrences manually
    ag_count <- 0
    ct_count <- 0
    for (i in 1:(nchar(seq) - 1)) {
        if (substr(seq, i, i + 1) == kmer) {
            ag_count <- ag_count + 1
        }
        if (substr(seq, i, i + 1) == rc_kmer) {
            ct_count <- ct_count + 1
        }
    }

    # Calculate positions where a kmer could start
    # For a sequence of length n and kmer of length k, there are (n-k+1) possible positions
    possible_positions <- nchar(seq) - nchar(kmer) + 1

    # Expected fractions
    expected_frac_fwd <- ag_count / possible_positions
    expected_frac_rev <- ct_count / possible_positions

    # For both strands, the denominator is doubled
    expected_frac_both <- (ag_count + ct_count) / (possible_positions * 2)

    # Test that fractions match expected values
    expect_equal(scores$frac_fwd, expected_frac_fwd)
    expect_equal(scores$frac_rev, expected_frac_rev)
    expect_equal(scores$frac_rc_fwd, expected_frac_rev) # Should match rev strand
    expect_equal(scores$frac_both, expected_frac_both)

    # Verify that both strands calculation is consistent
    expect_equal(scores$frac_both * 2 * possible_positions, (ag_count + ct_count))

    # Should equal the sum of counts divided by double the positions
    expect_equal(scores$frac_both, (ag_count + ct_count) / (2 * possible_positions))
})

test_that("kmer.frac handles edge cases with strand=0 correctly", {
    remove_all_vtracks()

    # Test with very short sequences
    short_interval <- gintervals(1, 200, 202) # 2bp sequence

    # Create kmer tracks for testing
    gvtrack.create("frac_1bp_both", NULL, "kmer.frac", list(kmer = "A", strand = 0))
    gvtrack.create("frac_2bp_both", NULL, "kmer.frac", list(kmer = "AG", strand = 0))

    # Extract scores
    scores_short <- gextract(
        c("frac_1bp_both", "frac_2bp_both"),
        short_interval,
        iterator = short_interval
    )

    # For 1bp kmer in 2bp sequence, should have 4 possible positions (2 positions * 2 strands)
    # For 2bp kmer in 2bp sequence, should have 2 possible positions (1 position * 2 strands)
    # We can't be absolutely precise about the values without knowing the exact sequence,
    # but we can verify they're in range
    expect_true(scores_short$frac_1bp_both >= 0 && scores_short$frac_1bp_both <= 1)
    expect_true(scores_short$frac_2bp_both >= 0 && scores_short$frac_2bp_both <= 1)

    # Test with palindromic kmers (where forward = reverse complement)
    palindromic_interval <- gintervals(1, 200, 220)

    # "AT" is its own reverse complement
    gvtrack.create("frac_at_fwd", NULL, "kmer.frac", list(kmer = "AT", strand = 1))
    expect_warning(gvtrack.create("frac_at_both", NULL, "kmer.frac", list(kmer = "AT", strand = 0)))

    # Extract scores
    scores_palindrome <- gextract(
        c("frac_at_fwd", "frac_at_both"),
        palindromic_interval,
        iterator = palindromic_interval
    )

    # For palindromic kmers, the count should be the same for either strand
    # But the denominator is still doubled for strand=0, so the fraction should be half
    expect_equal(scores_palindrome$frac_at_both, scores_palindrome$frac_at_fwd / 2)

    # Test with multiple sequence lengths and kmer lengths to ensure denominator calculation is correct
    for (seq_len in c(10, 20, 50)) {
        for (k_len in c(1, 2, 4)) {
            if (k_len <= seq_len) {
                test_interval <- gintervals(1, 200, 200 + seq_len)
                kmer <- paste0(rep("A", k_len), collapse = "")

                vtrack_name <- paste0("frac_", k_len, "mer_", seq_len, "_both")
                gvtrack.create(vtrack_name, NULL, "kmer.frac", list(kmer = kmer, strand = 0))

                scores <- gextract(vtrack_name, test_interval, iterator = test_interval)

                # Calculate expected number of positions
                possible_positions <- seq_len - k_len + 1

                # The total possible positions should be doubled for strand=0
                # We can verify this by forcing a sequence with all As
                # which would yield a count of possible_positions and fraction of 0.5
                expect_true(scores[[vtrack_name]] <= 1 && scores[[vtrack_name]] >= 0)
            }
        }
    }
})

test_that("kmer.frac with strand=0 handles boundary extension correctly", {
    remove_all_vtracks()

    # Test with extension enabled (default) and disabled
    test_interval <- gintervals(1, 200, 210)
    kmer <- "AAG" # 3-letter kmer

    gvtrack.create(
        "frac_both_extend", NULL, "kmer.frac",
        list(kmer = kmer, strand = 0, extend = TRUE)
    )
    gvtrack.create(
        "frac_both_no_extend", NULL, "kmer.frac",
        list(kmer = kmer, strand = 0, extend = FALSE)
    )

    # Extract scores
    scores <- gextract(
        c("frac_both_extend", "frac_both_no_extend"),
        test_interval,
        iterator = test_interval
    )

    # With extension, kmers that start at the last two positions should be counted
    # Without extension, they should not
    expect_true(scores$frac_both_extend >= scores$frac_both_no_extend)

    # Create very small interval where extension matters significantly
    small_interval <- gintervals(1, 200, 203) # 3bp

    scores_small <- gextract(
        c("frac_both_extend", "frac_both_no_extend"),
        small_interval,
        iterator = small_interval
    )

    # For a 3-letter kmer in a 3bp interval:
    # - Without extension: 1 possible position on each strand (2 total)
    # - With extension: 3 possible positions on each strand (6 total)
    # Again, we can't predict exact values without knowing the sequence
    expect_true(scores_small$frac_both_extend >= scores_small$frac_both_no_extend)

    # With non-extending 3bp kmer in a 3bp interval, the denominator should be 2 (1 per strand)
    expect_true(scores_small$frac_both_no_extend >= 0 && scores_small$frac_both_no_extend <= 1)
})

test_that("sum of base fractions equals 1 when strand=0", {
    remove_all_vtracks()

    test_interval <- gintervals(1, 200, 240)

    # Create tracks for all bases with strand=0
    gvtrack.create("frac_a", NULL, "kmer.frac", list(kmer = "A", strand = 0))
    gvtrack.create("frac_c", NULL, "kmer.frac", list(kmer = "C", strand = 0))
    gvtrack.create("frac_g", NULL, "kmer.frac", list(kmer = "G", strand = 0))
    gvtrack.create("frac_t", NULL, "kmer.frac", list(kmer = "T", strand = 0))

    # Extract sum of all fractions
    scores <- gextract("frac_a + frac_c + frac_g + frac_t",
        test_interval,
        iterator = test_interval,
        colnames = "sum"
    )

    # Sum should be 1 for single bases even with strand=0
    # because each position is counted only once per strand
    expect_equal(scores$sum, 1, tolerance = 1e-5)

    # For comparison, create tracks with strand=1
    gvtrack.create("frac_a_fwd", NULL, "kmer.frac", list(kmer = "A", strand = 1))
    gvtrack.create("frac_c_fwd", NULL, "kmer.frac", list(kmer = "C", strand = 1))
    gvtrack.create("frac_g_fwd", NULL, "kmer.frac", list(kmer = "G", strand = 1))
    gvtrack.create("frac_t_fwd", NULL, "kmer.frac", list(kmer = "T", strand = 1))

    # Extract and sum forward strand fractions
    scores_fwd <- gextract("frac_a_fwd + frac_c_fwd + frac_g_fwd + frac_t_fwd",
        test_interval,
        iterator = test_interval,
        colnames = "sum_fwd"
    )

    # Sum should be 1 for forward strand too
    expect_equal(scores_fwd$sum_fwd, 1, tolerance = 1e-5)

    scores_compare <- gextract(
        c("frac_a", "frac_a_fwd", "frac_c", "frac_c_fwd", "frac_g", "frac_g_fwd", "frac_t", "frac_t_fwd"),
        test_interval,
        iterator = test_interval
    )

    seq <- toupper(gseq.extract(test_interval))
    A_fwd <- stringr::str_count(seq, "A")
    C_fwd <- stringr::str_count(seq, "C")
    G_fwd <- stringr::str_count(seq, "G")
    T_fwd <- stringr::str_count(seq, "T")
    A_rev <- stringr::str_count(grevcomp(seq), "A")
    C_rev <- stringr::str_count(grevcomp(seq), "C")
    G_rev <- stringr::str_count(grevcomp(seq), "G")
    T_rev <- stringr::str_count(grevcomp(seq), "T")
    seq_l <- nchar(seq)

    # Check that fractions are consistent with expected values
    expect_equal(scores_compare$frac_a_fwd, A_fwd / seq_l, tolerance = 1e-5)
    expect_equal(scores_compare$frac_a, (A_fwd + A_rev) / (2 * seq_l), tolerance = 1e-5)
    expect_equal(scores_compare$frac_c_fwd, C_fwd / seq_l, tolerance = 1e-5)
    expect_equal(scores_compare$frac_c, (C_fwd + C_rev) / (2 * seq_l), tolerance = 1e-5)
    expect_equal(scores_compare$frac_g_fwd, G_fwd / seq_l, tolerance = 1e-5)
    expect_equal(scores_compare$frac_g, (G_fwd + G_rev) / (2 * seq_l), tolerance = 1e-5)
    expect_equal(scores_compare$frac_t_fwd, T_fwd / seq_l, tolerance = 1e-5)
    expect_equal(scores_compare$frac_t, (T_fwd + T_rev) / (2 * seq_l), tolerance = 1e-5)
})

test_that("kmer honors iterator shifts across multiple magnitudes (extend=TRUE, fwd strand)", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 2000, 2040)
    k <- "TA" # palindromic is fine for forward-only tests

    shift_values <- c(0, 1, 5, 10, 25, 50, 100, 250, 1000)
    count_names <- sprintf("count_shift_%d", shift_values)
    frac_names <- sprintf("frac_shift_%d", shift_values)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        gvtrack.create(count_names[i], NULL, "kmer.count", kmer = k, strand = 1)
        gvtrack.iterator(count_names[i], sshift = -s, eshift = s)
        gvtrack.create(frac_names[i], NULL, "kmer.frac", kmer = k, strand = 1)
        gvtrack.iterator(frac_names[i], sshift = -s, eshift = s)
    }

    scores <- gextract(c(count_names, frac_names), test_intervals, iterator = test_intervals)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        # Original (per-vtrack) iterator after shift
        orig <- test_intervals
        orig$start <- pmax(0, orig$start - s)
        orig$end <- orig$end + s
        # Fetch window extends at end by k-1 for forward strand when extend=TRUE (default)
        fetch <- orig
        fetch$end <- fetch$end + (nchar(k) - 1)
        seq_ext <- toupper(gseq.extract(fetch))

        res <- count_kmers_manually(seq_ext, k, orig, fetch)

        # Counts should match exactly; fractions within small tolerance
        expect_equal(scores[[count_names[i]]], res$count)
        tol <- if (s >= 250) 1e-5 else 1e-6
        expect_equal(scores[[frac_names[i]]], res$fraction, tolerance = tol)
    }
})

test_that("kmer honors iterator shifts without extension (extend=FALSE, fwd strand)", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 2000, 2040)
    k <- "TA"

    shift_values <- c(0, 5, 10, 50, 200)
    count_names <- sprintf("count_noext_shift_%d", shift_values)
    frac_names <- sprintf("frac_noext_shift_%d", shift_values)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        gvtrack.create(count_names[i], NULL, "kmer.count", kmer = k, strand = 1, extend = FALSE)
        gvtrack.iterator(count_names[i], sshift = -s, eshift = s)
        gvtrack.create(frac_names[i], NULL, "kmer.frac", kmer = k, strand = 1, extend = FALSE)
        gvtrack.iterator(frac_names[i], sshift = -s, eshift = s)
    }

    scores <- gextract(c(count_names, frac_names), test_intervals, iterator = test_intervals)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        orig <- test_intervals
        orig$start <- pmax(0, orig$start - s)
        orig$end <- orig$end + s
        fetch <- orig # no extension at fetch
        seq_noext <- toupper(gseq.extract(fetch))

        res <- count_kmers_manually(seq_noext, k, orig, fetch)

        expect_equal(scores[[count_names[i]]], res$count)
        tol <- if (s >= 50) 1e-5 else 1e-6
        expect_equal(scores[[frac_names[i]]], res$fraction, tolerance = tol)
    }
})

test_that("kmer honors iterator shifts for reverse strand (extend=TRUE)", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 2000, 2040)
    k <- "AG" # non-palindromic; reverse complement is CT
    rc_k <- grevcomp(k)

    shift_values <- c(0, 10, 200)
    count_names <- sprintf("count_rev_shift_%d", shift_values)
    frac_names <- sprintf("frac_rev_shift_%d", shift_values)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        gvtrack.create(count_names[i], NULL, "kmer.count", kmer = k, strand = -1)
        gvtrack.iterator(count_names[i], sshift = -s, eshift = s)
        gvtrack.create(frac_names[i], NULL, "kmer.frac", kmer = k, strand = -1)
        gvtrack.iterator(frac_names[i], sshift = -s, eshift = s)
    }

    scores <- gextract(c(count_names, frac_names), test_intervals, iterator = test_intervals)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        # Original (per-vtrack) iterator after shift
        orig <- test_intervals
        orig$start <- pmax(0, orig$start - s)
        orig$end <- orig$end + s
        # For reverse strand with extend=TRUE, fetch extends at the start by k-1
        fetch <- orig
        fetch$start <- pmax(0, fetch$start - (nchar(k) - 1))
        seq_ext <- toupper(gseq.extract(fetch))

        # Count occurrences of reverse complement in forward sequence window
        res <- count_kmers_manually(seq_ext, rc_k, orig, fetch)

        expect_equal(scores[[count_names[i]]], res$count)
        expect_equal(scores[[frac_names[i]]], res$fraction, tolerance = 1e-6)
    }
})

test_that("kmer honors iterator shifts for both strands (extend=TRUE)", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 2000, 2040)
    k <- "AG" # non-palindromic
    rc_k <- grevcomp(k)

    shift_values <- c(0, 10, 200)
    count_names <- sprintf("count_both_shift_%d", shift_values)
    frac_names <- sprintf("frac_both_shift_%d", shift_values)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        gvtrack.create(count_names[i], NULL, "kmer.count", list(kmer = k, strand = 0))
        gvtrack.iterator(count_names[i], sshift = -s, eshift = s)
        gvtrack.create(frac_names[i], NULL, "kmer.frac", list(kmer = k, strand = 0))
        gvtrack.iterator(frac_names[i], sshift = -s, eshift = s)
    }

    scores <- gextract(c(count_names, frac_names), test_intervals, iterator = test_intervals)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        orig <- test_intervals
        orig$start <- pmax(0, orig$start - s)
        orig$end <- orig$end + s
        # Forward fetch extends at end; reverse fetch extends at start
        fetch_fwd <- orig
        fetch_fwd$end <- fetch_fwd$end + (nchar(k) - 1)
        seq_fwd <- toupper(gseq.extract(fetch_fwd))
        res_fwd <- count_kmers_manually(seq_fwd, k, orig, fetch_fwd)

        fetch_rev <- orig
        fetch_rev$start <- pmax(0, fetch_rev$start - (nchar(k) - 1))
        seq_rev <- toupper(gseq.extract(fetch_rev))
        res_rev <- count_kmers_manually(seq_rev, rc_k, orig, fetch_rev)

        expected_count <- res_fwd$count + res_rev$count
        possible_positions <- (orig$end - orig$start) - nchar(k) + 1
        possible_positions <- max(0, possible_positions)
        expected_frac <- if (possible_positions > 0) expected_count / (2 * possible_positions) else 0

        expect_equal(scores[[count_names[i]]], expected_count)
        expect_equal(scores[[frac_names[i]]], expected_frac, tolerance = 1e-5)
    }
})
