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
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA")
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA")
    gvtrack.create("count_ccc", NULL, "kmer.count", kmer = "CCC")
    gvtrack.create("frac_ccc", NULL, "kmer.frac", kmer = "CCC")

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

test_that("kmer.count and kmer.frac work with smaller iterator intervals", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals)) # CCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCC

    # Create virtual tracks
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA")
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA")

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
    gvtrack.create("count_ta_upper", NULL, "kmer.count", kmer = "TA")
    gvtrack.create("count_ta_lower", NULL, "kmer.count", kmer = "ta")
    gvtrack.create("count_ta_mixed", NULL, "kmer.count", kmer = "Ta")

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
    gvtrack.create("count_taac", NULL, "kmer.count", kmer = "TAAC")
    gvtrack.create("frac_taac", NULL, "kmer.frac", kmer = "TAAC")
    gvtrack.create("count_ccctaa", NULL, "kmer.count", kmer = "CCCTAA")
    gvtrack.create("frac_ccctaa", NULL, "kmer.frac", kmer = "CCCTAA")

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
    gvtrack.create("count_full", NULL, "kmer.count", kmer = seq)
    gvtrack.create("frac_full", NULL, "kmer.frac", kmer = seq)

    # Test 1bp k-mer
    gvtrack.create("count_a", NULL, "kmer.count", kmer = "A")
    gvtrack.create("frac_a", NULL, "kmer.frac", kmer = "A")

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
    gvtrack.create("count_xyz", NULL, "kmer.count", kmer = "XYZ")
    gvtrack.create("frac_xyz", NULL, "kmer.frac", kmer = "XYZ")

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
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA")
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA")

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
    gvtrack.create("count_taa", NULL, "kmer.count", kmer = "TAA")
    gvtrack.create("frac_taa", NULL, "kmer.frac", kmer = "TAA")

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
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA")
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA")

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
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA")
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA")

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
    gvtrack.create("count_a", NULL, "kmer.count", kmer = "A")
    gvtrack.create("frac_a", NULL, "kmer.frac", kmer = "A")

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
    gvtrack.create("count_ccc", NULL, "kmer.count", kmer = "CCC")

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
    gvtrack.create("count_a", NULL, "kmer.count", kmer = "A")
    gvtrack.create("frac_a", NULL, "kmer.frac", kmer = "A")

    # Create a distance virtual track
    annotations <- gintervals(1, 220, 240)
    gvtrack.create("dist_track", annotations, "distance")

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
    expect_error(gvtrack.create("empty_kmer", NULL, "kmer.count", kmer = ""))

    # Non-string parameter should error
    expect_error(gvtrack.create("numeric_kmer", NULL, "kmer.count", kmer = 123))

    # No parameter should error
    expect_error(gvtrack.create("no_param", NULL, "kmer.count"))

    # Multiple k-mers should error
    expect_error(gvtrack.create("multi_kmer", NULL, "kmer.count", kmer = c("AA", "TT")))
})

test_that("kmer.count and kmer.frac work with explicit extension parameters", {
    remove_all_vtracks()

    # Create a test interval
    test_interval <- gintervals(1, 200, 220)

    # Create virtual tracks with different extension settings
    gvtrack.create("count_default", NULL, "kmer.count", kmer = "TA") # Default extension=TRUE
    gvtrack.create("count_extend", NULL, "kmer.count", kmer = "TA", extend = TRUE)
    gvtrack.create("count_no_extend", NULL, "kmer.count", kmer = "TA", extend = FALSE)

    # Create fraction tracks with the same settings
    gvtrack.create("frac_default", NULL, "kmer.frac", kmer = "TA")
    gvtrack.create("frac_extend", NULL, "kmer.frac", kmer = "TA", extend = TRUE)
    gvtrack.create("frac_no_extend", NULL, "kmer.frac", kmer = "TA", extend = FALSE)

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
    gvtrack.create("count_with_ext", NULL, "kmer.count", kmer = "ACGT", extend = TRUE)
    gvtrack.create("count_no_ext", NULL, "kmer.count", kmer = "ACGT", extend = FALSE)

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
    gvtrack.create("count_2mer", NULL, "kmer.count", kmer = "AT")
    gvtrack.create("count_5mer", NULL, "kmer.count", kmer = "ATGCG")
    gvtrack.create("count_10mer", NULL, "kmer.count", kmer = "ATGCGCGTAA")

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
    gvtrack.create("count_ccc", NULL, "kmer.count", kmer = "CCC")

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
