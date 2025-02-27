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

    # Manual count of "TA" occurrences - manual counting to handle potential overlapping matches
    ta_count <- 0
    for (i in 1:(nchar(seq) - nchar("TA") + 1)) {
        if (substr(seq, i, i + nchar("TA") - 1) == "TA") {
            ta_count <- ta_count + 1
        }
    }
    ta_possible_positions <- nchar(seq) - nchar("TA") + 1

    # Manual count of "CCC" occurrences - need to count overlapping occurrences
    ccc_count <- 0
    for (i in 1:(nchar(seq) - nchar("CCC") + 1)) {
        if (substr(seq, i, i + nchar("CCC") - 1) == "CCC") {
            ccc_count <- ccc_count + 1
        }
    }
    ccc_possible_positions <- nchar(seq) - nchar("CCC") + 1

    # Test that counts match expected values
    expect_equal(scores$count_ta, ta_count)
    expect_equal(scores$frac_ta, ta_count / ta_possible_positions, tolerance = 1e-5)
    expect_equal(scores$count_ccc, ccc_count)
    expect_equal(scores$frac_ccc, ccc_count / ccc_possible_positions, tolerance = 1e-5)
})

test_that("kmer.count and kmer.frac work with smaller iterator intervals", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals))

    # Create virtual tracks
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA")
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA")

    # Extract scores with 10bp iterator
    scores_10bp <- gextract(c("count_ta", "frac_ta"), test_intervals, iterator = 10)

    # Calculate expected counts and fractions for each 10bp bin
    expected_count <- numeric(4)
    expected_frac <- numeric(4)

    for (i in 1:4) {
        bin_start <- (i - 1) * 10 + 1
        bin_end <- min(i * 10, nchar(seq))
        bin_seq <- substr(seq, bin_start, bin_end)

        # Count occurrences of "TA" in this bin
        bin_count <- stringr::str_count(bin_seq, "TA")
        bin_possible_positions <- nchar(bin_seq) - nchar("TA") + 1

        expected_count[i] <- bin_count
        expected_frac[i] <- if (bin_possible_positions > 0) bin_count / bin_possible_positions else 0
    }

    # Test count values for each bin
    expect_equal(scores_10bp$count_ta, expected_count)
    expect_equal(scores_10bp$frac_ta, expected_frac, tolerance = 1e-5)
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

    # All should give the same results (case-insensitive matching)
    ta_count <- stringr::str_count(seq, "TA")

    expect_equal(scores$count_ta_upper, ta_count)
    expect_equal(scores$count_ta_lower, ta_count)
    expect_equal(scores$count_ta_mixed, ta_count)
})

test_that("kmer.count and kmer.frac handle longer k-mers correctly", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals))

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

    # Manual counts - using explicit loop to correctly count overlapping occurrences
    taac_count <- 0
    for (i in 1:(nchar(seq) - nchar("TAAC") + 1)) {
        if (substr(seq, i, i + nchar("TAAC") - 1) == "TAAC") {
            taac_count <- taac_count + 1
        }
    }
    taac_possible_positions <- nchar(seq) - nchar("TAAC") + 1

    ccctaa_count <- 0
    for (i in 1:(nchar(seq) - nchar("CCCTAA") + 1)) {
        if (substr(seq, i, i + nchar("CCCTAA") - 1) == "CCCTAA") {
            ccctaa_count <- ccctaa_count + 1
        }
    }
    ccctaa_possible_positions <- nchar(seq) - nchar("CCCTAA") + 1

    # Test that counts match expected values
    expect_equal(scores$count_taac, taac_count)
    expect_equal(scores$frac_taac, taac_count / taac_possible_positions, tolerance = 1e-5)
    expect_equal(scores$count_ccctaa, ccctaa_count)
    expect_equal(scores$frac_ccctaa, ccctaa_count / ccctaa_possible_positions, tolerance = 1e-5)
})

test_that("kmer.count and kmer.frac handle edge cases correctly", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)

    # Test k-mer length equal to sequence length
    gvtrack.create("count_full", NULL, "kmer.count", kmer = toupper(gseq.extract(test_intervals)))
    gvtrack.create("frac_full", NULL, "kmer.frac", kmer = toupper(gseq.extract(test_intervals)))

    # Test 1bp k-mer
    gvtrack.create("count_a", NULL, "kmer.count", kmer = "A")
    gvtrack.create("frac_a", NULL, "kmer.frac", kmer = "A")

    # Extract scores
    scores <- gextract(c("count_full", "frac_full", "count_a", "frac_a"),
        test_intervals,
        iterator = test_intervals
    )

    seq <- toupper(gseq.extract(test_intervals))

    # For full-length k-mer
    expect_equal(scores$count_full, 1)
    expect_equal(scores$frac_full, 1)

    # For single base k-mer - count manually for single characters too to maintain consistency
    a_count <- 0
    for (i in 1:nchar(seq)) {
        if (substr(seq, i, i) == "A") {
            a_count <- a_count + 1
        }
    }
    expect_equal(scores$count_a, a_count)
    expect_equal(scores$frac_a, a_count / nchar(seq), tolerance = 1e-5)
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

    # Total counts should sum to the same value
    expect_equal(sum(scores_5bp$count_ta), sum(scores_10bp$count_ta))
    expect_equal(sum(scores_10bp$count_ta), sum(scores_20bp$count_ta))

    # The average of fractions should be similar (but not exactly equal due to different bin sizes)
    ta_count <- sum(scores_5bp$count_ta)
    seq <- toupper(gseq.extract(test_intervals))
    possible_positions <- nchar(seq) - nchar("TA") + 1
    expected_avg_frac <- ta_count / possible_positions

    # Check that average fractions are reasonably close to expected
    expect_true(abs(mean(scores_5bp$frac_ta) - expected_avg_frac) < 0.1)
    expect_true(abs(mean(scores_10bp$frac_ta) - expected_avg_frac) < 0.1)
    expect_true(abs(mean(scores_20bp$frac_ta) - expected_avg_frac) < 0.1)
})

test_that("kmer.count and kmer.frac work with overlapping k-mers", {
    remove_all_vtracks()

    # Create a test interval with repeating sequence
    test_intervals <- gintervals(1, 200, 220)
    seq <- toupper(gseq.extract(test_intervals))

    # Create virtual tracks for overlapping k-mer
    gvtrack.create("count_taa", NULL, "kmer.count", kmer = "TAA")
    gvtrack.create("frac_taa", NULL, "kmer.frac", kmer = "TAA")

    # Extract scores
    scores <- gextract(c("count_taa", "frac_taa"), test_intervals, iterator = test_intervals)

    # Manual count (consider overlapping occurrences)
    taa_count <- 0
    for (i in 1:(nchar(seq) - 2)) {
        if (substr(seq, i, i + 2) == "TAA") {
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

    # The sum of sub-interval counts should equal the full interval count
    # (Only if no k-mer spans the boundary between sub-intervals)
    expect_equal(full_scores$count_ta, sub1_scores$count_ta + sub2_scores$count_ta)
})

test_that("kmer.count and kmer.frac return correct values for small intervals", {
    remove_all_vtracks()

    # Create a test interval with known content
    test_interval <- gintervals(1, 200, 203) # Just 3bp
    seq <- toupper(gseq.extract(test_interval))

    # Create virtual tracks
    gvtrack.create("count_2mer", NULL, "kmer.count", kmer = substr(seq, 1, 2))
    gvtrack.create("frac_2mer", NULL, "kmer.frac", kmer = substr(seq, 1, 2))
    gvtrack.create("count_3mer", NULL, "kmer.count", kmer = seq)
    gvtrack.create("frac_3mer", NULL, "kmer.frac", kmer = seq)
    gvtrack.create("count_4mer", NULL, "kmer.count", kmer = paste0(seq, "A"))
    gvtrack.create("frac_4mer", NULL, "kmer.frac", kmer = paste0(seq, "A"))

    # Extract scores
    scores <- gextract(c("count_2mer", "frac_2mer", "count_3mer", "frac_3mer", "count_4mer", "frac_4mer"),
        test_interval,
        iterator = test_interval
    )

    # 2-mer should be found once
    expect_equal(scores$count_2mer, 1)
    expect_equal(scores$frac_2mer, 1 / (nchar(seq) - 1), tolerance = 1e-5)

    # 3-mer should be found once
    expect_equal(scores$count_3mer, 1)
    expect_equal(scores$frac_3mer, 1)

    # 4-mer can't be found in 3bp
    expect_equal(scores$count_4mer, 0)
    expect_equal(scores$frac_4mer, 0)
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

test_that("kmer.count and kmer.frac handle negatively stranded intervals", {
    remove_all_vtracks()

    # Create normal and negative strand intervals
    pos_interval <- gintervals(1, 200, 220)
    neg_interval <- gintervals(1, 200, 220)
    neg_interval$strand <- -1

    # Get sequences (one forward, one reverse complement)
    seq_pos <- toupper(gseq.extract(pos_interval))
    seq_neg <- toupper(gseq.extract(neg_interval))

    # Create kmer count tracks
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA")
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA")

    # Extract scores for both intervals
    scores_pos <- gextract(c("count_ta", "frac_ta"), pos_interval, iterator = pos_interval)
    scores_neg <- gextract(c("count_ta", "frac_ta"), neg_interval, iterator = neg_interval)

    # Calculate expected values for reverse complement sequence
    ta_count_neg <- 0
    for (i in 1:(nchar(seq_neg) - nchar("TA") + 1)) {
        if (substr(seq_neg, i, i + nchar("TA") - 1) == "TA") {
            ta_count_neg <- ta_count_neg + 1
        }
    }
    ta_possible_positions_neg <- nchar(seq_neg) - nchar("TA") + 1

    # Test negative strand scores
    expect_equal(scores_neg$count_ta, ta_count_neg)
    expect_equal(scores_neg$frac_ta, ta_count_neg / ta_possible_positions_neg, tolerance = 1e-5)

    # Also verify that "TA" on negative strand = "AT" on positive strand
    gvtrack.create("count_at", NULL, "kmer.count", kmer = "AT")
    scores_at <- gextract("count_at", pos_interval, iterator = pos_interval)
    expect_equal(scores_neg$count_ta, scores_at$count_at)
})

test_that("kmer.count and kmer.frac handle repeated k-mers", {
    remove_all_vtracks()

    # Create test interval with repeating pattern
    test_interval <- gintervals(1, 200, 220)
    seq <- toupper(gseq.extract(test_interval))

    # Track for repeating k-mer
    gvtrack.create("count_rep", NULL, "kmer.count", kmer = "TAACCCT")
    gvtrack.create("frac_rep", NULL, "kmer.frac", kmer = "TAACCCT")

    # Extract scores
    scores <- gextract(c("count_rep", "frac_rep"), test_interval, iterator = test_interval)

    # Manual count of repeating k-mer
    kmer <- "TAACCCT"
    rep_count <- 0
    for (i in 1:(nchar(seq) - nchar(kmer) + 1)) {
        if (substr(seq, i, i + nchar(kmer) - 1) == kmer) {
            rep_count <- rep_count + 1
        }
    }
    rep_possible_positions <- nchar(seq) - nchar(kmer) + 1

    # Test counts
    expect_equal(scores$count_rep, rep_count)
    expect_equal(scores$frac_rep, rep_count / rep_possible_positions, tolerance = 1e-5)
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
    scores <- gextract(c("chrom", "start", "end", "count_ccc"),
        test_interval,
        iterator = 10
    )

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

test_that("kmer.count and kmer.frac performance is acceptable", {
    remove_all_vtracks()

    # Skip on CRAN
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
