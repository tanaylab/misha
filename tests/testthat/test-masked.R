create_isolated_test_db()

# Helper function to manually count masked base pairs
count_masked_manually <- function(sequence) {
    if (nchar(sequence) == 0) {
        return(list(count = 0, fraction = 0))
    }

    # Count lowercase characters
    chars <- strsplit(sequence, "")[[1]]
    masked_count <- sum(sapply(chars, function(c) {
        c %in% c("a", "c", "g", "t", "n")
    }))

    total_length <- nchar(sequence)
    fraction <- if (total_length > 0) masked_count / total_length else 0

    return(list(
        count = masked_count,
        fraction = fraction
    ))
}

test_that("masked.count and masked.frac work with basic inputs", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 300)
    seq <- gseq.extract(test_intervals)

    # Create virtual tracks
    gvtrack.create("masked_count", NULL, "masked.count")
    gvtrack.create("masked_frac", NULL, "masked.frac")

    # Extract scores
    scores <- gextract(c("masked_count", "masked_frac"), test_intervals, iterator = test_intervals)

    # Verify using manual counting
    manual_results <- count_masked_manually(seq)

    expect_equal(scores$masked_count, manual_results$count)
    expect_equal(scores$masked_frac, manual_results$fraction, tolerance = 1e-6)
})

test_that("masked functions handle all-uppercase sequences", {
    remove_all_vtracks()

    # Find interval with all uppercase sequence
    test_intervals <- gintervals(1, 10000000, 10000100)
    seq <- gseq.extract(test_intervals)

    gvtrack.create("masked_count", NULL, "masked.count")
    gvtrack.create("masked_frac", NULL, "masked.frac")

    scores <- gextract(c("masked_count", "masked_frac"), test_intervals, iterator = test_intervals)

    # Verify with manual counting
    manual_results <- count_masked_manually(seq)
    expect_equal(scores$masked_count, manual_results$count)
    expect_equal(scores$masked_frac, manual_results$fraction, tolerance = 1e-6)
})

test_that("masked functions handle sequences with lowercase bases", {
    remove_all_vtracks()

    # Test on various intervals
    test_intervals <- gintervals(1, 500, 600)
    seq <- gseq.extract(test_intervals)

    gvtrack.create("masked_count", NULL, "masked.count")
    gvtrack.create("masked_frac", NULL, "masked.frac")

    scores <- gextract(c("masked_count", "masked_frac"), test_intervals, iterator = test_intervals)

    manual_results <- count_masked_manually(seq)
    expect_equal(scores$masked_count, manual_results$count)
    expect_equal(scores$masked_frac, manual_results$fraction, tolerance = 1e-6)
})

test_that("masked functions work with iterator shifts", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 1000, 1100)

    # Create vtracks with different iterator shifts
    gvtrack.create("masked_full", NULL, "masked.frac")
    gvtrack.create("masked_shift", NULL, "masked.frac")
    gvtrack.iterator("masked_shift", sshift = 10, eshift = -10)

    scores <- gextract(c("masked_full", "masked_shift"), test_intervals, iterator = test_intervals)

    # Shifted version should evaluate on [1010, 1090)
    full_seq <- gseq.extract(test_intervals)
    shift_seq <- gseq.extract(gintervals(1, 1010, 1090))

    full_manual <- count_masked_manually(full_seq)
    shift_manual <- count_masked_manually(shift_seq)

    expect_equal(scores$masked_full, full_manual$fraction, tolerance = 1e-6)
    expect_equal(scores$masked_shift, shift_manual$fraction, tolerance = 1e-6)
})

test_that("masked functions work with multiple intervals", {
    remove_all_vtracks()

    test_intervals <- gintervals(
        1,
        c(200, 500, 1000),
        c(300, 600, 1100)
    )

    gvtrack.create("masked_count", NULL, "masked.count")
    gvtrack.create("masked_frac", NULL, "masked.frac")

    scores <- gextract(c("masked_count", "masked_frac"), test_intervals, iterator = test_intervals)

    expect_equal(nrow(scores), 3)

    # Verify each interval separately
    for (i in 1:3) {
        interval <- gintervals(test_intervals$chrom[i], test_intervals$start[i], test_intervals$end[i])
        seq <- gseq.extract(interval)
        manual <- count_masked_manually(seq)

        expect_equal(scores$masked_count[i], manual$count)
        expect_equal(scores$masked_frac[i], manual$fraction, tolerance = 1e-6)
    }
})

test_that("masked functions work with different chromosomes", {
    remove_all_vtracks()

    test_intervals <- rbind(
        gintervals(1, 200, 300),
        gintervals(2, 200, 300)
    )

    gvtrack.create("masked_frac", NULL, "masked.frac")

    scores <- gextract("masked_frac", test_intervals, iterator = test_intervals)

    expect_equal(nrow(scores), 2)

    # Verify each chromosome separately
    for (i in 1:2) {
        interval <- gintervals(test_intervals$chrom[i], test_intervals$start[i], test_intervals$end[i])
        seq <- gseq.extract(interval)
        manual <- count_masked_manually(seq)
        expect_equal(scores$masked_frac[i], manual$fraction, tolerance = 1e-6)
    }
})

test_that("masked.count returns integer-like values", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 300)

    gvtrack.create("masked_count", NULL, "masked.count")
    scores <- gextract("masked_count", test_intervals, iterator = test_intervals)

    # Count should be a whole number
    expect_equal(scores$masked_count, round(scores$masked_count))
})

test_that("masked.frac returns values between 0 and 1", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, c(200, 500, 1000), c(300, 600, 1100))

    gvtrack.create("masked_frac", NULL, "masked.frac")
    scores <- gextract("masked_frac", test_intervals, iterator = test_intervals)

    # Fraction should be between 0 and 1
    expect_true(all(scores$masked_frac >= 0))
    expect_true(all(scores$masked_frac <= 1))
})

test_that("masked functions work in expressions", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 400)

    gvtrack.create("masked_frac", NULL, "masked.frac")
    gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")

    # Test arithmetic with other vtracks
    scores <- gextract("masked_frac + g_frac", test_intervals,
        iterator = test_intervals,
        colnames = "combined"
    )

    expect_true("combined" %in% names(scores))
    expect_true(is.numeric(scores$combined))
})

test_that("masked functions ignore extra parameters with warning", {
    remove_all_vtracks()

    # Should warn about extra parameters
    expect_warning(
        gvtrack.create("masked_test", NULL, "masked.count", extra_param = 123),
        "do not accept parameters"
    )
})
