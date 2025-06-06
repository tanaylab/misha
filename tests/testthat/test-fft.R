# Helper function to manually calculate FFT on kmer occurrence signal
manual_kmer_fft <- function(sequence, kmer, window = "hann", freq = NULL, original_len = -1, seq_offset = 0) {
    kmer <- toupper(kmer)
    sequence <- toupper(sequence)

    signal_len <- if (original_len == -1) nchar(sequence) else original_len

    # Create binary signal for kmer occurrences
    signal <- numeric(signal_len)
    kmer_length <- nchar(kmer)

    if (kmer_length > 0 && signal_len > 0) {
        # Mark positions where kmer occurs
        for (i in 1:signal_len) {
            seq_idx <- i + seq_offset
            if (seq_idx + kmer_length - 1 <= nchar(sequence)) {
                if (substr(sequence, seq_idx, seq_idx + kmer_length - 1) == kmer) {
                    signal[i] <- 1
                }
            }
        }
    }

    # Apply window function
    n <- length(signal)
    if (n < 2) {
        if (!is.null(freq)) {
            return(0)
        } else {
            return(list(
                frequencies = numeric(0),
                power_spectrum = numeric(0),
                peak_freq = NaN,
                peak_power = NaN
            ))
        }
    }

    if (window == "hann") {
        window_func <- 0.5 * (1 - cos(2 * pi * (0:(n - 1)) / (n - 1)))
    } else if (window == "hamming") {
        window_func <- 0.54 - 0.46 * cos(2 * pi * (0:(n - 1)) / (n - 1))
    } else if (window == "blackman") {
        window_func <- 0.42 - 0.5 * cos(2 * pi * (0:(n - 1)) / (n - 1)) + 0.08 * cos(4 * pi * (0:(n - 1)) / (n - 1))
    } else { # "none"
        window_func <- rep(1, n)
    }

    windowed_signal <- signal * window_func

    # Calculate FFT
    fft_result <- fft(windowed_signal)
    power_spectrum <- abs(fft_result)^2 / n

    # Frequency bins (cycles per base)
    frequencies <- (0:(n - 1)) / n

    if (!is.null(freq)) {
        # Find closest frequency bin
        freq_idx <- which.min(abs(frequencies - freq))
        return(power_spectrum[freq_idx])
    } else {
        # Return peak frequency and power
        # Exclude DC component (frequency 0) and only check positive frequencies
        n_positive <- floor(n / 2) + 1
        if (n_positive < 2) { # Only DC component
            peak_idx <- 1
        } else {
            peak_idx <- which.max(power_spectrum[2:n_positive]) + 1
        }
        peak_freq <- frequencies[peak_idx]
        peak_power <- power_spectrum[peak_idx]

        return(list(
            frequencies = frequencies,
            power_spectrum = power_spectrum,
            peak_freq = peak_freq,
            peak_power = peak_power
        ))
    }
}

test_that("kmer.fft basic functionality works", {
    remove_all_vtracks()

    # Create test interval with known content
    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals))

    # Create FFT virtual tracks
    gvtrack.create("fft_cg", NULL, "kmer.fft",
        params = list(kmer = "CG", freq = 0.1, extend = FALSE, window = "none")
    )
    gvtrack.create("fft_peak", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", extend = FALSE, window = "none")
    )
    gvtrack.create("fft_peak_power", NULL, "kmer.fft.peak.power",
        params = list(kmer = "CG", extend = FALSE, window = "none")
    )

    # Extract scores
    scores <- gextract(c("fft_cg", "fft_peak", "fft_peak_power"),
        test_intervals,
        iterator = test_intervals
    )

    # Manual calculation
    manual_result <- manual_kmer_fft(seq, "CG", window = "none")
    manual_freq_power <- manual_kmer_fft(seq, "CG", window = "none", freq = 0.1)

    # Test basic functionality
    expect_true(is.numeric(scores$fft_cg))
    expect_true(is.numeric(scores$fft_peak))
    expect_true(is.numeric(scores$fft_peak_power))

    # Test ranges
    expect_true(scores$fft_peak >= 0 && scores$fft_peak <= 0.5) # Nyquist limit
    expect_true(scores$fft_peak_power >= 0)
    expect_true(scores$fft_cg >= 0)

    # Test correctness - compare with manual calculations
    expect_equal(scores$fft_cg, manual_freq_power, tolerance = 1e-6)
    expect_equal(scores$fft_peak, manual_result$peak_freq, tolerance = 1e-6)
    expect_equal(scores$fft_peak_power, manual_result$peak_power, tolerance = 1e-6)
})

test_that("kmer.fft works with different window functions", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)

    # Create FFT tracks with different window functions
    gvtrack.create("fft_none", NULL, "kmer.fft",
        params = list(kmer = "CG", freq = 0.1, window = "none")
    )
    gvtrack.create("fft_hann", NULL, "kmer.fft",
        params = list(kmer = "CG", freq = 0.1, window = "hann")
    )
    gvtrack.create("fft_hamming", NULL, "kmer.fft",
        params = list(kmer = "CG", freq = 0.1, window = "hamming")
    )
    gvtrack.create("fft_blackman", NULL, "kmer.fft",
        params = list(kmer = "CG", freq = 0.1, window = "blackman")
    )

    # Extract scores
    scores <- gextract(c("fft_none", "fft_hann", "fft_hamming", "fft_blackman"),
        test_intervals,
        iterator = test_intervals
    )

    # All should return valid numeric results
    expect_true(all(sapply(scores[2:5], is.numeric)))
    expect_true(all(sapply(scores[2:5], function(x) x >= 0)))

    # Different windows should potentially give different results
    # (but we don't enforce they must be different as it depends on the signal)
    expect_true(is.numeric(scores$fft_none))
    expect_true(is.numeric(scores$fft_hann))
    expect_true(is.numeric(scores$fft_hamming))
    expect_true(is.numeric(scores$fft_blackman))
})

test_that("kmer.fft works with different frequencies", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)

    # Test different frequencies
    frequencies <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

    for (i in seq_along(frequencies)) {
        vtrack_name <- paste0("fft_freq_", i)
        gvtrack.create(vtrack_name, NULL, "kmer.fft",
            params = list(kmer = "CG", freq = frequencies[i], window = "hann")
        )
    }

    # Extract all scores
    vtrack_names <- paste0("fft_freq_", seq_along(frequencies))
    scores <- gextract(vtrack_names, test_intervals, iterator = test_intervals)

    # All should return valid results
    for (name in vtrack_names) {
        expect_true(is.numeric(scores[[name]]))
        expect_true(scores[[name]] >= 0)
    }

    # Test frequency at exactly 0.0 and 0.5 (boundary cases)
    gvtrack.create("fft_dc", NULL, "kmer.fft",
        params = list(kmer = "CG", freq = 0.0, window = "none")
    )
    gvtrack.create("fft_nyquist", NULL, "kmer.fft",
        params = list(kmer = "CG", freq = 0.5, window = "none")
    )

    boundary_scores <- gextract(c("fft_dc", "fft_nyquist"), test_intervals, iterator = test_intervals)
    expect_true(is.numeric(boundary_scores$fft_dc))
    expect_true(is.numeric(boundary_scores$fft_nyquist))
})

test_that("kmer.fft correctness with different frequencies", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 260)
    kmer <- "CG"
    kmer_len <- nchar(kmer)
    extended_interval <- gintervals(test_intervals$chrom, test_intervals$start - (kmer_len - 1), test_intervals$end + (kmer_len - 1))
    extended_seq <- toupper(gseq.extract(extended_interval))
    original_len <- test_intervals$end - test_intervals$start
    seq_offset <- test_intervals$start - extended_interval$start

    # Test multiple frequencies
    frequencies <- c(0.0, 0.05, 0.1, 0.25, 0.4, 0.5)

    for (i in seq_along(frequencies)) {
        vtrack_name <- paste0("fft_", i)
        gvtrack.create(vtrack_name, NULL, "kmer.fft",
            params = list(kmer = kmer, freq = frequencies[i], window = "none")
        )
    }

    # Extract all scores
    vtrack_names <- paste0("fft_", seq_along(frequencies))
    scores <- gextract(vtrack_names, test_intervals, iterator = test_intervals)

    # Compare with manual calculations
    for (i in seq_along(frequencies)) {
        manual_power <- manual_kmer_fft(extended_seq, kmer,
            window = "none", freq = frequencies[i],
            original_len = original_len, seq_offset = seq_offset
        )
        expect_equal(scores[[vtrack_names[i]]], manual_power,
            tolerance = 1e-6,
            info = paste("Failed for frequency", frequencies[i])
        )
    }
})

test_that("kmer.fft correctness with different window functions", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 260)
    kmer <- "CG"
    kmer_len <- nchar(kmer)
    extended_interval <- gintervals(test_intervals$chrom, test_intervals$start - (kmer_len - 1), test_intervals$end + (kmer_len - 1))
    extended_seq <- toupper(gseq.extract(extended_interval))
    original_len <- test_intervals$end - test_intervals$start
    seq_offset <- test_intervals$start - extended_interval$start

    windows <- c("none", "hann", "hamming", "blackman")
    test_freq <- 0.1

    for (window in windows) {
        vtrack_name <- paste0("fft_", window)
        gvtrack.create(vtrack_name, NULL, "kmer.fft",
            params = list(kmer = kmer, freq = test_freq, window = window)
        )
    }

    # Extract scores
    vtrack_names <- paste0("fft_", windows)
    scores <- gextract(vtrack_names, test_intervals, iterator = test_intervals)

    # Compare with manual calculations
    for (window in windows) {
        vtrack_name <- paste0("fft_", window)
        manual_power <- manual_kmer_fft(extended_seq, kmer,
            window = window, freq = test_freq,
            original_len = original_len, seq_offset = seq_offset
        )
        expect_equal(scores[[vtrack_name]], manual_power,
            tolerance = 1e-6,
            info = paste("Failed for window", window)
        )
    }
})

test_that("kmer.fft correctness with numeric iterators", {
    remove_all_vtracks()

    # Test with numeric iterator - these are relative to chromosome start, not interval start
    test_intervals <- gintervals(1, 200, 280)
    kmer <- "CG"
    kmer_len <- nchar(kmer)

    gvtrack.create("fft_cg", NULL, "kmer.fft",
        params = list(kmer = kmer, freq = 0.1, window = "none")
    )
    gvtrack.create("fft_peak", NULL, "kmer.fft.peak",
        params = list(kmer = kmer, window = "none")
    )

    # Extract with 20bp bins - iterator positions are relative to chromosome start
    scores <- gextract(c("fft_cg", "fft_peak"), test_intervals, iterator = 20)

    expect_equal(nrow(scores), 4) # 80bp interval / 20bp bins = 4 rows

    # Verify each bin by manual calculation
    for (i in 1:nrow(scores)) {
        bin_interval <- gintervals(scores$chrom[i], scores$start[i], scores$end[i])
        extended_interval <- gintervals(scores$chrom[i], scores$start[i] - (kmer_len - 1), scores$end[i] + (kmer_len - 1))
        bin_seq <- toupper(gseq.extract(extended_interval))
        original_len <- bin_interval$end - bin_interval$start
        seq_offset <- bin_interval$start - extended_interval$start

        # Manual calculations for this bin
        manual_result <- manual_kmer_fft(bin_seq, kmer, window = "none", original_len = original_len, seq_offset = seq_offset)
        manual_freq_power <- manual_kmer_fft(bin_seq, kmer,
            window = "none", freq = 0.1,
            original_len = original_len, seq_offset = seq_offset
        )

        # Compare
        expect_equal(scores$fft_cg[i], manual_freq_power,
            tolerance = 1e-6,
            info = paste("Failed for bin", i, "at position", scores$start[i], "-", scores$end[i])
        )
        expect_equal(scores$fft_peak[i], manual_result$peak_freq,
            tolerance = 1e-6,
            info = paste("Failed for peak frequency in bin", i)
        )
    }
})

test_that("kmer.fft.peak and kmer.fft.peak.power are consistent", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)

    # Create peak detection tracks
    gvtrack.create("peak_freq", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", window = "hann")
    )
    gvtrack.create("peak_power", NULL, "kmer.fft.peak.power",
        params = list(kmer = "CG", window = "hann")
    )

    # Extract scores
    scores <- gextract(c("peak_freq", "peak_power"), test_intervals, iterator = test_intervals)

    # Now create a kmer.fft track at the detected peak frequency
    gvtrack.create("fft_at_peak", NULL, "kmer.fft",
        params = list(kmer = "CG", freq = scores$peak_freq, window = "hann")
    )

    peak_validation <- gextract("fft_at_peak", test_intervals, iterator = test_intervals)

    # The power at the peak frequency should match the peak power
    expect_equal(peak_validation$fft_at_peak, scores$peak_power, tolerance = 1e-6)

    # Peak frequency should be in valid range
    expect_true(scores$peak_freq >= 0 && scores$peak_freq <= 0.5)
    expect_true(scores$peak_power >= 0)
})

test_that("kmer.fft correctness with different kmers", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 280)
    test_freq <- 0.1

    # Test with different kmer lengths and sequences
    kmers <- c("A", "CG", "TAA", "ACGT")

    for (i in seq_along(kmers)) {
        vtrack_name <- paste0("fft_", i)
        gvtrack.create(vtrack_name, NULL, "kmer.fft",
            params = list(kmer = kmers[i], freq = test_freq, window = "none")
        )

        peak_name <- paste0("peak_", i)
        gvtrack.create(peak_name, NULL, "kmer.fft.peak",
            params = list(kmer = kmers[i], window = "none")
        )
    }

    vtrack_names <- paste0("fft_", seq_along(kmers))
    peak_names <- paste0("peak_", seq_along(kmers))
    all_names <- c(vtrack_names, peak_names)

    scores <- gextract(all_names, test_intervals, iterator = test_intervals)

    # Compare each kmer with manual calculations
    for (i in seq_along(kmers)) {
        kmer <- kmers[i]
        kmer_len <- nchar(kmer)

        # Account for extension by kmer_len - 1
        extended_interval <- gintervals(
            test_intervals$chrom,
            test_intervals$start - (kmer_len - 1),
            test_intervals$end + (kmer_len - 1)
        )
        extended_seq <- toupper(gseq.extract(extended_interval))
        original_len <- test_intervals$end - test_intervals$start
        seq_offset <- test_intervals$start - extended_interval$start

        manual_result <- manual_kmer_fft(extended_seq, kmer, window = "none", original_len = original_len, seq_offset = seq_offset)
        manual_freq_power <- manual_kmer_fft(extended_seq, kmer, window = "none", freq = test_freq, original_len = original_len, seq_offset = seq_offset)

        expect_equal(scores[[vtrack_names[i]]], manual_freq_power,
            tolerance = 1e-6,
            info = paste("Failed for kmer", kmers[i], "at frequency", test_freq)
        )
        expect_equal(scores[[peak_names[i]]], manual_result$peak_freq,
            tolerance = 1e-6,
            info = paste("Failed for kmer", kmers[i], "peak frequency")
        )
    }
})

test_that("kmer.fft handles case sensitivity correctly", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)

    # Create tracks with different case patterns
    gvtrack.create("fft_upper", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", window = "hann")
    )
    gvtrack.create("fft_lower", NULL, "kmer.fft.peak",
        params = list(kmer = "cg", window = "hann")
    )
    gvtrack.create("fft_mixed", NULL, "kmer.fft.peak",
        params = list(kmer = "Cg", window = "hann")
    )

    scores <- gextract(c("fft_upper", "fft_lower", "fft_mixed"),
        test_intervals,
        iterator = test_intervals
    )

    # Case should not matter - all should give same results
    expect_equal(scores$fft_upper, scores$fft_lower, tolerance = 1e-10)
    expect_equal(scores$fft_lower, scores$fft_mixed, tolerance = 1e-10)
})

test_that("kmer.fft handles extension parameter correctly", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 210) # Small interval

    # Create tracks with and without extension
    gvtrack.create("fft_extend", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", extend = TRUE, window = "none")
    )
    gvtrack.create("fft_no_extend", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", extend = FALSE, window = "none")
    )

    scores <- gextract(c("fft_extend", "fft_no_extend"),
        test_intervals,
        iterator = test_intervals
    )

    # Both should return valid results
    expect_true(is.numeric(scores$fft_extend))
    expect_true(is.numeric(scores$fft_no_extend))
    expect_true(scores$fft_extend >= 0 && scores$fft_extend <= 0.5)
    expect_true(scores$fft_no_extend >= 0 && scores$fft_no_extend <= 0.5)

    # Extension might affect results depending on the sequence
    # We just verify both are valid, not necessarily different
})

test_that("kmer.fft works with small intervals", {
    remove_all_vtracks()

    # Test with very small intervals
    small_intervals <- gintervals(1, 200, 205) # 5bp

    gvtrack.create("fft_small", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", window = "hann")
    )

    scores <- gextract("fft_small", small_intervals, iterator = small_intervals)

    # Should handle small intervals gracefully
    expect_true(is.numeric(scores$fft_small))
    expect_true(scores$fft_small >= 0 && scores$fft_small <= 0.5)
})

test_that("kmer.fft works with multiple intervals", {
    remove_all_vtracks()

    # Test with iterator creating multiple bins
    test_intervals <- gintervals(1, 200, 240)

    gvtrack.create("fft_peak", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", window = "hann")
    )
    gvtrack.create("fft_power", NULL, "kmer.fft.peak.power",
        params = list(kmer = "CG", window = "hann")
    )

    # Extract with 10bp bins
    scores <- gextract(c("fft_peak", "fft_power"), test_intervals, iterator = 10)

    # Should have 4 rows
    expect_equal(nrow(scores), 4)

    # All values should be valid
    expect_true(all(scores$fft_peak >= 0 & scores$fft_peak <= 0.5))
    expect_true(all(scores$fft_power >= 0))
    expect_true(all(is.finite(scores$fft_peak)))
    expect_true(all(is.finite(scores$fft_power)))
})

test_that("kmer.fft parameter validation works", {
    remove_all_vtracks()

    # Missing kmer parameter
    expect_error(gvtrack.create("bad1", NULL, "kmer.fft",
        params = list(freq = 0.1)
    ))

    # Missing freq parameter for kmer.fft
    expect_error(gvtrack.create("bad3", NULL, "kmer.fft",
        params = list(kmer = "CG")
    ))

    # Invalid window function
    expect_error(gvtrack.create("bad4", NULL, "kmer.fft",
        params = list(kmer = "CG", freq = 0.1, window = "invalid")
    ))

    # Invalid frequency (negative)
    expect_error(gvtrack.create("bad5", NULL, "kmer.fft",
        params = list(kmer = "CG", freq = -0.1)
    ))

    # Invalid frequency (> 0.5)
    expect_error(gvtrack.create("bad6", NULL, "kmer.fft",
        params = list(kmer = "CG", freq = 0.6)
    ))

    # Multiple kmers should error
    expect_error(gvtrack.create("bad7", NULL, "kmer.fft",
        params = list(kmer = c("CG", "AT"), freq = 0.1)
    ))

    # Non-string kmer
    expect_error(gvtrack.create("bad8", NULL, "kmer.fft",
        params = list(kmer = 123, freq = 0.1)
    ))
})

test_that("kmer.fft works with different iterator sizes", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 1, 80) # 80bp interval

    gvtrack.create("fft_peak", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", window = "hann")
    )

    # Test different iterator sizes
    scores_10bp <- gextract("fft_peak", test_intervals, iterator = 10)
    scores_20bp <- gextract("fft_peak", test_intervals, iterator = 20)
    scores_40bp <- gextract("fft_peak", test_intervals, iterator = 40)
    scores_80bp <- gextract("fft_peak", test_intervals, iterator = 80)

    # Check expected number of bins
    expect_equal(nrow(scores_10bp), 8)
    expect_equal(nrow(scores_20bp), 4)
    expect_equal(nrow(scores_40bp), 2)
    expect_equal(nrow(scores_80bp), 1)

    # All should return valid frequencies
    expect_true(all(scores_10bp$fft_peak >= 0 & scores_10bp$fft_peak <= 0.5))
    expect_true(all(scores_20bp$fft_peak >= 0 & scores_20bp$fft_peak <= 0.5))
    expect_true(all(scores_40bp$fft_peak >= 0 & scores_40bp$fft_peak <= 0.5))
    expect_true(all(scores_80bp$fft_peak >= 0 & scores_80bp$fft_peak <= 0.5))
})

test_that("kmer.fft works with gdist", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 300)

    gvtrack.create("fft_peak", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", window = "hann")
    )
    gvtrack.create("fft_power", NULL, "kmer.fft.peak.power",
        params = list(kmer = "CG", window = "hann")
    )

    # Use gdist to create distributions
    peak_dist <- gdist("fft_peak",
        breaks = seq(0, 0.5, by = 0.1),
        intervals = test_intervals, iterator = 10
    )
    power_dist <- gdist("fft_power",
        breaks = c(0, 1, 10, 100, 1000, Inf),
        intervals = test_intervals, iterator = 10
    )

    # Verify distributions
    expect_true(all(peak_dist >= 0))
    expect_true(all(power_dist >= 0))
    expect_equal(length(peak_dist), 5) # 5 bins for 6 breaks
    expect_equal(length(power_dist), 5) # 5 bins for 6 breaks
})

test_that("kmer.fft preserves position information", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 240)

    gvtrack.create("fft_peak", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", window = "hann")
    )

    scores <- gextract("fft_peak", test_intervals, iterator = 10) %>%
        arrange(intervalID)

    # Verify position columns are preserved
    expect_true("chrom" %in% colnames(scores))
    expect_true("start" %in% colnames(scores))
    expect_true("end" %in% colnames(scores))
    expect_true("fft_peak" %in% colnames(scores))

    # Check positions make sense
    expect_equal(scores$start[1], 200)
    expect_equal(scores$end[4], 240)
    expect_equal(nrow(scores), 4)
})

test_that("kmer.fft works with other virtual track functions", {
    remove_all_vtracks()

    test_intervals <- gintervals(1, 200, 300)

    # Create different types of virtual tracks
    gvtrack.create("fft_peak", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", window = "hann")
    )
    gvtrack.create("count_cg", NULL, "kmer.count", kmer = "CG", strand = 1)
    gvtrack.create("frac_cg", NULL, "kmer.frac", kmer = "CG", strand = 1)

    # Extract all together
    scores <- gextract(c("fft_peak", "count_cg", "frac_cg"),
        test_intervals,
        iterator = 20
    )

    # Verify all tracks return results
    expect_true(all(c("fft_peak", "count_cg", "frac_cg") %in% colnames(scores)))
    expect_equal(nrow(scores), 5)

    # All should be valid numeric values
    expect_true(all(is.numeric(scores$fft_peak)))
    expect_true(all(is.numeric(scores$count_cg)))
    expect_true(all(is.numeric(scores$frac_cg)))
})

test_that("kmer.fft functions work across chromosomes", {
    remove_all_vtracks()

    # Test on different chromosomes
    interval_chr1 <- gintervals(1, 200, 240)
    interval_chr2 <- gintervals(2, 200, 240)

    gvtrack.create("fft_peak", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", window = "hann")
    )

    scores_chr1 <- gextract("fft_peak", interval_chr1, iterator = interval_chr1)
    scores_chr2 <- gextract("fft_peak", interval_chr2, iterator = interval_chr2)

    # Both should return valid results
    expect_true(is.numeric(scores_chr1$fft_peak))
    expect_true(is.numeric(scores_chr2$fft_peak))
    expect_true(scores_chr1$fft_peak >= 0 && scores_chr1$fft_peak <= 0.5)
    expect_true(scores_chr2$fft_peak >= 0 && scores_chr2$fft_peak <= 0.5)
})

test_that("kmer.fft correctness with extend parameter", {
    remove_all_vtracks()

    # Test with small interval to see extend effect
    test_intervals <- gintervals(1, 200, 210)
    kmer <- "CG"
    kmer_len <- nchar(kmer)

    gvtrack.create("fft_extend", NULL, "kmer.fft",
        params = list(kmer = kmer, freq = 0.1, extend = TRUE, window = "none")
    )
    gvtrack.create("fft_no_extend", NULL, "kmer.fft",
        params = list(kmer = kmer, freq = 0.1, extend = FALSE, window = "none")
    )
    gvtrack.create("fft_peak_extend", NULL, "kmer.fft.peak",
        params = list(kmer = kmer, extend = TRUE, window = "none")
    )
    gvtrack.create("fft_peak_no_extend", NULL, "kmer.fft.peak",
        params = list(kmer = kmer, extend = FALSE, window = "none")
    )

    scores <- gextract(c("fft_extend", "fft_no_extend", "fft_peak_extend", "fft_peak_no_extend"),
        test_intervals,
        iterator = test_intervals
    )

    # Manual calculations
    extended_interval <- gintervals(test_intervals$chrom, test_intervals$start - (kmer_len - 1), test_intervals$end + (kmer_len - 1))
    seq_extend <- toupper(gseq.extract(extended_interval))
    seq_no_extend <- toupper(gseq.extract(test_intervals))
    original_len <- test_intervals$end - test_intervals$start
    seq_offset <- test_intervals$start - extended_interval$start

    manual_extend <- manual_kmer_fft(seq_extend, kmer,
        window = "none", freq = 0.1,
        original_len = original_len, seq_offset = seq_offset
    )
    manual_no_extend <- manual_kmer_fft(seq_no_extend, kmer, window = "none", freq = 0.1)

    manual_peak_extend <- manual_kmer_fft(seq_extend, kmer,
        window = "none",
        original_len = original_len, seq_offset = seq_offset
    )
    manual_peak_no_extend <- manual_kmer_fft(seq_no_extend, kmer, window = "none")

    expect_equal(scores$fft_extend, manual_extend, tolerance = 1e-6)
    expect_equal(scores$fft_no_extend, manual_no_extend, tolerance = 1e-6)
    expect_equal(scores$fft_peak_extend, manual_peak_extend$peak_freq, tolerance = 1e-6)
    expect_equal(scores$fft_peak_no_extend, manual_peak_no_extend$peak_freq, tolerance = 1e-6)
})

test_that("kmer.fft edge cases at chromosome boundaries", {
    remove_all_vtracks()

    # Test at chromosome start and end
    start_interval <- gintervals(1, 0, 20)
    end_interval <- gintervals(1, 9980, 10000)

    gvtrack.create("fft_peak_ext", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", extend = TRUE, window = "hann")
    )
    gvtrack.create("fft_peak_no_ext", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", extend = FALSE, window = "hann")
    )

    start_scores <- gextract(c("fft_peak_ext", "fft_peak_no_ext"),
        start_interval,
        iterator = start_interval
    )
    end_scores <- gextract(c("fft_peak_ext", "fft_peak_no_ext"),
        end_interval,
        iterator = end_interval
    )

    # Should handle boundaries gracefully
    expect_true(!is.na(start_scores$fft_peak_ext))
    expect_true(!is.na(start_scores$fft_peak_no_ext))
    expect_true(!is.na(end_scores$fft_peak_ext))
    expect_true(!is.na(end_scores$fft_peak_no_ext))

    # All should be in valid frequency range
    expect_true(start_scores$fft_peak_ext >= 0 && start_scores$fft_peak_ext <= 0.5)
    expect_true(end_scores$fft_peak_ext >= 0 && end_scores$fft_peak_ext <= 0.5)
})

test_that("kmer.fft performance is acceptable", {
    skip_on_cran()

    # Test with large interval
    large_interval <- gintervals(1, 1, 10000)

    gvtrack.create("fft_perf", NULL, "kmer.fft.peak",
        params = list(kmer = "CG", window = "hann")
    )

    # Measure timing
    timing <- system.time(
        gextract("fft_perf", large_interval, iterator = 1000)
    )

    # Should complete in reasonable time
    expect_true(timing["elapsed"] < 30) # 30 seconds max
})
