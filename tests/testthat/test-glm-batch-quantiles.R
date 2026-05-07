local_mm10_db <- Sys.getenv("MISHA_MM10_MOTIF_DB", "")
mm10_trackdb <- Sys.getenv("MISHA_MM10_TRACKDB", "")

test_that("glm_batch_quantiles matches gquantiles for multiple tracks", {
    # Requires local mm10 database with motif energy tracks
    skip_if(!dir.exists(local_mm10_db), "Local mm10 misha db not available")
    gdb.init(local_mm10_db)
    gdataset.load(mm10_trackdb)
    options(gmax.data.size = 1e9, gmultitasking = FALSE)

    # Pick 3 motifs for testing
    motifs <- c("JASPAR_Zic2", "HOMER_Unknown_ESC_element", "JOLMA_ZIC3_mono_full")
    track_dir <- "th_epi.all_db_motifs"
    track_names <- paste0(track_dir, ".", motifs)

    # Test with single percentile (returns named vector)
    percentile <- 0.9999
    sshift <- -140L
    eshift <- 140L
    iterator <- 20L

    # ---- R reference: gquantiles with virtual tracks ----
    ref_q <- sapply(seq_along(motifs), function(i) {
        vt_name <- paste0("bq", i)
        gvtrack.create(vtrack = vt_name, src = paste0(track_dir, ".", motifs[i]), func = "lse")
        gvtrack.iterator(vt_name, sshift = sshift, eshift = eshift)
        q <- gquantiles(vt_name, percentiles = percentile, iterator = iterator)
        gvtrack.rm(vt_name)
        q
    })
    names(ref_q) <- track_names

    # ---- C++ batch quantiles ----
    batch_q <- glm_batch_quantiles(
        track_names = track_names,
        percentiles = percentile,
        iterator = iterator,
        sshift = sshift,
        eshift = eshift,
        n_threads = 3L
    )

    cat("Reference quantiles (gquantiles w/ StreamPercentiler):\n")
    print(ref_q)
    cat("Batch quantiles (exact nth_element):\n")
    print(batch_q)

    # Results should be close but not bit-identical because gquantiles uses
    # StreamPercentiler (approximate) while glm_batch_quantiles uses exact
    # nth_element. Tolerance of 0.1 is generous; typical differences are ~0.02.
    expect_equal(length(batch_q), length(ref_q))
    expect_equal(names(batch_q), names(ref_q))

    for (i in seq_along(track_names)) {
        diff <- abs(batch_q[i] - ref_q[i])
        cat(sprintf(
            "  %s: ref=%.6f, batch=%.6f, diff=%.4e\n",
            track_names[i], ref_q[i], batch_q[i], diff
        ))
        expect_lt(diff, 0.1,
            label = paste("quantile difference for", track_names[i])
        )
    }

    cat("Single-percentile comparison PASSED\n")
})


test_that("glm_batch_quantiles works with multiple percentiles", {
    skip_if(!dir.exists(local_mm10_db), "Local mm10 misha db not available")
    gdb.init(local_mm10_db)
    gdataset.load(mm10_trackdb)
    options(gmax.data.size = 1e9, gmultitasking = FALSE)

    track_names <- c(
        "th_epi.all_db_motifs.JASPAR_Zic2",
        "th_epi.all_db_motifs.HOMER_Unknown_ESC_element"
    )
    percentiles <- c(0.5, 0.99, 0.9999)

    result <- glm_batch_quantiles(
        track_names = track_names,
        percentiles = percentiles,
        iterator = 20L,
        sshift = -140L,
        eshift = 140L,
        n_threads = 2L
    )

    # Should return a matrix
    expect_true(is.matrix(result))
    expect_equal(nrow(result), length(track_names))
    expect_equal(ncol(result), length(percentiles))
    expect_equal(rownames(result), track_names)

    # Values should be monotonically increasing across percentiles for each track
    for (i in seq_len(nrow(result))) {
        expect_true(all(diff(result[i, ]) >= 0),
            label = paste("monotonicity for", track_names[i])
        )
    }

    cat("Multiple-percentile test PASSED\n")
    cat("Result matrix:\n")
    print(result)
})


test_that("glm_batch_quantiles handles single track", {
    skip_if(!dir.exists(local_mm10_db), "Local mm10 misha db not available")
    gdb.init(local_mm10_db)
    gdataset.load(mm10_trackdb)
    options(gmax.data.size = 1e9, gmultitasking = FALSE)

    result <- glm_batch_quantiles(
        track_names = "th_epi.all_db_motifs.JASPAR_Zic2",
        percentiles = 0.9999,
        iterator = 20L,
        sshift = -140L,
        eshift = 140L,
        n_threads = 1L
    )
    expect_true(is.numeric(result))
    expect_equal(length(result), 1)
    expect_true(is.finite(result[1]))

    cat("Single track test PASSED, value:", result[1], "\n")
})


test_that("glm_batch_quantiles is self-consistent across thread counts", {
    skip_if(!dir.exists(local_mm10_db), "Local mm10 misha db not available")
    gdb.init(local_mm10_db)
    gdataset.load(mm10_trackdb)
    options(gmax.data.size = 1e9, gmultitasking = FALSE)

    tracks <- c(
        "th_epi.all_db_motifs.JASPAR_Zic2",
        "th_epi.all_db_motifs.HOMER_Unknown_ESC_element",
        "th_epi.all_db_motifs.JOLMA_ZIC3_mono_full"
    )

    # Run with 1 thread (sequential)
    q1 <- glm_batch_quantiles(tracks, 0.9999, n_threads = 1L)

    # Run with 3 threads (parallel)
    q3 <- glm_batch_quantiles(tracks, 0.9999, n_threads = 3L)

    # Should be identical (same code path, just concurrent)
    expect_equal(q1, q3, tolerance = 0)

    cat("Thread consistency test PASSED\n")
})
