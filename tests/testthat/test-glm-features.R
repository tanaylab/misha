local_mm10_db <- Sys.getenv("MISHA_MM10_MOTIF_DB", "")
mm10_trackdb <- Sys.getenv("MISHA_MM10_TRACKDB", "")

test_that("glm_extract_features matches R reference pipeline", {
    # Requires local mm10 database with motif energy tracks
    skip_if(!dir.exists(local_mm10_db), "Local mm10 misha db not available")
    gdb.init(local_mm10_db)
    gdataset.load(mm10_trackdb)
    options(gmax.data.size = 1e9)

    # Pick 3 motifs and 50 peaks for a fast comparison
    motifs <- c("JASPAR_Zic2", "HOMER_Unknown_ESC_element", "JOLMA_ZIC3_mono_full")
    track_dir <- "th_epi.all_db_motifs"
    track_names <- paste0(track_dir, ".", motifs)

    # Create 50 peaks on chr1 at known positions (300bp each)
    set.seed(42)
    peak_starts <- sort(sample(seq(5e6, 50e6, by = 1000), 50))
    peaks <- data.frame(
        chrom = "chr1",
        start = peak_starts,
        end = peak_starts + 300L
    )
    n_peaks <- nrow(peaks)

    # Config matching the LASSO training pipeline
    tile_size <- 200L
    flank_size <- 350L
    scale_factor <- 10
    dist_from_max <- 10
    top_percentile <- 0.9999

    # ---- R reference pipeline ----

    # 1. Create virtual tracks (LSE aggregation)
    for (motif in motifs) {
        gvtrack.create(vtrack = motif, src = paste0(track_dir, ".", motif), func = "lse")
    }

    # 2. Compute genome-wide quantiles
    shift <- (300 - 20) / 2 # (peak_size - 20) / 2
    gw_max_q <- sapply(motifs, function(motif) {
        gvtrack.iterator(motif, sshift = -shift, eshift = shift)
        gquantiles(motif, percentiles = top_percentile, iterator = 20)
    })
    names(gw_max_q) <- motifs

    # 3. Tile peaks
    tile_peaks_r <- function(peaks, tile_size, flank_size) {
        peaks_ext <- peaks
        peaks_ext$start <- peaks$start - flank_size
        peaks_ext$end <- peaks$end + flank_size
        peaks_ext <- gintervals.force_range(peaks_ext)

        do.call(rbind, lapply(seq_len(nrow(peaks_ext)), function(i) {
            starts <- seq(peaks_ext$start[i], peaks_ext$end[i] - 1, by = tile_size)
            ends <- pmin(starts + tile_size, peaks_ext$end[i])
            peak_mid <- (peaks$start[i] + peaks$end[i]) / 2
            tile_mid <- starts + tile_size / 2
            data.frame(
                chrom = peaks_ext$chrom[i],
                start = starts,
                end = ends,
                peak_id = i,
                dist = tile_mid - peak_mid
            )
        }))
    }

    tiled <- tile_peaks_r(peaks, tile_size, flank_size)
    n_tiles <- length(unique(tiled$dist))

    # 4. Extract motif energies (R pipeline)
    # Reset virtual track iterators to no-shift for tile-level extraction
    for (motif in motifs) {
        gvtrack.iterator(motif, sshift = 0, eshift = 0)
    }

    motif_energy <- gextract(motifs, intervals = tiled, iterator = tiled) %>%
        dplyr::arrange(intervalID)
    e_mat <- as.matrix(motif_energy[, motifs])

    # 5. Scale
    m_names <- rep(motifs, each = 1) # motifs already match column order
    ceiled <- pmin(sweep(e_mat, 2, gw_max_q[motifs], "-"), 0)
    floored <- pmax(ceiled, -dist_from_max)
    scaled_r <- scale_factor * (floored + dist_from_max) / dist_from_max
    scaled_r[!is.finite(scaled_r)] <- 0

    # 6. Logistic transforms
    logist <- function(x, x0, L, k) L / (1 + exp(-k * (x - x0)))

    f1 <- logist(scaled_r, x0 = 0, L = 2, k = 0.5) - 1
    f2 <- logist(scaled_r, x0 = 10, L = 2, k = 0.50)
    f3 <- logist(scaled_r - 5, x0 = 0, L = 1, k = 1)
    f4 <- logist(scaled_r, x0 = 10, L = 2, k = 1)

    # 7. GC features
    gvtrack.create("gc_sum", src = "seq.G_or_C", func = "sum")
    gc_raw <- gextract("gc_sum",
        intervals = tiled, iterator = tiled,
        colnames = "gc"
    ) %>%
        dplyr::arrange(intervalID) %>%
        dplyr::pull(gc)
    gc_raw[is.na(gc_raw)] <- 0
    gc_scaled_r <- (gc_raw / tile_size) * scale_factor

    # 8. Reshape to wide format (direct matrix indexing)
    dist_vals <- sort(unique(tiled$dist))
    n_motifs <- length(motifs)

    # For each transform, create wide matrix
    ref_motif_block <- matrix(NA_real_, nrow = n_peaks, ncol = n_motifs * n_tiles * 4)
    ref_gc <- matrix(NA_real_, nrow = n_peaks, ncol = n_tiles)

    for (ti in seq_len(n_tiles)) {
        row_idx <- seq(ti, nrow(tiled), by = n_tiles)

        for (mi in seq_len(n_motifs)) {
            # Column in C++ output: (tile * n_motifs + motif) * n_transforms + transform
            # (0-based in C++, 1-based here)
            base_col <- ((ti - 1) * n_motifs + (mi - 1)) * 4

            ref_motif_block[, base_col + 1] <- f1[row_idx, mi]
            ref_motif_block[, base_col + 2] <- f2[row_idx, mi]
            ref_motif_block[, base_col + 3] <- f3[row_idx, mi]
            ref_motif_block[, base_col + 4] <- f4[row_idx, mi]
        }

        ref_gc[, ti] <- gc_scaled_r[row_idx]
    }

    # GC interactions
    n_gc_inter <- n_tiles * (n_tiles - 1) / 2
    ref_gc_inter <- matrix(NA_real_, nrow = n_peaks, ncol = n_gc_inter)
    idx <- 1
    for (a in seq_len(n_tiles - 1)) {
        for (b in (a + 1):n_tiles) {
            ref_gc_inter[, idx] <- ref_gc[, a] * ref_gc[, b] / scale_factor
            idx <- idx + 1
        }
    }

    ref_full <- cbind(ref_motif_block, ref_gc, ref_gc_inter)

    # ---- C++ pipeline ----
    max_cap <- gw_max_q
    names(max_cap) <- track_names

    cpp_result <- glm_extract_features(
        track_names = track_names,
        intervals = peaks,
        tile_size = tile_size,
        flank_size = flank_size,
        max_cap = max_cap,
        dis_from_cap = dist_from_max,
        scale_factor = scale_factor,
        gc_track = "seq.G_or_C",
        gc_scale_factor = scale_factor
    )

    # ---- Compare ----
    expect_equal(nrow(cpp_result), n_peaks)
    expect_equal(ncol(cpp_result), ncol(ref_full))

    # Compare values — should be bit-identical since we match misha's
    # float-precision lse_accumulate and use std::exp for transforms
    max_diff <- max(abs(cpp_result - ref_full), na.rm = TRUE)
    cat("Max absolute difference:", max_diff, "\n")
    expect_equal(max_diff, 0)

    # Check motif block specifically
    motif_cols <- seq_len(n_motifs * n_tiles * 4)
    max_motif_diff <- max(abs(cpp_result[, motif_cols] - ref_full[, motif_cols]), na.rm = TRUE)
    cat("Max motif block diff:", max_motif_diff, "\n")
    expect_equal(max_motif_diff, 0)

    # Check GC block
    gc_cols <- n_motifs * n_tiles * 4 + seq_len(n_tiles)
    max_gc_diff <- max(abs(cpp_result[, gc_cols] - ref_full[, gc_cols]), na.rm = TRUE)
    cat("Max GC diff:", max_gc_diff, "\n")
    expect_equal(max_gc_diff, 0)

    # Check GC interactions
    inter_cols <- n_motifs * n_tiles * 4 + n_tiles + seq_len(n_gc_inter)
    max_inter_diff <- max(abs(cpp_result[, inter_cols] - ref_full[, inter_cols]), na.rm = TRUE)
    cat("Max GC interaction diff:", max_inter_diff, "\n")
    expect_equal(max_inter_diff, 0)

    cat("All comparisons PASSED\n")
})
