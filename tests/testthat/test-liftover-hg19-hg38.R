load_test_db()

# These tests use real hg19/hg38 genomes and the UCSC hg19ToHg38 chain file
# to verify our liftover implementation matches Kent's liftOver binary exactly

test_that("gintervals.liftover matches Kent liftOver on hg19->hg38 random intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    # Check if required genomes and chain file exist
    hg19_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19"
    hg38_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg38"
    chain_file <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19ToHg38.over.chain"

    skip_if_not(dir.exists(hg19_path), "hg19 genome not available")
    skip_if_not(dir.exists(hg38_path), "hg38 genome not available")
    skip_if_not(file.exists(chain_file), "hg19ToHg38 chain file not available")

    local_db_state()

    # Set up hg19 as source
    gsetroot(hg19_path)
    gdb.reload()

    # Generate random intervals
    set.seed(60427)
    random_intervals <- gintervals.random(n = 5000, size = 500)
    random_intervals <- gintervals.canonic(random_intervals)
    random_intervals$id <- paste0("id_", seq_len(nrow(random_intervals)))

    # Create BED file for Kent's liftOver
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    write.table(
        data.frame(
            chrom = random_intervals$chrom,
            start = random_intervals$start,
            end = random_intervals$end,
            name = random_intervals$id,
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run Kent's liftOver
    system2("liftOver",
        args = c(bed_input, chain_file, bed_output, bed_unmapped),
        stdout = FALSE, stderr = FALSE
    )

    # Read Kent's output
    kent_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(kent_result) <- c("chrom", "start", "end", "name", "score", "strand")
    kent_result <- kent_result[order(kent_result$chrom, kent_result$start), c("chrom", "start", "end", "name")]

    # Switch to hg38 (target database) to load chain and do liftover
    # This validates both chain target coords and result coords against hg38
    gsetroot(hg38_path)
    gdb.reload()
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(random_intervals, chain)

    # Order results the same way
    misha_result <- misha_result[order(misha_result$chrom, misha_result$start), ]

    # Compare number of lifted intervals (allow small differences for edge cases)
    pct_diff <- abs(nrow(misha_result) - nrow(kent_result)) / nrow(kent_result) * 100
    expect_lt(pct_diff, 1.0,
        label = sprintf(
            "Misha lifted %d intervals, Kent lifted %d (%.2f%% difference)",
            nrow(misha_result), nrow(kent_result), pct_diff
        )
    )

    # Match intervals by intervalID/name to handle potential ordering differences
    # intervalID is numeric, but BED name is "id_N", so we need to match on the numeric part
    kent_result$id_numeric <- as.integer(sub("id_", "", kent_result$name))
    misha_merged <- merge(misha_result, kent_result,
        by.x = "intervalID", by.y = "id_numeric",
        suffixes = c("_misha", "_kent")
    )

    # At least 99% of Kent intervals should be found in misha
    pct_matched <- nrow(misha_merged) / nrow(kent_result) * 100
    expect_gt(pct_matched, 99.0,
        label = sprintf("%.2f%% of Kent intervals found in misha results", pct_matched)
    )

    # For matched intervals, coordinates should match exactly (or very closely for edge cases)
    if (nrow(misha_merged) > 0) {
        # Check chromosomes match
        chrom_match <- sum(as.character(misha_merged$chrom_misha) == misha_merged$chrom_kent)
        chrom_pct <- chrom_match / nrow(misha_merged) * 100
        expect_gt(chrom_pct, 99.9,
            label = sprintf("%.2f%% of chromosomes match", chrom_pct)
        )

        # Check coordinates - allow small differences for edge cases
        coord_match <- sum(misha_merged$start_misha == misha_merged$start_kent &
            misha_merged$end_misha == misha_merged$end_kent)
        coord_pct <- coord_match / nrow(misha_merged) * 100
        expect_gt(coord_pct, 98.0,
            label = sprintf("%.2f%% of coordinates match exactly", coord_pct)
        )
    }
})

test_that("gintervals.liftover matches Kent liftOver on hg19->hg38 with larger random set", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    # Check if required genomes and chain file exist
    hg19_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19"
    hg38_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg38"
    chain_file <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19ToHg38.over.chain"

    skip_if_not(dir.exists(hg19_path), "hg19 genome not available")
    skip_if_not(dir.exists(hg38_path), "hg38 genome not available")
    skip_if_not(file.exists(chain_file), "hg19ToHg38 chain file not available")

    local_db_state()

    # Set up hg19 as source
    gsetroot(hg19_path)
    gdb.reload()

    # Generate larger set of random intervals with varying sizes
    set.seed(12345)
    random_intervals <- gintervals.random(n = 5000, size = 1000)
    random_intervals <- gintervals.canonic(random_intervals)
    random_intervals <- gintervals.normalize(random_intervals, 1000)
    random_intervals$id <- paste0("id_", seq_len(nrow(random_intervals)))

    # Create BED file for Kent's liftOver
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    write.table(
        data.frame(
            chrom = random_intervals$chrom,
            start = random_intervals$start,
            end = random_intervals$end,
            name = random_intervals$id,
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run Kent's liftOver
    system2("liftOver",
        args = c(bed_input, chain_file, bed_output, bed_unmapped),
        stdout = FALSE, stderr = FALSE
    )

    # Read Kent's output
    kent_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(kent_result) <- c("chrom", "start", "end", "name", "score", "strand")
    kent_result <- kent_result[, c("chrom", "start", "end", "name")]

    # Switch to hg38 (target database) to load chain and do liftover
    # This validates both chain target coords and result coords against hg38
    gsetroot(hg38_path)
    gdb.reload()
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(random_intervals, chain)

    # Compare number of lifted intervals (allow small differences for edge cases)
    pct_diff <- abs(nrow(misha_result) - nrow(kent_result)) / nrow(kent_result) * 100
    expect_lt(pct_diff, 1.0,
        label = sprintf(
            "Misha lifted %d intervals, Kent lifted %d (%.2f%% difference)",
            nrow(misha_result), nrow(kent_result), pct_diff
        )
    )

    # Match intervals by intervalID/name
    # intervalID is numeric, but BED name is "id_N", so we need to match on the numeric part
    kent_result$id_numeric <- as.integer(sub("id_", "", kent_result$name))
    misha_merged <- merge(misha_result, kent_result,
        by.x = "intervalID", by.y = "id_numeric",
        suffixes = c("_misha", "_kent")
    )

    # At least 99% of Kent intervals should be found in misha
    pct_matched <- nrow(misha_merged) / nrow(kent_result) * 100
    expect_gt(pct_matched, 99.0,
        label = sprintf("%.2f%% of Kent intervals found in misha results", pct_matched)
    )

    # For matched intervals, coordinates should match exactly (or very closely for edge cases)
    if (nrow(misha_merged) > 0) {
        # Check chromosomes match
        chrom_match <- sum(as.character(misha_merged$chrom_misha) == misha_merged$chrom_kent)
        chrom_pct <- chrom_match / nrow(misha_merged) * 100
        expect_gt(chrom_pct, 99.9,
            label = sprintf("%.2f%% of chromosomes match", chrom_pct)
        )

        # Check coordinates - allow small differences for edge cases
        coord_match <- sum(misha_merged$start_misha == misha_merged$start_kent &
            misha_merged$end_misha == misha_merged$end_kent)
        coord_pct <- coord_match / nrow(misha_merged) * 100
        expect_gt(coord_pct, 98.0,
            label = sprintf("%.2f%% of coordinates match exactly", coord_pct)
        )
    }

    # Verify percentage of original intervals that were successfully lifted
    pct_lifted_misha <- length(unique(misha_result$intervalID)) / nrow(random_intervals) * 100
    pct_lifted_kent <- length(unique(kent_result$name)) / nrow(random_intervals) * 100

    # Should be within 1% of each other
    expect_lt(abs(pct_lifted_misha - pct_lifted_kent), 1.0,
        label = sprintf(
            "Misha lifted %.2f%%, Kent lifted %.2f%% of original intervals",
            pct_lifted_misha, pct_lifted_kent
        )
    )
})

test_that("gintervals.liftover matches Kent liftOver on hg19->hg38 small intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    # Check if required genomes and chain file exist
    hg19_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19"
    hg38_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg38"
    chain_file <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19ToHg38.over.chain"

    skip_if_not(dir.exists(hg19_path), "hg19 genome not available")
    skip_if_not(dir.exists(hg38_path), "hg38 genome not available")
    skip_if_not(file.exists(chain_file), "hg19ToHg38 chain file not available")

    local_db_state()

    # Set up hg19 as source
    gsetroot(hg19_path)
    gdb.reload()

    # Generate random intervals with small sizes (1-10 bp)
    set.seed(99999)
    random_intervals <- gintervals.random(n = 2000, size = 10)
    random_intervals <- gintervals.canonic(random_intervals)
    # Randomly adjust sizes to be between 1 and 10 bp
    random_intervals$end <- random_intervals$start + sample(1:10, nrow(random_intervals), replace = TRUE)
    random_intervals$id <- paste0("id_", seq_len(nrow(random_intervals)))

    # Create BED file for Kent's liftOver
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    write.table(
        data.frame(
            chrom = random_intervals$chrom,
            start = random_intervals$start,
            end = random_intervals$end,
            name = random_intervals$id,
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run Kent's liftOver
    system2("liftOver",
        args = c(bed_input, chain_file, bed_output, bed_unmapped),
        stdout = FALSE, stderr = FALSE
    )

    # Read Kent's output
    if (file.exists(bed_output) && file.info(bed_output)$size > 0) {
        kent_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
        colnames(kent_result) <- c("chrom", "start", "end", "name", "score", "strand")
        kent_result <- kent_result[, c("chrom", "start", "end", "name")]
    } else {
        kent_result <- data.frame(
            chrom = character(), start = numeric(),
            end = numeric(), name = character()
        )
    }

    # Switch to hg38 (target database) to load chain and do liftover
    # This validates both chain target coords and result coords against hg38
    gsetroot(hg38_path)
    gdb.reload()
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(random_intervals, chain)

    if (is.null(misha_result)) {
        misha_result <- data.frame(
            chrom = character(), start = numeric(),
            end = numeric(), intervalID = integer()
        )
    }

    # Compare number of lifted intervals (allow small differences for edge cases)
    if (nrow(kent_result) > 0) {
        pct_diff <- abs(nrow(misha_result) - nrow(kent_result)) / nrow(kent_result) * 100
        expect_lt(pct_diff, 5.0,
            label = sprintf(
                "Misha lifted %d intervals, Kent lifted %d (%.2f%% difference)",
                nrow(misha_result), nrow(kent_result), pct_diff
            )
        )

        # Match intervals by intervalID/name
        # intervalID is numeric, but BED name is "id_N", so we need to match on the numeric part
        kent_result$id_numeric <- as.integer(sub("id_", "", kent_result$name))
        misha_merged <- merge(misha_result, kent_result,
            by.x = "intervalID", by.y = "id_numeric",
            suffixes = c("_misha", "_kent")
        )

        # At least 95% of Kent intervals should be found in misha (small intervals may differ more)
        pct_matched <- nrow(misha_merged) / nrow(kent_result) * 100
        expect_gt(pct_matched, 95.0,
            label = sprintf("%.2f%% of Kent intervals found in misha results", pct_matched)
        )

        # For matched intervals, coordinates should match exactly
        if (nrow(misha_merged) > 0) {
            expect_equal(as.character(misha_merged$chrom_misha), misha_merged$chrom_kent)
            expect_equal(misha_merged$start_misha, misha_merged$start_kent)
            expect_equal(misha_merged$end_misha, misha_merged$end_kent)
        }
    } else {
        expect_equal(nrow(misha_result), 0)
    }
})

test_that("gintervals.liftover matches Kent liftOver on hg19->hg38 specific chromosomes", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    # Check if required genomes and chain file exist
    hg19_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19"
    hg38_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg38"
    chain_file <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19ToHg38.over.chain"

    skip_if_not(dir.exists(hg19_path), "hg19 genome not available")
    skip_if_not(dir.exists(hg38_path), "hg38 genome not available")
    skip_if_not(file.exists(chain_file), "hg19ToHg38 chain file not available")

    local_db_state()

    # Set up hg19 as source
    gsetroot(hg19_path)
    gdb.reload()

    # Generate intervals on specific chromosomes with known differences
    set.seed(55555)
    test_chroms <- c("chr1", "chr2", "chrX", "chrY", "chrM")

    random_intervals <- data.frame()
    for (chrom in test_chroms) {
        chrom_intervals <- gintervals.random(n = 200, size = 500)
        chrom_intervals <- chrom_intervals[chrom_intervals$chrom == chrom, ]
        if (nrow(chrom_intervals) > 0) {
            random_intervals <- rbind(random_intervals, chrom_intervals)
        }
    }

    if (nrow(random_intervals) == 0) {
        skip("No intervals generated for test chromosomes")
    }

    random_intervals <- gintervals.canonic(random_intervals)
    random_intervals$id <- paste0("id_", seq_len(nrow(random_intervals)))

    # Create BED file for Kent's liftOver
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    write.table(
        data.frame(
            chrom = random_intervals$chrom,
            start = random_intervals$start,
            end = random_intervals$end,
            name = random_intervals$id,
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run Kent's liftOver
    system2("liftOver",
        args = c(bed_input, chain_file, bed_output, bed_unmapped),
        stdout = FALSE, stderr = FALSE
    )

    # Read Kent's output
    if (file.exists(bed_output) && file.info(bed_output)$size > 0) {
        kent_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
        colnames(kent_result) <- c("chrom", "start", "end", "name", "score", "strand")
        kent_result <- kent_result[, c("chrom", "start", "end", "name")]
    } else {
        kent_result <- data.frame(
            chrom = character(), start = numeric(),
            end = numeric(), name = character()
        )
    }

    # Switch to hg38 (target database) to load chain and do liftover
    # This validates both chain target coords and result coords against hg38
    gsetroot(hg38_path)
    gdb.reload()
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(random_intervals, chain)

    if (is.null(misha_result)) {
        misha_result <- data.frame(
            chrom = character(), start = numeric(),
            end = numeric(), intervalID = integer()
        )
    }

    # Compare number of lifted intervals (allow small differences for edge cases)
    if (nrow(kent_result) > 0) {
        pct_diff <- abs(nrow(misha_result) - nrow(kent_result)) / nrow(kent_result) * 100
        expect_lt(pct_diff, 5.0,
            label = sprintf(
                "Misha lifted %d intervals, Kent lifted %d (%.2f%% difference)",
                nrow(misha_result), nrow(kent_result), pct_diff
            )
        )

        # Match intervals by intervalID/name
        # intervalID is numeric, but BED name is "id_N", so we need to match on the numeric part
        kent_result$id_numeric <- as.integer(sub("id_", "", kent_result$name))
        misha_merged <- merge(misha_result, kent_result,
            by.x = "intervalID", by.y = "id_numeric",
            suffixes = c("_misha", "_kent")
        )

        # At least 95% of Kent intervals should be found in misha (small intervals may differ more)
        pct_matched <- nrow(misha_merged) / nrow(kent_result) * 100
        expect_gt(pct_matched, 95.0,
            label = sprintf("%.2f%% of Kent intervals found in misha results", pct_matched)
        )

        # For matched intervals, coordinates should match exactly
        if (nrow(misha_merged) > 0) {
            expect_equal(as.character(misha_merged$chrom_misha), misha_merged$chrom_kent)
            expect_equal(misha_merged$start_misha, misha_merged$start_kent)
            expect_equal(misha_merged$end_misha, misha_merged$end_kent)
        }
    } else {
        expect_equal(nrow(misha_result), 0)
    }
})

test_that("gtrack.liftover matches Kent liftOver on hg19->hg38 random sparse track", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    # Check if required genomes and chain file exist
    hg19_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19"
    hg38_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg38"
    chain_file <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19ToHg38.over.chain"

    skip_if_not(dir.exists(hg19_path), "hg19 genome not available")
    skip_if_not(dir.exists(hg38_path), "hg38 genome not available")
    skip_if_not(file.exists(chain_file), "hg19ToHg38 chain file not available")

    local_db_state()

    # Set up hg19 as source
    gsetroot(hg19_path)
    gdb.reload()

    # Generate random intervals with values
    set.seed(77777)
    random_intervals <- gintervals.random(n = 1000, size = 500)
    random_intervals <- gintervals.canonic(random_intervals)
    random_intervals$id <- paste0("id_", seq_len(nrow(random_intervals)))
    random_intervals$value <- runif(nrow(random_intervals), 0, 100)

    # Create sparse track
    temp_track <- "test_temp_track_hg19_hg38"
    if (gtrack.exists(temp_track)) gtrack.rm(temp_track, force = TRUE)

    withr::defer({
        if (gtrack.exists(temp_track)) gtrack.rm(temp_track, force = TRUE)
    })

    gtrack.create_sparse(
        temp_track, "Random test intervals",
        random_intervals[, c("chrom", "start", "end")],
        random_intervals$value
    )

    src_track_dir <- file.path(hg19_path, "tracks", paste0(temp_track, ".track"))

    # Create BED file for Kent's liftOver
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    write.table(
        data.frame(
            chrom = random_intervals$chrom,
            start = random_intervals$start,
            end = random_intervals$end,
            name = random_intervals$id,
            score = random_intervals$value,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run Kent's liftOver
    system2("liftOver",
        args = c(bed_input, chain_file, bed_output, bed_unmapped),
        stdout = FALSE, stderr = FALSE
    )

    # Read Kent's output
    kent_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(kent_result) <- c("chrom", "start", "end", "name", "value", "strand")
    kent_result <- kent_result[order(kent_result$chrom, kent_result$start), ]

    # Switch to hg38 for misha liftover
    gsetroot(hg38_path)
    gdb.reload()

    # Run gtrack.liftover
    lifted_track <- "test_lifted_track_hg19_hg38"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted test track", src_track_dir, chain_file)

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    misha_result <- misha_result[order(misha_result$chrom, misha_result$start), ]
    colnames(misha_result)[colnames(misha_result) == lifted_track] <- "value"

    # Compare number of lifted intervals (allow small differences for edge cases)
    pct_diff <- abs(nrow(misha_result) - nrow(kent_result)) / nrow(kent_result) * 100
    expect_lt(pct_diff, 1.0,
        label = sprintf(
            "Misha lifted %d intervals, Kent lifted %d (%.2f%% difference)",
            nrow(misha_result), nrow(kent_result), pct_diff
        )
    )

    # Find common intervals by coordinates
    misha_coords <- paste(misha_result$chrom, misha_result$start, misha_result$end, sep = "_")
    kent_coords <- paste(kent_result$chrom, kent_result$start, kent_result$end, sep = "_")

    common_coords <- intersect(misha_coords, kent_coords)
    pct_common <- length(common_coords) / nrow(kent_result) * 100

    # At least 99% of intervals should match
    expect_gt(pct_common, 99.0,
        label = sprintf("%.2f%% of Kent intervals found in misha results", pct_common)
    )

    # For common intervals, values should match closely
    misha_common <- misha_result[misha_coords %in% common_coords, ]
    kent_common <- kent_result[kent_coords %in% common_coords, ]

    # Sort both by coordinates to ensure proper matching
    misha_common <- misha_common[order(misha_common$chrom, misha_common$start), ]
    kent_common <- kent_common[order(kent_common$chrom, kent_common$start), ]

    # Compare values with tolerance for floating point
    expect_equal(misha_common$value, kent_common$value, tolerance = 1e-5)
})
