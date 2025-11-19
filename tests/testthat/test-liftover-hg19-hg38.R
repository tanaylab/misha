load_test_db()

# These tests use real hg19/hg38 genomes and the UCSC hg19ToHg38 chain file
# to verify our liftover implementation against Kent's liftOver binary
#
# Key behavioral differences between misha and Kent:
#
# 1. Split handling: Kent rejects intervals that split into multiple non-contiguous
#    parts in the target genome (marked as "#Split in new" in unmapped file).
#    Misha lifts these and produces multiple fragments.
#
# 2. Partial deletion: Kent rejects intervals that are "#Partially deleted in new"
#    when less than minMatch (default 95%) of bases map. Misha may still lift what remains.
#
# 3. minMatch threshold: Kent has -minMatch=0.95 by default. Misha doesn't have this
#    parameter and may lift intervals that Kent rejects.
#
# 4. Overlap policy: Kent doesn't pre-filter chains by target overlap - it keeps all
#    chains and maps based on source coverage. Use tgt_overlap_policy="keep" to match.
#
# 5. Coordinate differences: When both tools lift an interval, coordinate
#    differences may occur due to:
#    - Different fragment merging strategies (small differences)
#    - Misha producing multiple fragments when Kent merges them (can be large)
#    - Different chain selection for overlapping source regions
#
# These tests verify:
# - Intervals lifted only by misha have documented Kent unmapped reasons
# - Intervals lifted by both are on the same chromosome
# - Small percentage of intervals may have coordinate differences (< 1%)

# Helper function to parse Kent's unmapped file
parse_kent_unmapped <- function(unmapped_file) {
    lines <- readLines(unmapped_file)
    if (length(lines) == 0) return(data.frame(id = integer(), reason = character()))

    # Unmapped format: reason line (starting with #), then BED line
    reasons <- character()
    ids <- integer()

    i <- 1
    while (i <= length(lines)) {
        if (startsWith(lines[i], "#")) {
            reason <- lines[i]
            if (i + 1 <= length(lines) && !startsWith(lines[i + 1], "#")) {
                # Parse the BED line to get the ID
                bed_parts <- strsplit(lines[i + 1], "\t")[[1]]
                if (length(bed_parts) >= 4) {
                    id <- as.integer(sub("id_", "", bed_parts[4]))
                    ids <- c(ids, id)
                    reasons <- c(reasons, reason)
                }
                i <- i + 2
            } else {
                i <- i + 1
            }
        } else {
            i <- i + 1
        }
    }

    data.frame(id = ids, reason = reasons, stringsAsFactors = FALSE)
}

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

    # Run Kent's liftOver with default parameters
    # Kent uses chain scores to resolve overlapping targets by default
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

    # Load chain with policies that match Kent's behavior:
    # - tgt_overlap_policy = "keep": don't pre-filter chains by target overlap (Kent keeps all)
    # - src_overlap_policy = "keep": allow overlapping source regions in chains
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Run misha liftover with canonic = TRUE to merge adjacent intervals (matches Kent output)
    misha_result <- gintervals.liftover(random_intervals, chain, canonic = TRUE)

    # Order results the same way
    misha_result <- misha_result %>% arrange(chrom, start)

    # Match intervals by intervalID/name
    kent_result$id_numeric <- as.integer(sub("id_", "", kent_result$name))

    # Parse Kent's unmapped file to understand rejections
    kent_unmapped <- parse_kent_unmapped(bed_unmapped)

    # Get unique source intervals lifted by each tool
    misha_ids <- unique(misha_result$intervalID)
    kent_ids <- unique(kent_result$id_numeric)

    # Categorize differences
    only_in_misha <- setdiff(misha_ids, kent_ids)
    only_in_kent <- setdiff(kent_ids, misha_ids)
    in_both <- intersect(misha_ids, kent_ids)

    # For intervals only in misha: these should be explained by Kent's unmapped reasons
    # (split, partial deletion, etc.) - this is expected behavior difference
    if (length(only_in_misha) > 0) {
        # Check that these are in Kent's unmapped with expected reasons
        expected_reasons <- c("#Split in new", "#Partially deleted in new", "#Deleted in new")

        # First check all misha-only intervals are in Kent's unmapped
        unexplained <- only_in_misha[!only_in_misha %in% kent_unmapped$id]
        expect_equal(length(unexplained), 0,
            label = sprintf("Intervals lifted by misha but not in Kent's unmapped: %s",
                paste(head(unexplained, 10), collapse = ", "))
        )

        # Then verify the reasons are the expected ones
        misha_only_reasons <- kent_unmapped$reason[kent_unmapped$id %in% only_in_misha]
        unexpected_reasons <- unique(misha_only_reasons[!misha_only_reasons %in% expected_reasons])
        expect_equal(length(unexpected_reasons), 0,
            label = sprintf("Unexpected Kent unmapped reasons for misha-lifted intervals: %s",
                paste(unexpected_reasons, collapse = ", "))
        )
    }

    # For intervals only in Kent: misha rejected something Kent accepted
    # This may happen due to different chain coverage handling
    # Track but allow small percentage
    if (length(only_in_kent) > 0) {
        message(sprintf("Misha rejected %d intervals that Kent lifted: %s",
            length(only_in_kent), paste(head(only_in_kent, 10), collapse = ", ")))
    }
    only_kent_pct <- length(only_in_kent) / nrow(random_intervals) * 100
    expect_lt(only_kent_pct, 0.5,
        label = sprintf("%.2f%% of intervals lifted only by Kent (%d of %d)",
            only_kent_pct, length(only_in_kent), nrow(random_intervals))
    )

    # For intervals in both, compare coordinates
    misha_merged <- merge(misha_result, kent_result,
        by.x = "intervalID", by.y = "id_numeric",
        suffixes = c("_misha", "_kent")
    )

    # Check for coordinate mismatches
    mismatches <- misha_merged %>%
        filter(
            as.character(chrom_misha) != chrom_kent |
                start_misha != start_kent |
                end_misha != end_kent
        )

    if (nrow(mismatches) > 0) {
        message(sprintf("Coordinate mismatches: %d intervals out of %d", nrow(mismatches), length(in_both)))

        # Verify mismatches are on the same chromosome (critical - ensures same region)
        # Coordinate differences can be large due to:
        # - Different chain alignment selection (both valid but different paths)
        # - Different fragment merging when intervals split
        same_chrom <- mismatches %>%
            filter(as.character(chrom_misha) == chrom_kent)

        # All mismatches should be on the same chromosome
        expect_equal(nrow(same_chrom), nrow(mismatches),
            label = sprintf("Chromosome mismatches found: %d intervals on different chromosomes",
                nrow(mismatches) - nrow(same_chrom))
        )

        # Track coordinate differences for informational purposes
        # Note: Large differences can occur when misha produces multiple fragments
        # that Kent merges into one (due to different chain selection logic)
        if (nrow(same_chrom) > 0) {
            start_diffs <- abs(same_chrom$start_misha - same_chrom$start_kent)
            end_diffs <- abs(same_chrom$end_misha - same_chrom$end_kent)
            max_diff <- max(c(start_diffs, end_diffs))
            message(sprintf("Max coordinate difference: %d bp", max_diff))
        }
    }

    # Coordinates should match for intervals lifted by both
    # Allow small number of edge cases due to different merging strategies
    mismatch_pct <- nrow(mismatches) / length(in_both) * 100
    expect_lt(mismatch_pct, 1.0,
        label = sprintf("%.2f%% coordinate mismatches (%d of %d)",
            mismatch_pct, nrow(mismatches), length(in_both))
    )
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

    # Run Kent's liftOver with default parameters
    system2("liftOver",
        args = c(bed_input, chain_file, bed_output, bed_unmapped),
        stdout = FALSE, stderr = FALSE
    )

    # Read Kent's output
    kent_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(kent_result) <- c("chrom", "start", "end", "name", "score", "strand")
    kent_result <- kent_result[, c("chrom", "start", "end", "name")]

    # Switch to hg38 (target database) to load chain and do liftover
    gsetroot(hg38_path)
    gdb.reload()

    # Load chain with Kent-matching policies
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Run misha liftover with canonic = TRUE to merge adjacent fragments (like Kent does)
    misha_result <- gintervals.liftover(random_intervals, chain, canonic = TRUE)

    # Match intervals by intervalID/name
    kent_result$id_numeric <- as.integer(sub("id_", "", kent_result$name))

    # Parse Kent's unmapped file
    kent_unmapped <- parse_kent_unmapped(bed_unmapped)

    # Get unique source intervals lifted by each tool
    misha_ids <- unique(misha_result$intervalID)
    kent_ids <- unique(kent_result$id_numeric)

    # Categorize differences
    only_in_misha <- setdiff(misha_ids, kent_ids)
    only_in_kent <- setdiff(kent_ids, misha_ids)
    in_both <- intersect(misha_ids, kent_ids)

    # For intervals only in misha: check Kent's unmapped reasons
    if (length(only_in_misha) > 0) {
        expected_reasons <- c("#Split in new", "#Partially deleted in new", "#Deleted in new")

        unexplained <- only_in_misha[!only_in_misha %in% kent_unmapped$id]
        expect_equal(length(unexplained), 0,
            label = sprintf("Unexplained misha-only intervals: %s",
                paste(head(unexplained, 10), collapse = ", "))
        )

        # Verify the reasons are the expected ones
        misha_only_reasons <- kent_unmapped$reason[kent_unmapped$id %in% only_in_misha]
        unexpected_reasons <- unique(misha_only_reasons[!misha_only_reasons %in% expected_reasons])
        expect_equal(length(unexpected_reasons), 0,
            label = sprintf("Unexpected Kent unmapped reasons: %s",
                paste(unexpected_reasons, collapse = ", "))
        )
    }

    # For intervals only in Kent: track but allow small percentage
    only_kent_pct <- length(only_in_kent) / nrow(random_intervals) * 100
    expect_lt(only_kent_pct, 0.5,
        label = sprintf("%.2f%% of intervals lifted only by Kent (%d of %d)",
            only_kent_pct, length(only_in_kent), nrow(random_intervals))
    )

    # For intervals in both, compare coordinates
    misha_merged <- merge(misha_result, kent_result,
        by.x = "intervalID", by.y = "id_numeric",
        suffixes = c("_misha", "_kent")
    )

    mismatches <- misha_merged %>%
        filter(
            as.character(chrom_misha) != chrom_kent |
                start_misha != start_kent |
                end_misha != end_kent
        )

    if (nrow(mismatches) > 0) {
        # Verify mismatches are on the same chromosome
        same_chrom <- mismatches %>%
            filter(as.character(chrom_misha) == chrom_kent)

        expect_equal(nrow(same_chrom), nrow(mismatches),
            label = sprintf("Chromosome mismatches: %d intervals", nrow(mismatches) - nrow(same_chrom))
        )
    }

    # Allow small percentage of coordinate mismatches
    mismatch_pct <- nrow(mismatches) / length(in_both) * 100
    expect_lt(mismatch_pct, 1.0,
        label = sprintf("%.2f%% coordinate mismatches (%d of %d)",
            mismatch_pct, nrow(mismatches), length(in_both))
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
    gsetroot(hg38_path)
    gdb.reload()

    # Load chain with Kent-matching policies
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Use canonic = TRUE to merge adjacent fragments (like Kent does)
    misha_result <- gintervals.liftover(random_intervals, chain, canonic = TRUE)

    if (is.null(misha_result)) {
        misha_result <- data.frame(
            chrom = character(), start = numeric(),
            end = numeric(), intervalID = integer()
        )
    }

    # Compare which intervals were lifted
    if (nrow(kent_result) > 0) {
        kent_result$id_numeric <- as.integer(sub("id_", "", kent_result$name))
        kent_unmapped <- parse_kent_unmapped(bed_unmapped)

        misha_ids <- unique(misha_result$intervalID)
        kent_ids <- unique(kent_result$id_numeric)

        only_in_misha <- setdiff(misha_ids, kent_ids)
        only_in_kent <- setdiff(kent_ids, misha_ids)
        in_both <- intersect(misha_ids, kent_ids)

        # For intervals only in misha: check Kent's unmapped reasons
        if (length(only_in_misha) > 0) {
            expected_reasons <- c("#Split in new", "#Partially deleted in new", "#Deleted in new")

            unexplained <- only_in_misha[!only_in_misha %in% kent_unmapped$id]
            expect_equal(length(unexplained), 0,
                label = sprintf("Unexplained misha-only intervals: %s",
                    paste(head(unexplained, 10), collapse = ", "))
            )

            misha_only_reasons <- kent_unmapped$reason[kent_unmapped$id %in% only_in_misha]
            unexpected_reasons <- unique(misha_only_reasons[!misha_only_reasons %in% expected_reasons])
            expect_equal(length(unexpected_reasons), 0,
                label = sprintf("Unexpected Kent unmapped reasons: %s",
                    paste(unexpected_reasons, collapse = ", "))
            )
        }

        # For intervals only in Kent: allow small percentage
        only_kent_pct <- length(only_in_kent) / nrow(random_intervals) * 100
        expect_lt(only_kent_pct, 0.5,
            label = sprintf("%.2f%% of intervals lifted only by Kent (%d of %d)",
                only_kent_pct, length(only_in_kent), nrow(random_intervals))
        )

        # For intervals in both, compare coordinates
        if (length(in_both) > 0) {
            misha_merged <- merge(misha_result, kent_result,
                by.x = "intervalID", by.y = "id_numeric",
                suffixes = c("_misha", "_kent")
            )

            mismatches <- misha_merged %>%
                filter(
                    as.character(chrom_misha) != chrom_kent |
                        start_misha != start_kent |
                        end_misha != end_kent
                )

            if (nrow(mismatches) > 0) {
                same_chrom <- mismatches %>%
                    filter(as.character(chrom_misha) == chrom_kent)

                expect_equal(nrow(same_chrom), nrow(mismatches),
                    label = sprintf("Chromosome mismatches: %d intervals", nrow(mismatches) - nrow(same_chrom))
                )
            }

            # Allow small percentage of coordinate mismatches
            mismatch_pct <- nrow(mismatches) / length(in_both) * 100
            expect_lt(mismatch_pct, 1.0,
                label = sprintf("%.2f%% coordinate mismatches (%d of %d)",
                    mismatch_pct, nrow(mismatches), length(in_both))
            )
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
    gsetroot(hg38_path)
    gdb.reload()

    # Load chain with Kent-matching policies
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Use canonic = TRUE to merge adjacent fragments (like Kent does)
    misha_result <- gintervals.liftover(random_intervals, chain, canonic = TRUE)

    if (is.null(misha_result)) {
        misha_result <- data.frame(
            chrom = character(), start = numeric(),
            end = numeric(), intervalID = integer()
        )
    }

    # Compare which intervals were lifted
    if (nrow(kent_result) > 0) {
        kent_result$id_numeric <- as.integer(sub("id_", "", kent_result$name))
        kent_unmapped <- parse_kent_unmapped(bed_unmapped)

        misha_ids <- unique(misha_result$intervalID)
        kent_ids <- unique(kent_result$id_numeric)

        only_in_misha <- setdiff(misha_ids, kent_ids)
        only_in_kent <- setdiff(kent_ids, misha_ids)
        in_both <- intersect(misha_ids, kent_ids)

        # For intervals only in misha: check Kent's unmapped reasons
        if (length(only_in_misha) > 0) {
            expected_reasons <- c("#Split in new", "#Partially deleted in new", "#Deleted in new")

            unexplained <- only_in_misha[!only_in_misha %in% kent_unmapped$id]
            expect_equal(length(unexplained), 0,
                label = sprintf("Unexplained misha-only intervals: %s",
                    paste(head(unexplained, 10), collapse = ", "))
            )

            misha_only_reasons <- kent_unmapped$reason[kent_unmapped$id %in% only_in_misha]
            unexpected_reasons <- unique(misha_only_reasons[!misha_only_reasons %in% expected_reasons])
            expect_equal(length(unexpected_reasons), 0,
                label = sprintf("Unexpected Kent unmapped reasons: %s",
                    paste(unexpected_reasons, collapse = ", "))
            )
        }

        # For intervals only in Kent: allow small percentage
        only_kent_pct <- length(only_in_kent) / nrow(random_intervals) * 100
        expect_lt(only_kent_pct, 1.0,
            label = sprintf("%.2f%% of intervals lifted only by Kent (%d of %d)",
                only_kent_pct, length(only_in_kent), nrow(random_intervals))
        )

        # For intervals in both, compare coordinates
        if (length(in_both) > 0) {
            misha_merged <- merge(misha_result, kent_result,
                by.x = "intervalID", by.y = "id_numeric",
                suffixes = c("_misha", "_kent")
            )

            mismatches <- misha_merged %>%
                filter(
                    as.character(chrom_misha) != chrom_kent |
                        start_misha != start_kent |
                        end_misha != end_kent
                )

            if (nrow(mismatches) > 0) {
                same_chrom <- mismatches %>%
                    filter(as.character(chrom_misha) == chrom_kent)

                expect_equal(nrow(same_chrom), nrow(mismatches),
                    label = sprintf("Chromosome mismatches: %d intervals", nrow(mismatches) - nrow(same_chrom))
                )
            }

            # Allow small percentage of coordinate mismatches
            mismatch_pct <- nrow(mismatches) / length(in_both) * 100
            expect_lt(mismatch_pct, 1.0,
                label = sprintf("%.2f%% coordinate mismatches (%d of %d)",
                    mismatch_pct, nrow(mismatches), length(in_both))
            )
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

    # Create sparse track in hg19 database
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
    # Kent uses chain scores to resolve overlapping targets by default
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

    # Run Kent's liftOver with default parameters
    # Kent uses chain scores to resolve overlapping targets by default
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

    # Run gtrack.liftover with Kent-matching parameters:
    # - src_overlap_policy = "keep": allow overlapping source regions in chains
    # - tgt_overlap_policy = "keep": don't pre-filter chains by target overlap (Kent keeps all)
    lifted_track <- "test_lifted_track_hg19_hg38"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(
        lifted_track,
        "Lifted test track",
        src_track_dir,
        chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "keep"
    )

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    misha_result <- misha_result[order(misha_result$chrom, misha_result$start), ]
    colnames(misha_result)[colnames(misha_result) == lifted_track] <- "value"

    # For track liftover, compare by coordinates since gtrack.liftover
    # doesn't preserve source interval IDs in the output
    misha_coords <- paste(misha_result$chrom, misha_result$start, misha_result$end, sep = "_")
    kent_coords <- paste(kent_result$chrom, kent_result$start, kent_result$end, sep = "_")

    # Categorize differences
    only_in_misha <- setdiff(misha_coords, kent_coords)
    only_in_kent <- setdiff(kent_coords, misha_coords)
    common_coords <- intersect(misha_coords, kent_coords)

    total_intervals <- nrow(random_intervals)

    if (length(only_in_misha) > 0 || length(only_in_kent) > 0) {
        message(sprintf(
            "Coordinate differences: %d only in misha, %d only in Kent, %d in common",
            length(only_in_misha), length(only_in_kent), length(common_coords)
        ))
    }

    # Allow small percentage of differences
    # (due to split handling, minMatch differences)
    only_misha_pct <- length(only_in_misha) / total_intervals * 100
    only_kent_pct <- length(only_in_kent) / total_intervals * 100

    expect_lt(only_misha_pct, 1.0,
        label = sprintf("%.2f%% of intervals only in misha (%d of %d)",
            only_misha_pct, length(only_in_misha), total_intervals)
    )
    expect_lt(only_kent_pct, 1.0,
        label = sprintf("%.2f%% of intervals only in Kent (%d of %d)",
            only_kent_pct, length(only_in_kent), total_intervals)
    )

    # For common intervals, values should match exactly
    misha_common <- misha_result[misha_coords %in% common_coords, ]
    kent_common <- kent_result[kent_coords %in% common_coords, ]

    # Sort both by coordinates to ensure proper matching
    misha_common <- misha_common[order(misha_common$chrom, misha_common$start), ]
    kent_common <- kent_common[order(kent_common$chrom, kent_common$start), ]

    # Compare values with tolerance for floating point
    expect_equal(misha_common$value, kent_common$value, tolerance = 1e-5)
})
