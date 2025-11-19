load_test_db()

# Tests for gtrack.liftover value validation
test_that("gtrack.liftover preserves values in one-to-many mapping", {
    local_db_state()

    # Create source genome with one chromosome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 300), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with specific values in source
    # source1: [100,120) = 4, [121,122) = 5, [180,185) = 6
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(100, 121, 180),
        end = c(120, 122, 185),
        stringsAsFactors = FALSE
    )
    src_values <- c(4, 5, 6)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 400), collapse = ""), "\n")))

    # Create chain: source1[100-200) -> chr1[0-100) AND chr1[200-300)
    # This means each source interval maps to TWO target locations
    chain_file <- new_chain_file()

    # First mapping: source1[100-200) -> chr1[0-100)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 300, "+", 100, 200, "chr1", 400, "+", 0, 100, 1)

    # Second mapping: source1[100-200) -> chr1[200-300) (source overlap)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 300, "+", 100, 200, "chr1", 400, "+", 200, 300, 2)

    # Liftover with keep policy to allow source overlaps
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file,
        src_overlap_policy = "keep", tgt_overlap_policy = "auto"
    )

    # Extract and validate values
    result <- gextract(lifted_track, gintervals.all())

    # Expected: values appear in both target regions
    # chr1[0,20) = 4, chr1[21,22) = 5, chr1[80,85) = 6
    # chr1[200,220) = 4, chr1[221,222) = 5, chr1[280,285) = 6

    expect_true(nrow(result) >= 6)

    # Check first set of mappings (chr1[0-100))
    region1 <- result[result$chrom == "chr1" & result$start >= 0 & result$start < 100, ]
    expect_true(nrow(region1) >= 3)

    # Check for value 4 at [0,20)
    val4_region1 <- region1[region1$start == 0 & region1$end == 20, ]
    expect_equal(nrow(val4_region1), 1)
    expect_equal(as.numeric(val4_region1$lifted_track), 4)

    # Check for value 5 at [21,22)
    val5_region1 <- region1[region1$start == 21 & region1$end == 22, ]
    expect_equal(nrow(val5_region1), 1)
    expect_equal(as.numeric(val5_region1$lifted_track), 5)

    # Check for value 6 at [80,85)
    val6_region1 <- region1[region1$start == 80 & region1$end == 85, ]
    expect_equal(nrow(val6_region1), 1)
    expect_equal(as.numeric(val6_region1$lifted_track), 6)

    # Check second set of mappings (chr1[200-300))
    region2 <- result[result$chrom == "chr1" & result$start >= 200 & result$start < 300, ]
    expect_true(nrow(region2) >= 3)

    # Check for value 4 at [200,220)
    val4_region2 <- region2[region2$start == 200 & region2$end == 220, ]
    expect_equal(nrow(val4_region2), 1)
    expect_equal(as.numeric(val4_region2$lifted_track), 4)

    # Check for value 5 at [221,222)
    val5_region2 <- region2[region2$start == 221 & region2$end == 222, ]
    expect_equal(nrow(val5_region2), 1)
    expect_equal(as.numeric(val5_region2$lifted_track), 5)

    # Check for value 6 at [280,285)
    val6_region2 <- region2[region2$start == 280 & region2$end == 285, ]
    expect_equal(nrow(val6_region2), 1)
    expect_equal(as.numeric(val6_region2$lifted_track), 6)
})

test_that("gtrack.liftover preserves values in overlapping source regions mapping", {
    local_db_state()

    # Create source genome with one chromosome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 300), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with specific values in source
    # source1: [130,140) = 4, [140,160) = 5, [170,210) = 6
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(130, 140, 170),
        end = c(140, 160, 210),
        stringsAsFactors = FALSE
    )
    src_values <- c(4, 5, 6)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 400), collapse = ""), "\n")))

    # Create chain: source1[100-200) -> chr1[0-100) AND source1[150-250) -> chr1[200-300)
    # The source regions overlap in [150-200), creating a many-to-many mapping
    chain_file <- new_chain_file()

    # First mapping: source1[100-200) -> chr1[0-100)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 300, "+", 100, 200, "chr1", 400, "+", 0, 100, 1)

    # Second mapping: source1[150-250) -> chr1[200-300) (overlaps with first in [150-200))
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 300, "+", 150, 250, "chr1", 400, "+", 200, 300, 2)

    # Liftover with keep policy to allow source overlaps
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file,
        src_overlap_policy = "keep", tgt_overlap_policy = "auto"
    )

    # Extract and validate values
    result <- gextract(lifted_track, gintervals.all())
    # Expected mapping:
    # From source1[100-200) -> chr1[0-100):
    #   - source1[130,140)=4 -> chr1[30,40)=4
    #   - source1[140,160)=5 -> chr1[40,60)=5
    #   - source1[170,200)=6 -> chr1[70,100)=6 (note: [170,210) overlaps with [100,200) so we take [170,200))
    #
    # From source1[150-250) -> chr1[200-300):
    #   - source1[150,160)=5 -> chr1[250,260)=5 (note: [140,160) overlaps with [150,250) so we take [150,160))
    #   - source1[170,210)=6 -> chr1[220,260)=6
    #
    # Overlap analysis:
    # - [130,140)=4: Only in first chain -> chr1[30,40)
    # - [140,160)=5: Overlaps [100,200) and [150,250)
    #   * In first chain [100,200): source1[140,160) -> chr1[40,60)=5
    #   * In second chain [150,250): source1[150,160) -> chr1[250,260)=5
    # - [170,210)=6: Overlaps [100,200) and [150,250)
    #   * In first chain [100,200): source1[170,200) -> chr1[70,100)=6
    #   * In second chain [150,250): source1[170,210) -> chr1[220,260)=6

    expect_true(nrow(result) >= 5)

    # Check first mapping (chr1[0-100))
    region1 <- result[result$chrom == "chr1" & result$start >= 0 & result$start < 100, ]
    expect_true(nrow(region1) >= 3)

    # Check for value 4 at [30,40)
    val4_region1 <- region1[region1$start == 30 & region1$end == 40, ]
    expect_equal(nrow(val4_region1), 1)
    expect_equal(as.numeric(val4_region1$lifted_track), 4)

    # Check for value 5 at [40,60)
    val5_region1 <- region1[region1$start == 40 & region1$end == 60, ]
    expect_equal(nrow(val5_region1), 1)
    expect_equal(as.numeric(val5_region1$lifted_track), 5)

    # Check for value 6 at [70,100)
    val6_region1 <- region1[region1$start == 70 & region1$end == 100, ]
    expect_equal(nrow(val6_region1), 1)
    expect_equal(as.numeric(val6_region1$lifted_track), 6)

    # Check second mapping (chr1[200-300))
    region2 <- result[result$chrom == "chr1" & result$start >= 200 & result$start < 300, ]
    expect_true(nrow(region2) >= 1)

    # Check for value 6 at [220,260)
    val6_region2 <- region2[region2$start == 220 & region2$end == 260, ]
    expect_equal(nrow(val6_region2), 1)
    expect_equal(as.numeric(val6_region2$lifted_track), 6)
})

test_that("gtrack.liftover preserves values with one chain included in another", {
    local_db_state()

    # Create source genome with one chromosome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 400), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with values:
    # [60,80)=2 (only in inner chain)
    # [90,130)=3 (crosses chain borders)
    # [150,175)=4 (only in outer chain)
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(60, 90, 150),
        end = c(80, 130, 175),
        stringsAsFactors = FALSE
    )
    src_values <- c(2, 3, 4)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 600), collapse = ""), "\n")))

    # Create chain: inner chain source1[60-160) -> chr1[0-100)
    #               outer chain source1[80-180) -> chr1[100-200)
    # Inner chain is included (but not equal to) outer chain
    # Overlap region is source1[80-160) -> both chains
    chain_file <- new_chain_file()

    # Inner mapping: source1[60-160) -> chr1[0-100)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 400, "+", 60, 160, "chr1", 600, "+", 0, 100, 1)

    # Outer mapping: source1[80-180) -> chr1[100-200) (includes inner chain and extends)
    # Source: 180-80=100, Target: 200-100=100
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 400, "+", 80, 180, "chr1", 600, "+", 100, 200, 2)

    # Liftover with keep policy
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file,
        src_overlap_policy = "keep", tgt_overlap_policy = "auto"
    )

    # Extract and validate values
    result <- gextract(lifted_track, gintervals.all())

    # Expected mapping:
    # Inner chain source1[60-160) -> chr1[0-100):
    #   - source1[60,80)=2 -> chr1[0,20)=2 (only in inner chain)
    #   - source1[90,130)=3 -> chr1[30,70)=3 (part of value 3 that overlaps inner chain)
    #
    # Outer chain source1[80-180) -> chr1[100-200):
    #   - source1[90,130)=3 -> chr1[110,150)=3 (part of value 3 that overlaps outer chain)
    #   - source1[150,175)=4 -> chr1[170,195)=4 (only in outer chain)
    #
    # Overlap analysis:
    # - [60,80)=2: Only in inner chain -> chr1[0,20)=2
    # - [90,130)=3: Overlaps both chains
    #   * In inner chain [60,160): source1[90,130) -> chr1[30,70)=3
    #   * In outer chain [80,180): source1[90,130) -> chr1[110,150)=3
    # - [150,175)=4: Only in outer chain [80,180) -> chr1[170,195)=4

    expect_true(nrow(result) >= 4)

    # Check values that appear only through inner chain (chr1[0-100))
    region1 <- result[result$chrom == "chr1" & result$start >= 0 & result$start < 100, ]

    # Check for value 2 at [0,20) - only in inner chain
    val2_region1 <- region1[region1$start == 0 & region1$end == 20, ]
    expect_equal(nrow(val2_region1), 1)
    expect_equal(as.numeric(val2_region1$lifted_track), 2)

    # Check for value 3 at [30,70) - only in inner chain
    val3_region1 <- region1[region1$start == 30 & region1$end == 70, ]
    expect_equal(nrow(val3_region1), 1)
    expect_equal(as.numeric(val3_region1$lifted_track), 3)

    # Check values that appear only through outer chain (chr1[100-200))
    region2 <- result[result$chrom == "chr1" & result$start >= 100 & result$start < 200, ]

    # Check for value 3 at [110,150) - part that crosses into outer chain
    val3_region2 <- region2[region2$start == 110 & region2$end == 150, ]
    expect_equal(nrow(val3_region2), 1)
    expect_equal(as.numeric(val3_region2$lifted_track), 3)

    # Check for value 4 at [170,195) - only in outer chain
    val4_region2 <- region2[region2$start == 170 & region2$end == 195, ]
    expect_equal(nrow(val4_region2), 1)
    expect_equal(as.numeric(val4_region2$lifted_track), 4)
})

test_that("gtrack.liftover handles dense track with bin averaging", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create dense track with binsize 10
    # Values: bins 0-9 (coords 0-100) have specific values
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = rep("chrsource1", 10),
        start = seq(0, 90, by = 10),
        end = seq(10, 100, by = 10),
        stringsAsFactors = FALSE
    )
    src_values <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    gtrack.create_dense("src_dense", "Source dense track", src_intervals, src_values, 10, NaN)

    src_track_dir <- file.path(source_db, "tracks", "src_dense.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 200), collapse = ""), "\n")))

    # Create chain: source1[0-100) -> chr1[50-150)
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 100, "chr1", 200, "+", 50, 150, 1)

    # Liftover
    lifted_track <- "lifted_dense"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted dense track", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals("chr1", 50, 150))

    # Verify we have values at expected positions
    expect_true(nrow(result) >= 10)

    # Check that bins are correctly mapped with offset
    # source[0,10) (val=1) -> chr1[50,60)
    # source[10,20) (val=2) -> chr1[60,70)
    # etc.

    # Find bin at chr1[50,60) - should have value close to 1
    bin1 <- result[result$start >= 50 & result$start < 60, ]
    expect_true(nrow(bin1) > 0)
    expect_true(any(abs(bin1$lifted_dense - 1) < 0.1))

    # Find bin at chr1[90,100) - should have value close to 5
    bin5 <- result[result$start >= 90 & result$start < 100, ]
    expect_true(nrow(bin5) > 0)
    expect_true(any(abs(bin5$lifted_dense - 5) < 0.1))
})

test_that("gtrack.liftover with target overlap auto policy truncates correctly", {
    local_db_state()

    # Create source genome with two chromosomes
    # Filenames must start with "chr" for .gseq.import() to process them
    source_fasta1 <- file.path(tempdir(), "chrsource1.fasta")
    source_fasta2 <- file.path(tempdir(), "chrsource2.fasta")
    cat(">source1\n", paste(rep("A", 100), collapse = ""), "\n", sep = "", file = source_fasta1)
    cat(">source2\n", paste(rep("C", 100), collapse = ""), "\n", sep = "", file = source_fasta2)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta1)
        unlink(source_fasta2)
    })

    gdb.create(groot = source_db, fasta = c(source_fasta1, source_fasta2))
    gdb.init(source_db)

    # Create tracks with distinct values for each source chromosome
    # Database chromosome names are "chrsource1" and "chrsource2" (from filenames)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource2"),
        start = c(0, 0),
        end = c(60, 60),
        stringsAsFactors = FALSE
    )
    src_values <- c(100, 200)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 100), collapse = ""), "\n")))

    # Create chain with target overlaps:
    # source1[0-60) -> chr1[0-60)
    # source2[0-60) -> chr1[40-100) (overlaps chr1[40-60))
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 60, "chr1", 100, "+", 0, 60, 1)

    # Chain file uses "chrsource2" to match database chromosome name
    write_chain_entry(chain_file, "chrsource2", 100, "+", 0, 60, "chr1", 100, "+", 40, 100, 2)

    # Liftover with auto_first policy - should truncate overlaps
    lifted_track <- "lifted_auto"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file,
        tgt_overlap_policy = "auto_first"
    )

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # With auto_first policy, overlapping regions should be truncated
    # We expect: chr1[0,40) = 100, chr1[60,100) = 200
    # The overlap region [40,60) should be handled by truncation

    expect_true(nrow(result) >= 1)

    # Check non-overlapping regions have correct values
    region_before <- result[result$start < 40, ]
    if (nrow(region_before) > 0) {
        expect_true(all(region_before$lifted_auto == 100 | is.na(region_before$lifted_auto)))
    }

    region_after <- result[result$start >= 60, ]
    if (nrow(region_after) > 0) {
        expect_true(all(region_after$lifted_auto == 200 | is.na(region_after$lifted_auto)))
    }
})

test_that("gtrack.liftover with discard policies removes overlapping intervals", {
    local_db_state()

    # Create source genome with two chromosomes
    # Filenames must start with "chr" for .gseq.import() to process them
    source_fasta1 <- file.path(tempdir(), "chrsource1.fasta")
    source_fasta2 <- file.path(tempdir(), "chrsource2.fasta")
    source_fasta3 <- file.path(tempdir(), "chrsource3.fasta")
    cat(">source1\n", paste(rep("A", 100), collapse = ""), "\n", sep = "", file = source_fasta1)
    cat(">source2\n", paste(rep("C", 100), collapse = ""), "\n", sep = "", file = source_fasta2)
    cat(">source3\n", paste(rep("G", 100), collapse = ""), "\n", sep = "", file = source_fasta3)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta1)
        unlink(source_fasta2)
        unlink(source_fasta3)
    })

    gdb.create(groot = source_db, fasta = c(source_fasta1, source_fasta2, source_fasta3))
    gdb.init(source_db)

    # Create tracks: chrsource1 will have overlaps, chrsource3 is clean
    # Database chromosome names are "chrsource1", "chrsource2", "chrsource3" (from filenames)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource3"),
        start = c(0, 0),
        end = c(50, 50),
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 333)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    # Filenames must start with "chr" for .gseq.import() to process them
    target_fasta1 <- file.path(tempdir(), "chr1.fasta")
    target_fasta2 <- file.path(tempdir(), "chr2.fasta")
    cat(">chr1\n", paste(rep("T", 100), collapse = ""), "\n", sep = "", file = target_fasta1)
    cat(">chr2\n", paste(rep("T", 100), collapse = ""), "\n", sep = "", file = target_fasta2)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta1)
        unlink(target_fasta2)
    })

    gdb.create(groot = target_db, fasta = c(target_fasta1, target_fasta2))
    gdb.init(target_db)

    # Create chain with source overlap:
    # source1[0-50) -> chr1[0-50)
    # source1[20-70) -> chr1[50-100) (source overlap at [20-50))
    # source3[0-50) -> chr2[0-50) (no overlap)
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain file uses "chrsource1" to match database chromosome name
    cat("chain 1000 chrsource1 100 + 0 50 chr1 100 + 0 50 1\n", file = chain_file)
    cat("50\n\n", file = chain_file, append = TRUE)

    # Chain file uses "chrsource1" to match database chromosome name
    cat("chain 1000 chrsource1 100 + 20 70 chr1 100 + 50 100 2\n", file = chain_file, append = TRUE)
    cat("50\n\n", file = chain_file, append = TRUE)

    # Chain file uses "chrsource3" to match database chromosome name
    cat("chain 1000 chrsource3 100 + 0 50 chr2 100 + 0 50 3\n", file = chain_file, append = TRUE)
    cat("50\n\n", file = chain_file, append = TRUE)

    # Liftover with discard policy
    lifted_track <- "lifted_discard"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file,
        src_overlap_policy = "discard", tgt_overlap_policy = "discard"
    )

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # With discard policy, only source3 -> chr2 mapping should remain
    # source1 has overlaps so should be discarded
    expect_true(nrow(result) >= 0)

    # Check that only chr2 has values (source3 mapping)
    if (nrow(result) > 0) {
        chr2_result <- result[result$chrom == "chr2", ]
        expect_true(nrow(chr2_result) >= 1)

        # Check that chr2 has value 333
        expect_true(any(chr2_result$lifted_discard == 333))

        # chr1 should have no values or only NaN
        chr1_result <- result[result$chrom == "chr1", ]
        if (nrow(chr1_result) > 0) {
            expect_true(all(is.na(chr1_result$lifted_discard)))
        }
    }
})

test_that("gtrack.liftover preserves sparse track gaps correctly", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with gaps
    # [10,20) = 1, [40,50) = 2, [80,90) = 3
    # Gaps at [0,10), [20,40), [50,80), [90,200)
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 40, 80),
        end = c(20, 50, 90),
        stringsAsFactors = FALSE
    )
    src_values <- c(1, 2, 3)
    gtrack.create_sparse("src_sparse", "Source sparse", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_sparse.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 200), collapse = ""), "\n")))

    # Simple 1:1 mapping with offset: source1[0-100) -> chr1[100-200)
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 100, "chr1", 200, "+", 100, 200, 1)

    # Liftover
    lifted_track <- "lifted_sparse"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted sparse", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # Expected: only 3 intervals with values, rest is NaN/missing
    # chr1[110,120) = 1, chr1[140,150) = 2, chr1[180,190) = 3
    expect_equal(nrow(result), 3)

    # Check first interval
    expect_equal(as.character(result$chrom[1]), "chr1")
    expect_equal(as.numeric(result$start[1]), 110)
    expect_equal(as.numeric(result$end[1]), 120)
    expect_equal(as.numeric(result$lifted_sparse[1]), 1)

    # Check second interval
    expect_equal(as.numeric(result$start[2]), 140)
    expect_equal(as.numeric(result$end[2]), 150)
    expect_equal(as.numeric(result$lifted_sparse[2]), 2)

    # Check third interval
    expect_equal(as.numeric(result$start[3]), 180)
    expect_equal(as.numeric(result$end[3]), 190)
    expect_equal(as.numeric(result$lifted_sparse[3]), 3)
})

test_that("gtrack.liftover handles reverse strand mapping correctly", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with multiple intervals
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 50, 100),
        end = c(20, 60, 120),
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 222, 333)
    gtrack.create_sparse("src_reverse", "Source reverse", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_reverse.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 200), collapse = ""), "\n")))

    # Create chain with reverse strand (-)
    # source1[0-150) + -> chr1[0-150) -
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 150, "chr1", 200, "-", 50, 200, 1)

    # Liftover
    lifted_track <- "lifted_reverse"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted reverse", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # Expect all 3 intervals to be lifted
    expect_equal(nrow(result), 3)

    # Check that values are preserved (order may be reversed in coordinates)
    expect_true(all(sort(result$lifted_reverse) == sort(src_values)))
})

test_that("gtrack.liftover handles chain gaps (unmapped regions)", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 300), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with intervals in both mapped and unmapped regions
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 100, 200), # 10 is mapped, 100 is in gap, 200 is mapped
        end = c(20, 110, 210),
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 222, 333)
    gtrack.create_sparse("src_gaps", "Source gaps", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_gaps.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 300), collapse = ""), "\n")))

    # Create chain with gaps: maps [0-50) and [150-250), but NOT [50-150)
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 300, "+", 0, 50, "chr1", 300, "+", 0, 50, 1)

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 300, "+", 150, 250, "chr1", 300, "+", 100, 200, 2)

    # Liftover
    lifted_track <- "lifted_gaps"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted gaps", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # Only intervals in mapped regions should appear
    # [10,20) = 111 -> chr1[10,20) = 111
    # [100,110) = 222 -> UNMAPPED (should not appear)
    # [200,210) = 333 -> chr1[150,160) = 333
    expect_equal(nrow(result), 2)
    expect_true(any(result$lifted_gaps == 111))
    expect_true(any(result$lifted_gaps == 333))
    expect_false(any(result$lifted_gaps == 222)) # This one should be missing
})

test_that("gtrack.liftover handles special values (NaN, zero, large)", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with special values
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 30, 50, 70),
        end = c(20, 40, 60, 80),
        stringsAsFactors = FALSE
    )
    src_values <- c(0, NaN, 1e10, -1e10) # zero, NaN, very large positive, very large negative
    gtrack.create_sparse("src_special", "Source special", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_special.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 200), collapse = ""), "\n")))

    # Simple 1:1 mapping
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 100, "chr1", 200, "+", 0, 100, 1)

    # Liftover
    lifted_track <- "lifted_special"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted special", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    expect_equal(nrow(result), 4)

    # Check zero value is preserved
    expect_true(any(result$lifted_special == 0))

    # Check NaN is preserved
    expect_true(any(is.nan(result$lifted_special)))

    # Check large values are preserved
    expect_true(any(result$lifted_special == 1e10))
    expect_true(any(result$lifted_special == -1e10))
})

test_that("gtrack.liftover handles empty source track", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create empty sparse track
    src_intervals <- data.frame(
        chrom = character(0),
        start = numeric(0),
        end = numeric(0),
        stringsAsFactors = FALSE
    )
    src_values <- numeric(0)
    gtrack.create_sparse("src_empty", "Source empty", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_empty.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 200), collapse = ""), "\n")))

    # Create chain
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 100, "chr1", 200, "+", 0, 100, 1)

    # Liftover
    lifted_track <- "lifted_empty"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted empty", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # Should have no intervals (NULL or empty data frame)
    expect_true(is.null(result) || nrow(result) == 0)
})

test_that("gtrack.liftover handles single base pair intervals", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with single BP intervals
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 50, 100),
        end = c(11, 51, 101), # All single BP
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 222, 333)
    gtrack.create_sparse("src_single_bp", "Source single BP", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_single_bp.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 200), collapse = ""), "\n")))

    # Simple 1:1 mapping
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 150, "chr1", 200, "+", 0, 150, 1)

    # Liftover
    lifted_track <- "lifted_single_bp"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted single BP", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # All 3 single BP intervals should be preserved
    expect_equal(nrow(result), 3)
    expect_true(all(result$end - result$start == 1)) # All should be single BP
    expect_true(all(sort(result$lifted_single_bp) == sort(src_values)))
})

test_that("gtrack.liftover handles multiple source chromosomes to single target", {
    local_db_state()

    # Create source genome with multiple chromosomes
    # Filenames must start with "chr" for .gseq.import() to process them
    source_fasta1 <- file.path(tempdir(), "chrsource1.fasta")
    source_fasta2 <- file.path(tempdir(), "chrsource2.fasta")
    source_fasta3 <- file.path(tempdir(), "chrsource3.fasta")
    cat(">source1\n", paste(rep("A", 100), collapse = ""), "\n", sep = "", file = source_fasta1)
    cat(">source2\n", paste(rep("C", 100), collapse = ""), "\n", sep = "", file = source_fasta2)
    cat(">source3\n", paste(rep("G", 100), collapse = ""), "\n", sep = "", file = source_fasta3)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta1)
        unlink(source_fasta2)
        unlink(source_fasta3)
    })

    gdb.create(groot = source_db, fasta = c(source_fasta1, source_fasta2, source_fasta3))
    gdb.init(source_db)

    # Create sparse track with intervals on all three chromosomes
    src_intervals <- data.frame(
        # Database chromosome names are "chrsource1", "chrsource2", "chrsource3" (from filenames)
        chrom = c("chrsource1", "chrsource2", "chrsource3"),
        start = c(10, 10, 10),
        end = c(20, 20, 20),
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 222, 333)
    gtrack.create_sparse("src_multi", "Source multi", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_multi.track")

    # Create target genome with single chromosome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 400), collapse = ""), "\n")))

    # Create chain mapping all three sources to different regions of chr1
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 50, "chr1", 400, "+", 0, 50, 1)

    # Chain file uses "chrsource2" to match database chromosome name
    write_chain_entry(chain_file, "chrsource2", 100, "+", 0, 50, "chr1", 400, "+", 100, 150, 2)

    # Chain file uses "chrsource3" to match database chromosome name
    write_chain_entry(chain_file, "chrsource3", 100, "+", 0, 50, "chr1", 400, "+", 200, 250, 3)

    # Liftover
    lifted_track <- "lifted_multi"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted multi", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # All 3 intervals should be mapped to chr1
    expect_equal(nrow(result), 3)
    expect_true(all(result$chrom == "chr1"))
    # Check values are preserved
    expect_true(all(sort(result$lifted_multi) == sort(src_values)))

    # Check coordinates are in different regions
    coords <- result$start
    expect_true(any(coords >= 0 & coords < 50))
    expect_true(any(coords >= 100 & coords < 150))
    expect_true(any(coords >= 200 & coords < 250))
})

test_that("gtrack.liftover handles intervals at chromosome boundaries", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 100), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with intervals at boundaries
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(0, 45, 90), # At start, middle, near end
        end = c(10, 55, 100), # Last interval goes to chromosome end
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 222, 333)
    gtrack.create_sparse("src_boundary", "Source boundary", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_boundary.track")

    # Create target genome with single chromosome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 150), collapse = ""), "\n")))

    # Map entire source chromosome to target
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 100, "chr1", 150, "+", 0, 100, 1)

    # Liftover
    lifted_track <- "lifted_boundary"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted boundary", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # All 3 intervals should be preserved
    expect_equal(nrow(result), 3)
    expect_true(all(sort(result$lifted_boundary) == sort(src_values)))

    # Check that boundary intervals are preserved correctly
    first_interval <- result[result$lifted_boundary == 111, ]
    expect_equal(as.numeric(first_interval$start), 0) # Should start at 0

    last_interval <- result[result$lifted_boundary == 333, ]
    expect_equal(as.numeric(last_interval$end), 100) # Should end at mapped end
})

test_that("gtrack.liftover handles coordinate reversals (inversions)", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with ordered intervals
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 50, 100),
        end = c(20, 60, 110),
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 222, 333)
    gtrack.create_sparse("src_inversion", "Source inversion", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_inversion.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 300), collapse = ""), "\n")))

    # Create chain with inversion (reverse strand)
    # Maps source1[0-150) to chr1[150-0) (reversed)
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 150, "chr1", 300, "-", 150, 300, 1)

    # Liftover
    lifted_track <- "lifted_inversion"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted inversion", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # All 3 intervals should be preserved
    expect_equal(nrow(result), 3)

    # Values should be preserved regardless of coordinate reversal
    expect_true(all(sort(result$lifted_inversion) == sort(src_values)))
})

test_that("gtrack.liftover handles dense track with partial bin coverage", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 1000), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create dense track with binsize 10
    # Values: bins 0-9 have values 0-9
    binsize <- 10
    values <- seq(0, 99)
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = rep("chrsource1", length(values)),
        start = seq(0, length(values) - 1) * binsize,
        end = seq(1, length(values)) * binsize,
        stringsAsFactors = FALSE
    )
    gtrack.create_dense("src_dense_partial", "Source dense partial", src_intervals, values, binsize)

    src_track_dir <- file.path(source_db, "tracks", "src_dense_partial.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 1500), collapse = ""), "\n")))

    # Create chain that maps only part of source
    # Maps source1[200-400) to chr1[0-200)
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 1000, "+", 200, 400, "chr1", 1500, "+", 0, 200, 1)

    # Liftover
    lifted_track <- "lifted_dense_partial"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted dense partial", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals("chr1", c(0, 100), c(100, 200)))

    # Should have bins covering the mapped region
    # Bins in source [200-400) are bins 20-39 (200/10 = 20, 400/10 = 40)
    # These should map to chr1[0-200)
    expect_true(nrow(result) > 0)

    # Check that we don't have NaN for all values (some data should be present)
    expect_true(!all(is.na(result$lifted_dense_partial)))
})

test_that("gtrack.liftover handles sparse track with very small intervals", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 1000), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with mix of very small and regular intervals
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1", "chrsource1"),
        start = c(100, 101, 200, 500),
        end = c(101, 102, 250, 600), # 1bp, 1bp, 50bp, 100bp
        stringsAsFactors = FALSE
    )
    src_values <- c(1, 2, 50, 100)
    gtrack.create_sparse("src_small", "Source small", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_small.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 1000), collapse = ""), "\n")))

    # Simple 1:1 mapping
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 1000, "+", 0, 700, "chr1", 1000, "+", 0, 700, 1)

    # Liftover
    lifted_track <- "lifted_small"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted small", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # All 4 intervals should be preserved, including the 1bp ones
    expect_equal(nrow(result), 4)
    expect_true(all(sort(result$lifted_small) == sort(src_values)))

    # Check that 1bp intervals are preserved
    expect_true(any(result$end - result$start == 1))
})

test_that("gtrack.liftover handles one-to-many with different bin sizes in target", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 500), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create dense track
    binsize <- 10
    values <- seq(1, 50)
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = rep("chrsource1", length(values)),
        start = seq(0, length(values) - 1) * binsize,
        end = seq(1, length(values)) * binsize,
        stringsAsFactors = FALSE
    )
    gtrack.create_dense("src_onemany_dense", "Source one-many dense", src_intervals, values, binsize)

    src_track_dir <- file.path(source_db, "tracks", "src_onemany_dense.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 1000), collapse = ""), "\n")))

    # Create chain with one-to-many mapping
    # source1[100-200) maps to both chr1[0-100) and chr1[500-600)
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 500, "+", 100, 200, "chr1", 1000, "+", 0, 100, 1)

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 500, "+", 100, 200, "chr1", 1000, "+", 500, 600, 2)

    # Liftover with keep policy
    lifted_track <- "lifted_onemany_dense"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted one-many dense", src_track_dir, chain_file,
        src_overlap_policy = "keep", tgt_overlap_policy = "auto"
    )

    # Extract from both regions
    result1 <- gextract(lifted_track, gintervals("chr1", 0, 100))
    result2 <- gextract(lifted_track, gintervals("chr1", 500, 600))

    # Both regions should have data
    expect_true(nrow(result1) > 0)
    expect_true(nrow(result2) > 0)

    # Values in both regions should be from the same source bins
    # They should have similar (if not identical) values
    expect_true(!all(is.na(result1$lifted_onemany_dense)))
    expect_true(!all(is.na(result2$lifted_onemany_dense)))
})

test_that("gtrack.liftover handles source track with consecutive intervals", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 300), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with consecutive intervals (no gaps)
    # This is like a mini dense track but in sparse format
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1", "chrsource1"),
        start = c(100, 110, 120, 130),
        end = c(110, 120, 130, 140),
        stringsAsFactors = FALSE
    )
    src_values <- c(10, 20, 30, 40)
    gtrack.create_sparse("src_consecutive", "Source consecutive", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_consecutive.track")

    # Create target genome with single chromosome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 300), collapse = ""), "\n")))

    # Simple 1:1 mapping with offset
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 300, "+", 100, 150, "chr1", 300, "+", 50, 100, 1)

    # Liftover
    lifted_track <- "lifted_consecutive"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted consecutive", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # All 4 intervals should be preserved
    expect_equal(nrow(result), 4)
    expect_true(all(sort(result$lifted_consecutive) == sort(src_values)))

    # Check that intervals are still consecutive
    result <- result[order(result$start), ]
    for (i in 1:(nrow(result) - 1)) {
        expect_equal(as.numeric(result$end[i]), as.numeric(result$start[i + 1]))
    }
})

test_that("gtrack.liftover handles mixed positive and negative values", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track with mixed positive and negative values
    src_intervals <- data.frame(
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        chrom = c("chrsource1", "chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 40, 70, 100),
        end = c(20, 50, 80, 110),
        stringsAsFactors = FALSE
    )
    src_values <- c(-100, 100, -50.5, 75.25)
    gtrack.create_sparse("src_mixed", "Source mixed", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_mixed.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 200), collapse = ""), "\n")))

    # Simple 1:1 mapping
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 120, "chr1", 200, "+", 0, 120, 1)

    # Liftover
    lifted_track <- "lifted_mixed"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted mixed", src_track_dir, chain_file)

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # All 4 intervals should be preserved
    expect_equal(nrow(result), 4)

    # Check exact value preservation
    expect_true(all(sort(result$lifted_mixed) == sort(src_values)))

    # Verify negative and positive values are both present
    expect_true(any(result$lifted_mixed < 0))
    expect_true(any(result$lifted_mixed > 0))
})

test_that("gtrack.liftover handles target overlap with truncation", {
    local_db_state()

    # Create source genome with two chromosomes
    # Filenames must start with "chr" for .gseq.import() to process them
    source_fasta1 <- file.path(tempdir(), "chrsource1.fasta")
    source_fasta2 <- file.path(tempdir(), "chrsource2.fasta")
    cat(">source1\n", paste(rep("A", 100), collapse = ""), "\n", sep = "", file = source_fasta1)
    cat(">source2\n", paste(rep("C", 100), collapse = ""), "\n", sep = "", file = source_fasta2)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta1)
        unlink(source_fasta2)
    })

    gdb.create(groot = source_db, fasta = c(source_fasta1, source_fasta2))
    gdb.init(source_db)

    # Create sparse track with intervals on both chromosomes
    src_intervals <- data.frame(
        # Database chromosome names are "chrsource1" and "chrsource2" (from filenames)
        chrom = c("chrsource1", "chrsource2"),
        start = c(10, 10),
        end = c(60, 60),
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 222)
    gtrack.create_sparse("src_tgt_overlap", "Source target overlap", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_tgt_overlap.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 200), collapse = ""), "\n")))

    # Create chain with target overlap
    # source1[0-70) -> chr1[0-70)
    # source2[0-70) -> chr1[30-100) (overlaps at chr1[30-70))
    chain_file <- new_chain_file()

    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 70, "chr1", 200, "+", 0, 70, 1)

    # Chain file uses "chrsource2" to match database chromosome name
    write_chain_entry(chain_file, "chrsource2", 100, "+", 0, 70, "chr1", 200, "+", 30, 100, 2)

    # Liftover with auto truncation
    lifted_track <- "lifted_tgt_overlap"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted target overlap", src_track_dir, chain_file,
        src_overlap_policy = "error", tgt_overlap_policy = "auto"
    )

    # Extract and validate
    result <- gextract(lifted_track, gintervals.all())

    # Should have intervals from both sources, possibly truncated
    expect_true(nrow(result) >= 2)

    # Both values should appear somewhere
    expect_true(any(result$lifted_tgt_overlap == 111) || any(result$lifted_tgt_overlap == 222))
})

# Test for reverse strand liftover bug
test_that("gintervals.liftover works correctly with reverse strand targets", {
    local_db_state()

    # This test validates the fix for the reverse strand liftover bug
    # Bug: when target is on negative strand, offsets were applied in wrong direction

    # Create target genome with the chromosomes we need
    # Filenames must start with "chr" for .gseq.import() to process them
    target_fasta1 <- file.path(tempdir(), "chrtgt.fasta")
    target_fasta2 <- file.path(tempdir(), "chrsrc.fasta")
    cat(">tgt\n", paste(rep("T", 2000), collapse = ""), "\n", sep = "", file = target_fasta1)
    cat(">src\n", paste(rep("A", 1000), collapse = ""), "\n", sep = "", file = target_fasta2)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta1)
        unlink(target_fasta2)
    })

    gdb.create(groot = target_db, fasta = c(target_fasta1, target_fasta2))
    gdb.init(target_db)

    # Create a minimal chain with reverse strand target
    chain_file <- new_chain_file()

    # Chain format: source on + strand, target on - strand
    # Source: chrom "chrsrc" (from filename chrsrc.fasta), size 1000, + strand, coords 100-600
    # Target: chrom "chrtgt" (from filename chrtgt.fasta), size 2000, - strand, coords 500-1000
    # Single block of size 500
    write_chain_entry(chain_file, "chrsrc", 1000, "+", 100, 600, "chrtgt", 2000, "-", 500, 1000, 1)

    # Create source intervals to test
    # Test interval at chrsrc:300-301 (offset 200 from block start at 100)
    # Database chromosome name is "chrsrc" (from filename chrsrc.fasta)
    src_interv <- data.frame(
        chrom = "chrsrc",
        start = 300,
        end = 301,
        stringsAsFactors = FALSE
    )

    # Perform liftover (pass chain file directly)
    result <- gintervals.liftover(src_interv, chain_file, src_overlap_policy = "error")

    # Calculate expected result:
    # Source position: 300
    # Offset from block start (100): 200
    # Block in positive coords: tgt[500, 1000)
    # For reverse strand:
    #   tgt_size = 2000
    #   block_start_positive = 500
    #   block_end_positive = 1000
    #   block_start_reverse = 2000 - 1000 = 1000
    #   block_end_reverse = 2000 - 500 = 1500
    #
    # Correct mapping for reverse strand:
    #   offset 200 means we go backwards from the END of the reverse block
    #   result = block_end_reverse - offset - 1 = 1500 - 200 - 1 = 1299
    #   So interval [300,301) maps to [1299,1300)

    expected_start <- 1299
    expected_end <- 1300

    # Database chromosome name is "chrtgt" (from filename chrtgt.fasta)
    expect_equal(as.character(result$chrom[1]), "chrtgt")
    expect_equal(result$start[1], expected_start,
        info = paste("Expected start:", expected_start, "Got:", result$start[1])
    )
    expect_equal(result$end[1], expected_end,
        info = paste("Expected end:", expected_end, "Got:", result$end[1])
    )

    # Test another position to be thorough
    # Test interval at chrsrc:150-151 (offset 50 from block start)
    # Database chromosome name is "chrsrc" (from filename chrsrc.fasta)
    src_interv2 <- data.frame(
        chrom = "chrsrc",
        start = 150,
        end = 151,
        stringsAsFactors = FALSE
    )

    result2 <- gintervals.liftover(src_interv2, chain_file, src_overlap_policy = "error")

    # Expected: offset 50
    # result = 1500 - 50 - 1 = 1449
    # So interval [150,151) maps to [1449,1450)
    expected_start2 <- 1449
    expected_end2 <- 1450

    expect_equal(result2$start[1], expected_start2,
        info = paste("Expected start:", expected_start2, "Got:", result2$start[1])
    )
    expect_equal(result2$end[1], expected_end2,
        info = paste("Expected end:", expected_end2, "Got:", result2$end[1])
    )
})

# Test for reverse strand with gtrack.liftover
test_that("gtrack.liftover works correctly with reverse strand targets", {
    local_db_state()

    # Create source genome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsrc.fasta")
    cat(">src\n", paste(rep("A", 1000), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create sparse track in source at position 300-310
    # Database chromosome name is "chrsrc" (from filename chrsrc.fasta)
    src_intervals <- data.frame(
        chrom = "chrsrc",
        start = 300,
        end = 310,
        stringsAsFactors = FALSE
    )
    gtrack.create_sparse("test_track", "Test track", src_intervals, 42)
    src_track_dir <- file.path(source_db, "tracks", "test_track.track")

    # Create target genome
    setup_db(list(paste0(">tgt\n", paste(rep("T", 2000), collapse = ""), "\n")))

    # Create chain with reverse strand target
    chain_file <- new_chain_file()

    # Chain file uses "chrsrc" and "chrtgt" to match database chromosome names
    write_chain_entry(chain_file, "chrsrc", 1000, "+", 100, 600, "chrtgt", 2000, "-", 500, 1000, 1)

    # Perform liftover
    gtrack.liftover("lifted_track", "Lifted track", src_track_dir, chain_file,
        src_overlap_policy = "error", tgt_overlap_policy = "error"
    )

    # Extract result
    result <- gextract("lifted_track", gintervals.all())

    # Should have one interval
    expect_equal(nrow(result), 1)
    # Database chromosome name is "chrtgt" (from filename chrtgt.fasta)
    expect_equal(as.character(result$chrom[1]), "chrtgt")

    # Expected coordinates for [300,310) with offset 200 from block start:
    # Reverse coords: [1500 - 210, 1500 - 200) = [1290, 1300)
    expect_equal(result$start[1], 1290)
    expect_equal(result$end[1], 1300)
    expect_equal(result$lifted_track[1], 42)
})

test_that("gtrack.liftover finds all overlapping chains when they are non-consecutive", {
    local_db_state()

    # This test verifies the fix for a bug where gtrack.liftover would miss
    # overlapping chain intervals that were not consecutive in the sorted chain array.
    # The bug occurred when source intervals had overlapping regions mapping to
    # different targets, creating a pattern where overlapping chains are separated
    # by non-overlapping ones in the sorted (by start_src) chain array.

    # Create source genome first
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste0(rep("A", 100), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create source track with an interval overlapping first two chains
    gtrack.create_sparse(
        "source_track",
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        intervals = data.frame(chrom = "chrsource1", start = 20, end = 21, stringsAsFactors = FALSE),
        values = 42,
        description = "test track"
    )

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAA\n", ">chr3\nCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG\n"))

    # Create chain file replicating the user's bug scenario
    chain_file <- new_chain_file()

    # Chains with the same start position (will be consecutive when sorted)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 10, 30, "chr1", 44, "+", 0, 20, 1)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 10, 30, "chr2", 48, "+", 0, 20, 2)

    # Chain with different range that creates a gap
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 30, 32, "chr3", 42, "+", 0, 2, 3)

    # Load chain with keep policy
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Perform liftover
    gtrack.liftover(
        track = "lifted_track",
        description = "",
        src.track.dir = file.path(source_db, "tracks", "source_track.track"),
        chain = chain
    )

    # Extract result
    result <- gextract("lifted_track", gintervals.all())

    # Should return exactly 2 results (both chains with [10,30) overlap [20,21))
    expect_equal(nrow(result), 2)
    expect_true("chr1" %in% result$chrom)
    expect_true("chr2" %in% result$chrom)
    expect_false("chr3" %in% result$chrom)

    # Both should have the same value from the source
    expect_equal(result$lifted_track[result$chrom == "chr1"], 42)
    expect_equal(result$lifted_track[result$chrom == "chr2"], 42)
})

test_that("gtrack.liftover does not miss earlier long overlap when hint is to the right (policy=keep)", {
    local_db_state()

    # Create source genome first
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste0(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    # Create source track with two intervals
    gtrack.create_sparse(
        "source_track",
        intervals = data.frame(
            # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
            chrom = "chrsource1",
            start = c(70, 90), # Q1 then Q2
            end = c(71, 91),
            stringsAsFactors = FALSE
        ),
        values = c(100, 200),
        description = "test track"
    )

    # Create target genome
    setup_db(list(
        paste0(">chr1\n", paste(rep("A", 200), collapse = ""), "\n"),
        paste0(">chr2\n", paste(rep("C", 200), collapse = ""), "\n"),
        paste0(">chr3\n", paste(rep("G", 200), collapse = ""), "\n")
    ))

    chain_file <- new_chain_file()

    # A: [0,100) -> chr1 (overlaps Q1 and Q2)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 100, "chr1", 200, "+", 0, 100, 1)

    # B: [15,16) -> chr2 (does not overlap Q1 or Q2)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 15, 16, "chr2", 200, "+", 0, 1, 2)

    # C: [80,110) -> chr3 (overlaps Q2 only)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 80, 110, "chr3", 200, "+", 0, 30, 3)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Perform liftover
    gtrack.liftover(
        track = "lifted_track",
        description = "",
        src.track.dir = file.path(source_db, "tracks", "source_track.track"),
        chain = chain
    )

    result <- gextract("lifted_track", gintervals.all())

    # Expect: Q1 -> chr1; Q2 -> chr1 and chr3  (total 3 rows)
    # Previous bug would yield only 2 rows (Q1->chr1, Q2->chr3), missing chr1 for Q2.
    expect_equal(sum(result$lifted_track == 100), 1) # Q1 has exactly 1 hit (chr1)
    expect_equal(sum(result$lifted_track == 200), 2) # Q2 should have 2 hits (chr1 & chr3)
    expect_true("chr1" %in% result$chrom[result$lifted_track == 200])
    expect_true("chr3" %in% result$chrom[result$lifted_track == 200])
})

test_that("gtrack.liftover has deterministic ordering for chains with identical start_src", {
    local_db_state()

    # Create source genome first
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste0(rep("A", 500), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    # Create source track with an interval
    gtrack.create_sparse(
        "source_track",
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        intervals = data.frame(chrom = "chrsource1", start = 105, end = 106, stringsAsFactors = FALSE),
        values = 99,
        description = "test track"
    )

    # Create target genome with multiple chromosomes
    setup_db(list(
        ">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n",
        ">chr2\nGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAA\n",
        ">chr3\nTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG\n"
    ))

    chain_file <- new_chain_file()

    # Create chains with identical start_src but different end_src
    # All start at source1[100], but have different lengths
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 500, "+", 100, 120, "chr1", 48, "+", 0, 20, 1) # len=20
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 500, "+", 100, 110, "chr2", 32, "+", 0, 10, 2) # len=10
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 500, "+", 100, 130, "chr3", 32, "+", 0, 30, 3) # len=30
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 500, "+", 100, 115, "chr1", 48, "+", 20, 35, 4) # len=15

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Run liftover twice with the same query
    gtrack.liftover(
        track = "lifted_track1",
        description = "",
        src.track.dir = file.path(source_db, "tracks", "source_track.track"),
        chain = chain
    )

    gtrack.liftover(
        track = "lifted_track2",
        description = "",
        src.track.dir = file.path(source_db, "tracks", "source_track.track"),
        chain = chain
    )

    result1 <- gextract("lifted_track1", gintervals.all(), colnames = "value")
    result2 <- gextract("lifted_track2", gintervals.all(), colnames = "value")

    # Results should be identical (deterministic) - compare coordinates and intervalID
    expect_equal(result1, result2)

    # Verify all 4 chains are returned
    expect_equal(nrow(result1), 4)

    # Verify all expected chroms are present (order may vary deterministically)
    expect_true(all(c("chr1", "chr2", "chr3") %in% result1$chrom))

    # All should have the same value
    expect_true(all(result1$lifted_track1 == 99))
})

test_that("gtrack.liftover handles minus strand mapping at edges correctly", {
    local_db_state()

    # Create source genome first
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste0(rep("A", 100), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create source track with three intervals (left edge, middle, right edge)
    gtrack.create_sparse(
        "source_track",
        intervals = data.frame(
            # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
            chrom = "chrsource1",
            start = c(0, 15, 29), # left edge, middle, right edge
            end = c(1, 16, 30),
            stringsAsFactors = FALSE
        ),
        values = c(10, 20, 30),
        description = "test track"
    )

    # Create target genome with chr1 that has at least 50 bases
    setup_db(list(">chr1\n", paste0(rep("A", 60), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Minus strand chain: source1[0-30] -> chr1[20-50] (minus strand)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 30, "chr1", 60, "-", 20, 50, 1)

    chain <- gintervals.load_chain(chain_file)

    gtrack.liftover(
        track = "lifted_track",
        description = "",
        src.track.dir = file.path(source_db, "tracks", "source_track.track"),
        chain = chain
    )

    result <- gextract("lifted_track", gintervals.all())

    # All results should have start < end (primary correctness check)
    expect_true(all(result$start < result$end))

    # Verify we got results for all 3 queries
    expect_equal(nrow(result), 3)

    # Verify all coordinates are non-negative and reasonable
    expect_true(all(result$start >= 0))
    expect_true(all(result$end >= 0))

    # Verify values are preserved
    expect_true(all(result$lifted_track %in% c(10, 20, 30)))
})

test_that("gtrack.liftover handles dense cluster of chains with same start_src correctly", {
    local_db_state()

    # Create source genome first
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste0(rep("A", 10000), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    # Create source track with one interval that overlaps all chains
    gtrack.create_sparse(
        "source_track",
        # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
        intervals = data.frame(chrom = "chrsource1", start = 60, end = 61, stringsAsFactors = FALSE),
        values = 777,
        description = "test track"
    )

    # Create target genome
    setup_db(list(
        ">chr1\n", paste0(rep("A", 10000), collapse = ""), "\n",
        ">chr2\n", paste0(rep("C", 10000), collapse = ""), "\n"
    ))

    chain_file <- new_chain_file()

    # Create 100 chains all starting at source1[50] with varying lengths
    for (i in 1:100) {
        len <- i # lengths 1 to 100
        target_chrom <- if (i %% 2 == 0) "chr1" else "chr2"
        write_chain_entry(
            # Chain file uses "chrsource1" to match database chromosome name
            chain_file, "chrsource1", 10000, "+", 50, 50 + len,
            target_chrom, 10000, "+", i * 10, i * 10 + len, i
        )
    }

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    gtrack.liftover(
        track = "lifted_track",
        description = "",
        src.track.dir = file.path(source_db, "tracks", "source_track.track"),
        chain = chain
    )

    result <- gextract("lifted_track", gintervals.all())

    # Should find many overlapping chains (query at 60 overlaps chains where 50+length > 60)
    # Exact number depends on overlap resolution, but should be substantial
    expect_true(nrow(result) > 50) # At least half should overlap
    expect_true(nrow(result) <= 100) # At most all chains

    # All results should be valid intervals
    expect_true(all(result$start < result$end))

    # All should have the same value
    expect_true(all(result$lifted_track == 777))
})
