load_test_db()

test_that("gintervals.liftover matches liftOver binary - basic intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("A", 100), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Create chain: source1[10-50] -> chr1[5-45]
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 100 + 10 50 chr1 100 + 5 45 1\n", file = chain_file)
    cat("40\n\n", file = chain_file, append = TRUE)

    # Test intervals
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1"),
        start = c(12, 20, 35),
        end = c(18, 30, 45),
        stringsAsFactors = FALSE
    )

    # Create BED file for liftOver binary
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
            chrom = src_intervals$chrom,
            start = src_intervals$start,
            end = src_intervals$end,
            name = paste0("int", 1:3),
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    binary_result <- binary_result[order(binary_result$V2), c(1, 2, 3)]
    colnames(binary_result) <- c("chrom", "start", "end")

    # Run gintervals.liftover
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$start), c("chrom", "start", "end")]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - multiple chromosomes", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome with multiple chromosomes
    target_fasta1 <- file.path(tempdir(), "chr1.fasta")
    target_fasta2 <- file.path(tempdir(), "chr2.fasta")
    cat(">chr1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = target_fasta1)
    cat(">chr2\n", paste(rep("C", 200), collapse = ""), "\n", sep = "", file = target_fasta2)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta1)
        unlink(target_fasta2)
    })

    gdb.create(groot = target_db, fasta = c(target_fasta1, target_fasta2))
    gdb.init(target_db)

    # Create chain mapping both source chromosomes
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # source1[0-60] -> chr1[10-70]
    cat("chain 1000 source1 100 + 0 60 chr1 200 + 10 70 1\n", file = chain_file)
    cat("60\n\n", file = chain_file, append = TRUE)

    # source2[10-70] -> chr2[20-80]
    cat("chain 1000 source2 100 + 10 70 chr2 200 + 20 80 2\n", file = chain_file, append = TRUE)
    cat("60\n\n", file = chain_file, append = TRUE)

    # Test intervals on both chromosomes
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source2", "source2"),
        start = c(10, 40, 15, 50),
        end = c(20, 50, 25, 60),
        stringsAsFactors = FALSE
    )

    # Create BED file
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
            chrom = src_intervals$chrom,
            start = src_intervals$start,
            end = src_intervals$end,
            name = paste0("int", 1:4),
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    binary_result <- binary_result[order(binary_result$V1, binary_result$V2), c(1, 2, 3)]
    colnames(binary_result) <- c("chrom", "start", "end")

    # Run gintervals.liftover
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$chrom, misha_result$start), c("chrom", "start", "end")]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - reverse strand", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Create chain with reverse strand
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # source1[0-150] + -> chr1[50-200] - (reverse strand)
    cat("chain 1000 source1 200 + 0 150 chr1 200 - 0 150 1\n", file = chain_file)
    cat("150\n\n", file = chain_file, append = TRUE)

    # Test intervals
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1"),
        start = c(10, 50, 100),
        end = c(20, 60, 120),
        stringsAsFactors = FALSE
    )

    # Create BED file
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
            chrom = src_intervals$chrom,
            start = src_intervals$start,
            end = src_intervals$end,
            name = paste0("int", 1:3),
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    binary_result <- binary_result[order(binary_result$V2), c(1, 2, 3)]
    colnames(binary_result) <- c("chrom", "start", "end")

    # Run gintervals.liftover
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$start), c("chrom", "start", "end")]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - chain with gaps", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("A", 300), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Create chain with gaps: maps [0-50] and [150-250], but NOT [50-150]
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 300 + 0 50 chr1 300 + 0 50 1\n", file = chain_file)
    cat("50\n\n", file = chain_file, append = TRUE)

    cat("chain 1000 source1 300 + 150 250 chr1 300 + 100 200 2\n", file = chain_file, append = TRUE)
    cat("100\n\n", file = chain_file, append = TRUE)

    # Test intervals - some in mapped regions, some in gap
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1"),
        start = c(10, 100, 200),
        end = c(20, 110, 210),
        stringsAsFactors = FALSE
    )

    # Create BED file
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
            chrom = src_intervals$chrom,
            start = src_intervals$start,
            end = src_intervals$end,
            name = paste0("int", 1:3),
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    if (file.exists(bed_output) && file.info(bed_output)$size > 0) {
        binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
        binary_result <- binary_result[order(binary_result$V2), c(1, 2, 3)]
        colnames(binary_result) <- c("chrom", "start", "end")
    } else {
        binary_result <- data.frame(chrom = character(), start = numeric(), end = numeric())
    }

    # Run gintervals.liftover
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(src_intervals, chain)
    if (!is.null(misha_result) && nrow(misha_result) > 0) {
        misha_result <- misha_result[order(misha_result$start), c("chrom", "start", "end")]
    }

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    if (nrow(binary_result) > 0) {
        expect_equal(as.character(misha_result$chrom), binary_result$chrom)
        expect_equal(as.numeric(misha_result$start), binary_result$start)
        expect_equal(as.numeric(misha_result$end), binary_result$end)
    }
})

test_that("gintervals.liftover matches liftOver binary - small intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Simple 1:1 mapping
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 1000 + 0 700 chr1 1000 + 0 700 1\n", file = chain_file)
    cat("700\n\n", file = chain_file, append = TRUE)

    # Test intervals - mix of very small and regular intervals
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1", "source1"),
        start = c(100, 101, 200, 500),
        end = c(101, 102, 250, 600),
        stringsAsFactors = FALSE
    )

    # Create BED file
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
            chrom = src_intervals$chrom,
            start = src_intervals$start,
            end = src_intervals$end,
            name = paste0("int", 1:4),
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    binary_result <- binary_result[order(binary_result$V2), c(1, 2, 3)]
    colnames(binary_result) <- c("chrom", "start", "end")

    # Run gintervals.liftover
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$start), c("chrom", "start", "end")]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - boundary intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("A", 150), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Map entire source chromosome to target
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 100 + 0 100 chr1 150 + 0 100 1\n", file = chain_file)
    cat("100\n\n", file = chain_file, append = TRUE)

    # Test intervals at boundaries
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1"),
        start = c(0, 45, 90),
        end = c(10, 55, 100),
        stringsAsFactors = FALSE
    )

    # Create BED file
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
            chrom = src_intervals$chrom,
            start = src_intervals$start,
            end = src_intervals$end,
            name = paste0("int", 1:3),
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    binary_result <- binary_result[order(binary_result$V2), c(1, 2, 3)]
    colnames(binary_result) <- c("chrom", "start", "end")

    # Run gintervals.liftover
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$start), c("chrom", "start", "end")]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - consecutive intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("A", 300), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Simple 1:1 mapping with offset
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 300 + 100 150 chr1 300 + 50 100 1\n", file = chain_file)
    cat("50\n\n", file = chain_file, append = TRUE)

    # Test consecutive non-overlapping intervals
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1", "source1"),
        start = c(100, 110, 120, 130),
        end = c(110, 120, 130, 140),
        stringsAsFactors = FALSE
    )

    # Create BED file
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
            chrom = src_intervals$chrom,
            start = src_intervals$start,
            end = src_intervals$end,
            name = paste0("int", 1:4),
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    binary_result <- binary_result[order(binary_result$V2), c(1, 2, 3)]
    colnames(binary_result) <- c("chrom", "start", "end")

    # Run gintervals.liftover
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$start), c("chrom", "start", "end")]

    # Compare results - all 4 consecutive intervals should be preserved
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - one-to-many with keep policy", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("A", 400), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Create chain with source overlap: source1[100-200] maps to TWO target locations
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # First mapping: source1[100-200] -> chr1[0-100]
    cat("chain 1000 source1 300 + 100 200 chr1 400 + 0 100 1\n", file = chain_file)
    cat("100\n\n", file = chain_file, append = TRUE)

    # Second mapping: source1[100-200] -> chr1[200-300] (source overlap)
    cat("chain 1000 source1 300 + 100 200 chr1 400 + 200 300 2\n", file = chain_file, append = TRUE)
    cat("100\n\n", file = chain_file, append = TRUE)

    # Test intervals in the overlapping region
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1"),
        start = c(100, 121, 180),
        end = c(120, 122, 185),
        stringsAsFactors = FALSE
    )

    # Create BED file
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
            chrom = src_intervals$chrom,
            start = src_intervals$start,
            end = src_intervals$end,
            name = paste0("int", 1:3),
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run liftOver binary with -multiple and -noSerial flags
    system2("liftOver",
        args = c("-multiple", "-noSerial", bed_input, chain_file, bed_output, bed_unmapped),
        stdout = FALSE, stderr = FALSE
    )

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    binary_result <- binary_result[order(binary_result$V2), c(1, 2, 3)]
    colnames(binary_result) <- c("chrom", "start", "end")

    # Run gintervals.liftover with keep policy
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto")
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$start), c("chrom", "start", "end")]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - unmapped intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("A", 300), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Create chain with gaps: maps [0-50] and [150-250], but NOT [50-150]
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 300 + 0 50 chr1 300 + 0 50 1\n", file = chain_file)
    cat("50\n\n", file = chain_file, append = TRUE)

    cat("chain 1000 source1 300 + 150 250 chr1 300 + 100 200 2\n", file = chain_file, append = TRUE)
    cat("100\n\n", file = chain_file, append = TRUE)

    # Test intervals - some mapped, some unmapped
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1"),
        start = c(10, 100, 200),
        end = c(20, 110, 210),
        stringsAsFactors = FALSE
    )

    # Create BED file
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
            chrom = src_intervals$chrom,
            start = src_intervals$start,
            end = src_intervals$end,
            name = paste0("int", 1:3),
            score = 0,
            strand = "+"
        ),
        file = bed_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    if (file.exists(bed_output) && file.info(bed_output)$size > 0) {
        binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
        binary_result <- binary_result[order(binary_result$V2), c(1, 2, 3)]
        colnames(binary_result) <- c("chrom", "start", "end")
    } else {
        binary_result <- data.frame(chrom = character(), start = numeric(), end = numeric())
    }

    # Check unmapped - should contain int2
    unmapped_lines <- readLines(bed_unmapped)
    unmapped_intervals <- unmapped_lines[!grepl("^#", unmapped_lines)]
    expect_true(length(unmapped_intervals) > 0)
    expect_true(any(grepl("int2", unmapped_intervals)))

    # Run gintervals.liftover
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(src_intervals, chain)
    if (!is.null(misha_result) && nrow(misha_result) > 0) {
        misha_result <- misha_result[order(misha_result$start), c("chrom", "start", "end")]
    }

    # Compare results - should only have the 2 mapped intervals
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(nrow(misha_result), 2) # Only 2 intervals mapped
    if (nrow(binary_result) > 0) {
        expect_equal(as.character(misha_result$chrom), binary_result$chrom)
        expect_equal(as.numeric(misha_result$start), binary_result$start)
        expect_equal(as.numeric(misha_result$end), binary_result$end)
    }
})
