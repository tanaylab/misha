load_test_db()

# Helper function to check if liftOver binary is available
has_liftover_binary <- function() {
    result <- tryCatch(
        {
            system2("which", "liftOver", stdout = TRUE, stderr = FALSE)
            TRUE
        },
        error = function(e) FALSE,
        warning = function(w) FALSE
    )
    if (!isTRUE(result)) {
        result <- tryCatch(
            {
                system2("liftOver", stdout = FALSE, stderr = FALSE)
                TRUE
            },
            error = function(e) FALSE,
            warning = function(w) FALSE
        )
    }
    return(isTRUE(result))
}

test_that("gtrack.liftover matches liftOver binary - basic sparse track", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

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

    # Create sparse track with specific values
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(12, 20, 35),
        end = c(18, 30, 45),
        stringsAsFactors = FALSE
    )
    src_values <- c(100, 200, 300)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    # Filename must start with "chr" for .gseq.import() to process it
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("T", 100), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Create chain: chrsource1[10-50] -> chr1[5-45]
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 chrsource1 100 + 10 50 chr1 100 + 5 45 1\n", file = chain_file)
    cat("40\n\n", file = chain_file, append = TRUE)

    # Create BED file with values for liftOver binary
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # BED format: chrom start end name score
    # Using score field to store the value
    # Chromosome name must match chain file: "chrsource1" (from database)
    cat("chrsource1\t12\t18\tint1\t100\n", file = bed_input)
    cat("chrsource1\t20\t30\tint2\t200\n", file = bed_input, append = TRUE)
    cat("chrsource1\t35\t45\tint3\t300\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(binary_result) <- c("chrom", "start", "end", "name", "value")
    binary_result <- binary_result[order(binary_result$start), ]

    # Run gtrack.liftover
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file)

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    misha_result <- misha_result[order(misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
    expect_equal(as.numeric(misha_result[[lifted_track]]), binary_result$value)
})

test_that("gtrack.liftover matches liftOver binary - one-to-many mapping with keep policy", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

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

    # Create sparse track
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(100, 121, 180),
        end = c(120, 122, 185),
        stringsAsFactors = FALSE
    )
    src_values <- c(4, 5, 6)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    # Filename must start with "chr" for .gseq.import() to process it
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("T", 400), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Create chain with source overlap: chrsource1[100-200) maps to TWO target locations
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # First mapping: chrsource1[100-200) -> chr1[0-100)
    cat("chain 1000 chrsource1 300 + 100 200 chr1 400 + 0 100 1\n", file = chain_file)
    cat("100\n\n", file = chain_file, append = TRUE)

    # Second mapping: chrsource1[100-200) -> chr1[200-300) (source overlap)
    cat("chain 1000 chrsource1 300 + 100 200 chr1 400 + 200 300 2\n", file = chain_file, append = TRUE)
    cat("100\n\n", file = chain_file, append = TRUE)

    # Create BED file for liftOver binary
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # Chromosome name must match chain file: "chrsource1" (from database)
    cat("chrsource1\t100\t120\tint1\t4\n", file = bed_input)
    cat("chrsource1\t121\t122\tint2\t5\n", file = bed_input, append = TRUE)
    cat("chrsource1\t180\t185\tint3\t6\n", file = bed_input, append = TRUE)

    # Run liftOver binary with -multiple and -noSerial flags
    # -noSerial prevents serial numbers from being written to the value field
    system2("liftOver", args = c("-multiple", "-noSerial", bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(binary_result) <- c("chrom", "start", "end", "name", "value")
    binary_result <- binary_result[order(binary_result$start), ]

    # Run gtrack.liftover with keep policy
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file,
        src_overlap_policy = "keep", tgt_overlap_policy = "auto"
    )

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    misha_result <- misha_result[order(misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
    expect_equal(as.numeric(misha_result[[lifted_track]]), binary_result$value)
})

test_that("gtrack.liftover matches liftOver binary - multiple chromosomes", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create source genome with multiple chromosomes
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
    # Database chromosome names are "chrsource1" and "chrsource2" (from filenames)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource2", "chrsource2"),
        start = c(10, 40, 15, 50),
        end = c(20, 50, 25, 60),
        stringsAsFactors = FALSE
    )
    src_values <- c(10, 20, 30, 40)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    # Filenames must start with "chr" for .gseq.import() to process them
    target_fasta1 <- file.path(tempdir(), "chr1.fasta")
    target_fasta2 <- file.path(tempdir(), "chr2.fasta")
    cat(">chr1\n", paste(rep("T", 200), collapse = ""), "\n", sep = "", file = target_fasta1)
    cat(">chr2\n", paste(rep("G", 200), collapse = ""), "\n", sep = "", file = target_fasta2)

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

    # chrsource1[0-60] -> chr1[10-70]
    cat("chain 1000 chrsource1 100 + 0 60 chr1 200 + 10 70 1\n", file = chain_file)
    cat("60\n\n", file = chain_file, append = TRUE)

    # chrsource2[10-70] -> chr2[20-80]
    cat("chain 1000 chrsource2 100 + 10 70 chr2 200 + 20 80 2\n", file = chain_file, append = TRUE)
    cat("60\n\n", file = chain_file, append = TRUE)

    # Create BED file
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # Chromosome names must match chain file: "chrsource1" and "chrsource2" (from database)
    cat("chrsource1\t10\t20\tint1\t10\n", file = bed_input)
    cat("chrsource1\t40\t50\tint2\t20\n", file = bed_input, append = TRUE)
    cat("chrsource2\t15\t25\tint3\t30\n", file = bed_input, append = TRUE)
    cat("chrsource2\t50\t60\tint4\t40\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(binary_result) <- c("chrom", "start", "end", "name", "value")
    binary_result <- binary_result[order(binary_result$chrom, binary_result$start), ]

    # Run gtrack.liftover
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file)

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    misha_result <- misha_result[order(misha_result$chrom, misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
    expect_equal(as.numeric(misha_result[[lifted_track]]), binary_result$value)
})

test_that("gtrack.liftover matches liftOver binary - reverse strand", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

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

    # Create sparse track
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 50, 100),
        end = c(20, 60, 120),
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 222, 333)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    # Filename must start with "chr" for .gseq.import() to process it
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("T", 200), collapse = ""), "\n", sep = "", file = target_fasta)

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

    # chrsource1[0-150) + -> chr1[50-200) - (reverse strand)
    cat("chain 1000 chrsource1 200 + 0 150 chr1 200 - 0 150 1\n", file = chain_file)
    cat("150\n\n", file = chain_file, append = TRUE)

    # Create BED file
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # Chromosome name must match chain file: "chrsource1" (from database)
    cat("chrsource1\t10\t20\tint1\t111\n", file = bed_input)
    cat("chrsource1\t50\t60\tint2\t222\n", file = bed_input, append = TRUE)
    cat("chrsource1\t100\t120\tint3\t333\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(binary_result) <- c("chrom", "start", "end", "name", "value")
    binary_result <- binary_result[order(binary_result$start), ]

    # Run gtrack.liftover
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file)

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    misha_result <- misha_result[order(misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    # For reverse strand, check that values are preserved
    expect_true(all(sort(as.numeric(misha_result[[lifted_track]])) == sort(binary_result$value)))
})

test_that("gtrack.liftover matches liftOver binary - chain with gaps", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create source genome - needs to be at least 250bp for the chain
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

    # Create sparse track
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 100, 200),
        end = c(20, 110, 210),
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 222, 333)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    # Filename must start with "chr" for .gseq.import() to process it
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("T", 300), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Create chain with gaps: maps [0-50) and [150-250), but NOT [50-150)
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 chrsource1 300 + 0 50 chr1 300 + 0 50 1\n", file = chain_file)
    cat("50\n\n", file = chain_file, append = TRUE)

    cat("chain 1000 chrsource1 300 + 150 250 chr1 300 + 100 200 2\n", file = chain_file, append = TRUE)
    cat("100\n\n", file = chain_file, append = TRUE)

    # Create BED file
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # Chromosome name must match chain file: "chrsource1" (from database)
    cat("chrsource1\t10\t20\tint1\t111\n", file = bed_input)
    cat("chrsource1\t100\t110\tint2\t222\n", file = bed_input, append = TRUE) # In gap, won't be mapped
    cat("chrsource1\t200\t210\tint3\t333\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    if (file.exists(bed_output) && file.info(bed_output)$size > 0) {
        binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
        colnames(binary_result) <- c("chrom", "start", "end", "name", "value")
        binary_result <- binary_result[order(binary_result$start), ]
    } else {
        binary_result <- data.frame(chrom = character(), start = numeric(), end = numeric(), name = character(), value = numeric())
    }

    # Run gtrack.liftover
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file)

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    if (nrow(misha_result) > 0) {
        misha_result <- misha_result[order(misha_result$start), ]
    }

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    if (nrow(binary_result) > 0) {
        expect_equal(as.character(misha_result$chrom), binary_result$chrom)
        expect_equal(as.numeric(misha_result$start), binary_result$start)
        expect_equal(as.numeric(misha_result$end), binary_result$end)
        expect_equal(as.numeric(misha_result[[lifted_track]]), binary_result$value)
    }
})

test_that("gtrack.liftover matches liftOver binary - small intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

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
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource1", "chrsource1"),
        start = c(100, 101, 200, 500),
        end = c(101, 102, 250, 600),
        stringsAsFactors = FALSE
    )
    src_values <- c(1, 2, 50, 100)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    # Filename must start with "chr" for .gseq.import() to process it
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("T", 1000), collapse = ""), "\n", sep = "", file = target_fasta)

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

    cat("chain 1000 chrsource1 1000 + 0 700 chr1 1000 + 0 700 1\n", file = chain_file)
    cat("700\n\n", file = chain_file, append = TRUE)

    # Create BED file
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # Chromosome name must match chain file: "chrsource1" (from database)
    cat("chrsource1\t100\t101\tint1\t1\n", file = bed_input)
    cat("chrsource1\t101\t102\tint2\t2\n", file = bed_input, append = TRUE)
    cat("chrsource1\t200\t250\tint3\t50\n", file = bed_input, append = TRUE)
    cat("chrsource1\t500\t600\tint4\t100\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(binary_result) <- c("chrom", "start", "end", "name", "value")
    binary_result <- binary_result[order(binary_result$start), ]

    # Run gtrack.liftover
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file)

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    misha_result <- misha_result[order(misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
    expect_equal(as.numeric(misha_result[[lifted_track]]), binary_result$value)
})

test_that("gtrack.liftover matches liftOver binary - boundary intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

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
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(0, 45, 90),
        end = c(10, 55, 100),
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 222, 333)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    # Filename must start with "chr" for .gseq.import() to process it
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("T", 150), collapse = ""), "\n", sep = "", file = target_fasta)

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

    cat("chain 1000 chrsource1 100 + 0 100 chr1 150 + 0 100 1\n", file = chain_file)
    cat("100\n\n", file = chain_file, append = TRUE)

    # Create BED file
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # Chromosome name must match chain file: "chrsource1" (from database)
    cat("chrsource1\t0\t10\tint1\t111\n", file = bed_input)
    cat("chrsource1\t45\t55\tint2\t222\n", file = bed_input, append = TRUE)
    cat("chrsource1\t90\t100\tint3\t333\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(binary_result) <- c("chrom", "start", "end", "name", "value")
    binary_result <- binary_result[order(binary_result$start), ]

    # Run gtrack.liftover
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file)

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    misha_result <- misha_result[order(misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
    expect_equal(as.numeric(misha_result[[lifted_track]]), binary_result$value)
})

test_that("gtrack.liftover matches liftOver binary - special values", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

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

    # Create sparse track with special values (zero, large positive, large negative)
    # Note: NaN cannot be represented in BED files easily, so we skip it
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 50, 70),
        end = c(20, 60, 80),
        stringsAsFactors = FALSE
    )
    src_values <- c(0, 1e10, -1e10)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    # Filename must start with "chr" for .gseq.import() to process it
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("T", 200), collapse = ""), "\n", sep = "", file = target_fasta)

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

    cat("chain 1000 chrsource1 200 + 0 100 chr1 200 + 0 100 1\n", file = chain_file)
    cat("100\n\n", file = chain_file, append = TRUE)

    # Create BED file
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # Chromosome name must match chain file: "chrsource1" (from database)
    cat("chrsource1\t10\t20\tint1\t0\n", file = bed_input)
    cat("chrsource1\t50\t60\tint2\t1e10\n", file = bed_input, append = TRUE)
    cat("chrsource1\t70\t80\tint3\t-1e10\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(binary_result) <- c("chrom", "start", "end", "name", "value")
    binary_result <- binary_result[order(binary_result$start), ]

    # Run gtrack.liftover
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file)

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    misha_result <- misha_result[order(misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
    # Check special values
    expect_true(any(misha_result[[lifted_track]] == 0))
    expect_true(any(misha_result[[lifted_track]] == 1e10))
    expect_true(any(misha_result[[lifted_track]] == -1e10))
})

test_that("gtrack.liftover matches liftOver binary - unmapped intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

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
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 100, 200), # 10 is mapped, 100 is in gap, 200 is mapped
        end = c(20, 110, 210),
        stringsAsFactors = FALSE
    )
    src_values <- c(111, 222, 333)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    # Filename must start with "chr" for .gseq.import() to process it
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("T", 300), collapse = ""), "\n", sep = "", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta)
    gdb.init(target_db)

    # Create chain with gaps: maps [0-50) and [150-250), but NOT [50-150)
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 chrsource1 300 + 0 50 chr1 300 + 0 50 1\n", file = chain_file)
    cat("50\n\n", file = chain_file, append = TRUE)

    cat("chain 1000 chrsource1 300 + 150 250 chr1 300 + 100 200 2\n", file = chain_file, append = TRUE)
    cat("100\n\n", file = chain_file, append = TRUE)

    # Create BED file
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # Chromosome name must match chain file: "chrsource1" (from database)
    cat("chrsource1\t10\t20\tint1\t111\n", file = bed_input)
    cat("chrsource1\t100\t110\tint2\t222\n", file = bed_input, append = TRUE) # In gap
    cat("chrsource1\t200\t210\tint3\t333\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    if (file.exists(bed_output) && file.info(bed_output)$size > 0) {
        binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
        colnames(binary_result) <- c("chrom", "start", "end", "name", "value")
        binary_result <- binary_result[order(binary_result$start), ]
    } else {
        binary_result <- data.frame(chrom = character(), start = numeric(), end = numeric(), name = character(), value = numeric())
    }

    # Check unmapped - should contain int2 (value 222)
    unmapped_lines <- readLines(bed_unmapped)
    unmapped_intervals <- unmapped_lines[!grepl("^#", unmapped_lines)]
    expect_true(length(unmapped_intervals) > 0)
    expect_true(any(grepl("int2", unmapped_intervals)))

    # Run gtrack.liftover
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file)

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    if (nrow(misha_result) > 0) {
        misha_result <- misha_result[order(misha_result$start), ]
    }

    # Compare results - should only have the 2 mapped intervals (111 and 333)
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(nrow(misha_result), 2) # Only 2 intervals mapped
    if (nrow(binary_result) > 0) {
        expect_equal(as.character(misha_result$chrom), binary_result$chrom)
        expect_equal(as.numeric(misha_result$start), binary_result$start)
        expect_equal(as.numeric(misha_result$end), binary_result$end)
        expect_equal(as.numeric(misha_result[[lifted_track]]), binary_result$value)
        # Verify the unmapped interval (222) is NOT in the results
        expect_false(any(misha_result[[lifted_track]] == 222))
    }
})

test_that("gtrack.liftover matches liftOver binary - consecutive intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

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

    # Create sparse track with consecutive non-overlapping intervals
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource1", "chrsource1"),
        start = c(100, 110, 120, 130),
        end = c(110, 120, 130, 140),
        stringsAsFactors = FALSE
    )
    src_values <- c(10, 20, 30, 40)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    # Filename must start with "chr" for .gseq.import() to process it
    target_fasta <- file.path(tempdir(), "chr1.fasta")
    cat(">chr1\n", paste(rep("T", 300), collapse = ""), "\n", sep = "", file = target_fasta)

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

    cat("chain 1000 chrsource1 300 + 100 150 chr1 300 + 50 100 1\n", file = chain_file)
    cat("50\n\n", file = chain_file, append = TRUE)

    # Create BED file
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # Chromosome name must match chain file: "chrsource1" (from database)
    cat("chrsource1\t100\t110\tint1\t10\n", file = bed_input)
    cat("chrsource1\t110\t120\tint2\t20\n", file = bed_input, append = TRUE)
    cat("chrsource1\t120\t130\tint3\t30\n", file = bed_input, append = TRUE)
    cat("chrsource1\t130\t140\tint4\t40\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE)
    colnames(binary_result) <- c("chrom", "start", "end", "name", "value")
    binary_result <- binary_result[order(binary_result$start), ]

    # Run gtrack.liftover
    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file)

    # Extract results
    misha_result <- gextract(lifted_track, gintervals.all())
    misha_result <- misha_result[order(misha_result$start), ]

    # Compare results - all 4 consecutive intervals should be preserved
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
    expect_equal(as.numeric(misha_result[[lifted_track]]), binary_result$value)
})
