load_test_db()
test_that("gdb.convert_to_indexed converts per-chromosome database to indexed format", {
    local_db_state()

    # Create a small per-chromosome database with per-chromosome .seq files
    test_db <- tempfile()
    dir.create(test_db)
    dir.create(file.path(test_db, "seq"))

    withr::defer(unlink(test_db, recursive = TRUE))

    # Create chrom_sizes.txt
    chrom_sizes <- data.frame(
        chrom = c("chr1", "chr2"),
        size = c(12, 8)
    )
    write.table(chrom_sizes, file.path(test_db, "chrom_sizes.txt"),
        row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )

    # Create per-chromosome .seq files (raw sequence data)
    writeBin(charToRaw("ACTGACTGACTG"), file.path(test_db, "seq", "chr1.seq"))
    writeBin(charToRaw("GGGGCCCC"), file.path(test_db, "seq", "chr2.seq"))

    # Convert to indexed format (skip validation for minimal test DB)
    expect_message(
        gdb.convert_to_indexed(groot = test_db, force = TRUE, validate = FALSE, verbose = TRUE),
        "Converting database to indexed format"
    )

    # Check that indexed files were created
    expect_true(file.exists(file.path(test_db, "seq", "genome.idx")))
    expect_true(file.exists(file.path(test_db, "seq", "genome.seq")))

    # Verify the index file has reasonable content
    idx_size <- file.info(file.path(test_db, "seq", "genome.idx"))$size
    expect_true(idx_size > 0)

    # Verify the sequence file has the expected size (12 + 8 = 20 bytes)
    seq_size <- file.info(file.path(test_db, "seq", "genome.seq"))$size
    expect_equal(seq_size, 20)
})

test_that("gdb.convert_to_indexed preserves chrom_sizes.txt order (non-alphabetical)", {
    local_db_state()

    # Create a small per-chromosome database with NON-alphabetical order
    test_db <- tempfile()
    dir.create(test_db)
    dir.create(file.path(test_db, "seq"))

    withr::defer(unlink(test_db, recursive = TRUE))

    # Create chrom_sizes.txt in NON-alphabetical order
    # chr15 comes before chr10 alphabetically, but we want to preserve original order
    original_chrom_sizes <- data.frame(
        chrom = c("chr15", "chr10", "chr17_random", "chr1"),
        size = c(15, 10, 17, 20)
    )
    write.table(original_chrom_sizes, file.path(test_db, "chrom_sizes.txt"),
        row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )

    # Create per-chromosome .seq files
    writeBin(charToRaw(paste(rep("A", 15), collapse = "")), file.path(test_db, "seq", "chr15.seq"))
    writeBin(charToRaw(paste(rep("C", 10), collapse = "")), file.path(test_db, "seq", "chr10.seq"))
    writeBin(charToRaw(paste(rep("G", 17), collapse = "")), file.path(test_db, "seq", "chr17_random.seq"))
    writeBin(charToRaw(paste(rep("T", 20), collapse = "")), file.path(test_db, "seq", "chr1.seq"))

    # Convert to indexed format
    suppressMessages(
        gdb.convert_to_indexed(groot = test_db, force = TRUE, validate = FALSE, verbose = FALSE)
    )

    # Read the converted chrom_sizes.txt
    converted_chrom_sizes <- read.table(
        file.path(test_db, "chrom_sizes.txt"),
        header = FALSE, stringsAsFactors = FALSE, sep = "\t"
    )
    colnames(converted_chrom_sizes) <- c("chrom", "size")

    # Verify order is preserved (should match original order, not alphabetically sorted)
    expect_equal(converted_chrom_sizes$chrom, original_chrom_sizes$chrom)
    expect_equal(converted_chrom_sizes$size, original_chrom_sizes$size)

    # Verify it's NOT sorted alphabetically
    expect_false(identical(converted_chrom_sizes$chrom, sort(converted_chrom_sizes$chrom)))
})

test_that("gdb.convert_to_indexed preserves chrom_sizes.txt order when adding chr prefix", {
    local_db_state()

    # Create a database with chrom_sizes.txt WITHOUT chr prefix
    # but .seq files WITH chr prefix
    test_db <- tempfile()
    dir.create(test_db)
    dir.create(file.path(test_db, "seq"))

    withr::defer(unlink(test_db, recursive = TRUE))

    # Original chrom_sizes.txt WITHOUT chr prefix, in non-alphabetical order
    original_chrom_sizes <- data.frame(
        chrom = c("15", "10", "17_random", "1"),
        size = c(15, 10, 17, 20)
    )
    write.table(original_chrom_sizes, file.path(test_db, "chrom_sizes.txt"),
        row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )

    # Create .seq files WITH chr prefix
    writeBin(charToRaw(paste(rep("A", 15), collapse = "")), file.path(test_db, "seq", "chr15.seq"))
    writeBin(charToRaw(paste(rep("C", 10), collapse = "")), file.path(test_db, "seq", "chr10.seq"))
    writeBin(charToRaw(paste(rep("G", 17), collapse = "")), file.path(test_db, "seq", "chr17_random.seq"))
    writeBin(charToRaw(paste(rep("T", 20), collapse = "")), file.path(test_db, "seq", "chr1.seq"))

    # Convert to indexed format
    suppressMessages(
        gdb.convert_to_indexed(groot = test_db, force = TRUE, validate = FALSE, verbose = FALSE)
    )

    # Read the converted chrom_sizes.txt
    converted_chrom_sizes <- read.table(
        file.path(test_db, "chrom_sizes.txt"),
        header = FALSE, stringsAsFactors = FALSE, sep = "\t"
    )
    colnames(converted_chrom_sizes) <- c("chrom", "size")

    # Verify order is preserved (should match original order)
    # Note: chromosome names should now have chr prefix, but order should be preserved
    expected_order <- c("chr15", "chr10", "chr17_random", "chr1")
    expect_equal(converted_chrom_sizes$chrom, expected_order)
    expect_equal(converted_chrom_sizes$size, original_chrom_sizes$size)

    # Verify it's NOT sorted alphabetically
    expect_false(identical(converted_chrom_sizes$chrom, sort(converted_chrom_sizes$chrom)))
})

test_that("gdb.convert_to_indexed preserves chrom_sizes.txt order when removing chr prefix", {
    local_db_state()

    # Create a database with chrom_sizes.txt WITH chr prefix
    # but .seq files WITHOUT chr prefix
    test_db <- tempfile()
    dir.create(test_db)
    dir.create(file.path(test_db, "seq"))

    withr::defer(unlink(test_db, recursive = TRUE))

    # Original chrom_sizes.txt WITH chr prefix, in non-alphabetical order
    original_chrom_sizes <- data.frame(
        chrom = c("chr15", "chr10", "chr17_random", "chr1"),
        size = c(15, 10, 17, 20)
    )
    write.table(original_chrom_sizes, file.path(test_db, "chrom_sizes.txt"),
        row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )

    # Create .seq files WITHOUT chr prefix
    writeBin(charToRaw(paste(rep("A", 15), collapse = "")), file.path(test_db, "seq", "15.seq"))
    writeBin(charToRaw(paste(rep("C", 10), collapse = "")), file.path(test_db, "seq", "10.seq"))
    writeBin(charToRaw(paste(rep("G", 17), collapse = "")), file.path(test_db, "seq", "17_random.seq"))
    writeBin(charToRaw(paste(rep("T", 20), collapse = "")), file.path(test_db, "seq", "1.seq"))

    # Convert to indexed format
    suppressMessages(
        gdb.convert_to_indexed(groot = test_db, force = TRUE, validate = FALSE, verbose = FALSE)
    )

    # Read the converted chrom_sizes.txt
    converted_chrom_sizes <- read.table(
        file.path(test_db, "chrom_sizes.txt"),
        header = FALSE, stringsAsFactors = FALSE, sep = "\t"
    )
    colnames(converted_chrom_sizes) <- c("chrom", "size")

    # Verify order is preserved (should match original order)
    # Note: chromosome names should now be without chr prefix, but order should be preserved
    expected_order <- c("15", "10", "17_random", "1")
    expect_equal(converted_chrom_sizes$chrom, expected_order)
    expect_equal(converted_chrom_sizes$size, original_chrom_sizes$size)

    # Verify it's NOT sorted alphabetically
    expect_false(identical(converted_chrom_sizes$chrom, sort(converted_chrom_sizes$chrom)))
})

test_that("gdb.convert_to_indexed preserves chrom_sizes.txt order with mixed prefix patterns", {
    local_db_state()

    # Create a database with mixed naming patterns
    test_db <- tempfile()
    dir.create(test_db)
    dir.create(file.path(test_db, "seq"))

    withr::defer(unlink(test_db, recursive = TRUE))

    # Original chrom_sizes.txt in non-alphabetical order
    original_chrom_sizes <- data.frame(
        chrom = c("chrZ", "chrA", "chrM"),
        size = c(100, 200, 50)
    )
    write.table(original_chrom_sizes, file.path(test_db, "chrom_sizes.txt"),
        row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )

    # Create .seq files
    writeBin(charToRaw(paste(rep("Z", 100), collapse = "")), file.path(test_db, "seq", "chrZ.seq"))
    writeBin(charToRaw(paste(rep("A", 200), collapse = "")), file.path(test_db, "seq", "chrA.seq"))
    writeBin(charToRaw(paste(rep("M", 50), collapse = "")), file.path(test_db, "seq", "chrM.seq"))

    # Convert to indexed format
    suppressMessages(
        gdb.convert_to_indexed(groot = test_db, force = TRUE, validate = FALSE, verbose = FALSE)
    )

    # Read the converted chrom_sizes.txt
    converted_chrom_sizes <- read.table(
        file.path(test_db, "chrom_sizes.txt"),
        header = FALSE, stringsAsFactors = FALSE, sep = "\t"
    )
    colnames(converted_chrom_sizes) <- c("chrom", "size")

    # Verify order is preserved (chrZ, chrA, chrM should remain in that order)
    expect_equal(converted_chrom_sizes$chrom, original_chrom_sizes$chrom)
    expect_equal(converted_chrom_sizes$size, original_chrom_sizes$size)

    # Verify it's NOT sorted alphabetically (chrA would come first if sorted)
    expect_false(identical(converted_chrom_sizes$chrom, sort(converted_chrom_sizes$chrom)))
    expect_equal(converted_chrom_sizes$chrom[1], "chrZ")
})

test_that("gdb.convert_to_indexed preserves chrom_sizes.txt order and works with gdb.init", {
    local_db_state()

    # Create a database with non-alphabetical order
    test_db <- tempfile()
    dir.create(test_db)
    dir.create(file.path(test_db, "seq"))
    dir.create(file.path(test_db, "tracks"))

    withr::defer(unlink(test_db, recursive = TRUE))

    # Original chrom_sizes.txt in non-alphabetical order
    original_chrom_sizes <- data.frame(
        chrom = c("chr15", "chr10", "chr1"),
        size = c(15, 10, 20)
    )
    write.table(original_chrom_sizes, file.path(test_db, "chrom_sizes.txt"),
        row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )

    # Create .seq files
    writeBin(charToRaw(paste(rep("A", 15), collapse = "")), file.path(test_db, "seq", "chr15.seq"))
    writeBin(charToRaw(paste(rep("C", 10), collapse = "")), file.path(test_db, "seq", "chr10.seq"))
    writeBin(charToRaw(paste(rep("T", 20), collapse = "")), file.path(test_db, "seq", "chr1.seq"))

    # Convert to indexed format
    suppressMessages(
        gdb.convert_to_indexed(groot = test_db, force = TRUE, validate = FALSE, verbose = FALSE)
    )

    # Verify order is preserved after conversion
    converted_chrom_sizes <- read.table(
        file.path(test_db, "chrom_sizes.txt"),
        header = FALSE, stringsAsFactors = FALSE, sep = "\t"
    )
    colnames(converted_chrom_sizes) <- c("chrom", "size")
    expect_equal(converted_chrom_sizes$chrom, original_chrom_sizes$chrom)

    # Initialize the database and verify ALLGENOME preserves order
    suppressMessages(gdb.init(test_db))

    # ALLGENOME should match chrom_sizes.txt order (not sorted)
    allgenome <- .misha$ALLGENOME[[1]]
    expect_equal(as.character(allgenome$chrom), converted_chrom_sizes$chrom)

    # Verify chromosome sizes are correct
    expect_equal(allgenome$end[allgenome$chrom == "chr15"], 15)
    expect_equal(allgenome$end[allgenome$chrom == "chr10"], 10)
    expect_equal(allgenome$end[allgenome$chrom == "chr1"], 20)
})

test_that("gdb.convert_to_indexed with remove_old_files removes per-chromosome .seq files", {
    local_db_state()

    # Create a small per-chromosome database
    test_db <- tempfile()
    dir.create(test_db)
    dir.create(file.path(test_db, "seq"))

    withr::defer(unlink(test_db, recursive = TRUE))

    # Create chrom_sizes.txt
    chrom_sizes <- data.frame(
        chrom = c("test"),
        size = c(8)
    )
    write.table(chrom_sizes, file.path(test_db, "chrom_sizes.txt"),
        row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )

    # Create per-chromosome .seq file
    seq_file <- file.path(test_db, "seq", "test.seq")
    writeBin(charToRaw("ACTGACTG"), seq_file)

    # Convert with removal (skip validation for minimal per-chromosome test DB)
    gdb.convert_to_indexed(groot = test_db, remove_old_files = TRUE, force = TRUE, validate = FALSE)

    # Check that old file was removed
    expect_false(file.exists(seq_file))

    # Check that new files exist
    expect_true(file.exists(file.path(test_db, "seq", "genome.idx")))
    expect_true(file.exists(file.path(test_db, "seq", "genome.seq")))
})

test_that("gdb.convert_to_indexed skips if already converted", {
    local_db_state()

    # Create a database in indexed format
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
        # Try to convert - should skip
        expect_message(
            gdb.convert_to_indexed(groot = test_db, force = TRUE, verbose = TRUE),
            "already in indexed format"
        )
    })
})

test_that("gdb.convert_to_indexed validates converted sequences", {
    local_db_state()

    # Create a proper database first, then simulate per-chromosome format
    test_fasta <- tempfile(fileext = ".fasta")
    seq_data <- paste(rep("ACTG", 50), collapse = "")
    cat(">chr1\n", seq_data, "\n", sep = "", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    # Create a proper database with per-chromosome format
    withr::with_options(list(gmulticontig.indexed_format = FALSE), {
        suppressMessages(gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE))
    })

    # Now conversion with validation enabled
    expect_message(
        gdb.convert_to_indexed(groot = test_db, force = TRUE, validate = TRUE, verbose = TRUE),
        "Validating conversion"
    )

    # Check that validation passed
    expect_message(
        expect_message(
            gdb.convert_to_indexed(groot = test_db, force = TRUE, validate = TRUE, verbose = TRUE),
            "already in indexed format"
        ),
        NA # Should not get validation message since already converted
    )
})
