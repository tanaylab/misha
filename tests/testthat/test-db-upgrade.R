test_that("gdb.upgrade upgrades legacy database to indexed format", {
    # Create a small legacy database with per-chromosome .seq files
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

    # Create legacy .seq files (raw sequence data)
    writeBin(charToRaw("ACTGACTGACTG"), file.path(test_db, "seq", "chr1.seq"))
    writeBin(charToRaw("GGGGCCCC"), file.path(test_db, "seq", "chr2.seq"))

    # Upgrade to indexed format (skip validation for minimal test DB)
    expect_message(
        gdb.upgrade(groot = test_db, interactive = FALSE, validate = FALSE),
        "Upgrading database to indexed format"
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

test_that("gdb.upgrade with remove_old_files removes legacy .seq files", {
    # Create a small legacy database
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

    # Create legacy .seq file
    seq_file <- file.path(test_db, "seq", "test.seq")
    writeBin(charToRaw("ACTGACTG"), seq_file)

    # Upgrade with removal (skip validation for minimal test DB)
    gdb.upgrade(groot = test_db, remove_old_files = TRUE, interactive = FALSE, validate = FALSE)

    # Check that old file was removed
    expect_false(file.exists(seq_file))

    # Check that new files exist
    expect_true(file.exists(file.path(test_db, "seq", "genome.idx")))
    expect_true(file.exists(file.path(test_db, "seq", "genome.seq")))
})

test_that("gdb.upgrade skips if already upgraded", {
    # Create a database in indexed format
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta)

    # Try to upgrade - should skip
    expect_message(
        gdb.upgrade(groot = test_db, interactive = FALSE),
        "already in indexed format"
    )
})

test_that("gdb.upgrade validates upgraded sequences", {
    # Create a proper database first, then simulate legacy format
    test_fasta <- tempfile(fileext = ".fasta")
    seq_data <- paste(rep("ACTG", 50), collapse = "")
    cat(">chr1\n", seq_data, "\n", sep = "", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    # Create a proper database with legacy format
    options(gmulticontig.indexed_format = FALSE)
    suppressMessages(gdb.create(groot = test_db, fasta = test_fasta))

    # Now upgrade with validation enabled
    expect_message(
        gdb.upgrade(groot = test_db, interactive = FALSE, validate = TRUE),
        "Validating upgrade"
    )

    # Check that validation passed
    expect_message(
        expect_message(
            gdb.upgrade(groot = test_db, interactive = FALSE, validate = TRUE),
            "already in indexed format"
        ),
        NA # Should not get validation message since already upgraded
    )
})

# Restore the test database after all upgrade tests
suppressMessages(gdb.init("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"))
