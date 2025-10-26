test_that("gdb.convert_to_indexed converts per-chromosome database to indexed format", {
    withr::defer(gdb.init("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"))

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

test_that("gdb.convert_to_indexed with remove_old_files removes per-chromosome .seq files", {
    withr::defer(gdb.init("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"))

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
    withr::defer(gdb.init("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"))

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
    withr::defer(gdb.init("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"))

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
    options(gmulticontig.indexed_format = FALSE)
    suppressMessages(gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE))

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
