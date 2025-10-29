load_test_db()
# Integration tests for indexed genome format

test_that("indexed format works with gseq.extract on all strands", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTGACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
        gdb.init(test_db)

        # Forward strand
        fwd <- gseq.extract(gintervals("test", 0, 8, 1))
        expect_equal(fwd, "ACTGACTG")

        # Reverse strand
        rev <- gseq.extract(gintervals("test", 0, 8, -1))
        expect_equal(rev, "CAGTCAGT")

        # No strand specified (should default to forward)
        nostrand <- gseq.extract(gintervals("test", 0, 8))
        expect_equal(nostrand, "ACTGACTG")
    })
})

test_that("indexed format works with multiple interval extraction", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nAAAAAAAA\n>chr2\nCCCCCCCC\n>chr3\nGGGGGGGG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
        gdb.init(test_db)

        # Extract from multiple chromosomes
        intervals <- data.frame(
            chrom = c("chr1", "chr2", "chr3"),
            start = c(0, 0, 0),
            end = c(4, 4, 4)
        )

        seqs <- gseq.extract(intervals)
        expect_equal(seqs, c("AAAA", "CCCC", "GGGG"))
    })
})

test_that("indexed format works with gintervals.all()", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">a\nACTG\n>b\nGGGG\n>c\nCCCC\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::local_options(list(gmulticontig.indexed_format = TRUE))

    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)
    all_intervals <- gintervals.all()
    expect_equal(nrow(all_intervals), 3)
    expect_equal(as.character(all_intervals$chrom), c("a", "b", "c"))
    expect_equal(all_intervals$start, c(0, 0, 0))
    expect_equal(all_intervals$end, c(4, 4, 4))
})

test_that("indexed format works with gintervals.2d.all() for small genomes", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTG\n>chr2\nGGGG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    # Set threshold high so 2D is materialized
    withr::with_options(list(gmulticontig.indexed_format = TRUE, gmulticontig.2d.threshold = 100), {
        gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
        gdb.init(test_db)

        genome_2d <- gintervals.2d.all()

        expect_false(is.null(genome_2d))
        expect_equal(nrow(genome_2d), 4) # 2x2 combinations
    })
})

test_that("indexed format defers 2D for large genomes", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Create 10 contigs
    for (i in 1:10) {
        cat(sprintf(">contig%d\nACTG\n", i), file = test_fasta, append = (i > 1))
    }

    test_db <- tempfile()
    withr::defer(unlink(test_db, recursive = TRUE))

    # Set threshold low to force deferral
    withr::with_options(list(gmulticontig.indexed_format = TRUE, gmulticontig.2d.threshold = 5), {
        expect_message(
            gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE),
            "Deferring 2D genome generation"
        )
        gdb.init(test_db)

        # gintervals.2d.all() should still work (generates on demand)
        genome_2d <- gintervals.2d.all()
        expect_false(is.null(genome_2d))
        # With deferred generation, we get chromosome names only, not all pairs
        expect_true(nrow(genome_2d) >= 10) # At least 10 rows
    })
})

test_that("indexed format works with chromosome name lookups", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTG\n>chr2\nGGGG\n>chrX\nCCCC\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::local_options(list(gmulticontig.indexed_format = TRUE))
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Test that all chromosome names are accessible
    expect_silent(gseq.extract(gintervals("chr1", 0, 4)))
    expect_silent(gseq.extract(gintervals("chr2", 0, 4)))
    expect_silent(gseq.extract(gintervals("chrX", 0, 4)))

    # Test invalid chromosome name (error happens in gintervals, not gseq.extract)
    expect_error(
        gseq.extract(gintervals("chrY", 0, 4))
    )
})

test_that("indexed format handles boundary conditions", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTGACTGACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::local_options(list(gmulticontig.indexed_format = TRUE))
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Extract from start
    start <- gseq.extract(gintervals("test", 0, 1))
    expect_equal(start, "A")

    # Extract to end
    end <- gseq.extract(gintervals("test", 11, 12))
    expect_equal(end, "G")

    # Extract full sequence
    full <- gseq.extract(gintervals("test", 0, 12))
    expect_equal(full, "ACTGACTGACTG")

    # 0-length intervals are not allowed in misha
    # (intervals must have start < end)
})

test_that("indexed format persists across database reloads", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTGACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::local_options(list(gmulticontig.indexed_format = TRUE))
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Extract sequence
    seq1 <- gseq.extract(gintervals("test", 0, 8))

    # Reload database
    suppressMessages(gdb.init(test_db))

    # Extract again - should be identical
    seq2 <- gseq.extract(gintervals("test", 0, 8))

    expect_equal(seq1, seq2)
    expect_equal(seq1, "ACTGACTG")
})

test_that("indexed format is compatible with gdb.reload()", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::local_options(list(gmulticontig.indexed_format = TRUE))
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Reload should work
    expect_silent(suppressMessages(gdb.reload()))

    # Sequence extraction should still work
    seq <- gseq.extract(gintervals("test", 0, 4))
    expect_equal(seq, "ACTG")
})

test_that("indexed format validates on init", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::local_options(list(gmulticontig.indexed_format = TRUE))
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)

    # Corrupt the index file
    idx_file <- file.path(test_db, "seq", "genome.idx")
    idx_data <- readBin(idx_file, "raw", n = 1000)
    if (length(idx_data) > 50) {
        # Corrupt a byte
        idx_data[50] <- as.raw(bitwXor(as.integer(idx_data[50]), 0xFF))
        writeBin(idx_data, idx_file)

        # Init should fail with validation error
        expect_error(
            gdb.init(test_db),
            "checksum|corrupt|validation"
        )
    }
})

# Restore the test database after all integration tests
suppressMessages(gdb.init("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"))
