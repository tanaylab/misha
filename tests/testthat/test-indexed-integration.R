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

test_that("indexed dense track works correctly across many chromosomes", {
    # Regression test: indexed tracks previously re-mmapped the entire track.dat
    # for every chromosome during iterator init and chromosome transitions,
    # making them unusable on genomes with many contigs.
    n_chroms <- 50
    test_fasta <- tempfile(fileext = ".fasta")
    for (i in seq_len(n_chroms)) {
        cat(sprintf(">chr%d\n%s\n", i, paste(rep("ACGT", 25), collapse = "")),
            file = test_fasta, append = (i > 1))
    }

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::local_options(list(
        gmulticontig.indexed_format = TRUE,
        gmultitasking = FALSE
    ))
    gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE)
    gdb.init(test_db)

    # Create a dense track with known values on all chromosomes
    all_intervs <- gintervals.all()
    gdir.create("testidx")

    all_track_intervs <- do.call(rbind, lapply(seq_len(nrow(all_intervs)), function(i) {
        chrom <- as.character(all_intervs$chrom[i])
        chrom_end <- all_intervs$end[i]
        data.frame(
            chrom = chrom,
            start = seq(0, chrom_end - 10, by = 10),
            end = seq(10, chrom_end, by = 10)
        )
    }))
    all_vals <- rep(1.0, nrow(all_track_intervs))
    gtrack.create_dense("testidx.dense1", "Test dense track", all_track_intervs, all_vals, 10, 0)

    # Convert to indexed format
    gtrack.convert_to_indexed("testidx.dense1")

    # Verify indexed format was applied
    info <- gtrack.info("testidx.dense1")
    expect_equal(info$format, "indexed")

    # Single-chromosome extraction should work and return correct values
    first_chrom <- as.character(all_intervs$chrom[1])
    result <- gextract("testidx.dense1", gintervals(first_chrom, 0, 20))
    expect_equal(nrow(result), 2)
    expect_true(all(result$testidx.dense1 == 1.0))

    # Multi-chromosome extraction across all chromosomes should work
    # This exercises the start_chrom code path for every chromosome
    result_all <- gextract("testidx.dense1", gintervals.all())
    expect_equal(nrow(result_all), sum(all_intervs$end) / 10)

    # Virtual track with function should work across chromosome transitions
    gvtrack.create("v_max", "testidx.dense1", func = "max")
    result_vtrack <- gextract("v_max", gintervals(first_chrom, 0, 100))
    expect_true(all(result_vtrack$v_max == 1.0))
    gvtrack.rm("v_max")

    # Clean up
    gtrack.rm("testidx.dense1", force = TRUE)
    gdir.rm("testidx", force = TRUE, recursive = TRUE)
})

# Restore the test database after all integration tests
suppressMessages(gdb.init("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"))
