test_that("indexed database preserves chrom_sizes.txt order for chromid mapping", {
    # Regression test for bug where chromosome identifiers resolved to wrong
    # chromosomes due to sorting mismatch between chrom_sizes.txt and ALLGENOME
    #
    # BUG SCENARIO: When chrom_sizes.txt is unsorted but genome.idx is sorted,
    # chromosome lookups would use the wrong chromid

    withr::defer(gdb.init("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"))

    # Create a FASTA file with chromosomes in NON-alphabetical order
    test_fasta <- tempfile(fileext = ".fasta")
    cat(
        ">chr15\n", "AAAAAAAAAA", "\n",
        ">chr10\n", "CCCCCCCCCC", "\n",
        ">chr17_random\n", "GGGGGGGGGG", "\n",
        ">chr1\n", "TTTTTTTTTT", "\n",
        sep = "",
        file = test_fasta
    )

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    # Create database in indexed format
    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        suppressMessages(gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE))
    })

    # Initialize the database
    suppressMessages(gdb.init(test_db))

    # Check that chrom_sizes.txt is sorted (C++ import sorts it)
    chrom_sizes <- read.table(file.path(test_db, "chrom_sizes.txt"),
                              header = FALSE, stringsAsFactors = FALSE)
    colnames(chrom_sizes) <- c("chrom", "size")

    # After C++ import, chromosomes should be sorted alphabetically
    expect_equal(chrom_sizes$chrom, sort(chrom_sizes$chrom))

    # Get ALLGENOME - should match chrom_sizes.txt order (NO SORTING in R for indexed)
    allgenome <- .misha$ALLGENOME[[1]]
    expect_equal(as.character(allgenome$chrom), chrom_sizes$chrom)

    # KEY TEST: Extract sequences from each chromosome
    # Before the fix, chromid mismatches would cause wrong sequences to be returned

    seq_chr1 <- gseq.extract(gintervals("chr1", 0, 10))
    seq_chr10 <- gseq.extract(gintervals("chr10", 0, 10))
    seq_chr15 <- gseq.extract(gintervals("chr15", 0, 10))
    seq_chr17 <- gseq.extract(gintervals("chr17_random", 0, 10))

    # Each chromosome should have its expected sequence
    expect_equal(seq_chr1, "TTTTTTTTTT")
    expect_equal(seq_chr10, "CCCCCCCCCC")
    expect_equal(seq_chr15, "AAAAAAAAAA")
    expect_equal(seq_chr17, "GGGGGGGGGG")

    # Verify ALLGENOME lengths are correct
    expect_equal(allgenome$end[allgenome$chrom == "chr1"], 10)
    expect_equal(allgenome$end[allgenome$chrom == "chr10"], 10)
    expect_equal(allgenome$end[allgenome$chrom == "chr15"], 10)
    expect_equal(allgenome$end[allgenome$chrom == "chr17_random"], 10)
})

