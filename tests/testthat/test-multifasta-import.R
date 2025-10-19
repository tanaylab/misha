test_that("multi-FASTA import creates indexed format", {
    # Create a small test multi-FASTA file
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">contig1\nACTGACTGACTG\n>contig2\nGGGGCCCC\n>contig3\nTATATA\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    # Create database
    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)

    # Check indexed format files exist
    expect_true(file.exists(file.path(test_db, "seq", "genome.idx")))
    expect_true(file.exists(file.path(test_db, "seq", "genome.seq")))

    # Check legacy files don't exist
    expect_false(file.exists(file.path(test_db, "seq", "chrcontig1.seq")))
    expect_false(file.exists(file.path(test_db, "seq", "contig1.seq")))
})

test_that("multi-FASTA import extracts sequences correctly", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">seq1\nACTGACTGACTG\n>seq2\nGGGGCCCC\n>seq3\nTATATA\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Extract and verify sequences
    seq1 <- gseq.extract(gintervals("seq1", 0, 12))
    expect_equal(seq1, "ACTGACTGACTG")

    seq2 <- gseq.extract(gintervals("seq2", 0, 8))
    expect_equal(seq2, "GGGGCCCC")

    seq3 <- gseq.extract(gintervals("seq3", 0, 6))
    expect_equal(seq3, "TATATA")

    # Test partial extraction
    partial <- gseq.extract(gintervals("seq1", 4, 8))
    expect_equal(partial, "ACTG")
})

test_that("multi-FASTA import handles reverse complement", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTGACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Forward strand
    fwd <- gseq.extract(gintervals("test", 0, 8, 1))
    expect_equal(fwd, "ACTGACTG")

    # Reverse strand
    rev <- gseq.extract(gintervals("test", 0, 8, -1))
    expect_equal(rev, "CAGTCAGT")
})

test_that("multi-FASTA import sanitizes contig names", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Test various header formats
    cat(">lcl|scaffold_1 description here\nACTG\n", file = test_fasta)
    cat(">gi|12345|ref|NC_000001.1| Homo sapiens chromosome 1\nGGGG\n", file = test_fasta, append = TRUE)
    cat(">gnl|ASSEMBLY|ctg123\nCCCC\n", file = test_fasta, append = TRUE)
    cat(">simple\nTTTT\n", file = test_fasta, append = TRUE)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    chroms <- gintervals.all()

    # Check sanitized names (prefixes removed)
    expect_true("scaffold_1" %in% chroms$chrom)
    expect_true("NC_000001.1" %in% chroms$chrom)
    expect_true("ctg123" %in% chroms$chrom)
    expect_true("simple" %in% chroms$chrom)

    # Verify sequences are accessible
    expect_equal(gseq.extract(gintervals("scaffold_1", 0, 4)), "ACTG")
    expect_equal(gseq.extract(gintervals("simple", 0, 4)), "TTTT")
})

test_that("multi-FASTA import with small genome materializes 2D", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTG\n>chr2\nGGGG\n>chr3\nCCCC\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    # Set threshold high to force materialization
    options(gmulticontig.indexed_format = TRUE, gmulticontig.2d.threshold = 100)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # 2D genome should be materialized (3 contigs < 100)
    genome_2d <- gintervals.2d.all()
    expect_false(is.null(genome_2d))
    expect_equal(nrow(genome_2d), 9) # 3x3 combinations
    expect_true(all(c("chrom1", "chrom2") %in% colnames(genome_2d)))
})

test_that("multi-FASTA import with large genome defers 2D", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Create 10 small contigs
    for (i in 1:10) {
        cat(sprintf(">contig%d\n%s\n", i, paste(rep("ACTG", 3), collapse = "")),
            file = test_fasta, append = (i > 1)
        )
    }

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    # Set threshold low to force deferral
    options(gmulticontig.indexed_format = TRUE, gmulticontig.2d.threshold = 5)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # gintervals.2d.all() should generate on demand
    genome_2d <- gintervals.2d.all()
    expect_false(is.null(genome_2d))
    expect_true(nrow(genome_2d) > 0)
})

test_that("multi-FASTA import handles gzipped files", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test1\nACTGACTG\n>test2\nGGGGCCCC\n", file = test_fasta)

    # Gzip the file
    test_fasta_gz <- paste0(test_fasta, ".gz")
    system(sprintf("gzip -c %s > %s", test_fasta, test_fasta_gz))

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
        unlink(test_fasta_gz)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta_gz)
    gdb.init(test_db)

    # Verify sequences
    seq1 <- gseq.extract(gintervals("test1", 0, 8))
    expect_equal(seq1, "ACTGACTG")

    seq2 <- gseq.extract(gintervals("test2", 0, 8))
    expect_equal(seq2, "GGGGCCCC")
})

test_that("multi-FASTA import validates chromosome sizes", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTG\n>chr2\nGG\n>chr3\nTATATATATATATATA\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    chroms <- gintervals.all()

    # Check sizes
    expect_equal(chroms[chroms$chrom == "chr1", "end"], 16)
    expect_equal(chroms[chroms$chrom == "chr2", "end"], 2)
    expect_equal(chroms[chroms$chrom == "chr3", "end"], 16)
})

test_that("multi-FASTA import handles multi-line sequences", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Write sequence across multiple lines (common FASTA format)
    cat(">multiline\n", file = test_fasta)
    cat("ACTGACTGACTGACTG\n", file = test_fasta, append = TRUE)
    cat("GGGGCCCCAAAATTTT\n", file = test_fasta, append = TRUE)
    cat("TATATATATATATATA\n", file = test_fasta, append = TRUE)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Extract and verify concatenated sequence
    full_seq <- gseq.extract(gintervals("multiline", 0, 48))
    expected <- "ACTGACTGACTGACTGGGGGCCCCAAAATTTTTATATATATATATATA"
    expect_equal(full_seq, expected)

    # Test partial extraction across original line boundaries
    partial <- gseq.extract(gintervals("multiline", 14, 34))
    expect_equal(partial, "TGGGGGCCCCAAAATTTTTA")
})

test_that("multi-FASTA import handles N characters and gaps", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">withN\nACTGNNNNACTG\n>withGap\nACTG---ACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # N characters should be preserved (12 chars: 4 + 4 + 4)
    seq_n <- gseq.extract(gintervals("withN", 0, 12))
    expect_equal(seq_n, "ACTGNNNNACTG")

    # Gap characters (-) should be preserved (11 chars: 4 + 3 + 4)
    seq_gap <- gseq.extract(gintervals("withGap", 0, 11))
    expect_equal(seq_gap, "ACTG---ACTG")
})

test_that("multi-FASTA import rejects invalid characters", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Include invalid character (digit)
    cat(">invalid\nACTG1234ACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)

    # Should error on invalid character
    expect_error(
        gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE),
        "Invalid character"
    )
})

test_that("multi-FASTA import handles empty contigs gracefully", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Include contig with no sequence (should be skipped or have size 0)
    cat(">contig1\nACTG\n>empty\n>contig2\nGGGG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)

    # This may error or create a 0-length contig depending on implementation
    # At minimum, it shouldn't crash
    tryCatch(
        {
            gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
            gdb.init(test_db)
            chroms <- gintervals.all()
            # If empty contig is included, its size should be 0
            if ("empty" %in% chroms$chrom) {
                expect_equal(chroms[chroms$chrom == "empty", "end"], 0)
            }
        },
        error = function(e) {
            # Acceptable to error on empty contig
            expect_true(grepl("empty|invalid|zero", e$message, ignore.case = TRUE))
        }
    )
})

test_that("multi-FASTA import preserves exact sequence boundaries", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">boundary1\nAAAAAAAAAAAAAAAA\n>boundary2\nTTTTTTTTTTTTTTTT\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Extract from start
    start1 <- gseq.extract(gintervals("boundary1", 0, 4))
    expect_equal(start1, "AAAA")

    # Extract to end
    end1 <- gseq.extract(gintervals("boundary1", 12, 16))
    expect_equal(end1, "AAAA")

    # First base of second contig should be T, not A
    start2 <- gseq.extract(gintervals("boundary2", 0, 1))
    expect_equal(start2, "T")

    # Verify no cross-contamination
    full1 <- gseq.extract(gintervals("boundary1", 0, 16))
    full2 <- gseq.extract(gintervals("boundary2", 0, 16))
    expect_equal(full1, "AAAAAAAAAAAAAAAA")
    expect_equal(full2, "TTTTTTTTTTTTTTTT")
})

test_that("multi-FASTA import index checksum is valid", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTGACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)

    idx_path <- file.path(test_db, "seq", "genome.idx")
    expect_true(file.exists(idx_path))

    # Corrupt the index file
    idx_data <- readBin(idx_path, "raw", n = 1000)
    # Change a byte in the contig entry data - specifically in the offset/length fields
    # which won't break the read but will cause checksum mismatch
    # For a contig "test" (4 chars): header=24, chromid=4, name_len=2, name=4 bytes
    # So offset starts at byte 34 (24+4+2+4)
    if (length(idx_data) > 40) {
        # Corrupt byte in offset field (bytes 34-41)
        idx_data[35] <- as.raw(bitwXor(as.integer(idx_data[35]), 0xFF))
        writeBin(idx_data, idx_path)

        # Should fail to load with checksum error
        expect_error(gdb.init(test_db), "checksum|corrupt|mismatch")
    }
})

test_that("multi-FASTA import works with various contig name formats", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Test various naming conventions
    cat(">1\nACTG\n", file = test_fasta) # Numeric (Ensembl style)
    cat(">chr2\nGGGG\n", file = test_fasta, append = TRUE) # chr prefix
    cat(">MT\nCCCC\n", file = test_fasta, append = TRUE) # Mitochondrial
    cat(">scaffold_123\nTTTT\n", file = test_fasta, append = TRUE) # Scaffold
    cat(">contig.456\nAAAA\n", file = test_fasta, append = TRUE) # Dot notation
    cat(">X\nGGGG\n", file = test_fasta, append = TRUE) # Sex chromosome

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    chroms <- gintervals.all()

    # All names should be preserved
    expect_true("1" %in% chroms$chrom)
    expect_true("chr2" %in% chroms$chrom)
    expect_true("MT" %in% chroms$chrom)
    expect_true("scaffold_123" %in% chroms$chrom)
    expect_true("contig.456" %in% chroms$chrom)
    expect_true("X" %in% chroms$chrom)

    # All sequences should be accessible
    expect_equal(gseq.extract(gintervals("1", 0, 4)), "ACTG")
    expect_equal(gseq.extract(gintervals("chr2", 0, 4)), "GGGG")
    expect_equal(gseq.extract(gintervals("MT", 0, 4)), "CCCC")
})

test_that("multi-FASTA import handles large contigs efficiently", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Create a larger contig (1KB)
    large_seq <- paste(rep("ACTG", 250), collapse = "")
    cat(">large\n", file = test_fasta)
    # Write in chunks to simulate real FASTA
    for (i in 1:10) {
        cat(substr(large_seq, (i - 1) * 100 + 1, i * 100), "\n",
            file = test_fasta, append = TRUE
        )
    }

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Verify full sequence
    full <- gseq.extract(gintervals("large", 0, 1000))
    expect_equal(nchar(full), 1000)
    expect_equal(full, large_seq)

    # Test random access at different positions
    chunk1 <- gseq.extract(gintervals("large", 0, 100))
    expect_equal(chunk1, substr(large_seq, 1, 100))

    chunk2 <- gseq.extract(gintervals("large", 500, 600))
    expect_equal(chunk2, substr(large_seq, 501, 600))

    chunk3 <- gseq.extract(gintervals("large", 900, 1000))
    expect_equal(chunk3, substr(large_seq, 901, 1000))
})

test_that("legacy format option disables indexed format", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTGACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    # Disable indexed format
    options(gmulticontig.indexed_format = FALSE)

    # Should use per-chromosome format
    expect_message(
        gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE),
        "per-chromosome"
    )

    # Indexed files should NOT exist
    expect_false(file.exists(file.path(test_db, "seq", "genome.idx")))
    expect_false(file.exists(file.path(test_db, "seq", "genome.seq")))
})

test_that("multi-FASTA with multiple files falls back to legacy", {
    test_fasta1 <- tempfile(fileext = ".fasta")
    test_fasta2 <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTG\n", file = test_fasta1)
    cat(">chr2\nGGGG\n", file = test_fasta2)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta1)
        unlink(test_fasta2)
    })

    options(gmulticontig.indexed_format = TRUE)
    # Multiple files should trigger multi-chromosome mode
    expect_message(
        gdb.create(groot = test_db, fasta = c(test_fasta1, test_fasta2), verbose = TRUE),
        "per-chromosome"
    )
})

test_that("multi-FASTA import handles very long contig names", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Create a contig with a very long name (100+ characters)
    long_name <- paste(rep("contig", 30), collapse = "_")
    cat(sprintf(">%s\nACTG\n", long_name), file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Verify the contig was imported (name may be truncated)
    chroms <- gintervals.all()
    expect_true(nrow(chroms) > 0)
})

test_that("multi-FASTA import handles duplicate contig names", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Create FASTA with duplicate names (should error or handle gracefully)
    cat(">dup\nACTG\n>dup\nGGGG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)

    # Should either error or handle duplicates
    # Most implementations would error on duplicate names
    expect_error(
        gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE),
        "duplicate|unique"
    )
})

test_that("multi-FASTA import handles mixed case sequences", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nActGacTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Check that sequence is imported (case may or may not be preserved)
    seq <- gseq.extract(gintervals("test", 0, 8))
    expect_equal(toupper(seq), "ACTGACTG")
})

test_that("multi-FASTA import handles whitespace in headers", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">contig with spaces in header\nACTG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    chroms <- gintervals.all()
    # Name should be truncated at first space
    expect_true("contig" %in% chroms$chrom)
})

test_that("multi-FASTA import creates correct chromosome order", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Create contigs in non-alphabetical order
    cat(">zebra\nACTG\n>apple\nGGGG\n>middle\nCCCC\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Check that chrom_sizes.txt is alphabetically sorted
    chrom_sizes <- read.table(file.path(test_db, "chrom_sizes.txt"),
        header = FALSE, stringsAsFactors = FALSE
    )
    expect_equal(chrom_sizes$V1, c("apple", "middle", "zebra"))
})

test_that("multi-FASTA import handles ambiguous IUPAC codes", {
    test_fasta <- tempfile(fileext = ".fasta")
    # R=A/G, Y=C/T, W=A/T, S=G/C, K=G/T, M=A/C
    cat(">test\nRYWSKM\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    # Should preserve IUPAC codes
    seq <- gseq.extract(gintervals("test", 0, 6))
    expect_equal(seq, "RYWSKM")
})

test_that("multi-FASTA import works with read-only source file", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTG\n", file = test_fasta)

    # Make read-only
    Sys.chmod(test_fasta, mode = "0444")

    test_db <- tempfile()
    withr::defer({
        Sys.chmod(test_fasta, mode = "0644") # Restore for cleanup
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    # Should still work with read-only input (may produce messages)
    expect_no_error(gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE))
})

test_that("multi-FASTA import index file has correct structure", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTG\n>chr2\nGGGG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)

    # Check index file exists and has reasonable size
    idx_file <- file.path(test_db, "seq", "genome.idx")
    expect_true(file.exists(idx_file))

    # Index should be larger than header (24 bytes) + at least 2 entries
    idx_size <- file.info(idx_file)$size
    expect_true(idx_size > 50)
})

test_that("multi-FASTA import sequence file concatenates correctly", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">a\nAAAA\n>b\nCCCC\n>c\nGGGG\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)

    # Check sequence file has total length
    seq_file <- file.path(test_db, "seq", "genome.seq")
    seq_size <- file.info(seq_file)$size
    expect_equal(seq_size, 12) # 4 + 4 + 4 bytes
})

test_that("multi-FASTA import handles FASTA with trailing newlines", {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">test\nACTG\n\n\n", file = test_fasta)

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    seq <- gseq.extract(gintervals("test", 0, 4))
    expect_equal(seq, "ACTG")
    expect_equal(nrow(gintervals.all()), 1)
})

test_that("multi-FASTA import handles Windows line endings", {
    test_fasta <- tempfile(fileext = ".fasta")
    # Write with \r\n line endings
    writeLines(c(">test", "ACTG"), test_fasta, sep = "\r\n")

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    options(gmulticontig.indexed_format = TRUE)
    gdb.create(groot = test_db, fasta = test_fasta, verbose = TRUE)
    gdb.init(test_db)

    seq <- gseq.extract(gintervals("test", 0, 4))
    expect_equal(seq, "ACTG")
})

# Restore the test database after all multifasta-import tests
# This ensures subsequent test files have the correct database set
suppressMessages(gdb.init("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"))
