restore_groot_on_exit()

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

test_that("gdb.convert_to_indexed adds the chr prefix and keeps the chrom order gsetroot gives the database", {
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

    # This is a per-chromosome database (unprefixed chrom_sizes.txt, prefixed .seq files), so
    # gsetroot() names its chromosomes with the chr prefix and gives them chrom ids in the
    # order of the sorted names. The converted database keeps that order, not the
    # chrom_sizes.txt one, so that existing indexed tracks keep their chrom ids.
    expect_equal(converted_chrom_sizes$chrom, c("chr1", "chr10", "chr15", "chr17_random"))
    expect_equal(converted_chrom_sizes$size, c(20, 10, 15, 17))
})

test_that("gdb.convert_to_indexed keeps the chrom_sizes.txt names and order when the .seq files lack the chr prefix", {
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

    # The chromosomes keep the names gsetroot() gives them (those of chrom_sizes.txt), as they
    # do when the database is loaded; the .seq files are found without the prefix
    expect_equal(converted_chrom_sizes$chrom, original_chrom_sizes$chrom)
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
        gdb.convert_to_indexed(groot = test_db, force = TRUE, validate = TRUE, verbose = TRUE, convert_tracks = FALSE, convert_intervals = FALSE),
        "already in indexed format"
    )
})

test_that("gdb.convert_to_indexed keeps the chrom ids of a per-chromosome database whether or not it is loaded", {
    local_db_state()
    td <- tempfile("convert_order_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))

    # chrom_sizes.txt is unsorted, so gsetroot() gives the chromosomes the ids of the sorted
    # names. An indexed track made now is keyed by those ids and must stay readable.
    loaded_db <- create_db_with_unsorted_chrom_sizes(file.path(td, "loaded"))
    unloaded_db <- create_db_with_unsorted_chrom_sizes(file.path(td, "unloaded"))
    other_db <- create_test_db(file.path(td, "other"))
    snapshot <- function() {
        list(
            chroms = gintervals.all(),
            sp = gextract("sp", gintervals.all()),
            seq = gseq.extract(gintervals(gintervals.all()$chrom, 0, 5))
        )
    }
    before <- list()
    for (db in c(loaded_db, unloaded_db)) {
        gsetroot(db)
        cs <- read.table(file.path(db, "chrom_sizes.txt"), colClasses = "character")
        expect_false(identical(as.character(gintervals.all()$chrom), paste0("chr", cs$V1)))
        intervs <- gintervals.all()
        intervs$end <- 100
        gtrack.create_sparse("sp", "x", intervs, seq_len(nrow(intervs)))
        gtrack.convert_to_indexed("sp")
        before[[db]] <- snapshot()
    }

    gsetroot(loaded_db)
    suppressMessages(gdb.convert_to_indexed(force = TRUE))
    gsetroot(other_db)
    suppressMessages(gdb.convert_to_indexed(groot = unloaded_db, force = TRUE))

    for (db in c(loaded_db, unloaded_db)) {
        expect_true(misha:::.gdb.is_indexed_at(db))
        gsetroot(db)
        expect_equal(snapshot(), before[[db]])
    }
    # the same database whichever database was loaded
    for (f in c("chrom_sizes.txt", "seq/genome.idx", "seq/genome.seq")) {
        read_bytes <- function(db) readBin(file.path(db, f), "raw", n = file.size(file.path(db, f)))
        expect_identical(read_bytes(unloaded_db), read_bytes(loaded_db), info = f)
    }
})

test_that("a genome.idx left without genome.seq does not make a per-chromosome db count as indexed", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    ref <- create_db_with_unsorted_chrom_sizes(file.path(td, "ref"))
    file.create(file.path(db, "seq", "genome.idx"))
    cs <- utils::read.csv(file.path(db, "chrom_sizes.txt"),
        sep = "\t", header = FALSE, col.names = c("chrom", "size"), colClasses = c("character", "numeric")
    )
    expect_true(misha:::.gdb.chrom_order(db, cs)$per_chromosome)

    # it converts to the order a clean copy of the database converts to
    suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))
    suppressMessages(gdb.convert_to_indexed(groot = ref, force = TRUE, validate = FALSE))
    expect_equal(readLines(file.path(db, "chrom_sizes.txt")), readLines(file.path(ref, "chrom_sizes.txt")))
})

test_that("a conversion killed right after the sequence import leaves a database the next run converts", {
    skip_on_cran()
    skip_on_os("windows")
    skip_if_not_installed("callr")
    skip_if_not_installed("pkgload")
    root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
    skip_if_not(file.exists(file.path(root, "DESCRIPTION")), "needs the package source")
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))

    # SIGKILL right after gseq_multifasta_import, as a killed job: no R error handler runs
    try(callr::r(function(root, db) {
        pkgload::load_all(root, compile = FALSE, quiet = TRUE)
        gcall <- misha:::.gcall
        utils::assignInNamespace(".gcall", function(...) {
            res <- gcall(...)
            if (identical(..1, "gseq_multifasta_import")) {
                tools::pskill(Sys.getpid(), tools::SIGKILL)
            }
            res
        }, "misha")
        suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))
    }, args = list(root, db)), silent = TRUE)

    # not converted yet, so the next run converts it
    expect_false(file.exists(file.path(db, "seq", "genome.idx")))
    suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))
    expect_true(file.exists(file.path(db, "seq", "genome.idx")))
    gsetroot(db)
    expected <- c(chr1 = "G", chr1_KI270706v1_random = "N", chr10 = "C", chr2 = "A", chrX = "T")
    got <- vapply(names(expected), function(chrom) toupper(gseq.extract(gintervals(chrom, 0, 1))), character(1))
    expect_equal(got, expected)
})
