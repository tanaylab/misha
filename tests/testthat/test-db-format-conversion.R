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

# The bases of create_db_with_unsorted_chrom_sizes()'s chromosomes, by name
unsorted_db_bases <- c(chr1 = "G", chr1_KI270706v1_random = "N", chr10 = "C", chr2 = "A", chrX = "T")
first_bases <- function() {
    vapply(names(unsorted_db_bases), function(chrom) toupper(gseq.extract(gintervals(chrom, 0, 1))), character(1))
}

for (.kill_at in c("import", "rename")) {
    test_that(sprintf("a conversion killed at its %s leaves a database the next run converts", .kill_at), {
        skip_on_cran()
        skip_on_os("windows")
        local_db_state()
        td <- withr::local_tempdir()
        db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
        seq_dir <- file.path(db, "seq")
        original_chrom_sizes <- readLines(file.path(db, "chrom_sizes.txt"))
        # an index left by an earlier, interrupted conversion
        writeBin(as.raw(1:64), file.path(seq_dir, "genome.idx"))
        fastas <- list.files(tempdir(), pattern = "\\.fasta$", full.names = TRUE)

        # SIGKILL in a forked child, as a killed job: no R error handler runs
        kill_at <- .kill_at
        job <- parallel::mcparallel({
            if (kill_at == "import") {
                gcall <- misha:::.gcall
                utils::assignInNamespace(".gcall", function(...) {
                    res <- gcall(...)
                    if (identical(..1, "gseq_multifasta_import")) {
                        tools::pskill(Sys.getpid(), tools::SIGKILL)
                    }
                    res
                }, "misha")
            } else {
                # between replacing chrom_sizes.txt and moving genome.idx into place
                rename <- base::file.rename
                unlockBinding("file.rename", .BaseNamespaceEnv)
                assign("file.rename", function(from, to) {
                    if (basename(to) == "genome.idx") {
                        tools::pskill(Sys.getpid(), tools::SIGKILL)
                    }
                    rename(from, to)
                }, envir = .BaseNamespaceEnv)
            }
            suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))
        })
        # the child was killed: it delivers no result
        expect_warning(parallel::mccollect(job), "did not deliver a result")
        unlink(setdiff(list.files(tempdir(), pattern = "\\.fasta$", full.names = TRUE), fastas))

        # the child got as far as the import; genome.idx is not in place, the stale one included
        expect_true(all(c("genome.seq", "genome.idx.tmp") %in% list.files(seq_dir)))
        expect_false(file.exists(file.path(seq_dir, "genome.idx")))
        chrom_sizes <- readLines(file.path(db, "chrom_sizes.txt"))
        if (kill_at == "import") {
            expect_equal(chrom_sizes, original_chrom_sizes)
        } else {
            # chrom_sizes.txt is replaced before genome.idx is moved into place
            expect_equal(sub("\t.*", "", chrom_sizes), names(unsorted_db_bases))
        }

        suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))
        expect_true(file.exists(file.path(seq_dir, "genome.idx")))
        expect_false(file.exists(file.path(seq_dir, "genome.idx.tmp")))
        gsetroot(db)
        expect_equal(first_bases(), unsorted_db_bases)
    })
}

test_that("gdb.convert_to_indexed stops, leaving the database as it was, when the import does not match chrom_sizes.txt", {
    local_db_state()
    td <- withr::local_tempdir()
    state <- function(db) list(chrom_sizes = readLines(file.path(db, "chrom_sizes.txt")), seq = sort(list.files(file.path(db, "seq"))))

    # a .seq file shorter than chrom_sizes.txt says; its .seq files are kept although asked to go
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    writeBin(charToRaw(strrep("G", 900)), file.path(db, "seq", "chr1.seq"))
    before <- state(db)
    expect_error(
        suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE, remove_old_files = TRUE)),
        "seq/chr1.seq has 900 bytes and chrom_sizes.txt says 1000; they must agree before converting"
    )
    expect_equal(state(db), before)

    # a chromosome name the index cannot store as it is
    db2 <- file.path(td, "db2")
    dir.create(file.path(db2, "seq"), recursive = TRUE)
    dir.create(file.path(db2, "tracks"))
    writeBin(charToRaw(strrep("A", 100)), file.path(db2, "seq", "chr1.seq"))
    writeBin(charToRaw(strrep("C", 100)), file.path(db2, "seq", "HLA-A*01:01.seq"))
    writeLines(c("chr1\t100", "HLA-A*01:01\t100"), file.path(db2, "chrom_sizes.txt"))
    before <- state(db2)
    expect_error(suppressMessages(gdb.convert_to_indexed(groot = db2, force = TRUE, validate = FALSE)), "was imported as HLA-A_01_01")
    expect_equal(state(db2), before)
})

test_that("gdb.convert_to_indexed refuses a database whose seq/ is another database's", {
    local_db_state()
    td <- withr::local_tempdir()
    parent <- create_db_with_unsorted_chrom_sizes(file.path(td, "parent"))
    parent_state <- function() {
        list(seq = sort(list.files(file.path(parent, "seq"))), chrom_sizes = readLines(file.path(parent, "chrom_sizes.txt")))
    }
    before <- parent_state()
    refused <- "is in the database .*: convert that database instead"

    # a database made by gdb.create_linked(), its parent not loaded and loaded
    linked <- file.path(td, "linked")
    gdb.create_linked(linked, parent)
    expect_error(suppressMessages(gdb.convert_to_indexed(groot = linked, force = TRUE, validate = FALSE)), refused)
    gsetroot(parent)
    expect_error(suppressMessages(gdb.convert_to_indexed(groot = linked, force = TRUE, validate = FALSE)), refused)

    # a dataset saved without copy_seq
    gtrack.create_sparse("sp", "x", gintervals("chr1", 0, 10), 1)
    ds <- file.path(td, "ds")
    suppressMessages(gdataset.save(ds, "d", tracks = "sp"))
    expect_error(suppressMessages(gdb.convert_to_indexed(groot = ds, force = TRUE, validate = FALSE)), refused)
    expect_equal(parent_state(), before)

    # seq/ and chrom_sizes.txt linked to storage of their own elsewhere (a directory that is no
    # database: no tracks/, seq/ or chrom_sizes.txt but them) are converted there
    own <- create_db_with_unsorted_chrom_sizes(file.path(td, "own"))
    for (f in c("seq", "chrom_sizes.txt")) {
        storage <- file.path(td, paste0("storage_", f))
        dir.create(storage)
        expect_true(file.rename(file.path(own, f), file.path(storage, f)))
        expect_true(file.symlink(file.path(storage, f), file.path(own, f)))
    }
    suppressMessages(gdb.convert_to_indexed(groot = own, force = TRUE, validate = FALSE))
    expect_true(all(file.exists(file.path(td, "storage_seq", "seq", c("genome.idx", "genome.seq")))))
    expect_equal(readLines(file.path(td, "storage_chrom_sizes.txt", "chrom_sizes.txt"))[1], "chr1\t1000")
    gsetroot(own)
    expect_equal(first_bases(), unsorted_db_bases)
})

test_that("gdb.convert_to_indexed stops for a chromosome whose sequence file is seq/Genome.seq", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- file.path(td, "db")
    dir.create(file.path(db, "seq"), recursive = TRUE)
    dir.create(file.path(db, "tracks"))
    writeBin(charToRaw(strrep("A", 100)), file.path(db, "seq", "chr1.seq"))
    writeBin(charToRaw(strrep("C", 50)), file.path(db, "seq", "Genome.seq"))
    writeLines(c("chr1\t100", "Genome\t50"), file.path(db, "chrom_sizes.txt"))
    # seq/genome.seq on a file system that ignores case
    expect_error(
        suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE)),
        "seq/genome.seq, the name of the indexed format"
    )
    expect_false(file.exists(file.path(db, "seq", "genome.idx")))
})

test_that("gdb.convert_to_indexed stops for a chromosome whose sequence file is seq/genome.seq", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- file.path(td, "db")
    dir.create(file.path(db, "seq"), recursive = TRUE)
    dir.create(file.path(db, "tracks"))
    writeBin(charToRaw(strrep("A", 100)), file.path(db, "seq", "chr1.seq"))
    writeBin(charToRaw(strrep("C", 50)), file.path(db, "seq", "genome.seq"))
    writeLines(c("chr1\t100", "genome\t50"), file.path(db, "chrom_sizes.txt"))
    before <- readBin(file.path(db, "seq", "genome.seq"), "raw", 1000)
    expect_error(
        suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE, remove_old_files = TRUE)),
        "seq/genome.seq, the name of the indexed format"
    )
    expect_identical(readBin(file.path(db, "seq", "genome.seq"), "raw", 1000), before)
    expect_false(file.exists(file.path(db, "seq", "genome.idx")))
})

test_that("gdb.convert_to_indexed stops on a failed write, leaving the database as it was", {
    skip_if_not(file.exists("/dev/full"), "no /dev/full")
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    cs <- file.path(db, "chrom_sizes.txt")
    before <- list(chrom_sizes = readBin(cs, "raw", 10000), seq = sort(list.files(file.path(db, "seq"))))
    # writing genome.idx.tmp fails as on a full disk; only closing it shows it, and validate = FALSE
    # reads nothing back
    gcall <- misha:::.gcall
    local_mocked_bindings(.gcall = function(...) {
        if (identical(..1, "gseq_multifasta_import")) {
            file.symlink("/dev/full", ..4)
        }
        gcall(...)
    }, .package = "misha")
    expect_error(
        suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE)),
        "genome.idx.tmp: No space left on device"
    )
    expect_identical(list(chrom_sizes = readBin(cs, "raw", 10000), seq = sort(list.files(file.path(db, "seq")))), before)
})

test_that("gdb.convert_to_indexed stops when writing chrom_sizes.txt fails, leaving the database as it was", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    cs <- file.path(db, "chrom_sizes.txt")
    before <- list(chrom_sizes = readBin(cs, "raw", 10000), seq = sort(list.files(file.path(db, "seq"))))
    # R reports a failed write as a warning, as here from write.table() on a full disk
    trace("write.table", quote(if (grepl("chrom_sizes", file)) warning("No space left on device")), print = FALSE, where = asNamespace("misha"))
    withr::defer(untrace("write.table", where = asNamespace("misha")))
    expect_error(
        suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE)),
        "Conversion failed: No space left on device"
    )
    expect_identical(list(chrom_sizes = readBin(cs, "raw", 10000), seq = sort(list.files(file.path(db, "seq")))), before)
    expect_equal(list.files(db), c("chrom_sizes.txt", "seq", "tracks"))
})

test_that("gdb.convert_to_indexed stops, with chrom_sizes.txt as it was, when its backup copy is short", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    cs <- file.path(db, "chrom_sizes.txt")
    before <- readBin(cs, "raw", 10000)
    # hard links refused, and the copy comes out short
    trace("file.link", quote(to <- file.path(tempdir(), "no", "such", "dir", "x")), print = FALSE, where = baseenv())
    withr::defer(untrace("file.link", where = baseenv()))
    trace("file.copy", exit = quote(if (any(grepl("chrom_sizes", to))) writeBin(as.raw(1), to)), print = FALSE, where = baseenv())
    withr::defer(untrace("file.copy", where = baseenv()))
    expect_error(suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE)), "Failed to keep a copy of")
    expect_identical(readBin(cs, "raw", 10000), before)
    expect_equal(list.files(db), c("chrom_sizes.txt", "seq", "tracks"))
    expect_false(file.exists(file.path(db, "seq", "genome.idx")))
})

for (.interrupt_at in c("backup copy", "import", "chrom_sizes.txt rename")) {
    test_that(sprintf("gdb.convert_to_indexed interrupted at the %s leaves no temporary file next to chrom_sizes.txt", .interrupt_at), {
        local_db_state()
        td <- withr::local_tempdir()
        db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
        before <- readBin(file.path(db, "chrom_sizes.txt"), "raw", 10000)
        # Ctrl-C
        interrupt <- quote(stop(structure(class = c("interrupt", "condition"), list(message = "", call = NULL))))
        if (.interrupt_at == "backup copy") {
            # hard links refused, and the copy interrupted once made
            trace("file.link", quote(to <- file.path(tempdir(), "no", "such", "dir", "x")), print = FALSE, where = baseenv())
            withr::defer(untrace("file.link", where = baseenv()))
            trace("file.copy", exit = bquote(if (any(grepl("chrom_sizes", to))) .(interrupt)), print = FALSE, where = baseenv())
            withr::defer(untrace("file.copy", where = baseenv()))
        } else if (.interrupt_at == "import") {
            gcall <- misha:::.gcall
            local_mocked_bindings(.gcall = function(...) {
                res <- gcall(...)
                if (identical(..1, "gseq_multifasta_import")) eval(interrupt)
                res
            }, .package = "misha")
        } else {
            trace("file.rename", bquote(if (any(grepl("chrom_sizes.txt[.]", from))) .(interrupt)), print = FALSE, where = baseenv())
            withr::defer(untrace("file.rename", where = baseenv()))
        }
        expect_equal(
            tryCatch(suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE)), interrupt = function(i) "interrupted"),
            "interrupted"
        )
        expect_equal(list.files(db), c("chrom_sizes.txt", "seq", "tracks"))
        expect_identical(readBin(file.path(db, "chrom_sizes.txt"), "raw", 10000), before)
    })
}

test_that("a warning that is not from a write does not stop gdb.convert_to_indexed", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    gwith_umask <- misha:::.gwith_umask
    # as a finalizer's, raised while chrom_sizes.txt is written
    local_mocked_bindings(.gwith_umask = function(expr) {
        warning("closing unused connection 3")
        gwith_umask(expr)
    }, .package = "misha")
    expect_warning(
        suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE)),
        "closing unused connection"
    )
    expect_true(all(file.exists(file.path(db, "seq", c("genome.idx", "genome.seq")))))
    gsetroot(db)
    expect_equal(first_bases(), unsorted_db_bases)
})

test_that("gdb.convert_to_indexed stops when genome.seq is not the size of chrom_sizes.txt", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    cs <- file.path(db, "chrom_sizes.txt")
    before <- readBin(cs, "raw", 10000)
    gcall <- misha:::.gcall
    local_mocked_bindings(.gcall = function(...) {
        res <- gcall(...)
        if (identical(..1, "gseq_multifasta_import")) {
            # genome.seq cut short after the import
            genome_seq <- file.path(db, "seq", "genome.seq")
            writeBin(readBin(genome_seq, "raw", 100), genome_seq)
        }
        res
    }, .package = "misha")
    expect_error(suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE)), "genome.seq has 100 bytes")
    expect_identical(readBin(cs, "raw", 10000), before)
    expect_false(any(c("genome.seq", "genome.idx", "genome.idx.tmp") %in% list.files(file.path(db, "seq"))))
})

for (.backup in c("hard link", "copy")) {
    test_that(sprintf("a conversion that fails after replacing chrom_sizes.txt puts back the same bytes and mode (%s)", .backup), {
        local_db_state()
        td <- withr::local_tempdir()
        db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
        cs <- file.path(db, "chrom_sizes.txt")
        # no newline at the end, and group-writable
        writeBin(charToRaw("2\t2000\n10\t1500\n1\t1000\nX\t1200\n1_KI270706v1_random\t500"), cs)
        Sys.chmod(cs, "664", use_umask = FALSE)
        before <- list(bytes = readBin(cs, "raw", 10000), mode = file.info(cs)$mode)
        if (.backup == "copy") {
            # hard links refused, as for another user's file under fs.protected_hardlinks
            trace("file.link", quote(to <- file.path(tempdir(), "no", "such", "dir", "x")), print = FALSE, where = baseenv())
            withr::defer(untrace("file.link", where = baseenv()))
        }
        # genome.idx cannot be moved into place: a directory has its name
        dir.create(file.path(db, "seq", "genome.idx"))
        old_umask <- Sys.umask("077")
        withr::defer(Sys.umask(old_umask))
        expect_error(
            suppressWarnings(suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))),
            "genome.idx into place"
        )
        expect_identical(list(bytes = readBin(cs, "raw", 10000), mode = file.info(cs)$mode), before)
        expect_equal(list.files(db), c("chrom_sizes.txt", "seq", "tracks"))
    })
}

test_that("gdb.convert_to_indexed keeps the group of chrom_sizes.txt", {
    skip_on_os("windows")
    local_db_state()
    groups <- tryCatch(system2("id", "-G", stdout = TRUE), error = function(e) character(0), warning = function(w) character(0))
    primary <- tryCatch(system2("id", "-g", stdout = TRUE), error = function(e) "", warning = function(w) "")
    other <- setdiff(strsplit(paste(groups, collapse = " "), " ")[[1]], c(primary, ""))
    skip_if(length(other) == 0, "the test user belongs to no group besides the primary one")
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    cs <- file.path(db, "chrom_sizes.txt")
    skip_if(system2("chgrp", c(other[1], shQuote(cs))) != 0, "chgrp failed")
    gid <- file.info(cs, extra_cols = TRUE)$gid
    expect_equal(as.character(gid), other[1])

    # a failed conversion puts back the backup, here a copy (hard links refused)
    bytes <- readBin(cs, "raw", 10000)
    trace("file.link", quote(to <- file.path(tempdir(), "no", "such", "dir", "x")), print = FALSE, where = baseenv())
    dir.create(file.path(db, "seq", "genome.idx"))
    expect_error(suppressWarnings(suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))), "genome.idx into place")
    untrace("file.link", where = baseenv())
    unlink(file.path(db, "seq", "genome.idx"), recursive = TRUE)
    expect_equal(file.info(cs, extra_cols = TRUE)$gid, gid)
    expect_identical(readBin(cs, "raw", 10000), bytes)

    suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))
    expect_true(file.exists(file.path(db, "seq", "genome.idx")))
    expect_equal(file.info(cs, extra_cols = TRUE)$gid, gid)
})

test_that("gdb.convert_to_indexed leaves the user's files next to chrom_sizes.txt alone", {
    local_db_state()
    td <- withr::local_tempdir()
    mine <- c(chrom_sizes.txt.orig = "my backup", chrom_sizes.txt.tmp = "my notes")
    for (fail in c(TRUE, FALSE)) {
        db <- create_db_with_unsorted_chrom_sizes(file.path(td, paste0("db", fail)))
        for (f in names(mine)) writeLines(mine[[f]], file.path(db, f))
        if (fail) {
            dir.create(file.path(db, "seq", "genome.idx"))
        }
        r <- tryCatch(suppressWarnings(suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))), error = function(e) "failed")
        expect_equal(identical(r, "failed"), fail)
        for (f in names(mine)) expect_equal(readLines(file.path(db, f)), mine[[f]])
        expect_equal(list.files(db), sort(c("chrom_sizes.txt", names(mine), "seq", "tracks")))
    }
})

test_that("gsetroot stops when seq/genome.idx does not match chrom_sizes.txt, keeping the loaded database", {
    local_db_state()
    td <- withr::local_tempdir()
    parent <- create_db_with_unsorted_chrom_sizes(file.path(td, "parent"))
    # a database made from parent with its own chrom_sizes.txt and a link to its seq/
    child <- file.path(td, "child")
    dir.create(file.path(child, "tracks"), recursive = TRUE)
    expect_true(file.copy(file.path(parent, "chrom_sizes.txt"), child))
    expect_true(file.symlink(file.path(parent, "seq"), file.path(child, "seq")))

    # parent converted afterwards, with child loaded: its index follows its own, rewritten
    # chrom_sizes.txt, so child would not load again. Its session is kept, as its chrom ids are
    # still those of the index.
    gsetroot(child)
    expect_warning(
        suppressMessages(gdb.convert_to_indexed(groot = parent, force = TRUE)),
        "child is loaded as before the conversion, but gsetroot\\(\\) would not load it now: .*does not match chrom_sizes.txt.*copy .*parent/chrom_sizes.txt into .*child"
    )
    expect_equal(get("GROOT", envir = misha:::.misha), normalizePath(child))
    expect_equal(first_bases(), unsorted_db_bases)
    gtrack.create_sparse("written", "x", gintervals("chr1", 0, 10), 1)
    expect_true(dir.exists(file.path(child, "tracks", "written.track")))
    expect_false(dir.exists(file.path(parent, "tracks", "written.track")))

    expect_error(gsetroot(child), "does not match chrom_sizes.txt.*copy .*parent/chrom_sizes.txt into")
    expect_equal(get("GROOT", envir = misha:::.misha), normalizePath(child))

    # copying the converted parent's chrom_sizes.txt fixes it
    expect_true(file.copy(file.path(parent, "chrom_sizes.txt"), child, overwrite = TRUE))
    gsetroot(child)
    expect_equal(first_bases(), unsorted_db_bases)
})

test_that("a database that does not load after a conversion and no longer reads its own sequence is unloaded", {
    local_db_state()
    td <- withr::local_tempdir()
    parent <- create_db_with_unsorted_chrom_sizes(file.path(td, "parent"))
    # a database on parent's seq/ with prefixed names in chrom_sizes.txt order: its chrom ids are
    # not those of the sorted index the conversion writes
    child <- file.path(td, "child")
    dir.create(file.path(child, "tracks"), recursive = TRUE)
    writeLines(sub("^", "chr", readLines(file.path(parent, "chrom_sizes.txt"))), file.path(child, "chrom_sizes.txt"))
    expect_true(file.symlink(file.path(parent, "seq"), file.path(child, "seq")))
    gsetroot(child)
    expect_equal(as.character(gintervals.all()$chrom)[1], "chr2")
    expect_warning(
        suppressMessages(gdb.convert_to_indexed(groot = parent, force = TRUE)),
        "child, loaded before the conversion, no longer reads its own sequence: .*No database is loaded"
    )
    expect_null(get0("GROOT", envir = misha:::.misha))
})

for (.how in c("converted", "failed in the genome step", "failed in the track step", "converted under warn = 2")) {
    test_that(sprintf("gdb.convert_to_indexed of a loaded dataset (%s) keeps the rest of the session", .how), {
        local_db_state()
        td <- withr::local_tempdir()
        # X the root, Y a database loaded as a dataset of it (the same chrom_sizes.txt), Z another dataset
        x <- create_db_with_unsorted_chrom_sizes(file.path(td, "X"))
        y <- create_db_with_unsorted_chrom_sizes(file.path(td, "Y"))
        iv <- gintervals(c("chr1", "chr10", "chr2", "chrX"), 0, 10)
        gsetroot(y)
        gtrack.create_sparse("yt", "x", iv, c(1, 10, 2, 99))
        gsetroot(x)
        gtrack.create_sparse("zt", "x", gintervals("chr1", 0, 10), 7)
        z <- file.path(td, "Z")
        suppressMessages(gdataset.save(z, "z", tracks = "zt"))
        gtrack.rm("zt", force = TRUE)
        suppressMessages(gdataset.load(y))
        suppressMessages(gdataset.load(z))
        gvtrack.create("v", "zt", "max")
        expect_equal(gextract("yt", iv)$yt, c(1, 10, 2, 99))
        convert <- function() suppressMessages(gdb.convert_to_indexed(groot = y, force = TRUE, validate = FALSE, convert_tracks = TRUE))

        if (.how == "converted") {
            # Y's tracks are rewritten in its own chrom ids: Y is unloaded
            expect_warning(convert(), "Y was converted, so it is no longer loaded as a dataset\\.$")
        } else if (.how == "failed in the genome step") {
            # nothing changed: Y stays
            writeBin(charToRaw(strrep("A", 10)), file.path(y, "seq", "chrX.seq"))
            expect_error(convert(), "seq/chrX.seq has 10 bytes")
        } else if (.how == "failed in the track step") {
            local_mocked_bindings(.gdb.convert_to_indexed.tracks = function(...) stop("injected"), .package = "misha")
            expect_warning(expect_error(convert(), "injected"), "The conversion of .*Y failed after changing it, so it is no longer loaded as a dataset")
        } else {
            # the warning is an error: the session is put back before it
            withr::local_options(warn = 2)
            expect_error(convert(), "Y was converted, so it is no longer loaded as a dataset")
        }
        unloaded <- .how != "failed in the genome step"
        expect_equal(get("GROOT", envir = misha:::.misha), normalizePath(x))
        expect_equal(get("GDATASETS", envir = misha:::.misha), if (unloaded) normalizePath(z) else normalizePath(c(y, z)))
        expect_equal(gvtrack.ls(), "v")
        expect_equal(gextract("v", gintervals("chr1", 0, 10))$v, 7)
        if (unloaded) {
            expect_error(gextract("yt", iv))
        } else {
            expect_equal(gextract("yt", iv)$yt, c(1, 10, 2, 99))
        }
    })
}

test_that("gdb.convert_to_indexed of the loaded database leaves its session in the indexed format", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    gsetroot(db)
    expect_true(get("DB_IS_PER_CHROMOSOME", envir = misha:::.misha))
    suppressMessages(gdb.convert_to_indexed(force = TRUE, convert_tracks = TRUE))
    expect_false(get("DB_IS_PER_CHROMOSOME", envir = misha:::.misha))
    expect_equal(first_bases(), unsorted_db_bases)
})

for (.steps in c("validate", "tracks and interval sets")) {
    test_that(sprintf("gdb.convert_to_indexed of another database puts the session back as it was (%s)", .steps), {
        local_db_state()
        td <- withr::local_tempdir()
        db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
        other <- create_db_with_unsorted_chrom_sizes(file.path(td, "other"))
        gsetroot(other)
        gtrack.create_sparse("o", "x", gintervals("chr1", 0, 10), 1)
        gintervals.save("oi", gintervals("chr2", 0, 10))
        # a session with a dataset, another working directory and a virtual track
        gsetroot(db)
        gtrack.create_sparse("t", "x", gintervals("chr1", 0, 10), 1)
        ds <- file.path(td, "ds")
        suppressMessages(gdataset.save(ds, "d", tracks = "t", copy_seq = TRUE))
        gtrack.rm("t", force = TRUE)
        suppressMessages(gdataset.load(ds))
        gvtrack.create("v", "t", "avg")
        gdir.create("sub", showWarnings = FALSE)
        gdir.cd("sub")
        before <- mget(c("GROOT", "GWD", "GDATASETS", "GTRACKS", "GVTRACKS"), envir = misha:::.misha)
        if (.steps == "validate") {
            suppressMessages(gdb.convert_to_indexed(groot = other, force = TRUE))
        } else {
            suppressMessages(gdb.convert_to_indexed(groot = other, force = TRUE, validate = FALSE, convert_tracks = TRUE, convert_intervals = TRUE))
            expect_true(file.exists(file.path(other, "tracks", "o.track", "track.idx")))
        }
        expect_equal(mget(c("GROOT", "GWD", "GDATASETS", "GTRACKS", "GVTRACKS"), envir = misha:::.misha), before)
        gdir.cd("..")
        expect_equal(gextract("t", gintervals("chr1", 0, 10))$t, 1)
        expect_equal(first_bases(), unsorted_db_bases)
    })
}

test_that("gsetroot stops when seq/genome.idx has the chromosomes of chrom_sizes.txt in another order", {
    local_db_state()
    td <- withr::local_tempdir()
    # equal sizes: only the names tell the order
    parent <- file.path(td, "parent")
    dir.create(file.path(parent, "tracks"), recursive = TRUE)
    dir.create(file.path(parent, "seq"))
    for (i in 1:3) writeBin(charToRaw(strrep(c("A", "C", "G")[i], 1000)), file.path(parent, "seq", paste0("chr", c("2", "10", "1")[i], ".seq")))
    writeLines(paste(c("2", "10", "1"), 1000, sep = "\t"), file.path(parent, "chrom_sizes.txt"))
    gsetroot(parent)
    gtrack.create_sparse("s", "x", gintervals("chr1", 0, 10), 1)
    ds <- file.path(td, "ds")
    suppressMessages(gdataset.save(ds, "d", tracks = "s"))
    gdb.unload()
    suppressMessages(gdb.convert_to_indexed(groot = parent))
    expect_error(gsetroot(ds), "2 is chrom id 0 in chrom_sizes.txt and 2 in the index")
})

test_that("gsetroot loads an indexed database whose chrom_sizes.txt lists only the index's first contigs", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))
    cs <- file.path(db, "chrom_sizes.txt")
    full <- readLines(cs)
    # chr1 1000, chr1_KI270706v1_random 500, chr10 1500, chr2 2000, chrX 1200
    writeLines(full[1:3], cs)
    gsetroot(db)
    expect_equal(toupper(gseq.extract(gintervals(names(unsorted_db_bases)[1:3], 0, 1))), unname(unsorted_db_bases[1:3]))
    gdb.unload()

    stops <- list(
        "a contig left out of the middle" = full[c(1, 3, 4)],
        "a size changed" = c(full[1], "chr1_KI270706v1_random\t501"),
        "more contigs than the index" = c(full, "chrY\t100")
    )
    for (what in names(stops)) {
        writeLines(stops[[what]], cs)
        expect_error(gsetroot(db), "does not match chrom_sizes.txt", info = what)
    }
    # a name that differs among the first ones
    writeLines(c(full[1], "chrY\t500"), cs)
    expect_error(gsetroot(db), "lists 2 of the index's 5 contigs, which have to be its first ones")
})

test_that("gsetroot stops for an indexed database whose chrom_sizes.txt lists a chromosome twice", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))
    cs <- file.path(db, "chrom_sizes.txt")
    full <- readLines(cs)
    # chr1 1000, chr1_KI270706v1_random 500, chr10 1500, chr2 2000, chrX 1200
    # chr1 again without the prefix, in chrX's place
    writeLines(c(full[1:4], "1\t1200"), cs)
    expect_error(gsetroot(db), "lists the same chromosome twice: chr1 \\(chrom id 0\\) and 1 \\(chrom id 4\\)")
    writeLines(c(full[1:3], "chr10\t2000", full[5]), cs)
    expect_error(gsetroot(db), "lists the same chromosome twice: chr10 \\(chrom id 2\\) and chr10 \\(chrom id 3\\)")
    # a line with no size
    writeLines(c(full[1:4], "chrX"), cs)
    expect_error(gsetroot(db), "chrom id 4 is chrX \\(1200 bp\\) in the index and chrX \\(with no size\\) in chrom_sizes.txt")
})

test_that("gsetroot warns when only the names in seq/genome.idx differ from chrom_sizes.txt", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))
    # chrX renamed in chrom_sizes.txt only, as an interrupted rename leaves it; sizes and order agree
    cs <- file.path(db, "chrom_sizes.txt")
    writeLines(sub("^chrX\t", "chrW\t", readLines(cs)), cs)
    expect_warning(gsetroot(db), "names 1 of its 5 contigs differently")
    expect_equal(toupper(gseq.extract(gintervals("chrW", 0, 1))), "T")
})

test_that(".gdb.chrom_names_at and .gdb.is_indexed_at work with no database loaded", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    empty <- file.path(td, "empty")
    dir.create(file.path(empty, "seq"), recursive = TRUE)
    gdb.unload()
    expect_equal(.gdb.chrom_names_at(db), names(unsorted_db_bases))
    expect_false(.gdb.is_indexed_at(db))
    expect_false(.gdb.is_indexed_at(empty))
})
