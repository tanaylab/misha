test_that("gdb.export_fasta exports current database with wrapping and chunking", {
    local_db_state()

    test_fasta <- tempfile(fileext = ".fasta")
    test_db <- tempfile()
    out_fasta <- tempfile(fileext = ".fa")

    withr::defer({
        unlink(test_fasta)
        unlink(test_db, recursive = TRUE)
        unlink(out_fasta)
    })

    cat(">chrA\nACTGACTG\n>chrB\nTTAA\n", file = test_fasta)

    suppressMessages(gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE))
    suppressMessages(gdb.init(test_db))

    expect_no_error(gdb.export_fasta(out_fasta, line_width = 3, chunk_size = 2))
    expect_equal(
        readLines(out_fasta),
        c(">chrA", "ACT", "GAC", "TG", ">chrB", "TTA", "A")
    )
})

test_that("gdb.export_fasta can export explicit groot and restores previous root", {
    local_db_state()

    fasta_db1 <- tempfile(fileext = ".fasta")
    fasta_db2 <- tempfile(fileext = ".fasta")
    db1 <- tempfile()
    db2 <- tempfile()
    out_fasta <- tempfile(fileext = ".fa")

    withr::defer({
        unlink(fasta_db1)
        unlink(fasta_db2)
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
        unlink(out_fasta)
    })

    cat(">chr1\nAAAA\n", file = fasta_db1)
    cat(">chr2\nCCCC\n", file = fasta_db2)

    suppressMessages(gdb.create(groot = db1, fasta = fasta_db1, verbose = FALSE))
    suppressMessages(gdb.create(groot = db2, fasta = fasta_db2, verbose = FALSE))
    suppressMessages(gdb.init(db1))

    original_groot <- normalizePath(.misha$GROOT, mustWork = TRUE)

    expect_no_error(gdb.export_fasta(out_fasta, groot = db2, line_width = 2, chunk_size = 2))
    expect_equal(normalizePath(.misha$GROOT, mustWork = TRUE), original_groot)
    expect_equal(readLines(out_fasta), c(">chr2", "CC", "CC"))
})

test_that("gdb.export_fasta overwrite guard works", {
    local_db_state()

    test_fasta <- tempfile(fileext = ".fasta")
    test_db <- tempfile()
    out_fasta <- tempfile(fileext = ".fa")

    withr::defer({
        unlink(test_fasta)
        unlink(test_db, recursive = TRUE)
        unlink(out_fasta)
    })

    cat(">chr1\nACTG\n", file = test_fasta)

    suppressMessages(gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE))
    suppressMessages(gdb.init(test_db))

    writeLines("placeholder", out_fasta)
    expect_error(gdb.export_fasta(out_fasta), "already exists")

    expect_no_error(gdb.export_fasta(out_fasta, overwrite = TRUE))
    expect_equal(readLines(out_fasta), c(">chr1", "ACTG"))
})
