# Regression tests for the track import and creation findings of the 2026-11 review of 5.14.4.

restore_groot_on_exit()

# A small database of its own (chr1, chr2, chrX, chrM), in the given format, made the root for the
# calling test. Removed, and the previous root restored, when the test ends.
import_test_db <- function(format = "indexed", env = parent.frame()) {
    local_db_state(env)
    dir <- tempfile("import_review_")
    dir.create(dir)
    withr::defer(unlink(dir, recursive = TRUE), envir = env)
    sizes <- c(chr1 = 10000, chr2 = 10000, chrX = 5000, chrM = 500)
    records <- paste0(">", names(sizes), "\n", strrep("A", sizes), "\n")
    # an indexed database from one multi-FASTA, a per-chromosome one from a FASTA per chromosome
    fa <- if (format == "indexed") file.path(dir, "genome.fa") else file.path(dir, paste0(names(sizes), ".fa"))
    if (format == "indexed") {
        cat(records, sep = "", file = fa)
    } else {
        for (i in seq_along(sizes)) cat(records[i], file = fa[i])
    }
    db <- file.path(dir, "db")
    suppressMessages(gdb.create(groot = db, fasta = fa, format = format, verbose = FALSE))
    suppressMessages(gdb.init(db))
    db
}

write_file <- function(text, ext) {
    path <- tempfile(fileext = ext)
    cat(text, file = path)
    path
}

for (format in c("indexed", "per-chromosome")) {
    test_that(sprintf("a tab-delimited file without a trailing newline keeps its last line (%s)", format), {
        import_test_db(format)

        f <- write_file("chrom\tstart\tend\tvalue\nchr1\t0\t10\t1\nchr1\t20\t30\t2\nchrX\t5\t9\t3", ".tsv")
        gtrack.import("nonl", "", f, 0)
        res <- gextract("nonl", gintervals.all())
        expect_equal(as.character(res$chrom), c("chr1", "chr1", "chrX"))
        expect_equal(res$start, c(0, 20, 5))
        expect_equal(res$nonl, c(1, 2, 3))

        f <- write_file("chrom\tstart\tend\ta\tb\nchr1\t0\t10\t1\t2\nchr1\t20\t30\t3\t4", ".tsv")
        gtrack.array.import("arr_nonl", "", f)
        res <- gtrack.array.extract("arr_nonl", NULL, gintervals.all())
        expect_equal(res$start, c(0, 20))
        expect_equal(res$a, c(1, 3))
        expect_equal(res$b, c(2, 4))
    })
}

test_that("gtrack.2d.import_contacts reads a fends file not sorted by fend id", {
    import_test_db("indexed")

    fends <- write_file("fend\tchr\tcoord\n10\tchr1\t1000\n5\tchr1\t2000\n7\tchr2\t300\n", ".tsv")
    contacts <- write_file("fend1\tfend2\tcount\n10\t5\t3\n10\t7\t4\n5\t7\t1\n", ".tsv")
    gtrack.2d.import_contacts("unsorted_fends", "", contacts, fends)
    res <- gextract("unsorted_fends", gintervals.2d.all())
    got <- data.frame(
        chrom1 = as.character(res$chrom1), start1 = res$start1, chrom2 = as.character(res$chrom2),
        start2 = res$start2, value = res$unsorted_fends
    )
    expected <- data.frame(
        chrom1 = c("chr1", "chr1", "chr1", "chr1", "chr2", "chr2"),
        start1 = c(1000, 1000, 2000, 2000, 300, 300),
        chrom2 = c("chr1", "chr2", "chr1", "chr2", "chr1", "chr1"),
        start2 = c(2000, 300, 1000, 300, 1000, 2000),
        value = c(3, 4, 3, 1, 4, 1)
    )
    got <- got[do.call(order, got), ]
    rownames(got) <- NULL
    expect_equal(got, expected)

    neg <- write_file("fend\tchr\tcoord\n-1\tchr1\t1000\n2\tchr1\t2000\n", ".tsv")
    expect_error(gtrack.2d.import_contacts("neg_fend", "", contacts, neg), "fend is out of range")
    expect_false(gtrack.exists("neg_fend"))
})

test_that("gtrack.create_dense validates the intervals as gtrack.create_sparse does", {
    import_test_db("indexed")

    expect_error(
        gtrack.create_dense("neg", "", data.frame(chrom = "chr1", start = c(-10, 20), end = c(50, 60)), c(1, 3), 100),
        "start coordinate must be greater or equal than zero"
    )
    expect_error(
        gtrack.create_dense("rev", "", data.frame(chrom = "chr1", start = c(50, 20), end = c(10, 60)), c(1, 3), 100),
        "start coordinate must be lesser than end coordinate"
    )
    expect_error(
        gtrack.create_dense("na", "", data.frame(chrom = "chr1", start = c(NA, 20), end = c(10, 60)), c(1, 3), 100),
        "Invalid format of interval at index 1"
    )
    expect_error(
        gtrack.create_dense("past_end", "", data.frame(chrom = "chr1", start = 9990, end = 10010), 1, 100),
        "end coordinate exceeds chromosome boundaries"
    )
    expect_false(any(gtrack.exists(c("neg", "rev", "na", "past_end"))))

    # aliases still resolve
    gtrack.create_dense("alias", "", data.frame(chrom = c("1", "MT"), start = c(0, 0), end = c(50, 100)), c(2, 5), 100)
    expect_equal(gextract("alias", gintervals(c("chr1", "chrM"), 0, 100))$alias, c(2, 5))
})
