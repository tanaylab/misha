# Helper: build a tiny in-memory groot with two contigs.
make_test_groot <- function() {
    groot <- tempfile()
    dir.create(groot)
    dir.create(file.path(groot, "tracks"))
    dir.create(file.path(groot, "seq"))
    # Minimal chrom_sizes
    writeLines(c("chr1\t1000", "chr2\t1000"), file.path(groot, "chrom_sizes.txt"))
    # Empty .seq files (tests that don't touch sequence won't notice)
    file.create(file.path(groot, "seq", "chr1.seq"))
    file.create(file.path(groot, "seq", "chr2.seq"))
    groot
}

test_that(".save_intervals filters to ALLGENOME and writes the set", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    gdb.init(groot, rescan = TRUE)
    df <- data.frame(
        chrom = c("chr1", "chr1", "chrX"),
        start = c(0L, 100L, 0L),
        end = c(50L, 150L, 50L),
        stringsAsFactors = FALSE
    )
    .save_intervals("foo", df, overwrite = FALSE, verbose = FALSE)
    expect_true(gintervals.exists("foo"))
    expect_equal(nrow(gintervals.load("foo")), 2L)
})

test_that(".save_intervals errors on existing set when overwrite=FALSE", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    gdb.init(groot, rescan = TRUE)
    df <- data.frame(chrom = "chr1", start = 0L, end = 50L, stringsAsFactors = FALSE)
    .save_intervals("foo", df, overwrite = FALSE, verbose = FALSE)
    expect_error(
        .save_intervals("foo", df, overwrite = FALSE, verbose = FALSE),
        "already exists"
    )
})

test_that(".save_intervals overwrites when overwrite=TRUE", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    gdb.init(groot, rescan = TRUE)
    df1 <- data.frame(chrom = "chr1", start = 0L, end = 50L, stringsAsFactors = FALSE)
    df2 <- data.frame(
        chrom = "chr1",
        start = c(0L, 200L),
        end = c(50L, 300L),
        stringsAsFactors = FALSE
    )
    .save_intervals("foo", df1, overwrite = FALSE, verbose = FALSE)
    .save_intervals("foo", df2, overwrite = TRUE, verbose = FALSE)
    expect_equal(nrow(gintervals.load("foo")), 2L)
})

test_that(".save_intervals skips empty frame with message", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    gdb.init(groot, rescan = TRUE)
    df <- data.frame(
        chrom = character(0),
        start = integer(0),
        end = integer(0),
        stringsAsFactors = FALSE
    )
    expect_message(
        .save_intervals("empty", df, overwrite = FALSE, verbose = TRUE),
        "0 rows"
    )
    expect_false(gintervals.exists("empty"))
})
