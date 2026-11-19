# gintervals.import_bigwig reads the records of a bigWig over intervals. Fixtures as in
# test-gtrack.import-bigwig.R: bg.bw has three records per data block (-itemsPerSlot=3), so every
# chromosome spans many blocks and a query selects some of them; step.bw holds fixedStep and
# variableStep sections.

ensure_valid_groot()

bigwig_fixture <- function(name) normalizePath(testthat::test_path("fixtures", "bigwig", name))

setup_bigwig_db <- function() {
    create_test_db("bigwig_db") # chr1 and chr2, 10 kb each; scaffold_9 is not in it
    gdb.init("bigwig_db")
}

# the records of the text file, as gtrack.import reads them (clipped to the chromosomes)
text_records <- function(path) {
    gtrack.import("text_records", "", path, binsize = 0)
    on.exit(gtrack.rm("text_records", force = TRUE))
    x <- gextract("text_records", gintervals.all(), iterator = "text_records", colnames = "value")
    x$intervalID <- NULL
    x
}

# the records that overlap one of the intervals of query
overlapping <- function(records, query) {
    keep <- vapply(seq_len(nrow(records)), function(i) {
        any(query$chrom == records$chrom[i] & query$start < records$end[i] & query$end > records$start[i])
    }, logical(1))
    x <- records[keep, ]
    rownames(x) <- NULL
    x
}

random_query <- function(n, width) {
    start <- sample(0:(10000 - width), n, replace = TRUE)
    data.frame(chrom = sample(c("chr1", "chr2"), n, replace = TRUE), start = start, end = start + sample(1:width, n, replace = TRUE))
}

test_that("all the records of a bigWig are those of the file it was made from", {
    local_db_state()
    paths <- sapply(c("bg.bw", "bg_unc.bw", "bg_libbigwig.bw", "bg.bedGraph", "step.bw", "step.wig"), bigwig_fixture)
    withr::with_tempdir({
        setup_bigwig_db()
        suppressMessages({
            bg <- text_records(paths[["bg.bedGraph"]])
            step <- text_records(paths[["step.wig"]])
        })
        for (bw in c("bg.bw", "bg_unc.bw", "bg_libbigwig.bw")) {
            expect_message(x <- gintervals.import_bigwig(paths[[bw]]), "scaffold_9")
            expect_equal(x, bg, tolerance = 1e-7, label = bw)
        }
        expect_equal(gintervals.import_bigwig(paths[["step.bw"]]), step, tolerance = 1e-7)
    })
})

test_that("a query returns each overlapping record once and whole, reading only the blocks under it", {
    local_db_state()
    paths <- sapply(c("bg.bw", "bg_unc.bw", "step.bw"), bigwig_fixture)
    withr::with_tempdir({
        setup_bigwig_db()
        set.seed(17)
        for (bw in names(paths)) {
            all <- suppressMessages(gintervals.import_bigwig(paths[[bw]]))
            # overlapping and touching query intervals, spanning several blocks or within one
            for (width in c(5, 60, 700)) {
                query <- random_query(40, width)
                expect_equal(gintervals.import_bigwig(paths[[bw]], query), overlapping(all, query), label = sprintf("%s, width %d", bw, width))
            }
        }
    })
})

test_that("blocks away from the query are not read", {
    local_db_state()
    src <- bigwig_fixture("bg.bw")
    bytes <- readBin(src, "raw", file.size(src))
    data_offset <- sum(as.numeric(bytes[17:24]) * 256^(0:7))
    withr::with_tempdir({
        setup_bigwig_db()
        # a byte of the first data block (chr1 20-132, after the 8-byte block count) flipped: zlib's checksum fails
        f <- file.path(getwd(), "flipped.bw")
        writeBin(replace(bytes, data_offset + 19, as.raw(bitwXor(as.integer(bytes[data_offset + 19]), 0xff))), f)
        expect_error(gintervals.import_bigwig(f), "failed to decompress")
        expect_error(gintervals.import_bigwig(f, gintervals(1, 100, 110)), "failed to decompress")
        expect_equal(
            gintervals.import_bigwig(f, gintervals(c(1, 2), 5000, 5100)),
            gintervals.import_bigwig(src, gintervals(c(1, 2), 5000, 5100))
        )
    })
})

test_that("the records serve as a value-based virtual track", {
    local_db_state()
    bw <- bigwig_fixture("bg.bw")
    withr::with_tempdir({
        setup_bigwig_db()
        suppressMessages(gtrack.import("bg", "", bw, binsize = 0))
        peaks <- gintervals(c(1, 1, 2, 2), c(100, 132, 1000, 9000), c(400, 134, 3000, 9050))
        src <- gintervals.import_bigwig(bw, peaks)
        for (func in c("sum", "weighted.sum", "avg", "max", "size", "exists")) {
            gvtrack.create("sig", src = src, func = func)
            gvtrack.create("ref", src = "bg", func = func)
            x <- gextract(c("sig", "ref"), peaks, iterator = peaks)
            expect_equal(x$sig, x$ref, label = func)
            if (func %in% c("size", "exists")) {
                expect_equal(x$sig[2], 0) # no record in chr1 132-134
            } else {
                expect_true(is.nan(x$sig[2]), label = func)
            }
        }
    })
})

test_that("intervals no record overlaps give a result of no rows that serves as a virtual track", {
    local_db_state()
    bw <- bigwig_fixture("bg.bw")
    withr::with_tempdir({
        setup_bigwig_db()
        empty <- gintervals.import_bigwig(bw, gintervals(1, 132, 134)) # between the records 94-132 and 134-140
        expect_identical(empty, gintervals.import_bigwig(bw, gintervals(1, 100, 110))[0, ])
        gvtrack.create("sig", src = empty, func = "weighted.sum")
        x <- gextract("sig", gintervals(1, 0, 10000), iterator = gintervals(1, 0, 10000))
        expect_true(is.nan(x$sig))
    })
})

test_that("records past the chromosome end are clipped as gtrack.import clips them", {
    local_db_state()
    bw <- bigwig_fixture("bg.bw")
    withr::with_tempdir({
        # chr2 ends inside the record 9870-9910, before the record 9933-9981
        create_test_db("short_bg_db", chrom_sizes = data.frame(chrom = c("chr1", "chr2"), size = c(10000, 9900)))
        gdb.init("short_bg_db")
        msg <- "1 record\\(s\\) running past the end of their chromosome were clipped to it, and 1 record\\(s\\) starting at or past the chromosome end were dropped; the first: .*\\(chr2, 9870, 9910\\), chromosome size 9900"
        expect_warning(suppressMessages(x <- gintervals.import_bigwig(bw)), msg)
        expect_equal(tail(x$end, 1), 9900)
        # a query counts only the records it returns
        expect_warning(x <- gintervals.import_bigwig(bw, gintervals(2, 9880, 9890)), "1 record\\(s\\) running past .* clipped to it, and 0 record\\(s\\)")
        expect_equal(x$end, 9900)
        expect_no_warning(gintervals.import_bigwig(bw, gintervals(2, 9000, 9100)))
    })
})

test_that("chromosome names are matched through the database aliases", {
    local_db_state()
    bw <- bigwig_fixture("bg.bw")
    withr::with_tempdir({
        setup_bigwig_db()
        expected <- suppressMessages(gintervals.import_bigwig(bw, gintervals(2, 0, 500)))
        create_test_db("noprefix_db", chrom_sizes = data.frame(chrom = c("1", "2"), size = c(10000, 10000)))
        gdb.init("noprefix_db")
        x <- gintervals.import_bigwig(bw, gintervals(2, 0, 500))
        expect_equal(as.character(x$chrom), rep("2", nrow(expected)))
        expect_equal(x[, -1], expected[, -1])
    })
})

test_that("bad input is an error", {
    local_db_state()
    paths <- sapply(c("bg.bw", "bg.bedGraph"), bigwig_fixture)
    bw <- paths[["bg.bw"]]
    withr::with_tempdir({
        setup_bigwig_db()
        expect_error(gintervals.import_bigwig("no_such_file.bw"), "does not exist")
        expect_error(gintervals.import_bigwig(paths[["bg.bedGraph"]]), "is not a bigWig file")
        create_test_db("other_db", chrom_sizes = data.frame(chrom = c("contigA", "contigB"), size = c(10000, 10000)))
        gdb.init("other_db")
        expect_error(gintervals.import_bigwig(bw), "None of the chromosome names")
    })
})
