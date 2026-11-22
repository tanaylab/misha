# gintervals.import_bed reads bigBed files. The fixtures were made by UCSC's bedToBigBed with
# -blockSize=2 -itemsPerSlot=4, so the records span many data blocks:
#   typed.bb:  bedToBigBed -tab -type=bed9+7 -as=typed.as typed.bed chrom.sizes (extra fields of
#              several autoSql types; three records on scaffold_9, which the test database lacks;
#              chr1 5000-5100 before 5000-5050, so the file's order is not the sorted one;
#              chr2 9870-9910 is the last record of chr2)
#   bed12.bb:  bedToBigBed -type=bed12 bed12.bed chrom.sizes (UCSC's default BED autoSql)

ensure_valid_groot()

bigbed_dir <- normalizePath(testthat::test_path("fixtures", "bigbed"))
bigbed_fixture <- function(name) file.path(bigbed_dir, name)

setup_bigbed_db <- function() {
    create_test_db("bigbed_db") # chr1 and chr2, 10 kb each
    gdb.init("bigbed_db")
}

# the header field at byte pos, of size bytes (little-endian)
header_field <- function(bytes, pos, size) sum(as.numeric(bytes[pos + seq_len(size)]) * 256^(seq_len(size) - 1))

# a copy of a fixture with header bytes pos.. replaced by value (little-endian, size bytes)
patched <- function(name, pos, size, value) {
    bytes <- readBin(bigbed_fixture(name), "raw", file.size(bigbed_fixture(name)))
    bytes[pos + seq_len(size)] <- as.raw((value %/% 256^(seq_len(size) - 1)) %% 256)
    f <- tempfile(fileext = ".bb")
    writeBin(bytes, f)
    f
}

typed_text <- function() {
    x <- read.table(bigbed_fixture("typed.bed"), sep = "\t", quote = "", comment.char = "", colClasses = "character")
    x[x$V1 != "scaffold_9", ]
}

test_that("a bigBed's fields are named and typed by its autoSql, in the file's order", {
    local_db_state()
    bb <- bigbed_fixture("typed.bb")
    withr::with_tempdir({
        setup_bigbed_db()
        expect_message(x <- gintervals.import_bed(bb, keep_order = TRUE), "scaffold_9")
        expect_equal(names(x), c(
            "chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb",
            "count", "signal", "pval", "big", "note", "state", "triple"
        ))
        text <- typed_text()
        expect_equal(as.character(x$chrom), text$V1)
        expect_identical(x$start, as.numeric(text$V2))
        expect_identical(x$end, as.numeric(text$V3))
        expect_identical(x$strand, unname(c(`+` = 1, `-` = -1, `.` = 0)[text$V6]))
        numeric_cols <- c(score = 5, thickStart = 7, thickEnd = 8, count = 10, signal = 11, pval = 12, big = 13)
        for (col in names(numeric_cols)) {
            expect_identical(x[[col]], as.numeric(text[[numeric_cols[[col]]]]), label = col)
        }
        # itemRgb is "R,G,B" text although autoSql declares it a uint
        text_cols <- c(name = 4, itemRgb = 9, note = 14, state = 15, triple = 16)
        for (col in names(text_cols)) {
            expect_identical(x[[col]], text[[text_cols[[col]]]], label = col)
        }
        expect_true(any(grepl(",", x$itemRgb)))

        # sorted by default; the file has chr1 5000-5100 before 5000-5050
        by_pos <- suppressMessages(gintervals.import_bed(bb))
        expect_identical(by_pos, misha:::.gsort_intervals_df(x))
        same_start <- c("chr1_same_start_0", "chr1_same_start_1")
        expect_equal(x$end[x$name %in% same_start], c(5100, 5050))
        expect_equal(by_pos$end[by_pos$name %in% same_start], c(5050, 5100))
    })
})

test_that("name, score and strand govern the fields of those names", {
    local_db_state()
    bb <- bigbed_fixture("typed.bb")
    withr::with_tempdir({
        setup_bigbed_db()
        all <- suppressMessages(gintervals.import_bed(bb))
        x <- suppressMessages(gintervals.import_bed(bb, name = FALSE, score = FALSE, strand = FALSE))
        expect_identical(x, all[, setdiff(names(all), c("name", "score", "strand"))])
    })
})

test_that("a bigBed reads as the text BED it was made from", {
    local_db_state()
    paths <- sapply(c("bed12.bb", "bed12.bed"), bigbed_fixture)
    withr::with_tempdir({
        setup_bigbed_db()
        x <- gintervals.import_bed(paths[["bed12.bb"]])
        text <- gintervals.import_bed(paths[["bed12.bed"]])
        expect_equal(names(x), c(
            "chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "reserved",
            "blockCount", "blockSizes", "chromStarts"
        ))
        expect_identical(x[, names(text)], text)
    })
})

test_that("a bigBed without autoSql gets the BED standard field names", {
    local_db_state()
    withr::with_tempdir({
        setup_bigbed_db()
        with_as <- gintervals.import_bed(bigbed_fixture("bed12.bb"))
        # autoSqlOffset (8 bytes at 36) set to 0
        expect_identical(gintervals.import_bed(patched("bed12.bb", 36, 8, 0)), with_as)
        expect_error(
            suppressMessages(gintervals.import_bed(patched("typed.bb", 36, 8, 0))),
            "has 7 fields past its 9 BED standard ones, and no autoSql to name them"
        )
    })
})

test_that("format cannot be given for a bigBed, whatever its name", {
    local_db_state()
    withr::with_tempdir({
        setup_bigbed_db()
        expect_error(gintervals.import_bed(bigbed_fixture("bed12.bb"), format = "bed"), "a bigBed names its columns in its autoSql")
        f <- file.path(getwd(), "peaks.narrowPeak")
        file.copy(bigbed_fixture("bed12.bb"), f)
        expect_identical(gintervals.import_bed(f), gintervals.import_bed(bigbed_fixture("bed12.bb")))
    })
})

test_that("bigBed records past the chromosome end are clipped as gtrack.import clips them", {
    local_db_state()
    bb <- bigbed_fixture("bed12.bb")
    withr::with_tempdir({
        create_test_db("short_db", chrom_sizes = data.frame(chrom = c("chr1", "chr2"), size = c(10000, 9900)))
        gdb.init("short_db")
        expect_warning(
            x <- gintervals.import_bed(bb),
            "1 record\\(s\\) running past the end of their chromosome were clipped to it, and 0 record\\(s\\) .* the first: \\(chr2, 9870, 9910\\), chromosome size 9900"
        )
        expect_equal(tail(x$end, 1), 9900)
    })
})

test_that("a damaged bigBed is an error", {
    local_db_state()
    src <- bigbed_fixture("typed.bb")
    bytes <- readBin(src, "raw", file.size(src))
    data_offset <- header_field(bytes, 16, 8)
    num_records <- header_field(bytes, data_offset, 8)
    withr::with_tempdir({
        setup_bigbed_db()
        expect_error(gintervals.import_bed(patched("typed.bb", data_offset, 8, num_records + 1)), "records were read, but its header counts")
        f <- file.path(getwd(), "cut.bb")
        writeBin(bytes[seq_len(data_offset + 100)], f)
        expect_error(gintervals.import_bed(f), "bigBed")
        big_endian <- file.path(getwd(), "big_endian.bb")
        writeBin(c(rev(bytes[1:4]), bytes[-(1:4)]), big_endian)
        expect_error(gintervals.import_bed(big_endian), "big-endian")
    })
})
