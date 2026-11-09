# gtrack.import reads bigWig files natively. The fixtures were made from bg.bedGraph and
# step.wig by UCSC's bedGraphToBigWig (bg.bw with -blockSize=2 -itemsPerSlot=3, so both
# index trees have inner nodes and every chromosome spans many blocks; bg_unc.bw with -unc),
# wigToBigWig (step.bw: fixedStep with span < step, variableStep with span > 1) and libBigWig
# through pyBigWig (bg_libbigwig.bw, which also lists a chromosome that has no data).
# edge.bw (bedGraphToBigWig): chr1 100-110 12345.678, a zero-length record at 120, and
# 130-140 0.000123456789.

ensure_valid_groot()

setup_bigwig_db <- function() {
    create_test_db("bigwig_db") # chr1 and chr2, 10 kb each; scaffold_9 is not in it
    gdb.init("bigwig_db")
}

setup_short_db <- function() { # chr2 ends inside step.wig's last variableStep record (9982-9989)
    create_test_db("short_db", chrom_sizes = data.frame(chrom = c("chr1", "chr2"), size = c(10000, 9985)))
    gdb.init("short_db")
}

bigwig_fixture <- function(name) normalizePath(testthat::test_path("fixtures", "bigwig", name))

extract_all <- function(track, iterator) {
    x <- gextract(track, gintervals.all(), iterator = iterator, colnames = "v")
    x$intervalID <- NULL
    x
}

# header fields: fullDataOffset at byte 16, fullIndexOffset at byte 24
bigwig_offsets <- function(path) {
    h <- readBin(path, "raw", 64)
    at <- function(pos) sum(as.numeric(h[pos + 1:8]) * 256^(0:7))
    c(data = at(16), index = at(24))
}

test_that("a bigWig imports exactly like the text file it was made from", {
    local_db_state()
    sources <- c(bg.bw = "bg.bedGraph", bg_unc.bw = "bg.bedGraph", bg_libbigwig.bw = "bg.bedGraph", step.bw = "step.wig")
    paths <- sapply(c(names(sources), unique(sources)), bigwig_fixture)
    withr::with_tempdir({
        setup_bigwig_db()
        for (bw in names(sources)) {
            for (p in list(c(0, NaN), c(7, NaN), c(7, 0))) {
                suppressMessages({
                    gtrack.import("from_text", "", paths[[sources[[bw]]]], binsize = p[1], defval = p[2])
                    gtrack.import("from_bw", "", paths[[bw]], binsize = p[1], defval = p[2])
                })
                expect_identical(
                    extract_all("from_bw", if (p[1] == 0) "from_bw" else p[1]),
                    extract_all("from_text", if (p[1] == 0) "from_text" else p[1]),
                    label = sprintf("%s, binsize %g, defval %g", bw, p[1], p[2])
                )
                gtrack.rm("from_text", force = TRUE)
                gtrack.rm("from_bw", force = TRUE)
            }
        }
    })
})

test_that("step records are clipped at the database chromosome end, as in the text route", {
    local_db_state()
    paths <- sapply(c("step.bw", "step.wig"), bigwig_fixture)
    withr::with_tempdir({
        setup_short_db()
        for (binsize in c(0, 7)) {
            expect_warning(gtrack.import("from_text", "", paths[["step.wig"]], binsize = binsize), "1 record\\(s\\) running past the end of their chromosome were clipped to it")
            expect_warning(gtrack.import("from_bw", "", paths[["step.bw"]], binsize = binsize), "1 record\\(s\\) running past the end of their chromosome were clipped to it")
            bw <- extract_all("from_bw", if (binsize == 0) "from_bw" else binsize)
            expect_identical(bw, extract_all("from_text", if (binsize == 0) "from_text" else binsize))
            expect_equal(max(bw$end[bw$chrom == "chr2"]), 9985)
            gtrack.rm("from_text", force = TRUE)
            gtrack.rm("from_bw", force = TRUE)
        }
    })
})

test_that("bedGraph records past the database chromosome end are clipped or dropped, as in the text route", {
    local_db_state()
    paths <- sapply(c("bg.bw", "bg.bedGraph"), bigwig_fixture)
    withr::with_tempdir({
        # chr2 ends inside the record 9870-9910, before the record 9933-9981
        create_test_db("short_bg_db", chrom_sizes = data.frame(chrom = c("chr1", "chr2"), size = c(10000, 9900)))
        gdb.init("short_bg_db")
        msg <- "1 record\\(s\\) running past the end of their chromosome were clipped to it, and 1 record\\(s\\) starting at or past the chromosome end were dropped; the first: .*\\(chr2, 9870, 9910\\), chromosome size 9900"
        for (binsize in c(0, 7)) {
            suppressMessages({
                expect_warning(gtrack.import("from_text", "", paths[["bg.bedGraph"]], binsize = binsize), msg)
                expect_warning(gtrack.import("from_bw", "", paths[["bg.bw"]], binsize = binsize), msg)
            })
            bw <- extract_all("from_bw", if (binsize == 0) "from_bw" else binsize)
            expect_identical(bw, extract_all("from_text", if (binsize == 0) "from_text" else binsize))
            if (binsize == 0) {
                expect_equal(tail(bw$end[bw$chrom == "chr2"], 1), 9900)
            }
            gtrack.rm("from_text", force = TRUE)
            gtrack.rm("from_bw", force = TRUE)
        }
    })
})

test_that("the bigWig record past the chromosome end that the warning names is the first in the file", {
    local_db_state()
    paths <- sapply(c("bg.bw", "bg.bedGraph"), bigwig_fixture)
    withr::with_tempdir({
        # chr2 first in the database, chr1 first in the files; both end inside a record
        create_test_db("rev_db", chrom_sizes = data.frame(chrom = c("chr2", "chr1"), size = c(9900, 9900)))
        writeLines(c("chr2\t9900", "chr1\t9900"), file.path("rev_db", "chrom_sizes.txt"))
        gdb.init("rev_db")
        skip_if_not(identical(as.character(gintervals.all()$chrom), c("chr2", "chr1")), "the database sorts its chromosomes")
        for (p in paths) {
            suppressMessages(expect_warning(gtrack.import("first", "", p, binsize = 0), "the first: (line [0-9]+: )?\\(chr1, 9888, 9925\\)", label = basename(p)))
            gtrack.rm("first", force = TRUE)
        }
    })
})

test_that("bigWig values keep full float precision and zero-length records are skipped", {
    local_db_state()
    bw <- bigwig_fixture("edge.bw")
    withr::with_tempdir({
        setup_bigwig_db()
        gtrack.import("edge", "", bw, binsize = 0)
        x <- extract_all("edge", "edge")
        expect_equal(x$start, c(100, 130))
        expect_equal(x$end, c(110, 140))
        # bigWigToWig printed 6 significant digits: 12345.7, 0.000123457
        expect_equal(x$v, c(12345.678, 0.000123456789), tolerance = 1e-7)
    })
})

test_that("bigWig chromosomes the database lacks are reported, chromosomes without data are not", {
    local_db_state()
    bw <- bigwig_fixture("bg_libbigwig.bw")
    withr::with_tempdir({
        setup_bigwig_db()
        msgs <- character()
        expect_no_warning(withCallingHandlers(
            gtrack.import("skipped", "", bw, binsize = 10),
            message = function(m) {
                msgs <<- c(msgs, conditionMessage(m))
                invokeRestart("muffleMessage")
            }
        ))
        expect_true(any(grepl("scaffold_9", msgs)))
        expect_false(any(grepl("chrEmpty", msgs)))
        expect_true(gtrack.exists("skipped"))
    })
})

test_that("a truncated or corrupt bigWig fails and leaves no track", {
    local_db_state()
    src <- bigwig_fixture("bg.bw")
    off <- bigwig_offsets(src)
    bytes <- readBin(src, "raw", file.size(src))
    withr::with_tempdir({
        setup_bigwig_db()
        damaged <- list(
            # the full-resolution index comes after the data, so either cut loses it
            cut_in_data = bytes[seq_len(off[["data"]] + 20)],
            cut_in_index = bytes[seq_len(off[["index"]] + 60)],
            # inside the first compressed block (after the 8-byte block count): zlib's checksum fails
            flipped_byte = replace(bytes, off[["data"]] + 19, as.raw(bitwXor(as.integer(bytes[off[["data"]] + 19]), 0xff))),
            # the 8-byte block count at fullDataOffset no longer matches the index
            bad_block_count = replace(bytes, off[["data"]] + 1, as.raw((as.integer(bytes[off[["data"]] + 1]) + 1) %% 256)),
            # named .bw but empty, or text: not imported as a WIG
            empty = raw(0),
            text = charToRaw("chr1\t0\t100\t1\n")
        )
        for (what in names(damaged)) {
            f <- file.path(getwd(), paste0(what, ".bw"))
            writeBin(damaged[[what]], f)
            expect_error(gtrack.import("damaged", "", f, binsize = 10), "bigWig", label = what)
            expect_false(gtrack.exists("damaged"), label = what)
        }
    })
})
