# Paired-end fragments and SAM flag / MAPQ filtering in gtrack.import_mappedseq.

create_isolated_test_db()

sam_rec <- function(name, flag, pos, mapq, pnext, tlen, rname = "chr1") {
    paste(name, flag, rname, pos, mapq, "50M", "=", pnext, tlen, strrep("A", 50), "*", sep = "\t")
}

# Fragments (0-based, half-open):
#   p1, p2 [100, 300)    duplicates of each other
#   p3     [1000, 1150)  first mate on the reverse strand
#   p4     [2000, 2150)  first mate MAPQ 5
#   p5                   not a proper pair
#   p6                   secondary alignment
#   p7                   TLEN 5000
#   p8     [6000, 6050)  both mates start at the same position
#   u1                   unmapped mate placed at its partner's position
paired_sam_text <- function() {
    c(
        "@HD\tVN:1.6",
        "@SQ\tSN:chr1\tLN:247249719",
        sam_rec("p1", 99, 101, 40, 251, 200), sam_rec("p1", 147, 251, 40, 101, -200),
        sam_rec("p2", 99, 101, 40, 251, 200), sam_rec("p2", 147, 251, 40, 101, -200),
        sam_rec("p3", 163, 1001, 40, 1101, 150), sam_rec("p3", 83, 1101, 40, 1001, -150),
        sam_rec("p4", 99, 2001, 5, 2101, 150), sam_rec("p4", 147, 2101, 40, 2001, -150),
        sam_rec("p5", 65, 3001, 40, 3051, 100), sam_rec("p5", 129, 3051, 40, 3001, -100),
        sam_rec("p6", 355, 4001, 40, 4051, 100),
        sam_rec("p7", 99, 5001, 40, 9951, 5000), sam_rec("p7", 147, 9951, 40, 5001, -5000),
        sam_rec("p8", 99, 6001, 40, 6001, 50), sam_rec("p8", 147, 6001, 40, 6001, -50),
        paste("u1", 4, "chr1", 7001, 0, "*", "=", 7001, 0, "*", "*", sep = "\t")
    )
}

bins_at <- function(track, starts, binsize = 50) {
    r <- gextract(track, gintervals("chr1", starts, starts + binsize), iterator = binsize, colnames = "v")
    as.numeric(r$v)
}

import_tmp <- function(file, ...) {
    track <- random_track_name("test")
    gtrack.rm(track, force = TRUE)
    withr::defer(gtrack.rm(track, force = TRUE), envir = parent.frame())
    stats <- gtrack.import_mappedseq(track, "paired test", file, ...)
    list(track = track, stats = stats)
}

test_that("paired BAM: one fragment per proper pair, from the first mate", {
    bam <- make_test_bam(paired_sam_text())
    res <- import_tmp(bam, binsize = 50, paired = TRUE)

    expect_equal(bins_at(res$track, c(100, 150, 200, 250)), rep(1, 4))
    expect_equal(bins_at(res$track, c(1000, 1050, 1100)), rep(1, 3))
    expect_equal(bins_at(res$track, c(2000, 2050, 2100)), rep(1, 3))
    expect_equal(bins_at(res$track, 6000), 1)
    # outside the fragments: improper, secondary, too long, unmapped
    expect_equal(bins_at(res$track, c(50, 300, 3000, 3050, 4000, 5000, 7000, 9950)), rep(0, 8))

    expect_equal(res$stats[[1]][["total.mapped"]], 5)
    expect_equal(res$stats[[1]][["total.dups"]], 1)
    expect_equal(res$stats[[1]][["total.unmapped"]], 1)
    # first mates of p5 (improper), p6 (secondary) and p7 (too long); second mates are not counted
    expect_equal(res$stats[[1]][["total.filtered"]], 3)
})

test_that("paired import keeps duplicates when asked and filters by MAPQ and length", {
    bam <- make_test_bam(paired_sam_text())

    res <- import_tmp(bam, binsize = 50, paired = TRUE, remove.dups = FALSE)
    expect_equal(bins_at(res$track, c(100, 250)), c(2, 2))

    res <- import_tmp(bam, binsize = 50, paired = TRUE, min.mapq = 30)
    expect_equal(bins_at(res$track, c(2000, 2100)), c(0, 0))
    expect_equal(bins_at(res$track, 100), 1)
    expect_equal(res$stats[[1]][["total.filtered"]], 4)

    res <- import_tmp(bam, binsize = 50, paired = TRUE, max.fraglen = 6000)
    expect_equal(bins_at(res$track, c(5000, 9900)), c(1, 1))
    expect_equal(res$stats[[1]][["total.filtered"]], 2)
})

test_that("paired import reads a plain SAM file", {
    sam <- tempfile(fileext = ".sam")
    writeLines(paired_sam_text(), sam)
    res <- import_tmp(sam, binsize = 50, cols.order = NULL, paired = TRUE)
    expect_equal(bins_at(res$track, c(100, 1000, 6000)), c(1, 1, 1))
})

fragment_lines <- function() {
    c(
        "# comment line",
        "chr1\t100\t300\tAAAC-1\t3",
        "chr1\t100\t300\tAAAG-1\t1", # same span, another cell: counted again
        "chr1\t1000\t1025\tAAAC-1\t1",
        "chr1\t4000\t9000\tAAAC-1\t1", # longer than max.fraglen
        "chrUnknown\t1\t10\tAAAC-1\t1"
    )
}

test_that("paired import reads fragment files, plain, gzip and bgzip", {
    plain <- tempfile(fileext = ".tsv")
    writeLines(fragment_lines(), plain)
    gz <- tempfile(fileext = ".tsv.gz")
    con <- gzfile(gz, "w")
    writeLines(fragment_lines(), con)
    close(con)
    files <- list(plain = plain, gz = gz)
    if (nzchar(Sys.which("bgzip"))) {
        bgz <- tempfile(fileext = ".tsv")
        writeLines(fragment_lines(), bgz)
        system2("bgzip", bgz)
        files$bgzip <- paste0(bgz, ".gz")
    }

    for (f in names(files)) {
        res <- import_tmp(files[[f]], binsize = 50, paired = TRUE)
        expect_equal(bins_at(res$track, c(100, 250, 1000, 4000)), c(2, 2, 0.5, 0), info = f)
        expect_equal(res$stats[[1]][["total.mapped"]], 3, info = f)
        expect_equal(res$stats[[1]][["total.dups"]], 0, info = f)
        expect_equal(res$stats[[1]][["total.unmapped"]], 1, info = f)
        expect_equal(res$stats[[1]][["total.filtered"]], 1, info = f)
    }
})

test_that("paired import rejects arguments that do not apply", {
    frag <- tempfile(fileext = ".tsv")
    writeLines(fragment_lines(), frag)
    track <- random_track_name("test")
    withr::defer(gtrack.rm(track, force = TRUE))

    expect_error(gtrack.import_mappedseq(track, "x", frag, pileup = 100, binsize = 50, paired = TRUE), "pileup is not used")
    expect_error(gtrack.import_mappedseq(track, "x", frag, paired = TRUE), "binsize")
    expect_error(gtrack.import_mappedseq(track, "x", frag, binsize = 50, paired = TRUE, min.mapq = 30), "min.mapq requires SAM or BAM")
    expect_error(gtrack.import_mappedseq(track, "x", frag, binsize = 50, paired = TRUE, cols.order = c(1, 2, 3, 4)), "cols.order is not used")
})

test_that("single-end SAM import skips unmapped, secondary and low-MAPQ records", {
    sam <- tempfile(fileext = ".sam")
    writeLines(c(
        "@SQ\tSN:chr1\tLN:247249719",
        "r1\t0\tchr1\t100\t30\t10M\t*\t0\t0\tAAAAAAAAAA\t*",
        "r2\t16\tchr1\t200\t10\t10M\t*\t0\t0\tAAAAAAAAAA\t*",
        "r3\t4\tchr1\t300\t0\t*\t*\t0\t0\tAAAAAAAAAA\t*", # unmapped, placed on chr1
        "r4\t256\tchr1\t400\t30\t10M\t*\t0\t0\tAAAAAAAAAA\t*", # secondary
        "r5\t2048\tchr1\t500\t30\t10M\t*\t0\t0\tAAAAAAAAAA\t*" # supplementary
    ), sam)

    res <- import_tmp(sam, cols.order = NULL)
    expect_equal(res$stats[[1]][["total.mapped"]], 2)
    expect_equal(res$stats[[1]][["total.unmapped"]], 1)
    expect_equal(res$stats[[1]][["total.filtered"]], 2)

    res <- import_tmp(sam, cols.order = NULL, min.mapq = 20)
    expect_equal(res$stats[[1]][["total.mapped"]], 1)
    expect_equal(res$stats[[1]][["total.filtered"]], 3)
    expect_true(grepl("min.mapq=20", gtrack.attr.get(res$track, "created.by")))
})

test_that("paired import warns when a SAM file is read as fragments", {
    sam <- tempfile(fileext = ".sam")
    writeLines(paired_sam_text(), sam)
    expect_warning(import_tmp(sam, binsize = 50, paired = TRUE), "cols.order = NULL")
})

test_that("a bgzipped SAM is read as SAM with the default cols.order", {
    skip_if(!nzchar(Sys.which("bgzip")), "bgzip not on PATH")
    sam <- tempfile(fileext = ".sam")
    writeLines(default_sam_text(), sam)
    system2("bgzip", sam)
    res <- import_tmp(paste0(sam, ".gz"))
    expect_equal(res$stats[[1]][["total.mapped"]], 2)
    expect_equal(res$stats[[1]][["total.unmapped"]], 1)
})

test_that("a large max.fraglen is recorded without failing the import", {
    bam <- make_test_bam(paired_sam_text())
    res <- import_tmp(bam, binsize = 50, paired = TRUE, max.fraglen = 1e10)
    expect_equal(bins_at(res$track, c(5000, 9900)), c(1, 1))
    expect_true(grepl("max.fraglen=1e+10", gtrack.attr.get(res$track, "created.by"), fixed = TRUE))
})

test_that("fragment files with CRLF line endings are read", {
    frag <- tempfile(fileext = ".bed")
    writeBin(charToRaw("chr1\t100\t300\r\nchr1\t1000\t1050\r\n"), frag)
    res <- import_tmp(frag, binsize = 50, paired = TRUE)
    expect_equal(bins_at(res$track, c(100, 1000)), c(1, 1))
    expect_equal(res$stats[[1]][["total.unmapped"]], 0)
})

test_that("total counts each record once", {
    bam <- make_test_bam(paired_sam_text())
    s <- import_tmp(bam, binsize = 50, paired = TRUE)$stats[[1]]
    # first mates p1-p8 and the unmapped u1
    expect_equal(s[["total"]], 9)
    expect_equal(s[["total"]], s[["total.mapped"]] + s[["total.unmapped"]] + s[["total.filtered"]])
})

test_that("a reverse read past the chromosome end does not overflow the dense track", {
    chr1_end <- gintervals.all()$end[gintervals.all()$chrom == "chr1"]
    sam <- tempfile(fileext = ".sam")
    writeLines(paste("r1", 16, "chr1", chr1_end - 10, 30, "151M", "*", 0, 0, strrep("A", 151), "*", sep = "\t"), sam)
    res <- import_tmp(sam, cols.order = NULL, pileup = 100, binsize = 20)
    expect_equal(res$stats[[1]][["total.mapped"]], 1)
    expect_equal(sum(gextract(res$track, gintervals("chr1", chr1_end - 200, chr1_end), iterator = 20, colnames = "v")$v), 0)
})
