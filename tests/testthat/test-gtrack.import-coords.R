# Read placement in gtrack.import_mappedseq: SAM POS is 1-based, reverse reads end at their CIGAR span,
# and tab-delimited files are 0-based unless one.based = TRUE.

create_isolated_test_db()

import_points <- function(file, ...) {
    track <- random_track_name("test")
    gtrack.rm(track, force = TRUE)
    withr::defer(gtrack.rm(track, force = TRUE), envir = parent.frame())
    gtrack.import_mappedseq(track, "coords test", file, remove.dups = FALSE, ...)
    r <- gextract(track, gintervals("chr1", 0, 100000), colnames = "v")
    stats::setNames(r$v, r$start)
}

sam_line <- function(name, flag, pos, cigar, seq_len) {
    paste(name, flag, "chr1", pos, 30, cigar, "*", 0, 0, strrep("A", seq_len), "*", sep = "\t")
}

test_that("SAM reads are placed at their 0-based 5' base", {
    sam <- tempfile(fileext = ".sam")
    writeLines(c(
        sam_line("f1", 0, 100, "10M", 10), # 5' base 99
        sam_line("f2", 0, 1001, "3S7M", 10), # soft clip does not move POS: 1000
        sam_line("r1", 16, 200, "10M", 10), # 199 + 10 - 1 = 208
        sam_line("r2", 16, 300, "3S7M", 10), # 299 + 7 - 1 = 305
        sam_line("r3", 16, 400, "5M2D5M", 10), # 399 + 12 - 1 = 410
        sam_line("r4", 16, 500, "5M2I3M", 10), # 499 + 8 - 1 = 506
        sam_line("r5", 16, 600, "5M100N5M", 10), # 599 + 110 - 1 = 708
        sam_line("r6", 16, 900, "*", 10), # no CIGAR: unmapped, as in htslib
        sam_line("r7", 16, 1100, "2H3S5=2X1P2M", 10), # H, S and P take no reference: 1099 + 9 - 1 = 1107
        sam_line("f3", 0, 1201, "*", 10) # no CIGAR: unmapped
    ), sam)
    p <- import_points(sam, cols.order = NULL)
    expect_equal(as.numeric(names(p)), c(99, 208, 305, 410, 506, 708, 1000, 1107))
    expect_equal(as.numeric(p), rep(1, 8))
})

test_that("SAM dense pileup starts at the 5' base", {
    sam <- tempfile(fileext = ".sam")
    writeLines(c(sam_line("f1", 0, 101, "10M", 10), sam_line("r1", 16, 291, "10M", 10)), sam)
    track <- random_track_name("test")
    withr::defer(gtrack.rm(track, force = TRUE))
    gtrack.import_mappedseq(track, "dense", sam, cols.order = NULL, pileup = 50, binsize = 50)
    # forward [100, 150); reverse 5' base 299, so [250, 300)
    r <- gextract(track, gintervals("chr1", 0, 400), iterator = 50, colnames = "v")
    expect_equal(r$v, c(0, 0, 1, 0, 0, 1, 0, 0))
})

test_that("a read on the last base of a chromosome is imported", {
    chr1_end <- gintervals.all()$end[gintervals.all()$chrom == "chr1"]
    sam <- tempfile(fileext = ".sam")
    writeLines(sam_line("f1", 0, chr1_end, "1M", 1), sam)
    track <- random_track_name("test")
    withr::defer(gtrack.rm(track, force = TRUE))
    stats <- gtrack.import_mappedseq(track, "last base", sam, cols.order = NULL)
    expect_equal(stats[[1]][["total.mapped"]], 1)
    expect_equal(gextract(track, gintervals("chr1", chr1_end - 10, chr1_end))$start, chr1_end - 1)
})

test_that("tab-delimited coordinates are 0-based by default and 1-based with one.based", {
    tab <- tempfile(fileext = ".tsv")
    writeLines(c(
        paste(strrep("A", 10), "chr1", 100, "F", sep = "\t"),
        paste(strrep("A", 10), "chr1", 200, "R", sep = "\t")
    ), tab)
    # legacy placement: forward at the coordinate, reverse one past its 5' end (200 + 10)
    p <- import_points(tab, cols.order = 1:4)
    expect_equal(as.numeric(names(p)), c(100, 210))
    # 1-based: forward 99, reverse 199 + 10 - 1 = 208 (the same as SAM)
    p <- import_points(tab, cols.order = 1:4, one.based = TRUE)
    expect_equal(as.numeric(names(p)), c(99, 208))
})

test_that("one.based is rejected for fragment files", {
    frag <- tempfile(fileext = ".bed")
    writeLines("chr1\t100\t300", frag)
    track <- random_track_name("test")
    withr::defer(gtrack.rm(track, force = TRUE))
    expect_error(gtrack.import_mappedseq(track, "x", frag, binsize = 50, paired = TRUE, one.based = TRUE), "one.based is not used")
})

test_that("SAM records with a malformed or absurd CIGAR are not imported", {
    sam <- tempfile(fileext = ".sam")
    cigars <- c("10S", "5I5S", "0M", "10Q", "M10", "10M5", "10m", "9223372036854775807M")
    writeLines(c(
        mapply(sam_line, paste0("r", seq_along(cigars)), 16, 100, cigars, 10),
        sam_line("f1", 0, 200, "10Q", 10),
        sam_line("ok", 16, 300, "10M", 10)
    ), sam)
    track <- random_track_name("test")
    withr::defer(gtrack.rm(track, force = TRUE))
    stats <- gtrack.import_mappedseq(track, "bad cigar", sam, cols.order = NULL)
    expect_equal(stats[[1]][["total.mapped"]], 1)
    expect_equal(stats[[1]][["total.unmapped"]], length(cigars) + 1)
    expect_equal(gextract(track, gintervals("chr1", 0, 1000))$start, 308)
})

test_that("a tab-delimited reverse read ending at the chromosome end is not written past it", {
    chr1_end <- gintervals.all()$end[gintervals.all()$chrom == "chr1"]
    tab <- tempfile(fileext = ".tsv")
    writeLines(paste(strrep("A", 10), "chr1", chr1_end - 10, "R", sep = "\t"), tab)
    track <- random_track_name("test")
    withr::defer(gtrack.rm(track, force = TRUE))
    # legacy sparse: the point would be chr1_end itself
    expect_warning(stats <- gtrack.import_mappedseq(track, "end", tab, cols.order = 1:4), "No read was imported")
    expect_equal(stats[[1]][["total.mapped"]], 0)
    expect_null(gextract(track, gintervals("chr1", chr1_end - 100, chr1_end)))
})

test_that("reverse reads with the same 5' end are duplicates whatever their POS", {
    sam <- tempfile(fileext = ".sam")
    writeLines(c(
        sam_line("a", 16, 101, "10M", 10), # 5' base 109
        sam_line("b", 16, 103, "8M", 8), # 102 + 8 - 1 = 109
        sam_line("c", 16, 101, "8M2S", 10) # soft clip at the 5' end: 100 + 8 - 1 = 107
    ), sam)
    track <- random_track_name("test")
    withr::defer(gtrack.rm(track, force = TRUE))
    stats <- gtrack.import_mappedseq(track, "dups", sam, cols.order = NULL, remove.dups = TRUE)
    expect_equal(stats[[1]][["total.dups"]], 1)
    r <- gextract(track, gintervals("chr1", 0, 1000))
    expect_equal(r$start, c(107, 109))
})

test_that("SAM and BAM place reads the same way", {
    lines <- c(
        "@SQ\tSN:chr1\tLN:247249719",
        sam_line("f1", 0, 100, "10M", 10), sam_line("r1", 16, 200, "3S7M", 10),
        sam_line("r2", 16, 400, "5M2D5M", 10), sam_line("r3", 16, 600, "5M100N5M", 10)
    )
    bam <- make_test_bam(lines)
    sam <- tempfile(fileext = ".sam")
    writeLines(lines, sam)
    expect_equal(import_points(bam), import_points(sam, cols.order = NULL))
})

test_that("paired: a first mate on the last base is imported and POS 0 is unmapped", {
    chr1_end <- gintervals.all()$end[gintervals.all()$chrom == "chr1"]
    pair <- function(name, pos, pnext, tlen, flag1 = 83, flag2 = 163) {
        c(
            paste(name, flag1, "chr1", pos, 30, "1M", "=", pnext, -tlen, "A", "*", sep = "\t"),
            paste(name, flag2, "chr1", pnext, 30, "10M", "=", pos, tlen, strrep("A", 10), "*", sep = "\t")
        )
    }
    sam <- tempfile(fileext = ".sam")
    writeLines(c(
        pair("last", chr1_end, chr1_end - 49, 50), # fragment [chr1_end - 50, chr1_end)
        paste("zero", 67, "chr1", 0, 30, "10M", "=", 100, 0, strrep("A", 10), "*", sep = "\t")
    ), sam)
    track <- random_track_name("test")
    withr::defer(gtrack.rm(track, force = TRUE))
    stats <- gtrack.import_mappedseq(track, "last", sam, cols.order = NULL, paired = TRUE, binsize = 50)
    expect_equal(stats[[1]][["total.mapped"]], 1)
    expect_equal(stats[[1]][["total.unmapped"]], 1)
    last_bin <- (chr1_end - 50) %/% 50 * 50
    v <- gextract(track, gintervals("chr1", last_bin - 50, chr1_end), iterator = 50, colnames = "v")$v
    expect_equal(sum(v) * 50, 50, tolerance = 1e-5)
})
