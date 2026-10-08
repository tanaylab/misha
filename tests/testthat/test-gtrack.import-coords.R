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
        sam_line("r6", 16, 900, "*", 10) # no CIGAR: sequence length, 908
    ), sam)
    p <- import_points(sam, cols.order = NULL)
    expect_equal(as.numeric(names(p)), c(99, 208, 305, 410, 506, 708, 908, 1000))
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
