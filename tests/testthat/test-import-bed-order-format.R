create_isolated_test_db()

# gintervals.import_bed(keep_order =, format =): the file's row order, and the
# extra narrowPeak / broadPeak columns. Expected values come from a
# utils::read.table() parse of the same lines, not from the importer.

# a narrowPeak body: unsorted, with a -1 (no summit) row and a strand column
narrowpeak_lines <- function() {
    c(
        "chr2\t5000\t5600\tp1\t310\t.\t7.25\t20.5\t18.125\t260",
        "chr1\t9000\t9400\tp2\t1000\t+\t11.5\t31.75\t29.5\t180",
        "chr1\t2000\t2300\tp3\t55\t-\t3.5\t4.25\t2.125\t-1",
        "chr1\t2000\t2200\tp4\t70\t.\t5\t6\t1.5\t0"
    )
}

write_lines_to <- function(lines, ext, gz = FALSE) {
    f <- tempfile(fileext = ext)
    con <- if (gz) gzfile(f, "wt") else file(f, "wt")
    writeLines(lines, con)
    close(con)
    f
}

# the lines (header and comment lines skipped) parsed independently of misha
reference_table <- function(lines) {
    data_lines <- lines[!grepl("^(track|browser|#|$)", lines)]
    utils::read.table(text = data_lines, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
}

test_that("keep_order = TRUE returns the rows in file order, header lines skipped", {
    lines <- c(
        "track name=\"peaks\"",
        "# a comment",
        "chr2\t500\t600\tb\t2\t+",
        "browser position chr1:1-1000",
        "chr1\t300\t400\tc\t3\t-",
        "chr1\t100\t200\ta\t1\t."
    )
    ref <- reference_table(lines)
    f <- write_lines_to(lines, ".bed")
    withr::defer(unlink(f))

    for (fread_available in c(TRUE, FALSE)) {
        local({
            testthat::local_mocked_bindings(.gdata_table_available = function() fread_available, .package = "misha")
            out <- gintervals.import_bed(f, keep_order = TRUE)
            expect_equal(as.character(out$chrom), ref$V1)
            expect_equal(out$start, ref$V2)
            expect_equal(out$end, ref$V3)
            expect_equal(out$name, ref$V4)
            expect_equal(out$score, as.numeric(ref$V5))
            expect_equal(out$strand, c(1, -1, 0))
            expect_identical(rownames(out), as.character(seq_len(nrow(ref))))

            # the default is unchanged: sorted by chrom, then start
            sorted <- gintervals.import_bed(f)
            ord <- order(ref$V1, ref$V2)
            expect_equal(sorted$name, ref$V4[ord])
            expect_identical(rownames(sorted), as.character(seq_len(nrow(ref))))
        })
    }
})

test_that("keep_order must be TRUE or FALSE", {
    f <- write_lines_to("chr1\t100\t200", ".bed")
    withr::defer(unlink(f))
    expect_error(gintervals.import_bed(f, keep_order = NA), "keep_order must be TRUE or FALSE")
    expect_error(gintervals.import_bed(f, keep_order = "yes"), "keep_order must be TRUE or FALSE")
})

test_that("a .narrowPeak file gets signalValue, pValue, qValue and peak after the BED columns", {
    lines <- narrowpeak_lines()
    ref <- reference_table(lines)
    f <- write_lines_to(lines, ".narrowPeak")
    withr::defer(unlink(f))

    out <- gintervals.import_bed(f, keep_order = TRUE)
    expect_identical(
        colnames(out),
        c("chrom", "start", "end", "strand", "name", "score", "signalValue", "pValue", "qValue", "peak")
    )
    expect_equal(out$name, ref$V4)
    expect_equal(out$signalValue, ref$V7)
    expect_equal(out$pValue, ref$V8)
    expect_equal(out$qValue, ref$V9)
    expect_identical(out$peak, as.integer(ref$V10))
    # -1 (no summit) is kept as the file has it
    expect_identical(out$peak[ref$V4 == "p3"], -1L)

    # sorted, the format columns move with their rows
    sorted <- gintervals.import_bed(f)
    ord <- order(ref$V1, ref$V2, ref$V3)
    expect_equal(sorted$name, ref$V4[ord])
    expect_identical(sorted$peak, as.integer(ref$V10[ord]))
    expect_equal(sorted$qValue, ref$V9[ord])

    # name/score/strand still control columns 4-6 only
    bare <- gintervals.import_bed(f, name = FALSE, score = FALSE, strand = FALSE, keep_order = TRUE)
    expect_identical(colnames(bare), c("chrom", "start", "end", "signalValue", "pValue", "qValue", "peak"))
})

test_that("the format is guessed from .narrowPeak(.gz) / .broadPeak(.gz) in any case", {
    lines <- narrowpeak_lines()
    ref <- reference_table(lines)
    files <- c(
        np_gz = write_lines_to(lines, ".narrowPeak.gz", gz = TRUE),
        np_upper = write_lines_to(lines, ".NARROWPEAK"),
        bp = write_lines_to(vapply(strsplit(lines, "\t"), function(x) paste(x[1:9], collapse = "\t"), ""), ".broadPeak"),
        bp_gz = write_lines_to(vapply(strsplit(lines, "\t"), function(x) paste(x[1:9], collapse = "\t"), ""), ".broadPeak.gz", gz = TRUE),
        bed = write_lines_to(lines, ".bed"),
        bed3 = write_lines_to(c("chr1\t300\t400", "chr1\t100\t200"), ".bed")
    )
    withr::defer(unlink(files))

    # plain BED needs only 3 columns, sorted or not
    expect_equal(gintervals.import_bed(files[["bed3"]], keep_order = TRUE)$start, c(300, 100))
    expect_equal(gintervals.import_bed(files[["bed3"]])$start, c(100, 300))

    expect_identical(gintervals.import_bed(files[["np_gz"]], keep_order = TRUE)$peak, as.integer(ref$V10))
    expect_identical(gintervals.import_bed(files[["np_upper"]], keep_order = TRUE)$peak, as.integer(ref$V10))

    for (bp in files[c("bp", "bp_gz")]) {
        out <- gintervals.import_bed(bp, keep_order = TRUE)
        expect_identical(colnames(out), c("chrom", "start", "end", "strand", "name", "score", "signalValue", "pValue", "qValue"))
        expect_equal(out$signalValue, ref$V7)
        expect_equal(out$pValue, ref$V8)
        expect_equal(out$qValue, ref$V9)
    }

    # a .bed name reads BED6 even with 10 columns, unless format says otherwise
    expect_identical(colnames(gintervals.import_bed(files[["bed"]])), c("chrom", "start", "end", "strand", "name", "score"))
    expect_identical(gintervals.import_bed(files[["bed"]], format = "narrowPeak", keep_order = TRUE)$peak, as.integer(ref$V10))
    expect_identical(colnames(gintervals.import_bed(files[["np_upper"]], format = "bed")), c("chrom", "start", "end", "strand", "name", "score"))
})

test_that("format errors: an unknown value, too few columns, a peak that is not a whole number", {
    f6 <- write_lines_to("chr1\t100\t200\tn\t1\t+", ".narrowPeak")
    f9 <- write_lines_to("chr1\t100\t200\tn\t1\t+\t1\t2\t3", ".bed")
    frac <- write_lines_to(c("chr1\t100\t200\tn\t1\t+\t1\t2\t3\t10", "chr1\t300\t400\tm\t1\t+\t1\t2\t3\t12.5"), ".narrowPeak")
    withr::defer(unlink(c(f6, f9, frac)))

    expect_error(
        gintervals.import_bed(f9, format = "narrowpeak"),
        "format must be one of \"bed\", \"narrowPeak\", \"broadPeak\"",
        fixed = TRUE
    )
    expect_error(
        gintervals.import_bed(f6),
        sprintf("BED file %s has 6 columns, but format \"narrowPeak\" needs 10", f6),
        fixed = TRUE
    )
    expect_error(
        gintervals.import_bed(f9, format = "narrowPeak"),
        sprintf("BED file %s has 9 columns, but format \"narrowPeak\" needs 10", f9),
        fixed = TRUE
    )
    expect_error(
        gintervals.import_bed(f6, format = "broadPeak"),
        sprintf("BED file %s has 6 columns, but format \"broadPeak\" needs 9", f6),
        fixed = TRUE
    )
    expect_error(gintervals.import_bed(frac), "got \"12.5\" at data row 2", fixed = TRUE)
})

test_that("the documented summit window has the summit at (start + end) %/% 2", {
    lines <- narrowpeak_lines()
    ref <- reference_table(lines)
    f <- write_lines_to(lines, ".narrowPeak")
    withr::defer(unlink(f))

    peaks <- gintervals.import_bed(f, keep_order = TRUE)
    peaks <- peaks[peaks$peak >= 0, ]
    summit <- peaks$start + peaks$peak
    windows <- transform(peaks, start = summit - 100, end = summit + 100, peak = 100L)

    has_summit <- ref$V10 >= 0
    expect_equal((windows$start + windows$end) %/% 2, ref$V2[has_summit] + ref$V10[has_summit])
    expect_equal(windows$start + windows$peak, ref$V2[has_summit] + ref$V10[has_summit])
    expect_equal(windows$name, ref$V4[has_summit])
    # and the summit base is base 101 of the window's sequence
    summit_base <- gseq.extract(data.frame(chrom = ref$V1[has_summit], start = summit, end = summit + 1))
    expect_equal(substr(gseq.extract(windows[, c("chrom", "start", "end")]), 101, 101), summit_base)
})
