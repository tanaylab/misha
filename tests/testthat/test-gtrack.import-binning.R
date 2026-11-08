# Dense import bins checked against a mean computed here from the records, independently of
# the C++ binner (comparing two import routes that share the binner cannot catch its bugs):
# partly covered bins, defval fill, NaN records, records spanning bins, the last partial bin,
# and a binsize larger than the chromosome.

ensure_valid_groot()

expected_bins <- function(recs, sizes, binsize, defval) {
    do.call(rbind, lapply(names(sizes), function(chrom) {
        r <- recs[recs$chrom == chrom, ]
        starts <- seq(0, sizes[[chrom]] - 1, by = binsize)
        ends <- pmin(starts + binsize, sizes[[chrom]])
        v <- mapply(function(b, e) {
            ov <- pmax(0, pmin(r$end, e) - pmax(r$start, b))
            ok <- ov > 0 & !is.nan(r$value)
            s <- sum(ov[ok] * r$value[ok])
            n <- sum(ov[ok])
            if (n < e - b && !is.nan(defval)) {
                s <- s + defval * (e - b - n)
                n <- e - b
            }
            if (n) s / n else NaN
        }, starts, ends)
        data.frame(chrom = chrom, start = starts, end = ends, v = v)
    }))
}

check_import_bins <- function(file, recs, sizes) {
    for (binsize in c(1, 10, 100, 2000)) {
        for (defval in c(NaN, 0, 2.5)) {
            gtrack.import("binned", "", file, binsize = binsize, defval = defval)
            x <- gextract("binned", gintervals.all(), iterator = binsize, colnames = "v")
            want <- expected_bins(recs, sizes, binsize, defval)
            label <- sprintf("%s, binsize %d, defval %g", basename(file), binsize, defval)
            expect_equal(as.character(x$chrom), want$chrom, label = label)
            expect_equal(x$start, want$start, label = label)
            expect_equal(x$end, want$end, label = label)
            expect_equal(x$v, want$v, tolerance = 1e-6, label = label)
            gtrack.rm("binned", force = TRUE)
        }
    }
}

test_that("dense import bins are the defval-filled mean of the records", {
    local_db_state()
    withr::with_tempdir({
        sizes <- c(chr1 = 1003, chr2 = 25) # chr1 ends in a partial bin; chr2 is shorter than some binsizes
        create_test_db("binning_db", chrom_sizes = data.frame(chrom = names(sizes), size = unname(sizes)))
        gdb.init("binning_db")
        recs <- data.frame(
            chrom = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr2"),
            start = c(0, 5, 8, 40, 995, 3),
            end = c(5, 6, 23, 41, 1003, 7),
            value = c(1, NaN, 2.5, -1, 4, 1.5)
        )
        f <- file.path(getwd(), "recs.bedGraph")
        writeLines(sprintf("%s\t%d\t%d\t%s", recs$chrom, recs$start, recs$end, recs$value), f)
        check_import_bins(f, recs, sizes)
    })
})

test_that("bigWig import bins are the defval-filled mean of the records", {
    local_db_state()
    bw <- normalizePath(testthat::test_path("fixtures", "bigwig", "bg.bw"))
    bg <- read.delim(testthat::test_path("fixtures", "bigwig", "bg.bedGraph"),
        header = FALSE,
        col.names = c("chrom", "start", "end", "value")
    )
    withr::with_tempdir({
        sizes <- c(chr1 = 10000, chr2 = 10000)
        create_test_db("binning_db", chrom_sizes = data.frame(chrom = names(sizes), size = unname(sizes)))
        gdb.init("binning_db")
        suppressMessages(check_import_bins(bw, bg[bg$chrom %in% names(sizes), ], sizes))
    })
})

test_that("a value near the float limit survives binning", {
    local_db_state()
    withr::with_tempdir({
        create_test_db("binning_db", chrom_sizes = data.frame(chrom = "chr1", size = 1000))
        gdb.init("binning_db")
        f <- file.path(getwd(), "big.bedGraph")
        writeLines(c("chr1\t0\t32\t3.4e38", "chr1\t32\t64\t1e37", "chr1\t64\t96\t1"), f)
        # a float product value * 32 would overflow to Inf, which reads back as NaN
        float32 <- function(x) readBin(writeBin(x, raw(), size = 4), "numeric", size = 4, n = length(x))
        for (binsize in c(1, 32)) {
            gtrack.import("big", "", f, binsize = binsize)
            got <- gextract("big", gintervals("chr1", 0, 96), iterator = 32)$big
            expect_identical(got, float32(c(3.4e38, 1e37, 1)), label = paste("binsize", binsize))
            gtrack.rm("big", force = TRUE)
        }
    })
})

# Records a model writes: 32 bp each, end = start + 32 also at the chromosome end; chr1 has one
# more that starts past its end
past_end_records <- function(sizes) {
    recs <- do.call(rbind, lapply(names(sizes), function(chrom) {
        start <- seq(0, if (chrom == "chr1") 1024 else 576, by = 32)
        data.frame(chrom = chrom, start = start, end = start + 32)
    }))
    recs$value <- seq_len(nrow(recs)) / 2
    recs
}

# Captures the warnings of expr; returns them
collect_warnings <- function(expr) {
    warns <- character()
    withCallingHandlers(expr, warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    warns
}

test_that("records past the chromosome end are clipped, or dropped, with one warning, in every text format", {
    local_db_state()
    withr::with_tempdir({
        sizes <- c(chr1 = 1000, chr2 = 600)
        create_test_db("clip_db", chrom_sizes = data.frame(chrom = names(sizes), size = unname(sizes)))
        gdb.init("clip_db")
        recs <- past_end_records(sizes)
        # chr1 992-1024 and chr2 576-608 are clipped, chr1 1024-1056 is dropped
        expect_identical(sum(recs$end > sizes[recs$chrom] & recs$start < sizes[recs$chrom]), 2L)
        expect_identical(sum(recs$start >= sizes[recs$chrom]), 1L)
        kept <- recs[recs$start < sizes[recs$chrom], ]
        kept$end <- pmin(kept$end, sizes[kept$chrom])

        files <- list(
            bedGraph = sprintf("%s\t%d\t%d\t%s", recs$chrom, recs$start, recs$end, recs$value),
            tab = c("chrom\tstart\tend\tvalue", sprintf("%s\t%d\t%d\t%s", recs$chrom, recs$start, recs$end, recs$value)),
            fixedStep = unlist(lapply(names(sizes), function(chrom) {
                c(sprintf("fixedStep chrom=%s start=1 step=32 span=32", chrom), recs$value[recs$chrom == chrom])
            })),
            variableStep = unlist(lapply(names(sizes), function(chrom) {
                r <- recs[recs$chrom == chrom, ]
                c(sprintf("variableStep chrom=%s span=32", chrom), sprintf("%d\t%s", r$start + 1, r$value))
            }))
        )
        first_line <- c(bedGraph = 32, tab = 33, fixedStep = 33, variableStep = 33)
        for (fmt in names(files)) {
            f <- file.path(getwd(), paste0(fmt, ".txt"))
            writeLines(files[[fmt]], f)
            for (binsize in c(0, 32, 10)) {
                label <- paste(fmt, "binsize", binsize)
                warns <- collect_warnings(gtrack.import("clipped", "", f, binsize = binsize))
                expect_length(warns, 1)
                expect_match(warns, paste0(
                    "2 record\\(s\\) running past the end of their chromosome were clipped to it, and 1 record\\(s\\) ",
                    "starting at or past the chromosome end were dropped; the first: line ", first_line[[fmt]],
                    ": \\(chr1, 992, 1024\\), chromosome size 1000"
                ), label = label)
                if (binsize == 0) {
                    x <- gextract("clipped", gintervals.all(), iterator = "clipped", colnames = "v")
                    expect_equal(x$end, kept$end, label = label)
                    expect_equal(x$v, kept$value, label = label)
                } else {
                    x <- gextract("clipped", gintervals.all(), iterator = binsize, colnames = "v")
                    expect_equal(x$v, expected_bins(kept, sizes, binsize, NaN)$v, tolerance = 1e-6, label = label)
                }
                gtrack.rm("clipped", force = TRUE)
            }
        }
    })
})
