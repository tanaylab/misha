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
