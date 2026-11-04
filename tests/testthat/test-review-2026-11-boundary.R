# Regression tests for the R/C++ boundary findings of the 2026-11 review of 5.14.4:
# argument conversion and validation, and process state left after an error.

create_isolated_test_db()

scope <- gintervals(1, 0, 20000)

# The bin glookup() puts each value in: (b1, b2], ...; force.binning puts values outside the
# breaks into the first or the last bin.
lookup_bin <- function(vals, breaks) {
    pmin(pmax(findInterval(vals, breaks, left.open = TRUE), 1), length(breaks) - 1)
}

test_that("a chrom column that is not character or factor is rejected without leaving misha's state behind", {
    bad <- list(
        data.frame(chrom = 1L, start = 0, end = 100),
        data.frame(chrom = I(list("chr1")), start = 0, end = 100),
        c("chr1", "chr2")
    )
    for (intervs in bad) {
        expect_error(gextract("test.fixedbin", intervs), "Invalid format of intervals")
        expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    }
    iv2 <- gintervals.2d(1, 0, 100, 1, 0, 100)
    iv2$chrom2 <- 1L
    expect_error(gextract("test.rects", iv2), "Invalid format of intervals")
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
})
