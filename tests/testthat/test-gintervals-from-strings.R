test_that("gintervals.from_strings parses a single chrom:start-end string", {
    res <- gintervals.from_strings("chr1:100-200")
    expect_equal(nrow(res), 1)
    expect_equal(as.character(res$chrom), "chr1")
    expect_equal(res$start, 100)
    expect_equal(res$end, 200)
    expect_false("strand" %in% colnames(res))
})

test_that("gintervals.from_strings treats a chromosome-only string as the full extent", {
    res <- gintervals.from_strings("chr1")
    expect_equal(nrow(res), 1)
    expect_equal(as.character(res$chrom), "chr1")
    expect_equal(res$start, 0)
    all_iv <- gintervals.all()
    expect_equal(res$end, all_iv$end[as.character(all_iv$chrom) == "chr1"])
})

test_that("gintervals.from_strings accepts a dot separator", {
    res <- gintervals.from_strings("chr1:100..200")
    expect_equal(res$start, 100)
    expect_equal(res$end, 200)
})

test_that("gintervals.from_strings parses strand suffixes", {
    res <- gintervals.from_strings(c("chr1:100-200:+", "chr2:50-60:-"))
    expect_equal(nrow(res), 2)
    expect_true("strand" %in% colnames(res))
    expect_setequal(res$strand, c(1, -1))
})

test_that("gintervals.from_strings vectorizes over multiple strings", {
    res <- gintervals.from_strings(c("chr1:0-100", "chr2:0-200"))
    expect_equal(nrow(res), 2)
    expect_setequal(as.character(res$chrom), c("chr1", "chr2"))
})

test_that("gintervals.from_strings errors on a malformed coordinate string", {
    # Has a colon-structure but non-numeric coordinates -> cannot parse.
    expect_error(gintervals.from_strings("chr1:abc-def"), "Invalid interval string")
})
