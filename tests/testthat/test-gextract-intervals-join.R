create_isolated_test_db()

# Build a small 1D intervals data frame with metadata of every supported type.
make_1d_intervals_with_meta <- function() {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    n <- nrow(intervs)
    intervs$gene_id <- sprintf("g%04d", seq_len(n))
    intervs$score <- as.numeric(seq_len(n)) / 10
    intervs$rank <- as.integer(seq_len(n))
    intervs$is_top <- seq_len(n) %% 2L == 0L
    intervs$category <- factor(rep(c("a", "b", "c"), length.out = n))
    intervs
}

test_that("intervals_join defaults to 'id' (no behavior change)", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    a <- gextract("test.fixedbin", intervs)
    b <- gextract("test.fixedbin", intervs, intervals_join = "id")
    expect_equal(a, b)
    expect_true("intervalID" %in% colnames(a))
})

test_that("intervals_join='none' drops intervalID and preserves everything else", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    base <- gextract("test.fixedbin", intervs)
    none <- gextract("test.fixedbin", intervs, intervals_join = "none")
    expect_false("intervalID" %in% colnames(none))
    expect_equal(nrow(base), nrow(none))
    expect_equal(base[, setdiff(colnames(base), "intervalID")], none)
})

test_that("intervals_join='intervals' attaches input intervals' columns (1D, all metadata types)", {
    intervs <- make_1d_intervals_with_meta()
    res <- gextract("test.fixedbin", intervs, intervals_join = "intervals")

    expect_false("intervalID" %in% colnames(res))
    # iterator coords stay as-is, input coords get suffix '1'
    expect_true(all(c(
        "chrom", "start", "end", "test.fixedbin",
        "chrom1", "start1", "end1",
        "gene_id", "score", "rank", "is_top", "category"
    ) %in% colnames(res)))

    # Reference: do the join in R via positional indexing.
    base <- gextract("test.fixedbin", intervs)
    ref_intervs <- intervs
    conflict <- intersect(colnames(ref_intervs), colnames(base))
    colnames(ref_intervs)[match(conflict, colnames(ref_intervs))] <-
        paste0(conflict, "1")
    ref_intervs$intervalID <- seq_len(nrow(ref_intervs))
    ref <- base[order(base$intervalID), ]
    ref <- cbind(ref, ref_intervs[ref$intervalID, setdiff(colnames(ref_intervs), "intervalID"), drop = FALSE])
    ref$intervalID <- NULL
    rownames(ref) <- NULL

    res <- res[order(res$chrom, res$start), ]
    rownames(res) <- NULL
    ref <- ref[, colnames(res)]

    expect_equal(res$gene_id, ref$gene_id)
    expect_equal(res$score, ref$score)
    expect_equal(res$rank, ref$rank)
    expect_equal(res$is_top, ref$is_top)
    expect_equal(as.character(res$category), as.character(ref$category))
    expect_equal(as.character(res$chrom1), as.character(ref$chrom1))
    expect_equal(res$start1, ref$start1)
    expect_equal(res$end1, ref$end1)
})

test_that("intervals_join='intervals' works for 2D scope", {
    intervs <- gscreen("test.rects > 40", gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))
    intervs$tag <- sprintf("t%03d", seq_len(nrow(intervs)))
    intervs$weight <- as.numeric(seq_len(nrow(intervs)))

    res <- gextract("test.rects", intervs, intervals_join = "intervals")
    expect_false("intervalID" %in% colnames(res))
    # All 6 2D coord columns conflict: input chrom1/start1/... all get suffix '1'.
    expect_true(all(c(
        "chrom1", "start1", "end1", "chrom2", "start2", "end2",
        "chrom11", "start11", "end11", "chrom21", "start21", "end21",
        "tag", "weight"
    ) %in% colnames(res)))
})

test_that("multitask and serial produce identical results with intervals_join='intervals'", {
    intervs <- make_1d_intervals_with_meta()

    withr::local_options(gmultitasking = FALSE)
    res_serial <- gextract("test.fixedbin", intervs, intervals_join = "intervals")
    res_serial <- res_serial[order(res_serial$chrom, res_serial$start), ]
    rownames(res_serial) <- NULL

    options(gmultitasking = TRUE)
    res_mt <- gextract("test.fixedbin", intervs, intervals_join = "intervals")
    res_mt <- res_mt[order(res_mt$chrom, res_mt$start), ]
    rownames(res_mt) <- NULL

    expect_equal(colnames(res_serial), colnames(res_mt))
    expect_equal(res_serial, res_mt)
})

test_that("intervals_join='none' still allows the track-parallel strategy", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    withr::local_options(gmultitasking.strategy = "tracks", gmultitasking = TRUE)
    # 8 expressions is the threshold .gmultitasking_strategy uses for "tracks"; we
    # also force the option above to ensure the strategy is exercised.
    exprs <- rep("test.fixedbin", 9)
    cn <- paste0("v", seq_along(exprs))
    res <- gextract(exprs, intervals = intervs, intervals_join = "none", colnames = cn)
    expect_false("intervalID" %in% colnames(res))
    expect_true(all(c("chrom", "start", "end", cn) %in% colnames(res)))
})

test_that("intervals_join='intervals' errors on file output", {
    intervs <- make_1d_intervals_with_meta()
    tmp <- tempfile()
    expect_error(
        gextract("test.fixedbin", intervs, intervals_join = "intervals", file = tmp),
        regexp = "intervals_join"
    )
})

test_that("intervals_join='intervals' errors on intervals.set.out", {
    intervs <- make_1d_intervals_with_meta()
    expect_error(
        gextract("test.fixedbin", intervs,
            intervals_join = "intervals", intervals.set.out = "tmpset_join"
        ),
        regexp = "intervals_join"
    )
})

test_that("intervals_join='intervals' rejects unsupported column types", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    intervs$bad <- replicate(nrow(intervs), list(1L), simplify = FALSE)
    expect_error(
        gextract("test.fixedbin", intervs, intervals_join = "intervals"),
        regexp = "bad"
    )
})

test_that("intervals_join validates allowed values", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    expect_error(gextract("test.fixedbin", intervs, intervals_join = "bogus"))
})
