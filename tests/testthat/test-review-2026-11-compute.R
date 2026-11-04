# Regression tests for the compute and virtual-track findings of the 2026-11 review of 5.14.4.

create_isolated_test_db()

test_that("gcor is accurate at a large offset, single process and multitask", {
    iv <- gintervals(c(1, 2), 0, 2e5)
    x <- gextract("test.fixedbin", iv, iterator = "test.fixedbin")$test.fixedbin
    x <- x[!is.na(x)] + 1e8
    y <- 2 * x - 1e8
    for (mt in c(FALSE, TRUE)) {
        withr::local_options(list(gmultitasking = mt, gmin.scope4process = 1000, gmax.processes = 8))
        r <- gcor("test.fixedbin + 1e8", "test.fixedbin * 2 + 1e8", intervals = iv, iterator = "test.fixedbin", details = TRUE)
        expect_equal(r[1, "sd1"], sd(x), tolerance = 1e-6)
        expect_equal(r[1, "sd2"], sd(y), tolerance = 1e-6)
        expect_equal(r[1, "cov"], cov(x, y), tolerance = 1e-6)
        expect_equal(r[1, "cor"], 1, tolerance = 1e-9)
    }
})

test_that("summaries give a NaN standard deviation when a value is infinite, as sd() does", {
    expr <- "ifelse(test.fixedbin > 0.15, Inf, test.fixedbin)"
    iv <- gintervals(1, 0, 1000)
    x <- gextract(expr, iv, iterator = 100, colnames = "v")$v
    expect_true(any(is.infinite(x)))
    expect_true(is.nan(sd(x)))

    expect_true(is.nan(gsummary(expr, iv, iterator = 100)[["Std dev"]]))
    expect_true(is.nan(gintervals.summary(expr, iv, iterator = 100)[["Std dev"]]))
    expect_true(is.nan(gbins.summary("test.fixedbin", breaks = c(0, 0.5), expr = expr, intervals = iv, iterator = 100)[1, "Std dev"]))
    # a finite set still has a finite one
    expect_equal(gsummary("test.fixedbin", iv, iterator = 100)[["Std dev"]], sd(gextract("test.fixedbin", iv, iterator = 100)$test.fixedbin))
})

test_that("2D gintervals.annotate with tie_method and maxneighbors > 1 gives a clear error", {
    q <- gintervals.2d(1, 100, 200, 1, 300, 400)
    a <- gintervals.2d(1, c(500, 1000), c(600, 1100), 1, c(700, 1200), c(800, 1300))
    a$lab <- c("x", "y")
    expect_error(
        gintervals.annotate(q, a, annotation_columns = "lab", maxneighbors = 2, tie_method = "min.start"),
        "1D intervals only"
    )
    expect_equal(gintervals.annotate(q, a, annotation_columns = "lab", maxneighbors = 2)$lab, c("x", "y"))
})

test_that("first, last, sample, lse and their positions work on array tracks", {
    iv <- gintervals("chr1", c(100, 200, 300, 500, 700), c(150, 260, 400, 600, 750))
    a <- c(1, 2, 3, 4, 5)
    b <- c(11, 12, 13, 14, 15)
    vals <- (a + b) / 2
    f <- tempfile(fileext = ".tsv")
    on.exit(unlink(f), add = TRUE)
    write.table(cbind(iv, a = a, b = b), f, sep = "\t", quote = FALSE, row.names = FALSE)
    nm <- paste0("test.arr_fl_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(nm, force = TRUE))
    suppressMessages(gtrack.array.import(nm, "", f))
    withr::defer(remove_all_vtracks())

    q <- gintervals("chr1", c(0, 120, 230, 410), c(1000, 330, 240, 490))
    ext <- function(func) {
        gvtrack.create("arr_v", nm, func)
        gextract("arr_v", q, iterator = q, colnames = "v")$v
    }
    expect_equal(ext("first"), c(vals[1], vals[1], vals[2], NaN))
    expect_equal(ext("first.pos.abs"), c(100, 120, 230, NaN))
    expect_equal(ext("last"), c(vals[5], vals[3], vals[2], NaN))
    expect_equal(ext("last.pos.abs"), c(700, 300, 230, NaN))
    expect_equal(ext("last.pos.relative"), c(700, 180, 0, NaN))
    expect_equal(ext("lse"), c(log(sum(exp(vals))), log(sum(exp(vals[1:3]))), vals[2], NaN), tolerance = 1e-6)
    s <- ext("sample")
    expect_true(s[1] %in% vals && s[2] %in% vals[1:3])
    expect_equal(s[3:4], c(vals[2], NaN))
    sp <- ext("sample.pos.abs")
    expect_true(sp[1] %in% iv$start)
    expect_true(sp[2] %in% c(120, 200, 300))
    expect_equal(sp[3:4], c(230, NaN))
})
