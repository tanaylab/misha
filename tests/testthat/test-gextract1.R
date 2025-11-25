create_isolated_test_db()

test_that("gextract with fixedbin track works", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    expect_regression(gextract("test.fixedbin", intervs), "gextract.fixedbin")
})

test_that("gextract with sparse works", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    expect_regression(gextract("test.sparse", intervs), "gextract.sparse")
})

test_that("gextract with array works", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    expect_regression(gextract("test.array", intervs), "gextract.array")
})

test_that("gextract with computed2d works", {
    intervs <- gscreen("test.rects > 9")
    expect_regression(gextract("test.computed2d", intervs), "gextract.computed2d.1")
    intervs <- gscreen("test.computed2d > 9000000")
    expect_regression(gextract("test.computed2d", intervs), "gextract.computed2d.2")
})

test_that("gextract with .misha$ALLGENOME works", {
    withr::local_options(gmax.data.size = 1e9)
    expect_regression(gextract("test.fixedbin", .misha$ALLGENOME), "gextract.allgenome.fixedbin")
    expect_regression(gextract("test.sparse", .misha$ALLGENOME), "gextract.allgenome.sparse")
    expect_regression(gextract("test.array", .misha$ALLGENOME), "gextract.allgenome.array")
    expect_regression(gextract("test.rects", .misha$ALLGENOME), "gextract.allgenome.rects")
    expect_regression(gextract("test.computed2d", .misha$ALLGENOME), "gextract.allgenome.computed2d")
})

test_that("gextract with sparse track can save to a file", {
    tmp <- tempfile()
    gextract("test.sparse", gintervals(c(1, 2)), file = tmp)
    r <- readr::read_tsv(tmp, col_types = readr::cols(
        chrom = readr::col_character(),
        start = readr::col_double(),
        end = readr::col_double(),
        test.sparse = readr::col_double()
    ))
    r1 <- gextract("test.sparse", gintervals(c(1, 2))) %>%
        mutate(chrom = as.character(chrom)) %>%
        select(-intervalID)
    expect_equal(r, r1, ignore_attr = TRUE)
})


test_that("gextract with specific intervals works", {
    withr::local_options(gmax.data.size = 1e9)
    expect_regression(gextract("test.fixedbin", .misha$ALLGENOME), "gextract.allgenome.fixedbin")
    expect_regression(gextract("test.sparse", .misha$ALLGENOME), "gextract.allgenome.sparse")
    expect_regression(gextract("test.array", .misha$ALLGENOME), "gextract.allgenome.array")
    expect_regression(gextract("test.rects", .misha$ALLGENOME), "gextract.allgenome.rects")
    expect_regression(gextract("test.computed2d", .misha$ALLGENOME), "gextract.allgenome.computed2d")
})
