create_isolated_test_db()
test_that("gextract works with gscreen 1d", {
    intervs1 <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    intervs2 <- gscreen("test.fixedbin < 0.4", gintervals(c(1, 2)))
    expect_regression(gextract("test.fixedbin", intervs1, iterator = intervs2), "gextract.gscreen.1d.fixedbin")
    expect_regression(gextract("test.sparse", intervs1, iterator = intervs2), "gextract.gscreen.1d.sparse")
    expect_regression(gextract("test.array", intervs1, iterator = intervs2), "gextract.gscreen.1d.array")
    expect_error(gextract("test.rects", intervs1, iterator = intervs2))
    expect_error(gextract("test.computed2d", intervs1, iterator = intervs2))
})

test_that("gextract works with gscreen 2d", {
    intervs1 <- gscreen("test.rects > 40", gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))
    intervs2 <- gscreen("test.rects < 50", gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))
    expect_error(gextract("test.fixedbin", intervs1, iterator = intervs2))
    expect_error(gextract("test.sparse", intervs1, iterator = intervs2))
    expect_error(gextract("test.array", intervs1, iterator = intervs2))
    expect_regression(gextract("test.rects", intervs1, iterator = intervs2), "gextract.gscreen.2d.rects")
    expect_regression(gextract("test.computed2d", intervs1, iterator = intervs2), "gextract.gscreen.2d.computed2d")
})

test_that("gextract works with 2d .misha$ALLGENOME", {
    intervs <- rbind(
        gintervals.2d(1, 10, 100, 1, 10, 100),
        gintervals.2d(1, 400, 500, 1, 400, 500),
        gintervals.2d(2, 600, 700, 2, 600, 700),
        gintervals.2d(1, 200, 300, 2, 200, 300),
        gintervals.2d(1, 7000, 9100, "X", 7000, 9100),
        gintervals.2d(2, 9000, 18000, 2, 9000, 18000),
        gintervals.2d(1, 30000, 31000, 1, 30000, 31000),
        gintervals.2d(2, 1130, 15000, 1, 1130, 15000),
        gintervals.2d(1, 1100, 1120, 1, 1100, 1120),
        gintervals.2d(1, 1000, 1100, 2, 1000, 1100)
    )
    expect_regression(gextract("test.rects", intervals = .misha$ALLGENOME, iterator = intervs), "gextract.2d.ALLGENOME.rects")
})

test_that("gextract with giterator.intervals works", {
    expect_regression(gextract("test.generated_1d_1", intervals = giterator.intervals("test.generated_1d_2"), iterator = giterator.intervals("test.generated_1d_1")), "gextract.giterator.intervals.1")
    expect_regression(gextract("test.generated_1d_1", intervals = giterator.intervals("test.generated_1d_2"), iterator = "test.bigintervs_1d_1"), "gextract.giterator.intervals.2")
    expect_regression(gextract("test.generated_1d_1", intervals = giterator.intervals("test.generated_1d_2"), iterator = "test.generated_1d_1"), "gextract.giterator.intervals.3")
    expect_regression(gextract("test.generated_1d_1", intervals = "test.bigintervs_1d_2", iterator = giterator.intervals("test.generated_1d_1")), "gextract.giterator.intervals.4")
    expect_regression(gextract("test.generated_1d_1", intervals = "test.bigintervs_1d_2", iterator = "test.bigintervs_1d_1"), "gextract.giterator.intervals.5")
    expect_regression(gextract("test.generated_1d_1", intervals = "test.bigintervs_1d_2", iterator = "test.generated_1d_1"), "gextract.giterator.intervals.6")
    expect_regression(gextract("test.generated_1d_1", intervals = "test.generated_1d_2", iterator = giterator.intervals("test.generated_1d_1")), "gextract.giterator.intervals.7")
    expect_regression(gextract("test.generated_1d_1", intervals = "test.generated_1d_2", iterator = "test.bigintervs_1d_1"), "gextract.giterator.intervals.8")
    expect_regression(gextract("test.generated_1d_1", intervals = "test.generated_1d_2", iterator = "test.generated_1d_1"), "gextract.giterator.intervals.9")
    expect_regression(gextract("test.generated_2d_5", intervals = giterator.intervals("test.generated_2d_6"), iterator = giterator.intervals("test.generated_2d_5")), "gextract.giterator.intervals.10")
    expect_regression(gextract("test.generated_2d_5", intervals = giterator.intervals("test.generated_2d_6"), iterator = "test.bigintervs_2d_5"), "gextract.giterator.intervals.11")
    expect_regression(gextract("test.generated_2d_5", intervals = giterator.intervals("test.generated_2d_6"), iterator = "test.generated_2d_5"), "gextract.giterator.intervals.12")
    expect_regression(gextract("test.generated_2d_5", intervals = "test.bigintervs_2d_6", iterator = giterator.intervals("test.generated_2d_5")), "gextract.giterator.intervals.13")
    expect_regression(gextract("test.generated_2d_5", intervals = "test.bigintervs_2d_6", iterator = "test.bigintervs_2d_5"), "gextract.giterator.intervals.14")
    expect_regression(gextract("test.generated_2d_5", intervals = "test.bigintervs_2d_6", iterator = "test.generated_2d_5"), "gextract.giterator.intervals.15")
    expect_regression(gextract("test.generated_2d_5", intervals = "test.generated_2d_6", iterator = giterator.intervals("test.generated_2d_5")), "gextract.giterator.intervals.16")
    expect_regression(gextract("test.generated_2d_5", intervals = "test.generated_2d_6", iterator = "test.bigintervs_2d_5"), "gextract.giterator.intervals.17")
    expect_regression(gextract("test.generated_2d_5", intervals = "test.generated_2d_6", iterator = "test.generated_2d_5"), "gextract.giterator.intervals.18")
})
