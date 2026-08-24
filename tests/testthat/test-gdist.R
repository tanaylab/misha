# Re-roots into a database of its own; hand this process's overlay back to the next
# file in the parallel worker (helper-test_db.R).
restore_groot_on_exit()

create_isolated_test_db()

test_that("gdist works", {
    expect_regression(gdist("test.fixedbin", seq(0, 1, by = 0.001)), "gdist.1")
    expect_regression(gdist("test.fixedbin", seq(0.2, 1, by = 0.001)), "gdist.2")
    expect_regression(gdist("test.fixedbin", seq(0, 1, by = 0.01), "test.fixedbin+0.3", seq(0, 1, by = 0.01)), "gdist.3")
})

test_that("gscreen + gdist works (1)", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    expect_regression(gdist("test.fixedbin", seq(0.2, 1, by = 0.001), intervs), "gdist_gscreen.1")
})

test_that("gscreen + gdist 2d works (1)", {
    intervs <- gscreen("test.rects > 9")
    expect_regression(gdist("test.rects", seq(8, 10, by = 0.1), intervs), "gdist_gscreen_2d.1")
})

test_that("gscreen + gdist 2d works (2)", {
    intervs <- gscreen("test.computed2d > 9500000")
    expect_regression(gdist("test.computed2d", seq(9000000, 10000000, by = 100000), intervs), "gdist_gscreen_2d.2")
})

test_that("gdist works when dataframe = TRUE", {
    dst <- gdist("test.fixedbin", c(0, 0.2, 0.5, 1), dataframe = TRUE)
    dst_non_df <- gdist("test.fixedbin", c(0, 0.2, 0.5, 1), dataframe = FALSE)
    # Convert array to vector in same order as dataframe
    dst_non_df_vec <- as.vector(dst_non_df)
    expect_true(all(dst$n == dst_non_df_vec))
    expect_true(all(dst$test.fixedbin == dimnames(dst_non_df)[[1]]))
    expect_true(is.factor(dst$test.fixedbin))
    expect_true(all(levels(dst$test.fixedbin) == dimnames(dst_non_df)[[1]]))
})

test_that("gdist 2d works when dataframe = TRUE", {
    dst <- gdist("test.fixedbin", c(0, 0.2, 0.5, 1), "test.fixedbin+0.3", c(0, 0.2, 0.5, 1), iterator = 100, dataframe = TRUE)
    expect_equal(colnames(dst), c("test.fixedbin", "test.fixedbin+0.3", "n"))
    dst_non_df <- gdist("test.fixedbin", c(0, 0.2, 0.5, 1), "test.fixedbin+0.3", c(0, 0.2, 0.5, 1), iterator = 100, dataframe = FALSE)
    # Convert array to vector in same order as dataframe
    dst_non_df_vec <- as.vector(dst_non_df)
    expect_true(all(dst$n == dst_non_df_vec))
    expect_true(is.factor(dst$`test.fixedbin+0.3`))
    expect_true(is.factor(dst$test.fixedbin))
    expect_true(all(levels(dst$`test.fixedbin+0.3`) == colnames(dst_non_df)))
    expect_true(all(levels(dst$test.fixedbin) == rownames(dst_non_df)))
})

test_that("gdist dataframe = TRUE with names", {
    dst <- gdist("test.fixedbin", c(0, 0.2, 0.5, 1), "test.fixedbin+0.3", c(0, 0.2, 0.5, 1), iterator = 100, dataframe = TRUE, names = c("mytrack1", "mytrack2"))
    expect_equal(colnames(dst), c("mytrack1", "mytrack2", "n"))

    dst1 <- gdist("test.fixedbin", c(0, 0.2, 0.5, 1), "test.fixedbin+0.3", c(0, 0.2, 0.5, 1), iterator = 100, dataframe = TRUE)
    colnames(dst1) <- c("mytrack1", "mytrack2", "n")
    expect_equal(dst, dst1)
})

test_that("gdist errors when a positional scope collides with a named intervals=", {
    restore_groot_on_exit() # this test re-roots to the example db; put GROOT back for tests after it
    gdb.init_examples()
    br <- c(0, 0.2, 0.5, 1)
    small <- gintervals(1, 0, 50000)
    big <- gintervals(1, 0, 500000)

    # Non-colliding forms - verified numbers from the review's reproducer.
    dist_named_big <- sum(gdist("dense_track", br, intervals = big))
    dist_named_small <- sum(gdist("dense_track", br, intervals = small))
    expect_equal(dist_named_big, 6546)
    expect_equal(dist_named_small, 901)

    # Positional scope alone (no named intervals=) still works, and matches
    # the equivalent named call.
    expect_equal(sum(gdist("dense_track", br, small)), dist_named_small)

    # Colliding call: positional scope AND named intervals= must error
    # instead of silently discarding the named argument. Before the fix this
    # returned 901 (using `small`, discarding `intervals = big`) with no
    # message at all.
    expect_error(
        gdist("dense_track", br, small, intervals = big),
        "intervals"
    )
})
