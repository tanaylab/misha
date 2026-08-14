create_isolated_test_db()

test_that("gbins.quantiles works", {
    expect_regression(gbins.quantiles("test.fixedbin", c(0, 0.2, 0.3, 0.9, 1.2), "test.sparse", percentiles = c(0.2, 0.5, 0.6), iterator = 10), "gbins.quantiles.1")
    expect_regression(gbins.quantiles("test.fixedbin", c(0, 0.2, 0.3, 0.9, 1.2), "test.sparse", percentiles = c(0.2, 0.5, 0.6), iterator = 100), "gbins.quantiles.2")
})

test_that("gbins.summary works", {
    expect_regression(gbins.summary("test.fixedbin", c(0, 0.2, 0.3, 0.9, 1.2), "test.sparse", iterator = 100), "gbins.summary.1")
})

test_that("gbins.quantiles errors when a positional expr collides with a named expr=", {
    br <- c(0, 0.2, 0.3, 0.9, 1.2)

    # Non-colliding forms still return the same result either way.
    q_named <- gbins.quantiles("test.fixedbin", br, expr = "test.sparse", percentiles = c(0.2, 0.5, 0.6), iterator = 10)
    q_positional <- gbins.quantiles("test.fixedbin", br, "test.sparse", percentiles = c(0.2, 0.5, 0.6), iterator = 10)
    expect_equal(q_named, q_positional)

    # Colliding call: trailing positional expr AND named expr= must error
    # instead of silently discarding the named argument. Before the fix the
    # positional "test.sparse" silently won and "2 * test.sparse" was
    # dropped with no message.
    expect_error(
        gbins.quantiles("test.fixedbin", br, "test.sparse", expr = "2 * test.sparse", percentiles = c(0.2, 0.5, 0.6), iterator = 10),
        "expr"
    )
})

test_that("gbins.summary errors when a positional expr collides with a named expr=", {
    br <- c(0, 0.2, 0.3, 0.9, 1.2)

    s_named <- gbins.summary("test.fixedbin", br, expr = "test.sparse", iterator = 100)
    s_positional <- gbins.summary("test.fixedbin", br, "test.sparse", iterator = 100)
    expect_equal(s_named, s_positional)

    expect_error(
        gbins.summary("test.fixedbin", br, "test.sparse", expr = "2 * test.sparse", iterator = 100),
        "expr"
    )
})

# Multi-dimensional binning with a named expr=. `length(args)` is even in
# that case, and the pair loop used to be driven by `(length(args) - 1) / 2`,
# which yields one pair too few for every even length: a 2-pair call silently
# computed a 1-D binning and returned a plausible (wrong) array. `exprs` and
# `breaks` stayed length-consistent, so nothing downstream complained.
test_that("gbins.quantiles keeps every bin pair when expr is named (2 dimensions)", {
    br1 <- c(0, 0.2, 0.3, 0.9, 1.2)
    br2 <- c(0, 0.5, 1, 2)

    q_positional <- gbins.quantiles("test.fixedbin", br1, "test.sparse", br2, "test.fixedbin",
        percentiles = c(0.2, 0.5), iterator = 100
    )
    q_named <- gbins.quantiles("test.fixedbin", br1, "test.sparse", br2,
        expr = "test.fixedbin", percentiles = c(0.2, 0.5), iterator = 100
    )

    # 4 bins (br1) x 3 bins (br2) x 2 percentiles
    expect_equal(dim(q_positional), c(4L, 3L, 2L))
    expect_equal(dim(q_named), dim(q_positional))
    expect_equal(q_named, q_positional)
})

test_that("gbins.summary keeps every bin pair when expr is named (2 dimensions)", {
    br1 <- c(0, 0.2, 0.3, 0.9, 1.2)
    br2 <- c(0, 0.5, 1, 2)

    s_positional <- gbins.summary("test.fixedbin", br1, "test.sparse", br2, "test.fixedbin",
        iterator = 100
    )
    s_named <- gbins.summary("test.fixedbin", br1, "test.sparse", br2,
        expr = "test.fixedbin", iterator = 100
    )

    expect_equal(dim(s_positional)[1:2], c(4L, 3L))
    expect_equal(dim(s_named), dim(s_positional))
    expect_equal(s_named, s_positional)
})
