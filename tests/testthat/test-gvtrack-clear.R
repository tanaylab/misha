test_that("gvtrack.clear removes all virtual tracks for the current working directory", {
    gvtrack.create("vt_clear_1", "test.fixedbin", "max")
    gvtrack.create("vt_clear_2", "test.fixedbin", "avg")
    expect_true(all(c("vt_clear_1", "vt_clear_2") %in% gvtrack.ls()))

    gvtrack.clear()
    expect_null(gvtrack.ls())
})

test_that("gvtrack.clear is a no-op when there are no virtual tracks", {
    gvtrack.clear()
    expect_null(gvtrack.ls())
    expect_silent(gvtrack.clear())
})
