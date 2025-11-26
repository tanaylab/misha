create_isolated_test_db()

test_that("gtrack.var works", {
    expect_regression(gtrack.var.get("test.fixedbin", "pv.percentiles"), "test.fixedbin.pv.percentiles")
    expect_error(gtrack.var.get("test.fixedbin", "blablablabla"))
    expect_error(gtrack.var.get("test.blablablabla", "pv.percentiles"))
    expect_equal(gtrack.var.ls("test.fixedbin"), "pv.percentiles")
    expect_error(gtrack.var.ls("test.blablablabla"))
})

test_that("gtrack.var.set works", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create_sparse(tmptrack, "", gintervals(c(1, 2)), 1:2)
    gtrack.var.set(tmptrack, "testvar1", 1)
    gtrack.var.set(tmptrack, "testvar2", 1)
    r1 <- gtrack.var.ls(tmptrack)
    gtrack.var.rm(tmptrack, "testvar2")
    r2 <- gtrack.var.ls(tmptrack)
    expect_equal(list(r1, r2), list(c("testvar1", "testvar2"), "testvar1"))
})

test_that("gtrack.rm works", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    expect_warning(gtrack.var.rm("test.fixedbin", "blablablabla"))
    expect_error(gtrack.var.rm("test.blablablabla", "pv.percentiles"))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create_sparse(tmptrack, "", gintervals(c(1, 2)), 1:2)
    gtrack.var.set(tmptrack, "testvar1", 1:10)
    r1 <- gtrack.var.get(tmptrack, "testvar1")
    gtrack.var.set(tmptrack, "testvar1", 21:30)
    r2 <- gtrack.var.get(tmptrack, "testvar1")
    expect_equal(list(r1, r2), list(1:10, 21:30))

    expect_error(gtrack.var.set("test.blablablabla", "testvar1", 1:10))
})
