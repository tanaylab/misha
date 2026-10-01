create_isolated_test_db()

# test_that("gcluster.run works", {
#     v <- 17
#     r <- gcluster.run(2 + 3, gsummary("test.fixedbin+v+qwe"), 3 + 4, 4 + 5, gsummary("test.sparse-v"), gsummary("test.rects*v"))
#     expect_regression(list(r[[1]]$retv, r[[2]]$retv, r[[3]]$retv, r[[4]]$retv, r[[5]]$retv, r[[6]]$retv), "gcluster.run")
# })

# test_that("gcluster.run works (2)", {
#     v <- 17
#     r <- gcluster.run(2 + 3, gsummary("test.fixedbin+v+qwe"), 3 + 4, 4 + 5, gsummary("test.sparse-v"), gsummary("test.rects*v"), max.jobs = 2)
#     expect_regression(list(r[[1]]$retv, r[[2]]$retv, r[[3]]$retv, r[[4]]$retv, r[[5]]$retv, r[[6]]$retv), "gcluster.run.2")
# })

test_that("a gcluster.run job restores the caller's root, working dir, datasets and vtracks", {
    local_db_state()
    withr::with_tempdir({
        create_test_db("working_db")
        create_test_db("dataset_db")
        create_test_db("other_db", chrom_sizes = data.frame(chrom = "chrX", size = 5000))
        gsetroot("working_db")
        gdir.create("sub")
        gintervals.save("sub.w_set", gintervals(1, 0, 100))
        gdataset.load("dataset_db")
        gdir.cd("sub")
        gvtrack.create("v", gintervals(1, 0, 100), "distance")
        state <- function() {
            list(
                groot = .misha$GROOT, cwd = gdir.cwd(), datasets = gdataset.ls(),
                intervals = gintervals.ls(), vtracks = gvtrack.ls()
            )
        }
        expected <- state()
        expect_equal(expected$intervals, "w_set")
        save(.misha, file = "misha") # as gcluster.run does

        # a fresh job is rooted elsewhere (the example db) and has no vtracks
        gsetroot("other_db")
        rm("GVTRACKS", envir = .misha)
        .gcluster.restore_db("misha")

        expect_equal(state(), expected)
    })
})

test_that("a gcluster.run job restores a caller whose working dir is the db root", {
    local_db_state()
    withr::with_tempdir({
        create_test_db("working_db")
        create_test_db("other_db", chrom_sizes = data.frame(chrom = "chrX", size = 5000))
        gsetroot("working_db")
        gdir.cd("..")
        expected <- gdir.cwd()
        save(.misha, file = "misha")

        gsetroot("other_db")
        .gcluster.restore_db("misha")

        expect_equal(.misha$GROOT, normalizePath("working_db"))
        expect_equal(gdir.cwd(), expected)
    })
})

test_that("a gcluster.run job fails when the caller's db is not reachable", {
    local_db_state()
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        save(.misha, file = "misha")
        unlink("working_db", recursive = TRUE) # not mounted on the node
        gdb.init_examples(file.path(getwd(), "job"))
        expect_error(.gcluster.restore_db("misha"), "does not exist")
    })
})

test_that("a gcluster.run job keeps its own example db when the caller's is not reachable", {
    local_db_state()
    withr::with_tempdir({
        gdb.init_examples(file.path(getwd(), "caller"))
        save(.misha, file = "misha")
        unlink("caller", recursive = TRUE) # the caller's tempdir, as seen from another host
        gdb.init_examples(file.path(getwd(), "job"))
        job_root <- .misha$GROOT
        .gcluster.restore_db("misha")
        expect_equal(.misha$GROOT, job_root)
    })
})
