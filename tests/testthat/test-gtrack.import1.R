create_isolated_test_db()

test_that("import and extract from s_7_export.txt", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    intervs <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2)))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.rm(tmptrack, force = TRUE)
    gtrack.import_mappedseq(tmptrack, "", "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/input_files/s_7_export.txt", remove.dups = FALSE)
    r <- gextract(tmptrack, intervs, colnames = "test.tmptrack")
    expect_regression(r, "track.import_mappedseq.s_7_export")
})

test_that("import and extract from sample-small.sam", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    intervs <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2)))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.import_mappedseq(tmptrack, "", "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/input_files/sample-small.sam", cols.order = NULL, remove.dups = FALSE)
    r <- gextract(tmptrack, intervs, colnames = "test.tmptrack")
    expect_regression(r, "track.import_mappedseq.sample_small_sam")
})

test_that("import with pileup and binsize from s_7_export.txt", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    intervs <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2)))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.import_mappedseq(tmptrack, "", "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/input_files/s_7_export.txt", remove.dups = FALSE, pileup = 180, binsize = 50)
    r <- gextract(tmptrack, intervs, colnames = "test.tmptrack")
    expect_regression(r, "track.import_pileup_binsize")
})
