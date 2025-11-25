create_isolated_test_db()

test_that("gtrack.create_pwm_energy works", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    allgenome <- .misha$ALLGENOME
    .misha$ALLGENOME[[1]]$end <- 1000000
    .misha$ALLGENOME[[2]]$end1 <- 1000000
    .misha$ALLGENOME[[2]]$end2 <- 1000000
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create_pwm_energy(tmptrack, "", "misha_motifs", 0, 0.02, iterator = 50)
    r <- gextract(tmptrack, .misha$ALLGENOME, colnames = "test.tmptrack")
    .misha$ALLGENOME <- allgenome
    expect_regression(r, "gtrack.create_pwm_energy")
})
