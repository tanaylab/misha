load_test_db()
test_that("gtrack.create_pwm_energy works", {
    track_name <- random_track_name()
    allgenome <- .misha$ALLGENOME
    withr::defer({
        .misha$ALLGENOME <- allgenome
    })
    .misha$ALLGENOME[[1]]$end <- 1000000
    .misha$ALLGENOME[[2]]$end1 <- 1000000
    .misha$ALLGENOME[[2]]$end2 <- 1000000
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.create_pwm_energy(track_name, "", "misha_motifs", 0, 0.02, iterator = 50)
    r <- gextract(track_name, .misha$ALLGENOME, colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create_pwm_energy")
})
