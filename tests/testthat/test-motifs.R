{
    allgenome <- ALLGENOME
    ALLGENOME[[1]]$end <- 1000000
    ALLGENOME[[2]]$end1 <- 1000000
    ALLGENOME[[2]]$end2 <- 1000000
    gtrack.rm("test.tmptrack", force = T)
    gtrack.create_pwm_energy("test.tmptrack", "", "misha_motifs", 0, 0.02, iterator = 50)
    r <- gextract("test.tmptrack", ALLGENOME)
    ALLGENOME <- allgenome
    gtrack.rm("test.tmptrack", force = T)
    r
}