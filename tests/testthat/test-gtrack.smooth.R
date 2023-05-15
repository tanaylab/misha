{
    gtrack.rm("test.tmptrack", force=T)
    gtrack.smooth("test.tmptrack", "", "test.fixedbin", 10000, alg="LINEAR_RAMP")
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("test.tmptrack", force=T)
    r
}
{
    gtrack.rm("test.tmptrack", force=T)
    gtrack.smooth("test.tmptrack", "", "test.fixedbin", 10000, alg="MEAN")
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("test.tmptrack", force=T)
    r
}
{
    gtrack.rm("test.tmptrack", force=T)
    gtrack.smooth("test.tmptrack", "", "test.sparse", 10000, alg="LINEAR_RAMP")
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("test.tmptrack", force=T)
    r
}
{
    gtrack.rm("test.tmptrack", force=T)
    gtrack.smooth("test.tmptrack", "", "test.array", 10000, alg="LINEAR_RAMP")
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("test.tmptrack", force=T)
    r
}
{
    gtrack.rm("test.tmptrack", force=T)
    gtrack.smooth("test.tmptrack", "", "test.rects", 10000, alg="LINEAR_RAMP")
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("test.tmptrack", force=T)
    r
}