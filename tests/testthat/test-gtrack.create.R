
{
    gtrack.rm("aaaaaaaaaaaaa.bbbbbbbbbbb", force=T)
    gtrack.create("aaaaaaaaaaaaa.bbbbbbbbbbb", "", "test.fixedbin")
    r <- gextract("aaaaaaaaaaaaa.bbbbbbbbbbb", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("aaaaaaaaaaaaa.bbbbbbbbbbb", force=T)
    r
}
{
    gtrack.rm("test.tmptrack", force=T)
    gtrack.create("test.tmptrack", "", "test.fixedbin+1")
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("test.tmptrack", force=T)
    r
}
{
    gtrack.rm("test.tmptrack", force=T)
    gtrack.create("test.tmptrack", "", "test.fixedbin+1", iterator="test.sparse")
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("test.tmptrack", force=T)
    r
}
{
    gtrack.rm("test.tmptrack", force=T)
    gtrack.create("test.tmptrack", "", "test.fixedbin+1", iterator="test.array")
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("test.tmptrack", force=T)
    r
}
{
    gtrack.rm("test.tmptrack", force=T)
    gtrack.create("test.tmptrack", "", "test.rects+10")
    r <- gextract("test.tmptrack", gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))
    gtrack.rm("test.tmptrack", force=T)
    r
}
{
    intervs <- giterator.intervals("test.sparse", gintervals(c(1, 3, 4)))
    gtrack.rm("test.tmptrack", force=T)
    gtrack.create("test.tmptrack", "", "test.fixedbin+1", iterator=intervs)
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("test.tmptrack", force=T)
    r
}
{
    intervs <- giterator.intervals("test.array", gintervals(c(1, 3, 4)))
    gtrack.rm("test.tmptrack", force=T)
    gtrack.create("test.tmptrack", "", "test.fixedbin+1", iterator=intervs)
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("test.tmptrack", force=T)
    r
}
{
    intervs <- giterator.intervals("test.rects", gintervals.2d(chroms1 = c(2, 3, 5), chroms2 = c(2, 4, 7)))
    gtrack.rm("test.tmptrack", force=T)
    gtrack.create("test.tmptrack", "", "test.rects+10", iterator=intervs)
    r <- gextract("test.tmptrack", gintervals.2d(chroms1 = c(2, 3, 3), chroms2 = c(2, 3, 4)))
    gtrack.rm("test.tmptrack", force=T)
    r
}

{
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    gtrack.rm("test.tmptrack", force = T)
    gtrack.create_sparse("test.tmptrack", "", intervs, 1:dim(intervs)[1])
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000000))
    gtrack.rm("test.tmptrack", force = T)
    r
}
{
    gtrack.rm("test.tmptrack", force = T)
    gtrack.create_sparse("test.tmptrack", "", gintervals.2d(c(1, 2)), 1:2)
    gtrack.rm("test.tmptrack", force = T)
}


{
    intervs <- gscreen("test.rects > 80", gintervals.2d(c(1, 2)))
    gtrack.rm("test.tmptrack", force = T)
    gtrack.2d.create("test.tmptrack", "", intervs, 1:dim(intervs)[1])
    r <- gextract("test.tmptrackv", gintervals.2d(c(1, 2, 3)))
    gtrack.rm("test.tmptrack", force = T)
    r
}
{
    intervs <- gintervals(c(1, 2))
    gtrack.rm("test.tmptrack", force = T)
    gtrack.2d.create("test.tmptrack", "", intervs, 1:dim(intervs)[1])
    r <- gextract("test.tmptrack", gintervals.2d(c(1, 2, 3)))
    gtrack.rm("test.tmptrack", force = T)
    r
}
