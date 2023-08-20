gsegment("test.fixedbin", 10000, maxpval = 0.000001)
gsegment("test.sparse", 10000, maxpval = 0.000001)
gsegment("test.array", 10000, maxpval = 0.000001)
gsegment("test.rects", 10000, maxpval = 0.000001)
gsegment("test.computed2d", 10000, maxpval = 0.000001)
{
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    gsegment("test.fixedbin*2", 10000, maxpval = 0.000001, intervals = intervs)
}
{
    gintervals.rm("test.testintervs", force = T)
    intervs <- gscreen("test.fixedbin > 0.14", gintervals(c(1, 2), 0, -1))
    set.seed(60427)
    intervs <- intervs[sample(nrow(intervs)), ]
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 3200)
    try(gsegment("test.sparse", 10000, maxpval = 0.0001, iterator = 50, intervals.set.out = "test.testintervs"), silent = T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}
