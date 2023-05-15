
gwilcox("test.fixedbin", 100000, 1000, maxpval=0.000001, intervals=gintervals(c(1, 2), 0, -1))
gwilcox("test.sparse", 100000, 1000, maxpval=0.000001, intervals=gintervals(c(1, 2), 0, -1))
gwilcox("test.array", 100000, 1000, maxpval=0.000001, intervals=gintervals(c(1, 2), 0, -1))
gwilcox("test.rects", 100000, 1000, maxpval=0.000001, intervals=gintervals(c(1, 2), 0, -1))
gwilcox("test.computed2d", 100000, 1000, maxpval=0.000001, intervals=gintervals(c(1, 2), 0, -1))
{
    intervs <- gscreen("test.fixedbin < 0.2", gintervals(c(1, 2), 0, -1))
    gwilcox("test.fixedbin", 100000, 1000, maxpval=0.0001, intervals=intervs)
}
{
    gintervals.rm("test.testintervs", force=T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size=8700)
    try(gwilcox("test.fixedbin", 100000, 1000, maxpval=0.000001, intervals=gintervals(c(1, 2), 0, -1), intervals.set.out = "test.testintervs"), silent=T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force=T)
    r
}