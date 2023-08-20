{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = "test.sparse")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = "test.array")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "blabla")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = "test.sparse")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "blabla")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = "test.array")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "avg", 10)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "avg")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.sparse", func = "avg")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array", func = "avg")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects", func = "avg")
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d", func = "avg")
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "max", 10)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "max")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.sparse", func = "max")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array", func = "max")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects", func = "max")
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d", func = "max")
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "min", 10)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "min")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.sparse", func = "min")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array", func = "min")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects", func = "min")
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d", func = "min")
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "nearest", 10)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "nearest")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.sparse", func = "nearest")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array", func = "nearest")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects", func = "nearest")
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d", func = "nearest")
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "stddev", 10)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "stddev")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.sparse", func = "stddev")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array", func = "stddev")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects", func = "stddev")
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d", func = "stddev")
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "sum", 10)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "sum")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.sparse", func = "sum")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array", func = "sum")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects", func = "sum")
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d", func = "sum")
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "quantile")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "quantile", params = 0.5)
    gvtrack.create("v2", "test.fixedbin", func = "quantile", params = 0.9)
    r <- gextract("v1", "v2", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "quantile", params = 0.5)
    gvtrack.create("v2", "test.fixedbin", func = "quantile", params = 0.9)
    r <- gextract("v1", "v2", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "quantile", params = 0.5)
    gvtrack.create("v2", "test.fixedbin", func = "quantile", params = 0.9)
    r <- gextract("v1", "v2", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "quantile", params = 0.5)
    gvtrack.create("v2", "test.fixedbin", func = "quantile", params = 0.9)
    r <- gextract("v1", "v2", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "global.percentile", 10)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "global.percentile")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.sparse", func = "global.percentile")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array", func = "global.percentile")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects", func = "global.percentile")
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d", func = "global.percentile")
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "global.percentile.max", 10)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "global.percentile.max")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.sparse", func = "global.percentile.max")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array", func = "global.percentile.max")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects", func = "global.percentile.max")
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d", func = "global.percentile.max")
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "global.percentile.min", 10)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "global.percentile.min")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 233)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.sparse", func = "global.percentile.min")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array", func = "global.percentile.min")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 10000)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects", func = "global.percentile.min")
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d", func = "global.percentile.min")
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = c(2000000, 3000000))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    gvtrack.create("v1", intervs, "blabla")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 533)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    gvtrack.create("v1", intervs)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 533)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    intervs$strand <- 1
    gvtrack.create("v1", intervs, "distance")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 533)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    intervs$strand <- -1
    gvtrack.create("v1", intervs, "distance")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 533)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    gvtrack.create("v1", intervs, "distance")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = "test.sparse")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    gvtrack.create("v1", intervs, "distance")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = "test.array")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    gvtrack.create("v1", intervs, "distance")
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = "test.rects")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(6, 1, 5), 0, -1))
    gvtrack.create("v1", intervs, "distance")
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = "test.computed2d")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    intervs$strand <- 1
    gvtrack.create("v1", intervs, "distance", 200)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 533)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    intervs$strand <- -1
    gvtrack.create("v1", intervs, "distance", 200)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 533)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    intervs$strand <- 1
    gvtrack.create("v1", intervs, "distance.center", 100)
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 533)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    gvtrack.create("v1", intervs, "distance.center")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 533)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    intervs$strand <- 1
    gvtrack.create("v1", intervs, "distance.center")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 533)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    intervs$strand <- -1
    gvtrack.create("v1", intervs, "distance.center")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = 533)
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    gvtrack.create("v1", intervs, "distance.center")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = "test.sparse")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    gvtrack.create("v1", intervs, "distance.center")
    r <- gextract("v1", gintervals(c(1, 2)), iterator = "test.array")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(1, 3), 0, -1))
    gvtrack.create("v1", intervs, "distance.center")
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(1, 4), 3000000, -1), iterator = "test.rects")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    intervs <- gscreen("test.fixedbin > 0.5", gintervals(c(6, 1, 5), 0, -1))
    gvtrack.create("v1", intervs, "distance.center")
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = "test.computed2d")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    gvtrack.array.slice("v1")
    r <- gvtrack.ls()
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.sparse")
    gvtrack.array.slice("v1")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1")
    r <- gvtrack.ls()
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects")
    gvtrack.array.slice("v1")
    r <- gvtrack.ls()
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d")
    gvtrack.array.slice("v1")
    r <- gvtrack.ls()
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", -50)
    r <- gvtrack.ls()
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", 50)
    r <- gvtrack.ls()
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", "blabla")
    r <- gvtrack.ls()
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", c("col1", "col1", "col3", "col5"))
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", c("col1", "col3", "col5"))
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", c("col1", "col3", "col5"), "blabla")
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", c("col1", "col3", "col5"), "avg")
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", c("col1", "col3", "col5"), "avg", 25)
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", c("col1", "col3", "col5"), "min")
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", c("col1", "col3", "col5"), "max")
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", c("col1", "col3", "col5", "col6", "col8"), "stddev")
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", c("col1", "col3", "col5", "col6", "col8"), "sum")
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", c("col1", "col3", "col5", "col6", "col8"), "quantile")
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.array")
    gvtrack.array.slice("v1", c("col1", "col3", "col5", "col6", "col8"), "quantile", 0.4)
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin", func = "max")
    r <- gvtrack.info("v1")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    gvtrack.iterator("v1", dim = 1, 830, -724)
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    gvtrack.iterator("v1", dim = 0, sshift = -130, eshift = 224)
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    gvtrack.iterator("v1", sshift = -130, eshift = 224)
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects")
    gvtrack.iterator("v1", sshift = -130, eshift = 224)
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d")
    gvtrack.iterator("v1", sshift = -130, eshift = 224)
    r <- gextract("v1", gintervals(c(1, 2)))
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    gvtrack.iterator("v1", dim = 1)
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(2, 4), 3000000, -1), iterator = "test.rects")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    gvtrack.iterator("v1", dim = 2)
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(2, 4), 3000000, -1), iterator = "test.rects")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    gvtrack.iterator("v1", dim = 1, sshift = -130, eshift = 224)
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(2, 4), 3000000, -1), iterator = "test.rects")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    gvtrack.iterator("v1", dim = 1)
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = "test.computed2d")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    gvtrack.iterator("v1", dim = 2)
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = "test.computed2d")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    gvtrack.iterator("v1", dim = 1, sshift = -130, eshift = 224)
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = "test.computed2d")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.fixedbin")
    gvtrack.iterator.2d("v1")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects")
    gvtrack.iterator.2d("v1", sshift1 = -1000000, eshift1 = -500000, sshift2 = 2000000, eshift2 = 2800000)
    r <- gextract("v1", gintervals.2d(chroms1 = c(1, 3), 3000000, -1, chroms2 = c(2, 4), 3000000, -1), iterator = "test.rects")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.computed2d")
    gvtrack.iterator.2d("v1", sshift1 = -1000000, eshift1 = -500000, sshift2 = 2000000, eshift2 = 2800000)
    r <- gextract("v1", gintervals.2d(chroms1 = c(6, 1, 5), 3000000, -1, chroms2 = c(8, 1, 9), 3000000, -1), iterator = "test.computed2d")
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects")
    gvtrack.create("v2", "test.sparse")
    gvtrack.create("v3", "test.computed2d")
    r <- gvtrack.ls()
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v2", "test.sparse")
    gvtrack.rm("v1")
    r <- gvtrack.ls()
    remove.all.vtracks()
    r
}
{
    remove.all.vtracks()
    gvtrack.create("v1", "test.rects")
    gvtrack.create("v2", "test.sparse")
    r1 <- gvtrack.ls()
    gvtrack.rm("v1")
    r2 <- gvtrack.ls()
    remove.all.vtracks()
    list(r1, r2)
}
