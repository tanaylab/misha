{
    intervs <- gscreen("test.fixedbin > 0.14", gintervals(c(1, 2), 0, -1))
    set.seed(60427)
    intervs <- intervs[sample(nrow(intervs)), ]
    gpartition("test.fixedbin", seq(0, 1, by = 0.1), intervals = intervs)
}
gpartition("test.rects", seq(50, 100, by = 1), intervals = gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))
gpartition("test.computed2d", seq(5000000, 10000000, by = 1000000), intervals = gintervals.2d(chroms1 = c(6, 5), chroms2 = c(8, 9)))
{
    gintervals.rm("test.testintervs", force = T)
    intervs <- gscreen("test.fixedbin > 0.14", gintervals(c(1, 2), 0, -1))
    set.seed(60427)
    intervs <- intervs[sample(nrow(intervs)), ]
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 2000000)
    try(gpartition("test.fixedbin", seq(0, 1, by = 0.1), intervs, intervals.set.out = "test.testintervs"), silent = T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}
{
    gintervals.rm("test.testintervs", force = T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 18000)
    try(gpartition("test.rects", seq(0, 100, by = 1), gintervals.2d(chroms1 = c(6, 3), chroms2 = c(2, 4)), intervals.set.out = "test.testintervs"), silent = T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}
