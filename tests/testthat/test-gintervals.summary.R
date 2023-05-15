gintervals.summary("test.fixedbin", gintervals(c(1, 2), 0, -1))
gintervals.summary("test.sparse", gintervals(c(1, 2), 0, -1))
gintervals.summary("test.rects", gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))
gintervals.summary("test.computed2d", gintervals.2d(chroms1 = c(6,1,5), chroms2 = c(8,1,9)))
{
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
	set.seed(60427)
	intervs <- intervs[sample(nrow(intervs)), ]
    gintervals.summary("test.fixedbin", intervs)
}
{
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    gintervals.summary("test.sparse", intervs)
}
{
    intervs <- gscreen("test.rects > 40", gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))
	set.seed(60427)
	intervs <- intervs[sample(nrow(intervs)), ]
    gintervals.summary("test.rects", intervs)
}
{
    intervs <- gscreen("test.computed2d > 4000000", gintervals.2d(chroms1 = c(6,1,5), chroms2 = c(8,1,9)))
    gintervals.summary("test.computed2d", intervs)
}
gintervals.summary("test.generated_1d_1", intervals=giterator.intervals("test.generated_1d_2"), iterator=giterator.intervals("test.generated_1d_1"))
gintervals.summary("test.generated_1d_1", intervals=giterator.intervals("test.generated_1d_2"), iterator="test.bigintervs_1d_1")
gintervals.summary("test.generated_1d_1", intervals=giterator.intervals("test.generated_1d_2"), iterator="test.generated_1d_1")
gintervals.summary("test.generated_1d_1", intervals="test.bigintervs_1d_2", iterator=giterator.intervals("test.generated_1d_1"))
gintervals.summary("test.generated_1d_1", intervals="test.bigintervs_1d_2", iterator="test.bigintervs_1d_1")
gintervals.summary("test.generated_1d_1", intervals="test.bigintervs_1d_2", iterator="test.generated_1d_1")
gintervals.summary("test.generated_1d_1", intervals="test.generated_1d_2", iterator=giterator.intervals("test.generated_1d_1"))
gintervals.summary("test.generated_1d_1", intervals="test.generated_1d_2", iterator="test.bigintervs_1d_1")
gintervals.summary("test.generated_1d_1", intervals="test.generated_1d_2", iterator="test.generated_1d_1")
gintervals.summary("test.generated_2d_5", intervals=giterator.intervals("test.generated_2d_6"), iterator=giterator.intervals("test.generated_2d_5"))
gintervals.summary("test.generated_2d_5", intervals=giterator.intervals("test.generated_2d_6"), iterator="test.bigintervs_2d_5")
gintervals.summary("test.generated_2d_5", intervals=giterator.intervals("test.generated_2d_6"), iterator="test.generated_2d_5")
gintervals.summary("test.generated_2d_5", intervals="test.bigintervs_2d_6", iterator=giterator.intervals("test.generated_2d_5"))
gintervals.summary("test.generated_2d_5", intervals="test.bigintervs_2d_6", iterator="test.bigintervs_2d_5")
gintervals.summary("test.generated_2d_5", intervals="test.bigintervs_2d_6", iterator="test.generated_2d_5")
gintervals.summary("test.generated_2d_5", intervals="test.generated_2d_6", iterator=giterator.intervals("test.generated_2d_5"))
gintervals.summary("test.generated_2d_5", intervals="test.generated_2d_6", iterator="test.bigintervs_2d_5")
gintervals.summary("test.generated_2d_5", intervals="test.generated_2d_6", iterator="test.generated_2d_5")
{
    gintervals.rm("test.testintervs", force=T)
    intervs1 <- gscreen("test.fixedbin > 0.2 & test.fixedbin < 0.3", gintervals(c(1, 2, 3), 0, -1))
    intervs2 <- gscreen("test.fixedbin > 0.25 & test.fixedbin < 0.35", gintervals(c(1, 2), 0, -1))
    set.seed(60427)
    intervs2 <- intervs2[sample(nrow(intervs2)), ]
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = nrow(intervs2) - 100)
    try(gintervals.summary("test.fixedbin", intervals=intervs2, iterator=intervs1, intervals.set.out = "test.testintervs"), silent=T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force=T)
    r
}
{
    gintervals.rm("test.testintervs", force=T)
    intervs <- gscreen("test.rects > 40", gintervals.2d(c(1, 2, 5, 8), 0, -1))
    set.seed(60427)
    intervs <- intervs[sample(nrow(intervs)),]
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = nrow(intervs) - 100)
    try(gintervals.summary("test.rects", intervs, iterator = c(1, 1), intervals.set.out = "test.testintervs"), silent=T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force=T)
    r
}