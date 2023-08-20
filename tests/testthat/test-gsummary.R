gsummary("test.fixedbin")
gsummary("test.sparse")
gsummary("test.array")
gsummary("test.rects")
gsummary("test.computed2d")
{
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    gsummary("test.fixedbin", intervs)
}
{
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    gsummary("test.sparse", intervs)
}
{
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    gsummary("test.array", intervs)
}
{
    intervs <- gscreen("test.rects > 40", gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))
    gsummary("test.rects", intervs)
}
{
    intervs <- gscreen("test.computed2d > 4000000", gintervals.2d(chroms1 = c(6, 5), chroms2 = c(8, 9)))
    gsummary("test.computed2d", intervs)
}
{
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 150)
    r <- try(
        {
            gsummary("test.generated_1d_1", "test.generated_1d_2")
        },
        silent = T
    )
    options(gmax.data.size = max.data.size)
    r
}
{
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 15000)
    r <- try(
        {
            gsummary("test.generated_2d_6", "test.generated_2d_5")
        },
        silent = T
    )
    options(gmax.data.size = max.data.size)
    r
}
