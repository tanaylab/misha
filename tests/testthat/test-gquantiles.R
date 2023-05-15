{
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    gquantiles("test.fixedbin+0.2", percentile=c(0.5, 0.3, 0.2, 0.9), intervs)
}
gquantiles("test.rects", percentile=c(0.5, 0.3, 0.2, 0.9, 0.999), gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))
gquantiles("test.computed2d", percentile=c(0.5, 0.3, 0.2, 0.9, 0.999), gintervals.2d(chroms1 = c(6, 5), chroms2 = c(8, 9)))
{
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    gquantiles("test.fixedbin+0.2", percentile=c(0.5, 0.999))
}