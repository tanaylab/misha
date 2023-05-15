{
    m1 = matrix(1:15, nrow=5, ncol=3)
    glookup(m1, "test.fixedbin", seq(0.1, 0.2, length.out=6), "test.sparse", seq(0.25, 0.48, length.out=4), gintervals(c(1,2)), iterator="test.fixedbin")
}
{
    m1 = matrix(1:15, nrow=5, ncol=3)
    glookup(m1, "test.fixedbin", seq(0.1, 0.2, length.out=6), "test.sparse", seq(0.25, 0.48, length.out=4), gintervals(c(1,2)), force.binning=FALSE, iterator="test.fixedbin")
}
{
    m1 = matrix(1:15, nrow=5, ncol=3)
    glookup(m1, "test.rects", seq(50, 100, length.out=6), "test.rects / 2", seq(0, 40, length.out=4), gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)), force.binning=FALSE)
}
{
    m1 = matrix(1:15, nrow=5, ncol=3)
    glookup(m1, "test.computed2d", seq(5000000, 10000000, length.out=6), "test.computed2d / 2", seq(0, 4000000, length.out=4), gintervals.2d(chroms1 = c(6, 5), chroms2 = c(8, 9)), force.binning=FALSE)
}