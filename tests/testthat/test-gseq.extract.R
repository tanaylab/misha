{
    intervs <- gscreen("test.fixedbin > 0.6", gintervals(c(1, 2, 3)))
    gseq.extract(intervs)
}
{
    intervs <- gscreen("test.fixedbin > 0.6", gintervals(c(1, 2, 3)))
    intervs$strand <- -1
    gseq.extract(intervs)
}
gseq.extract(gintervals.2d(1, 10, 100, 2, 20, 300))
