gscreen("test.fixedbin")
gscreen("test.rects")
gscreen("2 * test.fixedbin+0.2 > 0.4")
gscreen("2 * test.rects+1 > 100")
gscreen("test.fixedbin < -1")
gscreen("test.rects < -1")
gscreen("test.fixedbin > 0.2", gintervals(1, c(0, 2000000, 4000000), c(1000000, 3000000, 5000000)))
gscreen("test.rects > 40", gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))