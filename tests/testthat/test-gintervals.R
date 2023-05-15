gintervals(c(1,2), c(0, 50, 2000, 50, 10000, 1500), c(100, 1300, 3000, 300, 11000, 2300))
gintervals(c(1,2), c(0, 50, 2000, 50, 10000, 1500), c(100, 1300, 3000, 300, 11000, 2300), "a")
gintervals(c(1,2), c(0, 50, 2000, 50, 10000, 1500), c(100, 1300, 3000, 300, 11000, 2300), 25)
gintervals(c(1,2), c(0, 50, 2000, 50, 10000, 1500), c(100, 1300, 3000, 300, 11000, 2300), 0.1)
gintervals(c(1,2), c(0, 50, 2000, 50, 10000, 1500), c(100, 1300, 3000, 300, 11000, 2300), -1)
gintervals.2d(c(1), c(0, 50), c(100, 200), c(3), c(0, 50), c(400, 600))
gintervals.2d(c(1,2), c(0, 1000, 2000, 50, 10000, 1500), c(100, 1300, 3000, 300, 11000, 2300), c(3,4), c(10, 1010, 2010, 60, 10010, 1510), c(110, 1310, 3010, 310, 11010, 2310))
gintervals.2d.all()
gintervals.all()
{
    intervs <- giterator.intervals("test.rects_big_rects", gintervals.2d(c(2,3)), iterator=c(123450, 97891))
    set.seed(60427)
    intervs <- intervs[sample(nrow(intervs)),]
    gintervals.2d.band_intersect(intervs, band=c(-198743, 23456))
}
{
    gintervals.rm("test.testintervs", force=T)
    intervs <- giterator.intervals("test.rects_big_rects", gintervals.2d(c(2,3)), iterator=c(123450, 97891))
    set.seed(60427)
    intervs <- intervs[sample(nrow(intervs)),]
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 16000)
    try(gintervals.2d.band_intersect(intervs, band=c(-198743, 23456), intervals.set.out = "test.testintervs"), silent=T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force=T)
    r
}
gintervals.chrom_sizes("bigintervs1d")
gintervals.chrom_sizes("bigintervs2d")
gintervals.chrom_sizes("test.tss")
gintervals.chrom_sizes("test.array")
gintervals.chrom_sizes("test.fixed_bin")
gintervals.chrom_sizes("test.sparse")
gintervals.chrom_sizes("test.rects")
{
    intervs1 <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    intervs2 <- gscreen("test.fixedbin > 0.4", gintervals(c(1, 2), 0, -1))
    gintervals.diff(intervs1, intervs2)
}
{
    intervs1 <- gscreen("test.fixedbin > 0.1 & test.fixedbin < 0.3", gintervals(c(1, 2)))
    intervs2 <- gscreen("test.fixedbin > 0.2 & test.fixedbin < 0.4", gintervals(c(1, 2)))
    intervs3 <- gscreen("(test.fixedbin > 0.25 & test.fixedbin < 0.32) | test.fixedbin > 0.35", gintervals(c(1, 2)))
    gintervals.diff(rbind(intervs1, intervs2), intervs3)
}
{
    gintervals.rm("test.testintervs", force = T)
    intervs1 <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2, 4, 8, 9), 0, -1))
    intervs2 <- gscreen("test.fixedbin > 0.4", gintervals(c(1, 2, 4, 7, 9), 0, -1))
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 1300000)
    try(gintervals.diff(intervs1, intervs2, intervals.set.out = "test.testintervs"), silent = T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}
gintervals.load("test.foodgene")
gintervals.load("bigintervs1d")
gintervals.load("bigintervs1d", chrom = 2)
gintervals.load("bigintervs2d")
gintervals.load("bigintervs2d", chrom1 = 2, chrom2 = 2)
gintervals.load("test.rects", chrom1 = 1, chrom2 = 2)
gintervals.load("test.generated_1d_1", chrom = 13)
gintervals.load("test.generated_1d_1", chrom = 12)
gintervals.load("test.generated_2d_5", chrom1 = 1, chrom2 = 2)
gintervals.load("test.generated_2d_5", chrom1 = 1, chrom2 = 3)
{
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 200000)
    r <- try(gintervals.load("test.generated_1d_2"), silent = T)
    options(gmax.data.size = max.data.size)
    r
}
{
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 100000)
    r <- try(gintervals.load("test.generated_1d_2", chrom = 1), silent = T)
    options(gmax.data.size = max.data.size)
    r
}
{
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 1000000)
    r <- try(gintervals.load("test.generated_2d_6"), silent = T)
    options(gmax.data.size = max.data.size)
    r
}
{
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 100)
    r <- try(gintervals.load("test.generated_2d_6", chrom1 = 1, chrom2 = 13), silent = T)
    options(gmax.data.size = max.data.size)
    r
}
gintervals.ls()
gintervals.ls("test")


{
    gintervals.rm("test.testintervs", force = T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 6000000)
    try(gextract("test.fixedbin", gintervals(c(1, 2)), intervals.set.out = "test.testintervs"), silent = T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}
{
    gintervals.rm("test.testintervs", force = T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 300000)
    try(gextract("test.rects", gintervals.2d(c(1, 2)), intervals.set.out = "test.testintervs"), silent = T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}

{
    gintervals.rm("test.testintervs", force = T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 6000000)
    try(glookup(m1, "test.fixedbin", seq(0.1, 0.2, length.out = 6), "test.sparse", seq(0.25, 0.48, length.out = 4), gintervals(c(1, 2)), iterator = "test.fixedbin", intervals.set.out = "test.testintervs"), silent = T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}
{
    gintervals.rm("test.testintervs", force = T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 100000)
    try(glookup(m1, "test.computed2d", seq(5000000, 10000000, length.out = 6), "test.computed2d / 2", seq(0, 4000000, length.out = 4), gintervals.2d(chroms1 = c(6, 5), chroms2 = c(8, 9)), force.binning = FALSE, intervals.set.out = "test.testintervs"), silent = T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}



{
    gintervals.rm("test.testintervs", force = T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 1000000)
    try(gscreen("2 * test.sparse+0.2 > 0.4", intervals.set.out = "test.testintervs"), silent = T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}
{
    gintervals.rm("test.testintervs", force = T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 130000)
    try(gscreen("test.rects > 40", intervals.set.out = "test.testintervs"), silent = T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}



gintervals.diff("test.bigintervs_1d_1", "test.bigintervs_1d_2")
gintervals.exists("test.tss")
gintervals.exists("test.blablablablabla")
gintervals.exists("blablablablabla.blablablablabla")
{
    gintervals.force_range(rbind(data.frame(chrom="chr1", start=10, end=100),
	                             data.frame(chrom="chr1", start=300, end=200),
								 data.frame(chrom="chr1", start=-100, end=50),
								 data.frame(chrom="chr1", start=-100, end=-30),
								 data.frame(chrom="chr1", start=-30, end=-100),
								 data.frame(chrom="chr1", start=100, end=1e+09),
								 data.frame(chrom="chr1", start=1e+09, end=10 + 1e+09),
								 data.frame(chrom="chr1", start=10 + 1e+09, end=1e+09)))
}
{
    gintervals.force_range(rbind(data.frame(chrom1="chr1", start1=10, end1=100, chrom2="chr2", start2=10, end2=100),
	                             data.frame(chrom1="chr1", start1=300, end1=200, chrom2="chr2", start2=300, end2=200),
								 data.frame(chrom1="chr1", start1=-100, end1=50, chrom2="chr2", start2=-100, end2=50),
								 data.frame(chrom1="chr1", start1=-100, end1=-30, chrom2="chr2", start2=-100, end2=-30),
								 data.frame(chrom1="chr1", start1=-30, end1=-100, chrom2="chr2", start2=-30, end2=-100),
								 data.frame(chrom1="chr1", start1=100, end1=1e+09, chrom2="chr2", start2=100, end2=1e+09),
								 data.frame(chrom1="chr1", start1=1e+09, end1=10 + 1e+09, chrom2="chr2", start2=1e+09, end2=10 + 1e+09),
								 data.frame(chrom1="chr1", start1=10 + 1e+09, end1=1e+09, chrom2="chr2", start2=10 + 1e+09, end2=1e+09)))
}

{
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = 100)
	try(r<-gintervals.is.bigset("bigintervs1d"), silent = T)
	options(gmax.data.size = max.data.size)
	r
}
{
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = 100)
	try(r<-gintervals.is.bigset("bigintervs2d"), silent = T)
	options(gmax.data.size = max.data.size)
	r
}
{
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = 100)
	try(r<-gintervals.is.bigset("test.tss"), silent = T)
	options(gmax.data.size = max.data.size)
	r
}


{
    gintervals.rm("test.testintervs", force = TRUE)
    intervs1 <- gextract("test.fixedbin", gintervals(c(1, 2), 1000, 4000))
    intervs2 <- gextract("test.fixedbin", gintervals(c(2, "X"), 2000, 5000))
    gintervals.save("test.testintervs", intervs2)
    r <- gintervals.rbind(intervs1, "test.testintervs")
    gintervals.rm("test.testintervs", force = TRUE)
    r
}
{
    gintervals.rm("test.testintervs", force = T)
    gintervals.save("test.testintervs", gintervals(c(1, 2)))
    r1 <- gintervals.ls()
    gintervals.rm("test.testintervs", force = T)
    r2 <- gintervals.ls()
    list(r1, r2)
}
{
    gintervals.rm("test.testintervs", force = T)
    gintervals.save("test.testintervs", gintervals(c(1, 2), 1000, 2000))
    gintervals.rm("test.testintervs", force = T)
    gextract("test.fixedbin", "test.testintervs")
}
gintervals.rm("test.aaaaaaaaaaaaaaaaaaa", force = T)
gintervals.rm("test.aaaaaaaaaaaaaaaaaaa")
{
    gintervals.rm("test.testintervs", force = T)
    r1 <- gintervals.ls()
    gintervals.save("test.testintervs", gintervals(c(1, 2), 1000, 2000))
    r2 <- gintervals.ls()
    gintervals.rm("test.testintervs", force = T)
    list(r1, r2)
}
{
    gintervals.rm("test.testintervs", force = T)
    gintervals.save("test.testintervs", gintervals(c(1, 2), 1000, 2000))
    r <- gextract("test.fixedbin", "test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}


{
    gintervals.rm("test.testintervs", force=T)
    intervs1 <- gscreen("test.fixedbin > 0.1 & test.fixedbin < 0.3", gintervals(c(1, 2, 4, 8, 9), 0, -1))
    intervs2 <- gscreen("test.fixedbin < 0.2", gintervals(c(1, 2, 4, 7, 9), 0, -1))
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 1000000)
    try(gintervals.union(intervs1, intervs2, intervals.set.out = "test.testintervs"), silent=T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force=T)
    r
}
{
    gintervals.rm("test.testintervs", force=T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100)
    gintervals.save("test.testintervs", "bigintervs1d")
    options(gmax.data.size = max.data.size)
    gintervals.update("test.testintervs", gintervals(1), chrom=1)
}
{
    gintervals.rm("test.testintervs", force=T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100)
    gintervals.save("test.testintervs", "bigintervs1d")
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs", chrom=2)
    gintervals.update("test.testintervs", r[c(2,3),], chrom1=1)
}
{
    gintervals.rm("test.testintervs", force=T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100)
    gintervals.save("test.testintervs", "bigintervs1d")
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs", chrom=2)
    gintervals.update("test.testintervs", r[c(2,3),], chrom=2)
    r <- list(gintervals.load("test.testintervs", chrom=2), gintervals.chrom_sizes("test.testintervs"))
    gintervals.rm("test.testintervs", force=T)
    r
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100)
	gintervals.save("test.testintervs", "bigintervs1d")
	options(gmax.data.size = max.data.size)
	gintervals.update("test.testintervs", NULL, chrom=2)
	r <- list(gintervals.load("test.testintervs", chrom=2), gintervals.chrom_sizes("test.testintervs"))
	gintervals.rm("test.testintervs", force=T)
	r
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100)
	gintervals.save("test.testintervs", "bigintervs2d")
	options(gmax.data.size = max.data.size)
	gintervals.update("test.testintervs", gintervals.2d(1), chrom1=1)
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100)
	gintervals.save("test.testintervs", "bigintervs2d")
	options(gmax.data.size = max.data.size)
	r <- gintervals.load("test.testintervs", chrom1=2, chrom2=2)
	gintervals.update("test.testintervs", r[c(2,3),], chrom=1)
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100)
	gintervals.save("test.testintervs", "bigintervs2d")
	options(gmax.data.size = max.data.size)
	r <- gintervals.load("test.testintervs", chrom1=2, chrom2=2)
	gintervals.update("test.testintervs", r[c(2,3),], chrom1=2, chrom2=2)
	r <- list(gintervals.load("test.testintervs", chrom1=2, chrom2=2), gintervals.chrom_sizes("test.testintervs"))
	gintervals.rm("test.testintervs", force=T)
	r
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100)
	gintervals.save("test.testintervs", "bigintervs2d")
	options(gmax.data.size = max.data.size)
	gintervals.update("test.testintervs", NULL, chrom1=2, chrom2=2)
	r <- list(gintervals.load("test.testintervs", chrom1=2, chrom2=2), gintervals.chrom_sizes("test.testintervs"))
	gintervals.rm("test.testintervs", force=T)
	r
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs1d")$size) - 100)
	gintervals.save("test.testintervs", "bigintervs1d")
	options(gmax.data.size = max.data.size)
	gintervals.update("test.testintervs", gintervals(1), chrom=1)
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs1d")$size) - 100)
	gintervals.save("test.testintervs", "bigintervs1d")
	options(gmax.data.size = max.data.size)
	r <- gintervals.load("test.testintervs", chrom=2)
	gintervals.update("test.testintervs", r[c(2,3),], chrom1=1)
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs1d")$size) - 100)
	gintervals.save("test.testintervs", "bigintervs1d")
	options(gmax.data.size = max.data.size)
	r <- gintervals.load("test.testintervs", chrom=2)
	gintervals.update("test.testintervs", r[c(2,3),], chrom=2)
	r <- list(gintervals.load("test.testintervs", chrom=2), gintervals.chrom_sizes("test.testintervs"))
	gintervals.rm("test.testintervs", force=T)
	r
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs1d")$size) - 100)
	gintervals.save("test.testintervs", "bigintervs1d")
	options(gmax.data.size = max.data.size)
	gintervals.update("test.testintervs", NULL, chrom=2)
	r <- list(gintervals.load("test.testintervs", chrom=2), gintervals.chrom_sizes("test.testintervs"))
	gintervals.rm("test.testintervs", force=T)
	r
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs2d")$size) - 100)
	gintervals.save("test.testintervs", "bigintervs2d")
	options(gmax.data.size = max.data.size)
	gintervals.update("test.testintervs", gintervals.2d(1), chrom1=1)
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs2d")$size) - 100)
	gintervals.save("test.testintervs", "bigintervs2d")
	options(gmax.data.size = max.data.size)
	r <- gintervals.load("test.testintervs", chrom1=2, chrom2=2)
	gintervals.update("test.testintervs", r[c(2,3),], chrom=1)
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs2d")$size) - 100)
	gintervals.save("test.testintervs", "bigintervs2d")
	options(gmax.data.size = max.data.size)
	r <- gintervals.load("test.testintervs", chrom1=2, chrom2=2)
	gintervals.update("test.testintervs", r[c(2,3),], chrom1=2, chrom2=2)
	r <- list(gintervals.load("test.testintervs", chrom1=2, chrom2=2), gintervals.chrom_sizes("test.testintervs"))
	gintervals.rm("test.testintervs", force=T)
	r
}
{
	gintervals.rm("test.testintervs", force=T)
	max.data.size <- getOption("gmax.data.size")
	options(gmax.data.size = sum(gintervals.chrom_sizes("bigintervs2d")$size) - 100)
	gintervals.save("test.testintervs", "bigintervs2d")
	options(gmax.data.size = max.data.size)
	gintervals.update("test.testintervs", NULL, chrom1=2, chrom2=2)
	r <- list(gintervals.load("test.testintervs", chrom1=2, chrom2=2), gintervals.chrom_sizes("test.testintervs"))
	gintervals.rm("test.testintervs", force=T)
	r
}

{
    gintervals.rm("test.testintervs", force=T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size=6000000)
    try(giterator.intervals("test.fixedbin", gintervals(c(1,2)), intervals.set.out = "test.testintervs"), silent=T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force=T)
    r
}
{
    gintervals.rm("test.testintervs", force=T)
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size=300000)
    try(giterator.intervals("test.rects", gintervals.2d(c(1,2)), intervals.set.out = "test.testintervs"), silent=T)
    options(gmax.data.size = max.data.size)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force=T)
    r
}