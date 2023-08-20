gtrack.array.extract("test.fixedbin", NULL, ALLGENOME)
gtrack.array.extract("test.array", NULL, gintervals(c(1, 2)))
gtrack.array.extract("test.array", c("col1", "col3", "col5"), gintervals(c(1, 2)))
{
    gtrack.array.extract("test.array", NULL, gintervals(c(1, 2)), file = "tmpresfile")
    r <- read.table("tmpresfile", sep = "	", nrows = 1000)
    unlink("tmpresfile")
    r
}
{
    gtrack.array.extract("test.array", c("col1", "col3", "col5"), gintervals(c(1, 2)), file = "tmpresfile")
    r <- read.table("tmpresfile", sep = "	", nrows = 1000)
    unlink("tmpresfile")
    r
}
{
    gintervals.rm("test.testintervs", force = T)
    intervs <- gscreen("test.fixedbin>0.2", gintervals(c(2, 4, 5, 10)))
    intervs <- intervs[sample(nrow(intervs)), ]
    try(gtrack.array.extract("test.array", c("col1", "col3", "col5"), intervals, intervals.set.out = "test.testintervs"), silent = T)
    r <- gintervals.load("test.testintervs")
    gintervals.rm("test.testintervs", force = T)
    r
}

gtrack.array.get_colnames("test.fixedbin")
gtrack.array.get_colnames("test.array")
gtrack.array.set_colnames("test.fixedbin", "col1")
gtrack.array.set_colnames("test.array", "col1")
{
    cols <- gtrack.array.get_colnames("test.array")
    gtrack.array.set_colnames("test.array", paste(cols, "blabla", sep = ""))
    r <- gtrack.array.get_colnames("test.array")
    gtrack.array.set_colnames("test.array", cols)
    r
}

{
    .ginteractive <- options(".ginteractive")[[1L]]
    options(.ginteractive = F)
    f1 <- tempfile()
    gextract("test.sparse", gintervals(c(1, 2)), file = f1)
    f2 <- tempfile()
    gtrack.array.extract("test.array", c("col2", "col3", "col4"), gintervals(c(1, 2)), file = f2)
    f3 <- tempfile()
    gtrack.array.extract("test.array", c("col1", "col3"), gintervals(c(1, 2)), file = f3)

    gtrack.array.import("test_track1", "", f1, f2)
    r1 <- gtrack.array.extract("test_track1", NULL, ALLGENOME)

    gtrack.array.import("test_track2", "", "test_track1", f3)
    r2 <- gtrack.array.extract("test_track2", NULL, ALLGENOME)
    r <- list()
    r$r1 <- r1
    r$r2 <- r2
    gtrack.rm("test_track1", TRUE)
    gtrack.rm("test_track2", TRUE)
    unlink(c(f1, f2, f3))
    options(.ginteractive = .ginteractive)
    r
}
