{
    gdir.cd("test")
    r <- gtrack.ls()
    gdir.cd("..")
    r
}
{
    gdir.cd("test")
    r <- gdir.cwd()
    gdir.cd("..")
    r
}
{
    gdir.create("testdir")
    r1 <- dir(gdir.cwd())
    gdir.rm("testdir", force = T)
    r2 <- dir(gdir.cwd())
    list(r1, r2)
}



gtrack.exists("aaaaaaaaaaa.nnnnnnnnnnnnnn")
gtrack.exists("test.rects")



gtrack.ls()
gtrack.ls("tes")
gtrack.ls(blalaattr = "bubu")
gtrack.ls(created.by = "import")
gtrack.ls("wig", created.by = "import")
{
    gtrack.rm("test.tmptrack", force = T)
    gtrack.create("test.tmptrack", "", "test.fixedbin")
    intervs <- gscreen("test.fixedbin > 0.17 | is.na(test.fixedbin)", gintervals(c(1, 7)))
    gtrack.modify("test.tmptrack", "test.fixedbin + test.fixedbin", intervs)
    r <- gextract("test.tmptrack", gintervals(c(1, 2)))
    gtrack.rm("test.tmptrack", force = T)
    r
}
{
    gtrack.rm("test.tmptrack", force = T)
    gtrack.create("test.tmptrack", "", "test.sparse")
    gtrack.modify("test.tmptrack", "test.fixedbin + test.fixedbin", gintervals(1, 1000, 2000))
    gtrack.rm("test.tmptrack", force = T)
}
{
    gtrack.rm("test.tmptrack", force = T)
    gtrack.create("test.tmptrack", "", "test.rects")
    gtrack.modify("test.tmptrack", "test.fixedbin + test.fixedbin", gintervals.2d(1, 1000, 2000, 2, 3000, 4000))
    gtrack.rm("test.tmptrack", force = T)
}
{
    gtrack.rm("test.tmptrack", force = T)
    gtrack.create_sparse("test.tmptrack", "", gintervals(c(1, 2), 100, 2000), c(100, 200))
    r1 <- gtrack.ls()
    gtrack.rm("test.tmptrack", force = T)
    r2 <- gtrack.ls()
    list(r1, r2)
}
