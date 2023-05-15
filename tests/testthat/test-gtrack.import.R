{
    intervs <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2)))
    gtrack.rm("test.tmptrack", force = T)
    gtrack.import_mappedseq("test.tmptrack", "", "db/test/s_7_export.txt", remove.dups = F)
    r <- gextract("test.tmptrack", intervs)
    gtrack.rm("test.tmptrack", force = T)
    r
}
{
    intervs <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2)))
    gtrack.rm("test.tmptrack", force = T)
    gtrack.import_mappedseq("test.tmptrack", "", "db/test/sample-small.sam", cols.order = NULL, remove.dups = F)
    r <- gextract("test.tmptrack", intervs)
    gtrack.rm("test.tmptrack", force = T)
    r
}
{
    intervs <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2)))
    gtrack.rm("test.tmptrack", force = T)
    gtrack.import_mappedseq("test.tmptrack", "", "db/test/s_7_export.txt", remove.dups = F, pileup = 180, binsize = 50)
    r <- gextract("test.tmptrack", intervs)
    gtrack.rm("test.tmptrack", force = T)
    r
}

{
    max.data.size <- getOption("gmax.data.size")
    options(gmax.data.size = 10000)
    gtrack.rm("test.tmptrack", T)
    try(gtrack.2d.import("test.tmptrack", "aaa7", c("gtrack.2d.import/f4")), silent = T)
    options(gmax.data.size = max.data.size)
    r <- gextract("test.tmptrack", ALLGENOME)
    gtrack.rm("test.tmptrack", T)
    r
}