gdb.get_readonly_attrs()
{
    old_ro_attrs <- gdb.get_readonly_attrs()
    r <- try(gdb.set_readonly_attrs(c(old_ro_attrs, "")), silent = T)
    gdb.set_readonly_attrs(old_ro_attrs)
    r
}
{
    old_ro_attrs <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(NULL)
    r <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(old_ro_attrs)
    r
}
{
    old_ro_attrs <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(c(old_ro_attrs, "testattr1", "testattr2"))
    r <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(old_ro_attrs)
    r
}

gtrack.attr.get("test.fixedbin", "created.by")
gtrack.attr.get("test.fixedbin", "blablabla")
{
    attrs <- gdb.get_readonly_attrs()
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)
    gtrack.attr.set("test.fixedbin", "testattr1", "value")
    gtrack.attr.get("test.fixedbin", "testattr1")
}
{
    attrs <- gdb.get_readonly_attrs()
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)
    gtrack.attr.set("test.fixedbin", "testattr1", "value")
    gtrack.attr.set("test.fixedbin", "testattr1", "")
    gtrack.attr.export("test.fixedbin")
}
gtrack.attr.export()
gtrack.attr.export("blablablatrack")
gtrack.attr.export(attrs = c("created.by"))
gtrack.attr.export(attrs = c("created.by", "created.date"))
gtrack.attr.export(c("test.fixedbin", "test.sparse"))
gtrack.attr.export(c("test.fixedbin", "test.sparse"), attrs = c("created.by", "created.date"))
{
    attrs <- gdb.get_readonly_attrs()
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)
    r <- gtrack.attr.export()
    r$testattr1 <- 1:dim(r)[1]
    gtrack.attr.import(r)
    gtrack.attr.export()
}
{
    attrs <- gdb.get_readonly_attrs()
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)
    r <- gtrack.attr.export()
    r$testattr1 <- 1:dim(r)[1]
    gtrack.attr.import(r)
    r$testattr1 <- NULL
    gtrack.attr.import(r)
    gtrack.attr.export()
}
{
    attrs <- gdb.get_readonly_attrs()
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)
    r <- gtrack.attr.export()
    r$testattr1 <- 1:dim(r)[1]
    gtrack.attr.import(r)
    r$testattr1 <- NULL
    gtrack.attr.import(r, replace = T)
    gtrack.attr.export()
}
{
    attrs <- gdb.get_readonly_attrs()
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)
    r <- gtrack.attr.export()
    r$testattr1 <- dim(r)[1]:1
    r <- r[c("test.fixedbin", "test.sparse"), ]
    gtrack.attr.import(r)
    gtrack.attr.export()
}
{
    attrs <- gdb.get_readonly_attrs()
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(c(attrs, "testattr1"))
    r <- gtrack.attr.export()
    r$testattr1 <- 10 + (dim(r)[1]:1)
    gtrack.attr.import(r)
}