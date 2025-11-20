library(dplyr, warn.conflicts = FALSE)
library(glue, warn.conflicts = FALSE)

if (getOption("gmulticontig.indexed_format", FALSE)) {
    testdb_path <- "/net/mraid20/ifs/wisdom/tanay_lab/tgdata/db/tgdb/misha_test_db_indexed/"
} else {
    testdb_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"
}

# Only set up the main testdb if it exists
if (dir.exists(testdb_path)) {
    gsetroot(testdb_path)
    gdb.reload()
    gdir.create("temp", showWarnings = FALSE)
    options(gmax.data.size = 1e9)
}
