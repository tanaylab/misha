library(dplyr, warn.conflicts = FALSE)
library(glue, warn.conflicts = FALSE)

# Only set up the main testdb if it exists
testdb_path <- "/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"
if (dir.exists(testdb_path)) {
    gsetroot(testdb_path)
    gdb.reload()
    options(gmax.data.size = 1e9)
}
