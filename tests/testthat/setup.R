library(dplyr, warn.conflicts = FALSE)
library(glue, warn.conflicts = FALSE)

if (getOption("gmulticontig.indexed_format", FALSE)) {
    gsetroot("/net/mraid20/ifs/wisdom/tanay_lab/tgdata/db/tgdb/misha_test_db_indexed/")
} else {
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/")
}
gdb.reload()
gdir.create("temp", showWarnings = FALSE)
options(gmax.data.size = 1e9)
