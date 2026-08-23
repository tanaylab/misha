library(dplyr, warn.conflicts = FALSE)
library(glue, warn.conflicts = FALSE)

# Every test process starts rooted at its own overlay of the shared test
# database, never at the shared database itself. Files that create tracks
# without asking for a database of their own therefore write into TMPDIR, so
# two suites can run at the same time without seeing each other's tracks.
# See build_test_db_overlay() in helper-test_db.R.
testdb_path <- test_db_root()

if (!is.null(testdb_path)) {
    gsetroot(testdb_path)
    gdb.reload()
    gdir.create("temp", showWarnings = FALSE)
    options(gmax.data.size = 1e9)
}
