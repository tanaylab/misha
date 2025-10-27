#' Tests if an object was changed since the last run.
#' If an rds file named \code{snapshot_dir/id.rds} exists its contents are compared with \{obj},
#' otherwise the file is created.
#'
#' @param obj an R object
#' @param id unique test id.
#' @param snapshot_dir directory with rds file containing snapshot of previous versions
expect_regression <- function(obj, id, snapshot_dir = "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot", tolerance = 1e-5) {
    regression_file <- file.path(snapshot_dir, glue::glue("{id}.rds"))

    if (!file.exists(regression_file)) {
        readr::write_rds(obj, regression_file)
        system(glue::glue("chmod a-w {regression_file}"))
    }

    # We need testthat to always find the `expect` statement (otherwise - the test would be skipped)
    old <- readr::read_rds(regression_file)
    expect_equal(old, obj, tolerance = tolerance)
}

load_regression_file <- function(id, snapshot_dir = "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot") {
    file_path <- file.path(snapshot_dir, glue::glue("{id}.rds"))
    if (!file.exists(file_path)) {
        stop(glue::glue("Regression file {id} not found"))
    }
    readr::read_rds(file_path)
}

#' Save and restore the current database state
#'
#' This is useful for tests that temporarily switch to a different database.
#' Uses withr-style automatic cleanup.
#'
#' @param env The environment to use for defer (defaults to parent frame)
local_db_state <- function(env = parent.frame()) {
    # Save current state
    original_groot <- if (exists("GROOT", envir = .misha, inherits = FALSE)) {
        get("GROOT", envir = .misha)
    } else {
        NULL
    }

    # Register cleanup
    withr::defer(
        {
            if (!is.null(original_groot)) {
                suppressMessages(gdb.init(original_groot))
            }
        },
        envir = env
    )
}
