library(testthat)
library(misha)

if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    # The suite never runs on CRAN, and its multitasking tests start more than the
    # two processes that R CMD check --as-cran allows (parallel's .check_ncores()).
    Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
    test_check("misha")
}
