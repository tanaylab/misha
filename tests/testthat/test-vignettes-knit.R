# The vignettes are executable documentation: every chunk that can run against
# the bundled examples database does run. This test is what keeps them that
# way. It knits each vignette source in a temp directory and fails on the first
# chunk that errors.
#
# knitr::knit() (not rmarkdown::render()) is used deliberately: it needs no
# pandoc, and it runs in this process, so the vignettes are evaluated against
# the misha that is under test rather than whatever copy happens to be
# installed. knitr's own default is error = TRUE (it records the error in the
# output and carries on), which would make this test pass on a vignette that
# hard-errors, so error = FALSE is forced below; the two Manual chunks that
# deliberately show an error set error = TRUE per chunk and still win.
#
# skip_on_cran: `R CMD build` already runs every vignette, so making CRAN pay
# for it twice buys nothing.

restore_groot_on_exit()

vignette_sources <- function() {
    # Source checkout / R CMD check on the unpacked tarball
    candidates <- c(
        file.path(testthat::test_path(), "..", "..", "vignettes"),
        system.file("doc", package = "misha")
    )
    for (dir in candidates) {
        if (!nzchar(dir) || !dir.exists(dir)) next
        files <- list.files(dir, pattern = "\\.Rmd$", full.names = TRUE)
        if (length(files)) {
            return(normalizePath(files))
        }
    }
    character(0)
}

test_that("every vignette knits without error", {
    skip_on_cran()
    skip_if_not_installed("knitr")

    vignettes <- vignette_sources()
    skip_if(!length(vignettes), "vignette sources are not available here")

    old_error <- knitr::opts_chunk$get("error")
    knitr::opts_chunk$set(error = FALSE)
    withr::defer(knitr::opts_chunk$set(error = old_error))

    for (vignette in vignettes) {
        wd <- withr::local_tempdir()
        stopifnot(file.copy(vignette, wd))
        err <- withr::with_dir(wd, {
            tryCatch(
                {
                    knitr::knit(basename(vignette), output = tempfile(fileext = ".md"), quiet = TRUE)
                    NULL
                },
                error = function(e) conditionMessage(e)
            )
        })
        expect_true(
            is.null(err),
            info = sprintf("%s failed to knit: %s", basename(vignette), err)
        )
    }
})

test_that("the runnable vignettes really do evaluate their chunks", {
    # Guards the other direction: a regression that re-disables evaluation
    # would leave the knit test above passing on an empty document.
    vignettes <- vignette_sources()
    skip_if(!length(vignettes), "vignette sources are not available here")

    runnable <- c("Misha-Basics.Rmd", "Manual.Rmd", "Database-Formats.Rmd")
    for (name in runnable) {
        f <- vignettes[basename(vignettes) == name]
        skip_if(!length(f), paste(name, "not found"))
        headers <- grep("^```\\{r", readLines(f[1]), value = TRUE)
        headers <- headers[!grepl("include\\s*=\\s*FALSE", headers)]
        disabled <- grepl("eval\\s*=\\s*FALSE", headers)
        expect_gt(sum(!disabled), 0)
        # Database-Formats keeps four chunks off (FASTA input, a multi-GB
        # download, and two destructive shell commands); the other two must be
        # fully live.
        limit <- if (identical(name, "Database-Formats.Rmd")) 4L else 0L
        expect_lte(sum(disabled), limit)
    }
})
