#' Tests if an object was changed since the last run.
#' If an rds file named \code{snapshot_dir/id.rds} exists its contents are compared with \{obj},
#' otherwise the file is created.
#'
#' @param obj an R object
#' @param id unique test id.
#' @param snapshot_dir directory with rds file containing snapshot of previous versions
expect_regression <- function(obj, id, snapshot_dir = "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot", tolerance = 1e-5, arrange_df = FALSE) {
    regression_file <- file.path(snapshot_dir, glue::glue("{id}.rds"))

    if (!file.exists(regression_file)) {
        readr::write_rds(obj, regression_file)
        system(glue::glue("chmod a-w {regression_file}"))
    }

    # We need testthat to always find the `expect` statement (otherwise - the test would be skipped)
    old <- readr::read_rds(regression_file)

    # Helper function to normalize chromosome factors for comparison
    # Converts factors to characters and sorts by genomic coordinates
    normalize_for_comparison <- function(df) {
        if (!is.data.frame(df)) {
            return(df)
        }

        # Handle 1D intervals
        if (all(c("chrom", "start", "end") %in% colnames(df))) {
            df$chrom <- as.character(df$chrom)
            df <- df[order(df$chrom, df$start, df$end), ]
            rownames(df) <- NULL
            return(df)
        }

        # Handle 2D intervals
        if (all(c("chrom1", "start1", "end1", "chrom2", "start2", "end2") %in% colnames(df))) {
            df$chrom1 <- as.character(df$chrom1)
            df$chrom2 <- as.character(df$chrom2)
            df <- df[order(df$chrom1, df$start1, df$end1, df$chrom2, df$start2, df$end2), ]
            rownames(df) <- NULL
            return(df)
        }

        return(df)
    }

    # Handle NULL comparisons - both should be NULL or both should be non-NULL
    if (is.null(obj) || is.null(old)) {
        expect_equal(old, obj, tolerance = tolerance)
        return(invisible())
    }

    expect_equal(old, obj, tolerance = tolerance)
}

load_regression_file <- function(id, snapshot_dir = "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot") {
    file_path <- file.path(snapshot_dir, glue::glue("{id}.rds"))
    if (!file.exists(file_path)) {
        stop(glue::glue("Regression file {id} not found"))
    }
    readr::read_rds(file_path)
}
