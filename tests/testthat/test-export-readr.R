test_that("gintervals.export_bed uses readr when available and option enabled", {
    tmp <- withr::local_tempfile(fileext = ".bed")
    intervals <- data.frame(chrom = "chr1", start = 0L, end = 10L)
    used_readr <- FALSE

    withr::local_options(misha.use_readr = TRUE)
    with_mocked_bindings(
        .misha_readr_available = function() TRUE,
        .misha_write_with_readr = function(bed_data, file, append, col_names, na) {
            used_readr <<- TRUE
            lines <- apply(bed_data, 1, paste, collapse = "\t")
            con <- file(file, if (isTRUE(append)) "a" else "w")
            on.exit(close(con), add = TRUE)
            writeLines(lines, con, sep = "\n", useBytes = TRUE)
        },
        .env = asNamespace("misha"),
        gintervals.export_bed(intervals, tmp)
    )

    expect_true(used_readr)
    expect_identical(readLines(tmp), "chr1\t0\t10")
})

test_that("gintervals.export_bed falls back to write.table when readr is unavailable", {
    tmp <- withr::local_tempfile(fileext = ".bed")
    intervals <- data.frame(chrom = "chr1", start = 0L, end = 10L)
    withr::local_options(misha.use_readr = FALSE)

    gintervals.export_bed(intervals, tmp)

    expect_identical(readLines(tmp), "chr1\t0\t10")
})

test_that("gtrack.export_bed uses readr path and NA handling when available", {
    tmp <- withr::local_tempfile(fileext = ".bed")
    intervals <- gintervals(c(1, 1), c(0, 10), c(10, 20))
    values <- c(NA_real_, 5)
    track <- paste0("track_readr_", sample(1:1e9, 1))
    gtrack.create_sparse(track, "tmp", intervals, values)
    withr::defer(gtrack.rm(track, force = TRUE))

    used_readr <- FALSE
    na_arg <- NULL

    withr::local_options(misha.use_readr = TRUE)
    with_mocked_bindings(
        .misha_readr_available = function() TRUE,
        .misha_write_with_readr = function(bed_data, file, append, col_names, na) {
            used_readr <<- TRUE
            na_arg <<- na
            bed_out <- bed_data
            bed_out[is.na(bed_out)] <- na
            lines <- vapply(seq_len(nrow(bed_out)), function(i) paste(bed_out[i, ], collapse = "\t"), character(1))
            con <- file(file, if (isTRUE(append)) "a" else "w")
            on.exit(close(con), add = TRUE)
            writeLines(lines, con, sep = "\n", useBytes = TRUE)
        },
        .env = asNamespace("misha"),
        gtrack.export_bed(track, tmp, intervals = intervals)
    )

    expect_true(used_readr)
    expect_identical(na_arg, ".")
    bed_lines <- readLines(tmp)
    expect_length(bed_lines, 2)
    expect_match(bed_lines[1], "^(chr)?1\\t0\\t10\\tinterval_1\\t\\.$")
    expect_match(bed_lines[2], "^(chr)?1\\t10\\t20\\tinterval_2\\t5$")
})
