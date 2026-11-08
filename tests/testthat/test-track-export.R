# Re-roots into a database of its own; hand this process's overlay back to the next
# file in the parallel worker (helper-test_db.R).
restore_groot_on_exit()

gdb.init_examples()

test_that("bedGraph export writes valid format for dense track", {
    outfile <- tempfile(fileext = ".bedgraph")
    withr::defer(unlink(outfile))

    gtrack.export_bedgraph("dense_track", outfile)

    expect_true(file.exists(outfile))

    lines <- readLines(outfile)
    expect_true(length(lines) > 1)

    # Check header line
    expect_match(lines[1], "^track type=bedGraph name=")

    # Check data lines have 4 tab-separated columns
    data_lines <- lines[-1]
    splits <- strsplit(data_lines[1:min(10, length(data_lines))], "\t")
    for (s in splits) {
        expect_equal(length(s), 4)
        # chrom should be a valid chromosome
        expect_match(s[1], "^chr")
        # start and end should be numeric
        expect_false(is.na(as.numeric(s[2])))
        expect_false(is.na(as.numeric(s[3])))
        # value should be numeric
        expect_false(is.na(as.numeric(s[4])))
    }
})

test_that("bedGraph export works with specific intervals", {
    outfile <- tempfile(fileext = ".bedgraph")
    withr::defer(unlink(outfile))

    intervs <- gintervals(1, 0, 1000)
    gtrack.export_bedgraph("dense_track", outfile, intervals = intervs)

    lines <- readLines(outfile)
    data_lines <- lines[-1]

    # All data lines should be from chr1
    chroms <- vapply(strsplit(data_lines, "\t"), `[`, character(1), 1)
    expect_true(all(chroms == "chr1"))

    # All start positions should be within the specified range
    starts <- as.numeric(vapply(strsplit(data_lines, "\t"), `[`, character(1), 2))
    expect_true(all(starts >= 0))

    ends <- as.numeric(vapply(strsplit(data_lines, "\t"), `[`, character(1), 3))
    expect_true(all(ends <= 1000))
})

test_that("bedGraph export works with iterator parameter", {
    outfile <- tempfile(fileext = ".bedgraph")
    withr::defer(unlink(outfile))

    intervs <- gintervals(1, 0, 1000)
    gtrack.export_bedgraph("dense_track", outfile,
        intervals = intervs,
        iterator = 200
    )

    lines <- readLines(outfile)
    data_lines <- lines[-1]

    # With iterator=200 and range 0-1000, we expect 5 bins
    expect_equal(length(data_lines), 5)

    # Check bin boundaries
    starts <- as.numeric(vapply(strsplit(data_lines, "\t"), `[`, character(1), 2))
    ends <- as.numeric(vapply(strsplit(data_lines, "\t"), `[`, character(1), 3))
    expect_equal(starts, c(0, 200, 400, 600, 800))
    expect_equal(ends, c(200, 400, 600, 800, 1000))
})

test_that("bedGraph export supports gzip compression", {
    outfile <- tempfile(fileext = ".bedgraph.gz")
    withr::defer(unlink(outfile))

    intervs <- gintervals(1, 0, 1000)
    gtrack.export_bedgraph("dense_track", outfile, intervals = intervs)

    expect_true(file.exists(outfile))

    # Read back from gzip
    lines <- readLines(gzfile(outfile))
    expect_true(length(lines) > 1)
    expect_match(lines[1], "^track type=bedGraph name=")
})

test_that("bedGraph export excludes NaN values", {
    outfile <- tempfile(fileext = ".bedgraph")
    withr::defer(unlink(outfile))

    # Sparse tracks have NaN for uncovered regions when using a fixed iterator
    gtrack.export_bedgraph("sparse_track", outfile,
        intervals = gintervals(1, 0, 10000),
        iterator = 100
    )

    lines <- readLines(outfile)
    data_lines <- lines[-1]

    # No data line should contain NaN
    for (line in data_lines) {
        parts <- strsplit(line, "\t")[[1]]
        val <- as.numeric(parts[4])
        expect_false(is.nan(val))
    }
})

test_that("bedGraph export uses custom name in header", {
    outfile <- tempfile(fileext = ".bedgraph")
    withr::defer(unlink(outfile))

    gtrack.export_bedgraph("dense_track", outfile,
        intervals = gintervals(1, 0, 500),
        name = "my_custom_name"
    )

    lines <- readLines(outfile)
    expect_match(lines[1], 'name="my_custom_name"')
})

test_that("bedGraph export errors on 2D track", {
    expect_error(
        gtrack.export_bedgraph("rects_track", tempfile()),
        "2D tracks are not supported"
    )
})

test_that("bedGraph export works with track expression", {
    outfile1 <- tempfile(fileext = ".bedgraph")
    outfile2 <- tempfile(fileext = ".bedgraph")
    withr::defer({
        unlink(outfile1)
        unlink(outfile2)
    })

    intervs <- gintervals(1, 0, 500)
    gtrack.export_bedgraph("dense_track", outfile1,
        intervals = intervs, iterator = 100
    )
    gtrack.export_bedgraph("dense_track * 2", outfile2,
        intervals = intervs, iterator = 100
    )

    lines1 <- readLines(outfile1)[-1]
    lines2 <- readLines(outfile2)[-1]

    vals1 <- as.numeric(vapply(strsplit(lines1, "\t"), `[`, character(1), 4))
    vals2 <- as.numeric(vapply(strsplit(lines2, "\t"), `[`, character(1), 4))

    expect_equal(vals2, vals1 * 2, tolerance = 1e-10)
})

test_that("bedGraph export sorts by chromosome order then start", {
    outfile <- tempfile(fileext = ".bedgraph")
    withr::defer(unlink(outfile))

    gtrack.export_bedgraph("dense_track", outfile, iterator = 10000)

    lines <- readLines(outfile)
    data_lines <- lines[-1]

    chroms <- vapply(strsplit(data_lines, "\t"), `[`, character(1), 1)
    starts <- as.numeric(vapply(strsplit(data_lines, "\t"), `[`, character(1), 2))

    # Get expected chromosome order from the genome
    genome_chroms <- gintervals.all()$chrom
    chrom_idx <- match(chroms, genome_chroms)

    # Should be sorted by chrom index then start
    expect_true(all(diff(chrom_idx) >= 0))
    for (ch in unique(chroms)) {
        ch_starts <- starts[chroms == ch]
        expect_true(all(diff(ch_starts) > 0))
    }
})

# The data lines gtrack.export_bedgraph wrote when R built them: NaN rows dropped, rows in genome order
r_bedgraph_lines <- function(track, intervals, iterator) {
    data <- gextract(track, intervals = intervals, iterator = iterator)
    data <- data[!is.nan(data[[track]]), , drop = FALSE]
    data <- data[order(match(data$chrom, gintervals.all()$chrom), data$start), , drop = FALSE]
    paste(data$chrom, sprintf("%.0f", data$start), sprintf("%.0f", data$end), data[[track]], sep = "\t")
}

test_that("bedGraph export writes the lines R's paste() writes, whatever scipen and OutDec", {
    outfile <- tempfile(fileext = ".bedgraph")
    withr::defer(unlink(outfile))
    intervs <- gintervals(c(1, 2), 0, 20000)
    # values of many magnitudes; a few exact zeros and round numbers
    exprs <- c("dense_track", "dense_track * 1e5 / 3", "dense_track / 7e8", "round(dense_track * 1000) * 100", "-dense_track * 1e300")
    for (expr in exprs) {
        for (opts in list(list(scipen = 0, OutDec = "."), list(scipen = 3, OutDec = "."), list(scipen = -4, OutDec = ","))) {
            withr::with_options(opts, {
                gtrack.export_bedgraph(expr, outfile, intervals = intervs, iterator = 50)
                expect_identical(readLines(outfile)[-1], r_bedgraph_lines(expr, intervs, 50),
                    label = paste(expr, opts$scipen, opts$OutDec)
                )
            })
        }
    }
})

test_that("C++ formats doubles as as.character() does", {
    set.seed(17)
    n <- 20000
    vals <- c(
        runif(n, -1, 1) * 10^runif(n, -310, 308),
        round(runif(n, -1e6, 1e6), sample(0:8, n, replace = TRUE)),
        sample(-1e5:1e5, n, replace = TRUE) * 10^sample(-20:20, n, replace = TRUE),
        10^(-20:20) - 10^(-20:20) * 1e-15, 10^(-20:20) + 10^(-20:20) * 1e-15,
        readBin(writeBin(runif(n, -100, 100), raw(), size = 4), "numeric", size = 4, n = n),
        0.5 + 1:1000, 1 / (1:1000), 2^(-1074:-1000), 2^(1000:1023), NA, Inf, -Inf, 0, -0
    )
    nv <- length(vals)
    for (opts in list(list(scipen = 0, OutDec = "."), list(scipen = 7, OutDec = "."), list(scipen = -2, OutDec = ","))) {
        withr::with_options(opts, {
            fmt <- .gbedgraph_format_options()
            got <- .gcall(
                "gtrack_bedgraph_lines", rep("c", nv), rep(0, nv), rep(1, nv), vals, NULL,
                1, 1e9, fmt$scipen, fmt$dec, .misha_env()
            )[[1]]
            expect_identical(strsplit(got, "\n")[[1]], paste("c", "0", "1", vals, sep = "\t"),
                label = paste(opts$scipen, opts$OutDec)
            )
        })
    }
})

test_that("C++ formats rows in threads as in one thread, in the order given", {
    withr::local_options(gmax.processes = 4)
    set.seed(3)
    n <- 300000
    vals <- c(runif(n / 2) * 1e3, rep(NaN, 1000), round(runif(n / 2 - 1000), 3))
    chrom <- factor(rep(c("chr1", "chr2"), each = n / 2))
    start <- seq_len(n) * 10
    ord <- sample(n)
    fmt <- .gbedgraph_format_options()
    got <- character()
    from <- 1
    # 250,000 rows in threads, then 50,000 in one
    while (from <= n) {
        res <- .gcall("gtrack_bedgraph_lines", chrom, start, start + 10, vals, ord, from, 250000, fmt$scipen, fmt$dec, .misha_env())
        got <- c(got, res[[1]])
        from <- res[[2]]
    }
    expect_length(got, 2)
    keep <- ord[!is.nan(vals[ord])]
    expected <- paste0(paste(chrom[keep], sprintf("%.0f", start[keep]), sprintf("%.0f", start[keep] + 10), vals[keep], sep = "\t"), "\n", collapse = "")
    expect_identical(paste(got, collapse = ""), expected)
})

test_that("bedGraph export over intervals with no values writes the header only", {
    outfile <- tempfile(fileext = ".bedgraph")
    withr::defer(unlink(outfile))
    # sparse_track has no values there (its first two intervals are chr1:0-50 and chr1:100-150)
    empty <- gintervals(1, 60, 90)
    expect_null(gextract("sparse_track", intervals = empty))
    expect_warning(gtrack.export_bedgraph("sparse_track", outfile, intervals = empty), "no values in the intervals")
    expect_identical(readLines(outfile), 'track type=bedGraph name="sparse_track"')
})

test_that("bedGraph export errors on missing arguments", {
    expect_error(gtrack.export_bedgraph(), "Usage")
})

test_that("bigwig export errors when converter not available", {
    # Skip if bedGraphToBigWig is available (we test success separately)
    converter <- Sys.which("bedGraphToBigWig")
    bundled <- system.file("bedGraphToBigWig", package = "misha")
    wig_converter <- Sys.which("wigToBigWig")
    bundled_wig <- system.file("wigToBigWig", package = "misha")

    has_converter <- nzchar(converter) ||
        (nzchar(bundled) && file.exists(bundled)) ||
        nzchar(wig_converter) ||
        (nzchar(bundled_wig) && file.exists(bundled_wig))

    skip_if(has_converter, "bedGraphToBigWig or wigToBigWig is available")

    expect_error(
        gtrack.export_bigwig("dense_track", tempfile(fileext = ".bw")),
        "bedGraphToBigWig or wigToBigWig not found"
    )
})

test_that("bigwig export works when converter is available", {
    converter <- Sys.which("bedGraphToBigWig")
    bundled <- system.file("bedGraphToBigWig", package = "misha")
    wig_converter <- Sys.which("wigToBigWig")
    bundled_wig <- system.file("wigToBigWig", package = "misha")

    has_converter <- nzchar(converter) ||
        (nzchar(bundled) && file.exists(bundled)) ||
        nzchar(wig_converter) ||
        (nzchar(bundled_wig) && file.exists(bundled_wig))

    skip_if_not(has_converter, "bedGraphToBigWig or wigToBigWig not available")

    outfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(outfile))

    gtrack.export_bigwig("dense_track", outfile,
        intervals = gintervals(1, 0, 10000)
    )

    expect_true(file.exists(outfile))
    expect_true(file.info(outfile)$size > 0)
})

test_that("bigwig export works when the genome order is not the C order of the chromosome names", {
    has_converter <- nzchar(Sys.which("bedGraphToBigWig")) || file.exists(system.file("bedGraphToBigWig", package = "misha"))
    skip_if_not(has_converter, "bedGraphToBigWig not available")
    local_db_state()
    withr::with_tempdir({
        sizes <- c("1" = 1000, "2" = 1000, "10" = 1000)
        create_test_db("num_db", chrom_sizes = data.frame(chrom = names(sizes), size = unname(sizes)))
        writeLines(paste0(names(sizes), "\t", sizes), file.path("num_db", "chrom_sizes.txt"))
        gdb.init("num_db")
        skip_if_not(identical(as.character(gintervals.all()$chrom), names(sizes)), "the database sorts its chromosomes")
        gtrack.create_dense("t", "", gintervals.all(), c(1, 2, 3), 100)
        # bedGraphToBigWig refuses "10" after "2"
        f <- file.path(getwd(), "t.bw")
        gtrack.export_bigwig("t", f)
        gtrack.import("back", "", f, binsize = 100)
        expect_identical(gextract("back", gintervals.all(), iterator = 100)$back, gextract("t", gintervals.all(), iterator = 100)$t)
    })
})

test_that("bigwig export errors on 2D track", {
    expect_error(
        gtrack.export_bigwig("rects_track", tempfile(fileext = ".bw")),
        "2D tracks are not supported"
    )
})

test_that("bigwig export errors on missing arguments", {
    expect_error(gtrack.export_bigwig(), "Usage")
})
