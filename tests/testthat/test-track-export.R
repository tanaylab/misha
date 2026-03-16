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

test_that("bigwig export errors on 2D track", {
    expect_error(
        gtrack.export_bigwig("rects_track", tempfile(fileext = ".bw")),
        "2D tracks are not supported"
    )
})

test_that("bigwig export errors on missing arguments", {
    expect_error(gtrack.export_bigwig(), "Usage")
})
