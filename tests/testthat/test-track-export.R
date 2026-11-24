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

test_that("bedGraph lines whose threads cannot start are formatted in the calling thread", {
    skip_on_cran()
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux"), "uses prlimit and ps -L")
    skip_if(!nzchar(Sys.which("prlimit")) || !nzchar(Sys.which("ps")), "prlimit or ps is not available")
    # in a child R that may not start more threads (RLIMIT_NPROC counts the threads of the user, and
    # the limit is set to half of those running): a thread that cannot start used to abort R
    script <- tempfile(fileext = ".R")
    out <- tempfile()
    withr::defer(unlink(c(script, out)))
    writeLines(c(
        sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
        "suppressMessages(library(misha))",
        "suppressMessages(gdb.init_examples())",
        "options(gmax.processes = 8)",
        "n <- 300000",
        "chrom <- factor(rep(c('chr1', 'chr2'), each = n / 2))",
        "start <- seq_len(n) * 10",
        "vals <- seq_len(n) / 7",
        "expected <- paste0(paste(chrom, sprintf('%.0f', start), sprintf('%.0f', start + 10), vals, sep = '\\t'), '\\n', collapse = '')",
        "num <- as.integer(system('ps -L -u $(id -un) --no-headers | wc -l', intern = TRUE)) %/% 2",
        "system(sprintf('prlimit --pid %d --nproc=%d:%d', Sys.getpid(), num, num))",
        "res <- misha:::.gcall('gtrack_bedgraph_lines', chrom, start, start + 10, vals, NULL, 1, n, 0L, '.', misha:::.misha_env())",
        sprintf("writeLines(as.character(identical(res[[1]], expected)), %s)", deparse(out))
    ), script)
    # system2() quotes the command itself
    status <- suppressWarnings(system2(file.path(R.home("bin"), "R"), c("--no-echo", "--no-restore", "-f", shQuote(script)),
        stdout = FALSE, stderr = FALSE, env = "R_TESTS="
    ))
    expect_identical(status, 0L)
    expect_identical(readLines(out), "TRUE")
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

# x rounded to float32, as a bigWig stores its values
f32 <- function(x) readBin(writeBin(as.numeric(x), raw(), size = 4), "double", n = length(x), size = 4)

# The records gextract returns, without NaN, values rounded to float32: what a bigWig export of the
# same expression holds
expected_records <- function(expr, intervals, iterator = NULL) {
    d <- gextract(expr, intervals = intervals, iterator = iterator, colnames = "value")
    d <- d[!is.nan(d$value), c("chrom", "start", "end", "value")]
    d$value <- f32(d$value)
    d
}

# The parts of a bigWig file that misha's reader does not read: the header, the total summary,
# the chromosome tree and the records of each zoom level
read_bigwig_meta <- function(f) {
    raw <- readBin(f, "raw", file.size(f))
    uint <- function(off, n) sum(as.numeric(raw[off + seq_len(n)]) * 256^(seq_len(n) - 1))
    real <- function(off, n, size) readBin(raw[off + seq_len(n * size)], "double", n = n, size = size)

    chroms <- data.frame(name = character(0), id = numeric(0), size = numeric(0))
    tree <- uint(8, 8)
    key_size <- uint(tree + 8, 4)
    read_node <- function(off) {
        count <- uint(off + 2, 2)
        for (i in seq_len(count) - 1) {
            item <- off + 4 + i * (key_size + 8)
            if (as.integer(raw[off + 1]) == 1) {
                key <- raw[item + seq_len(key_size)]
                chroms[nrow(chroms) + 1, ] <<- list(rawToChar(key[key != 0]), uint(item + key_size, 4), uint(item + key_size + 4, 4))
            } else {
                read_node(uint(item + key_size, 8))
            }
        }
    }
    read_node(tree + 32)

    # the blocks an R-tree indexes: offset and size of each
    rtree_leaves <- function(off) {
        res <- list()
        walk <- function(node) {
            count <- uint(node + 2, 2)
            for (i in seq_len(count) - 1) {
                if (as.integer(raw[node + 1]) == 1) {
                    item <- node + 4 + i * 32
                    res[[length(res) + 1]] <<- c(uint(item + 16, 8), uint(item + 24, 8))
                } else {
                    walk(uint(node + 4 + i * 24 + 16, 8))
                }
            }
        }
        walk(off + 48)
        res
    }

    num_zooms <- uint(6, 2)
    zooms <- lapply(seq_len(num_zooms) - 1, function(z) {
        hdr <- 64 + 24 * z
        recs <- do.call(rbind, lapply(rtree_leaves(uint(hdr + 16, 8)), function(block) {
            data <- memDecompress(raw[block[1] + seq_len(block[2])], "gzip")
            n <- length(data) / 32
            ints <- matrix(as.numeric(readBin(data, "integer", n = 8 * n, size = 4)) %% 2^32, nrow = 8)
            floats <- matrix(readBin(data, "double", n = 8 * n, size = 4), nrow = 8)
            data.frame(
                chromid = ints[1, ], start = ints[2, ], end = ints[3, ], count = ints[4, ],
                min = floats[5, ], max = floats[6, ], sum = floats[7, ], sumsq = floats[8, ]
            )
        }))
        list(reduction = uint(hdr, 4), records = recs)
    })

    summary_offset <- uint(44, 8)
    list(
        magic = uint(0, 4), version = uint(4, 2), chroms = chroms, zooms = zooms,
        summary = c(count = uint(summary_offset, 8), real(summary_offset + 8, 4, 8))
    )
}

test_that("bigwig export round-trips a dense track", {
    f <- tempfile(fileext = ".bw")
    withr::defer(unlink(f))
    gtrack.export_bigwig("dense_track", f)

    expect_identical(readBin(f, "integer", n = 1, size = 4), -2003829722L) # 0x888FFC26, the bigWig magic
    exp <- expected_records("dense_track", gintervals.all())
    got <- gintervals.import_bigwig(f, gintervals.all())
    expect_equal(as.character(got$chrom), as.character(exp$chrom))
    expect_identical(got$start, exp$start)
    expect_identical(got$end, exp$end)
    expect_identical(got$value, exp$value)

    withr::defer(gtrack.rm("bw_back_dense", force = TRUE))
    gtrack.import("bw_back_dense", "", f, binsize = 50)
    back <- gextract("bw_back_dense", gintervals.all())$bw_back_dense
    orig <- gextract("dense_track", gintervals.all())$dense_track
    expect_identical(is.nan(back), is.nan(orig))
    expect_identical(back[!is.nan(back)], f32(orig[!is.nan(orig)]))
})

test_that("bigwig export round-trips a sparse track", {
    f <- tempfile(fileext = ".bw")
    withr::defer(unlink(f))
    gtrack.export_bigwig("sparse_track", f)
    exp <- expected_records("sparse_track", gintervals.all())
    got <- gintervals.import_bigwig(f, gintervals.all())
    expect_identical(nrow(got), nrow(exp))
    expect_identical(got$start, exp$start)
    expect_identical(got$end, exp$end)
    expect_identical(got$value, exp$value)
})

test_that("bigwig export of an expression over overlapping intervals writes each base once", {
    f <- tempfile(fileext = ".bw")
    withr::defer(unlink(f))
    ivs <- rbind(gintervals(1, 1000, 50000), gintervals(1, 30000, 80000), gintervals(2, 0, 120000), gintervals("X", 5, 99999))
    gtrack.export_bigwig("dense_track * 2 + 1", f, intervals = ivs, iterator = 100)
    exp <- expected_records("dense_track * 2 + 1", gintervals.canonic(ivs, unify_touching_intervals = FALSE), iterator = 100)
    got <- gintervals.import_bigwig(f, gintervals.all())
    expect_equal(as.character(got$chrom), as.character(exp$chrom))
    expect_identical(got$start, exp$start)
    expect_identical(got$end, exp$end)
    expect_identical(got$value, exp$value)
    # the partial bins at the ends of the intervals are kept
    expect_identical(got$start[got$chrom == "chrX"][1:2], c(5, 100))
})

test_that("bigwig export of an array track and of a virtual track writes gextract's values", {
    f <- tempfile(fileext = ".bw")
    withr::defer(unlink(f))
    gtrack.export_bigwig("array_track", f)
    exp <- expected_records("array_track", gintervals.all())
    got <- gintervals.import_bigwig(f, gintervals.all())
    expect_identical(got$start, exp$start)
    expect_identical(got$value, exp$value)

    withr::defer(gvtrack.rm("bw_vt"))
    gvtrack.create("bw_vt", "dense_track", func = "max")
    gtrack.export_bigwig("bw_vt", f, intervals = gintervals(c(1, 2), 0, 33333), iterator = 37)
    exp <- expected_records("bw_vt", gintervals(c(1, 2), 0, 33333), iterator = 37)
    got <- gintervals.import_bigwig(f, gintervals.all())
    expect_identical(got$start, exp$start)
    expect_identical(got$end, exp$end)
    expect_identical(got$value, exp$value)
})

test_that("bigwig export over a big intervals set reads it one chromosome at a time", {
    withr::local_options(gbig.intervals.size = 2)
    set <- paste0("test_bw_bigset_", Sys.getpid())
    withr::defer(gintervals.rm(set, force = TRUE))
    ivs <- rbind(gintervals(1, c(0, 500, 400), c(1000, 2000, 450)), gintervals(2, 100, 900), gintervals("X", 0, 300))
    gintervals.save(set, ivs)
    expect_true(misha:::.gintervals.is_bigset(set, FALSE))
    f <- tempfile(fileext = ".bw")
    withr::defer(unlink(f))
    gtrack.export_bigwig("dense_track", f, intervals = set)
    exp <- expected_records("dense_track", gintervals.canonic(ivs, unify_touching_intervals = FALSE))
    got <- gintervals.import_bigwig(f, gintervals.all())
    expect_identical(got$start, exp$start)
    expect_identical(got$end, exp$end)
    expect_identical(got$value, exp$value)
})

test_that("bigwig export writes the summary and zoom levels of the values", {
    f <- tempfile(fileext = ".bw")
    withr::defer(unlink(f))
    gtrack.export_bigwig("dense_track", f)
    meta <- read_bigwig_meta(f)
    exp <- expected_records("dense_track", gintervals.all())
    w <- exp$end - exp$start

    expect_identical(meta$version, 4)
    expect_identical(meta$chroms$name, c("chr1", "chr2", "chrX"))
    expect_identical(meta$chroms$id, c(0, 1, 2))
    expect_identical(meta$chroms$size, gintervals.all()$end)
    expect_equal(unname(meta$summary), c(sum(w), min(exp$value), max(exp$value), sum(exp$value * w), sum(exp$value^2 * w)))

    # 16 times the mean record width (50 bp bins), then 4 times larger at each level
    reductions <- vapply(meta$zooms, `[[`, numeric(1), "reduction")
    expect_identical(reductions[1], 800)
    expect_identical(reductions, 800 * 4^(seq_along(reductions) - 1))
    for (zoom in meta$zooms) {
        r <- zoom$records
        expect_identical(sum(r$count), sum(w))
        expect_equal(sum(r$sum), sum(exp$value * w), tolerance = 1e-6)
        expect_identical(min(r$min), min(exp$value))
        expect_identical(max(r$max), max(exp$value))
        # each record lies inside one bin of the level's size
        expect_true(all(r$start %/% zoom$reduction == (r$end - 1) %/% zoom$reduction))
    }
    # the first level against records summarized here
    z <- meta$zooms[[1]]$records
    chrom1 <- exp[exp$chrom == "chr1" & exp$start < 800, ]
    expect_identical(z$count[1], sum(chrom1$end - chrom1$start))
    expect_equal(z$sum[1], sum(chrom1$value * (chrom1$end - chrom1$start)), tolerance = 1e-6)
})

test_that("bigwig export writes the chromosomes in name order when the database's order differs", {
    local_db_state()
    withr::with_tempdir({
        create_db_with_unsorted_chrom_sizes("unsorted_db")
        gdb.init("unsorted_db")
        all <- gintervals.all()
        bins <- giterator.intervals(intervals = all, iterator = 10)
        set.seed(1)
        v <- round(rnorm(nrow(bins)), 3)
        v[sample(length(v), length(v) %/% 4)] <- NaN # NaN gaps
        v[bins$chrom == "chr10"] <- NaN # a chromosome without values
        gtrack.create_dense("t", "", bins, v, 10, defval = NaN)

        f <- file.path(getwd(), "t.bw")
        gtrack.export_bigwig("t", f)
        meta <- read_bigwig_meta(f)
        names <- as.character(all$chrom)
        expect_identical(meta$chroms$name, sort(names, method = "radix"))
        expect_identical(meta$chroms$id, as.numeric(seq_along(names) - 1))
        expect_identical(meta$chroms$size, all$end[match(meta$chroms$name, names)])

        exp <- expected_records("t", all)
        got <- gintervals.import_bigwig(f, all)
        expect_identical(nrow(got), nrow(exp))
        expect_identical(got$start, exp$start)
        expect_identical(got$value, exp$value)
        expect_false(any(got$chrom == "chr10"))

        gtrack.import("back", "", f, binsize = 10)
        back <- gextract("back", all, iterator = 10)$back
        orig <- gextract("t", all, iterator = 10)$t
        expect_identical(is.nan(back), is.nan(orig))
        expect_identical(back[!is.nan(back)], f32(orig[!is.nan(orig)]))
    })
})

test_that("bigwig export writes coordinates past 2^31", {
    local_db_state()
    withr::with_tempdir({
        dir.create("large_db/tracks", recursive = TRUE)
        dir.create("large_db/seq")
        sizes <- c(4.2e9, 1e6)
        for (i in 1:2) {
            # a sparse file: only its size matters
            con <- file(file.path("large_db", "seq", paste0("chr", i, ".seq")), "wb")
            seek(con, sizes[i] - 1, rw = "write")
            writeBin(charToRaw("A"), con)
            close(con)
        }
        writeLines(paste(1:2, sprintf("%.0f", sizes), sep = "\t"), file.path("large_db", "chrom_sizes.txt"))
        gdb.init("large_db")
        starts <- c(10, 2^31 - 5, 2^31 + 100, 3.5e9, 4.2e9 - 100)
        ivs <- gintervals(1, starts, starts + c(50, 10, 1000, 1e6, 100))
        gtrack.create_sparse("s", "", ivs, c(1.5, -2, 3e10, 0.25, 7))

        f <- file.path(getwd(), "s.bw")
        gtrack.export_bigwig("s", f)
        got <- gintervals.import_bigwig(f, gintervals.all())
        expect_identical(got$start, ivs$start)
        expect_identical(got$end, ivs$end)
        expect_identical(got$value, f32(c(1.5, -2, 3e10, 0.25, 7)))

        # bins across a batch boundary of the numeric iterator (2^22 bins of 1000 bp)
        scope <- gintervals(1, 4.19e9, 4.2e9)
        batches <- misha:::.gbigwig_batches(scope, 1000, gintervals.all(), c("chr1", "chr2"))[[1]]()
        expect_identical(batches[[1]]$end[1], 1000 * 2^22)
        gtrack.export_bigwig("s", f, intervals = scope, iterator = 1000)
        exp <- expected_records("s", scope, iterator = 1000)
        got <- gintervals.import_bigwig(f, gintervals.all())
        expect_identical(got$start, exp$start)
        expect_identical(got$end, exp$end)
        expect_identical(got$value, exp$value)
    })
})

test_that("bigwig export over intervals with no values writes a file without data", {
    f <- tempfile(fileext = ".bw")
    withr::defer(unlink(f))
    expect_warning(gtrack.export_bigwig("sparse_track", f, intervals = gintervals(1, 60, 90)), "no values in the intervals")
    expect_identical(nrow(gintervals.import_bigwig(f, gintervals.all())), 0L)
    meta <- read_bigwig_meta(f)
    expect_identical(meta$chroms$name, c("chr1", "chr2", "chrX"))
    expect_length(meta$zooms, 0)
    expect_identical(meta$summary[["count"]], 0)
})

test_that("bigwig export leaves no file when it fails", {
    f <- tempfile(fileext = ".bw")
    expect_error(gtrack.export_bigwig("dense_track + bw_no_such_variable", f))
    expect_false(file.exists(f))
})

test_that("a bigwig export that cannot open its file keeps the file already there", {
    skip_on_os("windows")
    f <- tempfile(fileext = ".bw")
    writeLines("keep me", f)
    Sys.chmod(f, "0444")
    withr::defer({
        Sys.chmod(f, "0644")
        unlink(f)
    })
    skip_if(file.access(f, 2) == 0, "file is writable (running as root)")
    expect_error(gtrack.export_bigwig("dense_track", f))
    expect_equal(readLines(f), "keep me")
})

test_that("the bigWig writer refuses unsorted and overlapping records", {
    f <- tempfile(fileext = ".bw")
    withr::defer(unlink(f))
    chroms <- c("chr1", "chr2")
    state <- misha:::.gcall("gtrack_bigwig_create", f, chroms, c(1000, 1000), misha:::.misha_env())
    expect_error(
        misha:::.gcall("gtrack_bigwig_append", f, state, c(0L, 0L), c(0, 50), c(100, 200), c(1, 2), misha:::.misha_env()),
        "must not overlap"
    )
    expect_error(
        misha:::.gcall("gtrack_bigwig_append", f, state, c(1L, 0L), c(0, 0), c(10, 10), c(1, 2), misha:::.misha_env()),
        "sorted by chromosome"
    )
    expect_error(
        misha:::.gcall("gtrack_bigwig_append", f, state, 0L, 0, 1001, 1, misha:::.misha_env()),
        "1000 bases long"
    )
    expect_error(misha:::.gcall("gtrack_bigwig_create", f, c("chr2", "chr1"), c(1000, 1000), misha:::.misha_env()), "sorted by name")
})

test_that("bigwig export works when the genome order is not the C order of the chromosome names", {
    local_db_state()
    withr::with_tempdir({
        sizes <- c("1" = 1000, "2" = 1000, "10" = 1000)
        create_test_db("num_db", chrom_sizes = data.frame(chrom = names(sizes), size = unname(sizes)))
        writeLines(paste0(names(sizes), "\t", sizes), file.path("num_db", "chrom_sizes.txt"))
        gdb.init("num_db")
        skip_if_not(identical(as.character(gintervals.all()$chrom), names(sizes)), "the database sorts its chromosomes")
        gtrack.create_dense("t", "", gintervals.all(), c(1, 2, 3), 100)
        f <- file.path(getwd(), "t.bw")
        gtrack.export_bigwig("t", f)
        expect_identical(read_bigwig_meta(f)$chroms$name, c("1", "10", "2"))
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
