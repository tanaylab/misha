# Regression tests for the track import and creation findings of the 2026-11 review of 5.14.4.

restore_groot_on_exit()

# A small database of its own (chr1, chr2, chrX, chrM), in the given format, made the root for the
# calling test. Removed, and the previous root restored, when the test ends.
import_test_db <- function(format = "indexed", env = parent.frame()) {
    local_db_state(env)
    dir <- tempfile("import_review_")
    dir.create(dir)
    withr::defer(unlink(dir, recursive = TRUE), envir = env)
    sizes <- c(chr1 = 10000, chr2 = 10000, chrX = 5000, chrM = 500)
    records <- paste0(">", names(sizes), "\n", strrep("A", sizes), "\n")
    # an indexed database from one multi-FASTA, a per-chromosome one from a FASTA per chromosome
    fa <- if (format == "indexed") file.path(dir, "genome.fa") else file.path(dir, paste0(names(sizes), ".fa"))
    if (format == "indexed") {
        cat(records, sep = "", file = fa)
    } else {
        for (i in seq_along(sizes)) cat(records[i], file = fa[i])
    }
    db <- file.path(dir, "db")
    suppressMessages(gdb.create(groot = db, fasta = fa, format = format, verbose = FALSE))
    suppressMessages(gdb.init(db))
    db
}

write_file <- function(text, ext) {
    path <- tempfile(fileext = ext)
    cat(text, file = path)
    path
}

for (format in c("indexed", "per-chromosome")) {
    test_that(sprintf("a tab-delimited file without a trailing newline keeps its last line (%s)", format), {
        import_test_db(format)

        f <- write_file("chrom\tstart\tend\tvalue\nchr1\t0\t10\t1\nchr1\t20\t30\t2\nchrX\t5\t9\t3", ".tsv")
        gtrack.import("nonl", "", f, 0)
        res <- gextract("nonl", gintervals.all())
        expect_equal(as.character(res$chrom), c("chr1", "chr1", "chrX"))
        expect_equal(res$start, c(0, 20, 5))
        expect_equal(res$nonl, c(1, 2, 3))

        f <- write_file("chrom\tstart\tend\ta\tb\nchr1\t0\t10\t1\t2\nchr1\t20\t30\t3\t4", ".tsv")
        gtrack.array.import("arr_nonl", "", f)
        res <- gtrack.array.extract("arr_nonl", NULL, gintervals.all())
        expect_equal(res$start, c(0, 20))
        expect_equal(res$a, c(1, 3))
        expect_equal(res$b, c(2, 4))
    })
}

test_that("gtrack.2d.import_contacts reads a fends file not sorted by fend id", {
    import_test_db("indexed")

    fends <- write_file("fend\tchr\tcoord\n10\tchr1\t1000\n5\tchr1\t2000\n7\tchr2\t300\n", ".tsv")
    contacts <- write_file("fend1\tfend2\tcount\n10\t5\t3\n10\t7\t4\n5\t7\t1\n", ".tsv")
    gtrack.2d.import_contacts("unsorted_fends", "", contacts, fends)
    res <- gextract("unsorted_fends", gintervals.2d.all())
    got <- data.frame(
        chrom1 = as.character(res$chrom1), start1 = res$start1, chrom2 = as.character(res$chrom2),
        start2 = res$start2, value = res$unsorted_fends
    )
    expected <- data.frame(
        chrom1 = c("chr1", "chr1", "chr1", "chr1", "chr2", "chr2"),
        start1 = c(1000, 1000, 2000, 2000, 300, 300),
        chrom2 = c("chr1", "chr2", "chr1", "chr2", "chr1", "chr1"),
        start2 = c(2000, 300, 1000, 300, 1000, 2000),
        value = c(3, 4, 3, 1, 4, 1)
    )
    got <- got[do.call(order, got), ]
    rownames(got) <- NULL
    expect_equal(got, expected)

    neg <- write_file("fend\tchr\tcoord\n-1\tchr1\t1000\n2\tchr1\t2000\n", ".tsv")
    expect_error(gtrack.2d.import_contacts("neg_fend", "", contacts, neg), "fend -1 is out of range")
    expect_false(gtrack.exists("neg_fend"))
})

test_that("gtrack.create_dense validates the intervals as gtrack.create_sparse does", {
    import_test_db("indexed")

    expect_error(
        gtrack.create_dense("neg", "", data.frame(chrom = "chr1", start = c(-10, 20), end = c(50, 60)), c(1, 3), 100),
        "start coordinate must be greater or equal than zero"
    )
    expect_error(
        gtrack.create_dense("rev", "", data.frame(chrom = "chr1", start = c(50, 20), end = c(10, 60)), c(1, 3), 100),
        "start coordinate must be lesser than end coordinate"
    )
    expect_error(
        gtrack.create_dense("na", "", data.frame(chrom = "chr1", start = c(NA, 20), end = c(10, 60)), c(1, 3), 100),
        "Invalid format of interval at index 1"
    )
    expect_error(
        gtrack.create_dense("past_end", "", data.frame(chrom = "chr1", start = 9990, end = 10010), 1, 100),
        "end coordinate exceeds chromosome boundaries"
    )
    expect_false(any(gtrack.exists(c("neg", "rev", "na", "past_end"))))

    # aliases still resolve
    gtrack.create_dense("alias", "", data.frame(chrom = c("1", "MT"), start = c(0, 0), end = c(50, 100)), c(2, 5), 100)
    expect_equal(gextract("alias", gintervals(c("chr1", "chrM"), 0, 100))$alias, c(2, 5))
})

test_that("gtrack.create_dense counts the bases under overlapping intervals once for defval", {
    import_test_db("indexed")

    # two intervals over the same 50 bases of a 100 bp bin, the other 50 uncovered
    ov <- data.frame(chrom = "chr1", start = c(0, 0), end = c(50, 50))
    gtrack.create_dense("ov_mean", "", ov, c(1, 3), 100, defval = 0)
    expect_equal(gextract("ov_mean", gintervals("chr1", 0, 100))$ov_mean, (50 * 1 + 50 * 3 + 50 * 0) / 150, tolerance = 1e-6)
    gtrack.create_dense("ov_min", "", ov, c(1, 3), 100, defval = -1, func = "min")
    expect_equal(gextract("ov_min", gintervals("chr1", 0, 100))$ov_min, -1)

    # the same through a BED file
    bed <- write_file("chr1\t0\t50\tp1\t1\nchr1\t0\t50\tp2\t3\n", ".bed")
    suppressMessages(gtrack.import("ov_bed", "", bed, 100, defval = 0))
    expect_equal(gextract("ov_bed", gintervals("chr1", 0, 100))$ov_bed, 4 / 3, tolerance = 1e-6)

    # random overlapping intervals against a per-base computation
    set.seed(17)
    n <- 60
    starts <- sample(0:9800, n, replace = TRUE)
    ends <- pmin(starts + sample(1:500, n, replace = TRUE), 10000)
    vals <- round(runif(n, 0, 10), 2)
    defval <- 0.5
    gtrack.create_dense("ov_rand", "", data.frame(chrom = "chr1", start = starts, end = ends), vals, 100, defval = defval)
    bin_starts <- seq(0, 9900, by = 100)
    expected <- vapply(bin_starts, function(bs) {
        overlap <- pmax(0, pmin(ends, bs + 100) - pmax(starts, bs))
        covered <- sum(vapply(bs:(bs + 99), function(x) any(starts <= x & ends > x), logical(1)))
        (sum(vals * overlap) + defval * (100 - covered)) / (sum(overlap) + 100 - covered)
    }, numeric(1))
    expect_equal(gextract("ov_rand", gintervals("chr1", 0, 10000))$ov_rand, expected, tolerance = 1e-5)
})

test_that("gtrack.create_dense does not slow down when a long interval spans many short ones", {
    local_db_state()
    dir <- tempfile("import_review_big_")
    dir.create(dir)
    withr::defer(unlink(dir, recursive = TRUE))
    fa <- file.path(dir, "chr1.fa")
    cat(">chr1\n", strrep("A", 200000), "\n", sep = "", file = fa)
    suppressMessages(gdb.create(groot = file.path(dir, "db"), fasta = fa, format = "indexed", verbose = FALSE))
    suppressMessages(gdb.init(file.path(dir, "db")))

    n <- 100000
    short <- data.frame(chrom = "chr1", start = seq(0, by = 2, length.out = n), end = seq(1, by = 2, length.out = n))
    with_long <- rbind(data.frame(chrom = "chr1", start = 0, end = 200000), short)
    elapsed <- function(intervals) {
        min(replicate(3, {
            t <- system.time(gtrack.create_dense("timed", "", intervals, rep(1, nrow(intervals)), 10))[["elapsed"]]
            gtrack.rm("timed", force = TRUE)
            t
        }))
    }
    elapsed(short)
    # one more interval should cost about as much as one more interval, not a rescan per bin
    # (about 50x on 5.15.1)
    expect_lt(elapsed(with_long), 10 * max(elapsed(short), 0.005))
})

test_that("gtrack.import_mappedseq resolves chromosome aliases", {
    import_test_db("indexed")

    sam <- write_file(paste0(
        "@HD\tVN:1.6\n",
        "r1\t0\t1\t101\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n",
        "r2\t0\tMT\t11\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n",
        "r3\t0\tchr1\t501\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n",
        "r4\t0\tchrBogus\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n"
    ), ".sam")
    stats <- gtrack.import_mappedseq("ms_alias", "", sam, cols.order = NULL)
    expect_equal(stats[[1]][["total.mapped"]], 3)
    expect_equal(stats[[1]][["total.unmapped"]], 1)
    res <- gextract("ms_alias", gintervals.all())
    expect_equal(as.character(res$chrom), c("chr1", "chr1", "chrM"))
    expect_equal(res$start, c(100, 500, 10))
    expect_equal(res$ms_alias, c(1, 1, 1))
})

test_that("a WIG / bedGraph block that resumes after an unknown chromosome is an error, as after a known one", {
    import_test_db("indexed")

    bg <- write_file("chr1\t0\t10\t1\nchrZZ\t0\t10\t9\nchr1\t20\t30\t3\n", ".bedgraph")
    expect_error(gtrack.import("resume_bg", "", bg, 0), "not sorted by chromosomes")
    expect_false(gtrack.exists("resume_bg"))

    wig <- write_file("variableStep chrom=chr1\n1\t1\nvariableStep chrom=chrZZ\n1\t9\nvariableStep chrom=chr1\n21\t3\n", ".wig")
    expect_error(gtrack.import("resume_wig", "", wig, 0), "not sorted by chromosomes")

    # the sibling case: a known chromosome in the middle
    bg2 <- write_file("chr1\t0\t10\t1\nchr2\t0\t10\t9\nchr1\t20\t30\t3\n", ".bedgraph")
    expect_error(gtrack.import("resume_bg2", "", bg2, 0), "not sorted by chromosomes")

    # an unknown block between two different chromosomes is still skipped
    bg3 <- write_file("chr1\t0\t10\t1\nchrZZ\t0\t10\t9\nchr2\t20\t30\t3\n", ".bedgraph")
    expect_message(gtrack.import("skip_bg", "", bg3, 0), "chrZZ")
    res <- gextract("skip_bg", gintervals.all())
    expect_equal(as.character(res$chrom), c("chr1", "chr2"))
    expect_equal(res$skip_bg, c(1, 3))
})

test_that("bedGraph and WIG import skip UCSC browser lines", {
    import_test_db("indexed")

    bg <- write_file("browser position chr1:1-100\nbrowser hide all\ntrack type=bedGraph\nchr1\t0\t10\t1\nchr1\t20\t30\t2\n", ".bedgraph")
    gtrack.import("browser_bg", "", bg, 0)
    res <- gextract("browser_bg", gintervals.all())
    expect_equal(res$start, c(0, 20))
    expect_equal(res$browser_bg, c(1, 2))

    wig <- write_file("browser position chr1:1-100\ntrack type=wiggle_0\nfixedStep chrom=chr1 start=1 step=10 span=10\n4\n5\n", ".wig")
    gtrack.import("browser_wig", "", wig, 0)
    res <- gextract("browser_wig", gintervals.all())
    expect_equal(res$start, c(0, 10))
    expect_equal(res$browser_wig, c(4, 5))
})

test_that("BED import skips and reports chromosomes missing from the database", {
    import_test_db("indexed")

    bed <- write_file("chr1\t0\t10\nchrUn_zzz\t0\t10\nchrX\t5\t9\n", ".bed")
    expect_message(gtrack.import("bed_sparse", "", bed, 0), "skipped: chrUn_zzz")
    res <- gextract("bed_sparse", gintervals.all())
    expect_equal(as.character(res$chrom), c("chr1", "chrX"))
    expect_equal(res$start, c(0, 5))
    expect_equal(res$bed_sparse, c(1, 1))

    expect_message(gtrack.import("bed_dense", "", bed, 10), "skipped: chrUn_zzz")
    res <- gextract("bed_dense", gintervals(c("chr1", "chrX"), 0, 10))
    expect_equal(res$bed_dense, c(1, 1))

    # a primary-looking name gets the warning of the other formats
    bed7 <- write_file("chr1\t0\t10\nchr7\t0\t10\n", ".bed")
    expect_warning(suppressMessages(gtrack.import("bed_chr7", "", bed7, 0)), "primary chromosome")

    # nothing matched: an error, and no track
    none <- write_file("chrUn_a\t0\t10\nchrUn_b\t0\t10\n", ".bed")
    expect_error(suppressMessages(gtrack.import("bed_none", "", none, 0)), "None of the chromosome names")
    expect_false(gtrack.exists("bed_none"))
})

for (format in c("indexed", "per-chromosome")) {
    test_that(sprintf("2D imports skip and report unknown chromosomes, and fail when nothing matches (%s)", format), {
        import_test_db(format)
        header <- "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tvalue\n"

        f <- write_file(paste0(header, "1\t0\t10\tMT\t0\t10\t1\nchrBogus\t0\t10\tchr1\t0\t10\t2\nchr2\t0\t10\tX\t0\t10\t3\n"), ".tsv")
        expect_message(gtrack.2d.import("rects_unknown", "", f), "skipped: chrBogus")
        res <- gextract("rects_unknown", gintervals.2d.all())
        expect_equal(paste(res$chrom1, res$chrom2), c("chr1 chrM", "chr2 chrX"))
        expect_equal(res$rects_unknown, c(1, 3))

        # the names of every reading process are reported
        f2 <- write_file(paste0(header, "chr1\t20\t30\tchr2\t20\t30\t4\nchrUn_b\t0\t10\tchr1\t0\t10\t5\n"), ".tsv")
        withr::with_options(list(gmultitasking = TRUE, gmax.processes = 2), {
            expect_message(gtrack.2d.import("rects_two", "", c(f, f2)), "chrBogus, chrUn_b")
        })
        expect_equal(nrow(gextract("rects_two", gintervals.2d.all())), 3)

        none <- write_file(paste0(header, "chrBogus\t0\t10\tchrUn_q\t0\t10\t1\n"), ".tsv")
        expect_error(gtrack.2d.import("rects_none", "", none), "None of the chromosome names")
        expect_false(gtrack.exists("rects_none"))

        empty <- write_file(header, ".tsv")
        expect_error(gtrack.2d.import("rects_empty", "", empty), "No intervals were read")
        expect_false(gtrack.exists("rects_empty"))

        expect_message(gtrack.2d.import_contacts("contacts_unknown", "", f), "skipped: chrBogus")
        expect_equal(nrow(gextract("contacts_unknown", gintervals.2d.all())), 4)
        expect_error(gtrack.2d.import_contacts("contacts_none", "", none), "None of the chromosome names")
        expect_error(gtrack.2d.import_contacts("contacts_empty", "", empty), "No contacts were read")
        expect_false(any(gtrack.exists(c("contacts_none", "contacts_empty"))))

        # an unknown chromosome in a fends file: its fends are skipped
        fends <- write_file("fend\tchr\tcoord\n1\tchr1\t1000\n2\tchrUn_x\t5\n3\tchr2\t300\n", ".tsv")
        contacts <- write_file("fend1\tfend2\tcount\n1\t3\t4\n1\t2\t7\n", ".tsv")
        expect_message(gtrack.2d.import_contacts("fends_unknown", "", contacts, fends), "skipped: chrUn_x")
        res <- gextract("fends_unknown", gintervals.2d.all())
        expect_equal(paste(res$chrom1, res$start1, res$chrom2, res$start2), c("chr1 1000 chr2 300", "chr2 300 chr1 1000"))
        expect_equal(res$fends_unknown, c(4, 4))
    })
}

test_that("gtrack.array.import skips and reports unknown chromosomes, and fails when nothing matches", {
    import_test_db("indexed")

    f <- write_file("chrom\tstart\tend\ta\nchrBogus\t0\t10\t1\nchr1\t0\t10\t2\n", ".tsv")
    expect_message(gtrack.array.import("arr_unknown", "", f), "skipped: chrBogus")
    res <- gtrack.array.extract("arr_unknown", NULL, gintervals.all())
    expect_equal(as.character(res$chrom), "chr1")
    expect_equal(res$a, 2)

    none <- write_file("chrom\tstart\tend\ta\nchrBogus\t0\t10\t1\n", ".tsv")
    expect_error(gtrack.array.import("arr_none", "", none), "None of the chromosome names")
    expect_false(gtrack.exists("arr_none"))
})

test_that("gtrack.import_mappedseq and gtrack.array.import write the indexed format on an indexed database", {
    import_test_db("indexed")

    sam <- write_file("r1\t0\tchr1\t101\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\nr2\t0\tchr2\t11\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n", ".sam")
    gtrack.import_mappedseq("ms_indexed", "", sam, cols.order = NULL)
    expect_equal(gtrack.info("ms_indexed")$format, "indexed")
    expect_setequal(list.files(file.path(.misha$GWD, "ms_indexed.track")), c("track.dat", "track.idx"))
    res <- gextract("ms_indexed", gintervals.all())
    expect_equal(paste(res$chrom, res$start), c("chr1 100", "chr2 10"))

    f <- write_file("chrom\tstart\tend\ta\tb\nchr1\t0\t10\t1\t2\nchrX\t20\t30\t3\t4\n", ".tsv")
    gtrack.array.import("arr_indexed", "", f)
    expect_equal(gtrack.info("arr_indexed")$format, "indexed")
    expect_false(any(c("chr1", "chr2", "chrX", "chrM") %in% list.files(file.path(.misha$GWD, "arr_indexed.track"))))
    res <- rbind(
        gtrack.array.extract("arr_indexed", NULL, gintervals("chr1")),
        gtrack.array.extract("arr_indexed", NULL, gintervals("chrX"))
    )
    expect_equal(paste(res$chrom, res$start), c("chr1 0", "chrX 20"))
    expect_equal(res$a, c(1, 3))
    expect_equal(res$b, c(2, 4))
})

test_that("gtrack.export_bedgraph writes each base once when the intervals overlap", {
    import_test_db("indexed")

    gtrack.create_dense("bins", "", data.frame(chrom = "chr1", start = seq(0, 90, 10), end = seq(10, 100, 10)), 1:10, 10)
    out <- tempfile(fileext = ".bedgraph")
    gtrack.export_bedgraph("bins", out, intervals = gintervals("chr1", c(0, 30), c(50, 80)))
    res <- read.table(out, skip = 1, col.names = c("chrom", "start", "end", "value"))
    expect_equal(res$start, seq(0, 70, 10))
    expect_equal(res$end, seq(10, 80, 10))
    expect_equal(res$value, 1:8)
})

test_that("gtrack.2d.import_contacts rejects fend ids too large to index and fend ids listed twice", {
    import_test_db("indexed")
    contacts <- write_file("fend1\tfend2\tcount\n1\t2\t3\n", ".tsv")

    for (id in c("9223372036854775807", "99999999999999999999")) {
        fends <- write_file(sprintf("fend\tchr\tcoord\n%s\tchr1\t100\n2\tchr1\t200\n", id), ".tsv")
        expect_error(gtrack.2d.import_contacts("big_fend", "", contacts, fends), "fend .* is out of range")
        expect_false(gtrack.exists("big_fend"))
    }

    dup <- write_file("fend\tchr\tcoord\n1\tchr1\t100\n2\tchr1\t200\n1\tchr2\t300\n", ".tsv")
    expect_error(gtrack.2d.import_contacts("dup_fend", "", contacts, dup), "line 4: fend 1 appears more than once")
    expect_false(gtrack.exists("dup_fend"))
})

for (binsize in c(0, 10)) {
    test_that(sprintf("BED import clips records past the chromosome end and drops empty ones (binsize %g)", binsize), {
        import_test_db("indexed")

        # chrM is 500 bp long
        bed <- write_file("chr1\t0\t10\nchr1\t100\t100\nchrM\t490\t520\nchrM\t600\t700\nchrM\t500\t510\n", ".bed")
        expect_warning(
            suppressMessages(gtrack.import("bed_fit", "", bed, binsize)),
            "1 record\\(s\\) running past the end of their chromosome were clipped to it, and 3 empty record\\(s\\)"
        )
        res <- gextract("bed_fit", gintervals(c("chr1", "chrM"), c(0, 490), c(10, 500)))
        expect_equal(res$bed_fit, c(1, 1))
        if (binsize == 0) {
            res <- gextract("bed_fit", gintervals.all())
            expect_equal(paste(res$chrom, res$start, res$end), c("chr1 0 10", "chrM 490 500"))
        }

        # a BED file that needs neither: no warning
        ok <- write_file("chr1\t0\t10\n", ".bed")
        expect_no_warning(suppressMessages(gtrack.import("bed_ok", "", ok, binsize)))

        # the creators called directly stay strict
        expect_error(
            gtrack.create_sparse("strict", "", data.frame(chrom = "chrM", start = 490, end = 520), 1),
            "end coordinate exceeds chromosome boundaries"
        )
    })
}

test_that("gtrack.export_bedgraph exports over a big intervals set", {
    import_test_db("indexed")

    gtrack.create_dense("bins_big", "", data.frame(chrom = "chr1", start = seq(0, 90, 10), end = seq(10, 100, 10)), 1:10, 10)
    iv <- data.frame(chrom = c("chr1", "chr1", "chr2"), start = c(0, 30, 0), end = c(50, 80, 20))
    withr::with_options(list(gmax.data.size = 2), gintervals.save("big_export_set", iv))
    expect_true(.gintervals.is_bigset("big_export_set"))

    out <- tempfile(fileext = ".bedgraph")
    gtrack.export_bedgraph("bins_big", out, intervals = "big_export_set")
    res <- read.table(out, skip = 1, col.names = c("chrom", "start", "end", "value"))
    expect_equal(res$chrom, rep("chr1", 8))
    expect_equal(res$start, seq(0, 70, 10))
    expect_equal(res$value, 1:8)
})

test_that("2D imports say so when chromosome names matched but no record has both chromosomes in the database", {
    import_test_db("indexed")

    f <- write_file("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tvalue\nchrBogus\t0\t10\tchr1\t0\t10\t1\n", ".tsv")
    expect_error(gtrack.2d.import("one_sided", "", f), "No interval of .* has both of its chromosomes in the genome database")
    expect_error(gtrack.2d.import_contacts("one_sided", "", f), "No contact of .* has both of its chromosomes in the genome database")
    expect_false(gtrack.exists("one_sided"))

    # fends on known chromosomes, but no contact between two listed fends
    fends <- write_file("fend\tchr\tcoord\n1\tchr1\t100\n2\tchrUn_x\t5\n", ".tsv")
    contacts <- write_file("fend1\tfend2\tcount\n1\t2\t3\n1\t7\t1\n", ".tsv")
    expect_error(suppressMessages(gtrack.2d.import_contacts("one_sided", "", contacts, fends)), "No contact of .* could be placed")
    expect_false(gtrack.exists("one_sided"))
})

test_that("a forced 2D re-conversion refuses when only some pairs have a per-pair file", {
    import_test_db("indexed")

    gtrack.2d.create("partial", "", gintervals.2d(c("chr1", "chr2"), 0, 100, c("chr1", "chrX"), 0, 100), c(1, 2))
    gtrack.2d.convert_to_indexed("partial", remove.old = FALSE, force = TRUE)
    dir <- file.path(.misha$GWD, "partial.track")
    pair_files <- setdiff(list.files(dir), c("track.dat", "track.idx", ".attributes"))
    expect_length(pair_files, 2)
    unlink(file.path(dir, pair_files[1]))
    dat_before <- readBin(file.path(dir, "track.dat"), "raw", file.size(file.path(dir, "track.dat")))

    expect_error(gtrack.2d.convert_to_indexed("partial", remove.old = TRUE, force = TRUE), "per-pair files for only some")
    expect_identical(readBin(file.path(dir, "track.dat"), "raw", file.size(file.path(dir, "track.dat"))), dat_before)
    expect_true(file.exists(file.path(dir, pair_files[2])))
    expect_equal(gextract("partial", gintervals.2d.all())$partial, c(1, 2))
})
