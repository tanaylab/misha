# Regression tests for the interval-function findings of the 2026-11 review of 5.14.4.

create_isolated_test_db()

test_that("use_intervals1_strand reads a factor strand by its levels", {
    tss <- data.frame(chrom = "chr1", start = c(1000, 5000), end = c(1001, 5001), strand = factor(c("+", "-")))
    feat <- data.frame(chrom = "chr1", start = c(800, 4800), end = c(810, 4810))
    tss_chr <- transform(tss, strand = as.character(strand))

    expect_equal(gintervals.neighbors(tss, feat, use_intervals1_strand = TRUE)$dist, c(-190, 190))
    expect_equal(
        gintervals.neighbors(tss, feat, use_intervals1_strand = TRUE)$dist,
        gintervals.neighbors(tss_chr, feat, use_intervals1_strand = TRUE)$dist
    )
    up <- gintervals.neighbors.upstream(tss, feat, maxdist = 1000)
    expect_equal(nrow(up), 1)
    expect_equal(up$start, 1000)
    # a factor whose level codes differ from the strand order
    tss_rev <- transform(tss, strand = factor(c("+", "-"), levels = c("-", "+")))
    expect_equal(gintervals.neighbors(tss_rev, feat, use_intervals1_strand = TRUE)$dist, c(-190, 190))
})

test_that("gintervals.annotate takes the neighbor's columns when names collide with the query's", {
    q <- data.frame(chrom = "chr1", start = c(100, 5000), end = c(110, 5010), gene = c("Q1", "Q2"))
    ann <- data.frame(chrom = "chr1", start = c(200, 5100), end = c(210, 5110), gene = c("A1", "A2"))

    res <- gintervals.annotate(q, ann, annotation_columns = "gene", column_names = "ann_gene")
    expect_equal(res$ann_gene, c("A1", "A2"))
    expect_equal(res$gene, c("Q1", "Q2"))
    expect_equal(res$dist, c(90, 90))

    res <- gintervals.annotate(q, ann, annotation_columns = "gene", overwrite = TRUE)
    expect_equal(res$gene, c("A1", "A2"))

    # the query's own "dist" column is not the distance
    q2 <- data.frame(chrom = "chr1", start = c(100, 5000), end = c(110, 5010), dist = c(-7, -8))
    ann2 <- data.frame(chrom = "chr1", start = c(200, 5100), end = c(210, 5110), label = c("A1", "A2"))
    res <- gintervals.annotate(q2, ann2, annotation_columns = "label", dist_column = "distance")
    expect_equal(res$distance, c(90, 90))
    expect_equal(res$dist, c(-7, -8))
    res <- gintervals.annotate(q2, ann2, annotation_columns = "label", overwrite = TRUE)
    expect_equal(res$dist, c(90, 90))
    expect_equal(sum(names(res) == "dist"), 1)

    # max_dist blanks the annotation, not the query's column of the same name
    ann3 <- data.frame(chrom = "chr1", start = c(200, 5015), end = c(210, 5025), gene = c("A1", "A2"))
    res <- gintervals.annotate(q, ann3, annotation_columns = "gene", overwrite = TRUE, max_dist = 50)
    expect_equal(res$gene, c(NA, "A2"))
    expect_equal(res$dist, c(90, 5))
    res <- gintervals.annotate(q, ann3, annotation_columns = "gene", column_names = "ann_gene", max_dist = 50)
    expect_equal(res$gene, c("Q1", "Q2"))
    expect_equal(res$ann_gene, c(NA, "A2"))
})

test_that("gintervals.annotate tie_method orders neighbors that share a distance", {
    q <- data.frame(chrom = "chr1", start = 1000, end = 1100)
    ann <- data.frame(chrom = "chr1", start = c(1200, 800), end = c(1300, 900), label = c("right", "left"))
    res <- gintervals.annotate(q, ann, annotation_columns = "label", maxneighbors = 2, tie_method = "min.start")
    expect_equal(res$label, c("left", "right"))

    # also without keep_order (it errored), and no internal column leaks out
    res <- gintervals.annotate(q, ann, annotation_columns = "label", maxneighbors = 2, tie_method = "min.start", keep_order = FALSE)
    expect_equal(res$label, c("left", "right"))
    expect_false(".orig_order" %in% names(res))
})

test_that("giterator.intervals(interval_relative = TRUE) intervalID is the input row", {
    iv <- data.frame(chrom = c("chr2", "chr1", "chr1"), start = c(100, 500, 100), end = c(200, 600, 180))
    res <- giterator.intervals(NULL, iv, iterator = 50, interval_relative = TRUE)
    src <- iv[res$intervalID, ]
    expect_equal(as.character(res$chrom), src$chrom)
    expect_true(all(res$start >= src$start & res$end <= src$end))
    expect_equal(res$intervalID[as.character(res$chrom) == "chr2"], c(1L, 1L))
})

test_that("gintervals.as_chain: strand 0 is forward, and min_score is applied", {
    ch <- gintervals.as_chain(data.frame(
        chrom = "chr1", start = 1000, end = 2000, strand = 0,
        chromsrc = "chr1", startsrc = 5000, endsrc = 6000, strandsrc = 0, chain_id = 1L, score = 1000
    ))
    res <- gintervals.liftover(gintervals("chr1", 5100, 5200), ch)
    expect_equal(c(res$start, res$end), c(1100, 1200))

    two <- data.frame(
        chrom = "chr1", start = c(1000, 3000), end = c(2000, 4000), strand = 1,
        chromsrc = "chr2", startsrc = 100, endsrc = 1100, strandsrc = 1, chain_id = c(1L, 2L), score = c(10, 1000)
    )
    ch2 <- gintervals.as_chain(two, src_overlap_policy = "keep", tgt_overlap_policy = "keep", min_score = 100)
    res <- gintervals.liftover(gintervals("chr2", 200, 300), ch2, include_metadata = TRUE)
    expect_equal(res$chain_id, 2L)
    expect_equal(c(res$start, res$end), c(3100, 3200))

    # a chain without attributes takes min_score from the liftover call
    res <- gintervals.liftover(gintervals("chr2", 200, 300), two,
        src_overlap_policy = "keep", tgt_overlap_policy = "keep", min_score = 100, include_metadata = TRUE
    )
    expect_equal(res$chain_id, 2L)
})

test_that("an intervs_mat subset by row names keeps its intervals", {
    mat <- gintervals.to_mat(data.frame(
        chrom = c("chr1", "chr1", "chr2"), start = c(100L, 500L, 200L), end = c(200L, 700L, 400L),
        t1 = c(1.5, 2.5, 3.5), t2 = c(10, 20, 30)
    ))
    res <- gintervals.from_mat(mat[c("chr2:200-400", "chr1:100-200"), ])
    expect_equal(as.character(res$chrom), c("chr2", "chr1"))
    expect_equal(res$start, c(200, 100))
    expect_equal(res$t1, c(3.5, 1.5))
    res <- gintervals.from_mat(mat[c("chr2:200-400", "chr1:100-200"), "t2", drop = FALSE])
    expect_equal(res$start, c(200, 100))
    expect_equal(gintervals.from_mat(mat[-1, ])$start, c(500, 200))
})

test_that("gintervals.force_range resolves chromosome aliases", {
    res <- gintervals.force_range(data.frame(chrom = c("1", "chrZZ", "chr2"), start = c(-10, 0, 10), end = c(60, 10, 20)))
    expect_equal(as.character(res$chrom), c("chr1", "chr2"))
    expect_equal(res$start, c(0, 10))
    expect_equal(res$end, c(60, 20))

    nm <- paste0("test.force_range_alias_", sample(1:1e9, 1))
    on.exit(gintervals.rm(nm, force = TRUE), add = TRUE)
    gintervals.force_range(data.frame(chrom = "1", start = -10, end = 60), intervals.set.out = nm)
    expect_equal(nrow(gintervals.load(nm)), 1)

    res <- gintervals.force_range(data.frame(chrom1 = "1", start1 = -5, end1 = 100, chrom2 = "chr2", start2 = 0, end2 = 50))
    expect_equal(as.character(res$chrom1), "chr1")
    expect_equal(res$start1, 0)
})

test_that("gintervals.neighbors with intervals.set.out saves unique column names", {
    nm <- paste0("test.nbo_", sample(1:1e9, 1))
    on.exit(gintervals.rm(nm, force = TRUE), add = TRUE)
    gintervals.neighbors(gintervals("chr1", 100, 200), gintervals("chr1", 500, 600), intervals.set.out = nm)
    saved <- gintervals.load(nm)
    expect_equal(names(saved), c("chrom", "start", "end", "chrom1", "start1", "end1", "dist"))
    expect_equal(saved$start1, 500)
})

test_that("gsynth.train uses the whole scope when its start is not bin-aligned", {
    st <- 1000050
    en <- 1008000
    k <- 2L
    seq <- toupper(gseq.extract(gintervals("chr1", st, en)))
    windows <- substring(seq, 1:(nchar(seq) - k), (k + 1):nchar(seq))
    expected <- 2 * sum(!grepl("[^ACGT]", windows))
    skip_if(expected == 0, "no sequence in the test region")

    for (scope_start in c(st, st + 49)) {
        model <- suppressMessages(gsynth.train(intervals = gintervals("chr1", scope_start, en), iterator = 100, k = k))
        exp_here <- if (scope_start == st) expected else 2 * sum(!grepl("[^ACGT]", windows[-(1:49)]))
        expect_equal(model$total_kmers, exp_here)
    }
})

test_that("gintervals.random(filter =) can start anywhere in a gap between filters", {
    chr1_end <- gintervals.all()$end[as.character(gintervals.all()$chrom) == "chr1"]
    # leaves start positions 2000..2050 for intervals of size 100
    filter <- gintervals("chr1", c(0, 2150), c(2000, chr1_end))
    set.seed(17)
    res <- gintervals.random(100, 3000, dist_from_edge = 0, chromosomes = "chr1", filter = filter)
    expect_equal(sort(unique(res$start)), 2000:2050)
})

test_that("summary Std dev is accurate for values with a large offset", {
    scope <- gintervals(c("chr1", "chr2"), 0, 200000)
    expr <- "test.fixedbin + 1e8"
    vals <- gextract(expr, scope, iterator = "test.fixedbin", colnames = "v")$v
    vals <- vals[!is.na(vals)]

    withr::local_options(gmultitasking = FALSE)
    single <- gsummary(expr, scope, iterator = "test.fixedbin")
    expect_equal(unname(single["Std dev"]), sd(vals), tolerance = 1e-9)
    expect_equal(unname(single["Sum"]), sum(vals), tolerance = 1e-12)

    # whole genome, so that the multitask run forks kids and merges their summaries
    all_vals <- gextract(expr, gintervals.all(), iterator = 1000, colnames = "v")$v
    all_vals <- all_vals[!is.na(all_vals)]
    single_all <- gsummary(expr, gintervals.all(), iterator = 1000)
    withr::with_options(list(gmultitasking = TRUE, gmax.processes = 4, gmin.scope4process = 1), {
        multi <- gsummary(expr, gintervals.all(), iterator = 1000)
    })
    expect_equal(unname(single_all["Std dev"]), sd(all_vals), tolerance = 1e-9)
    expect_equal(unname(multi["Std dev"]), sd(all_vals), tolerance = 1e-9)
    expect_equal(multi[c("Total intervals", "NaN intervals", "Min", "Max")], single_all[c("Total intervals", "NaN intervals", "Min", "Max")])
    expect_equal(unname(multi["Sum"]), unname(single_all["Sum"]), tolerance = 1e-12)

    per_interval <- gintervals.summary(expr, scope, iterator = "test.fixedbin")
    per_vals <- gextract(expr, scope, iterator = "test.fixedbin", colnames = "v")
    expected_sd <- tapply(per_vals$v, per_vals$intervalID, sd, na.rm = TRUE)
    expect_equal(per_interval[["Std dev"]], unname(as.numeric(expected_sd)), tolerance = 1e-9)
})

test_that("gbins.summary and gintervals.summary(intervals.set.out) keep double precision", {
    scope <- gintervals("chr1", c(0, 100000), c(100000, 200000))
    expr <- "test.fixedbin + 100000.123"
    vals <- gextract(expr, scope, iterator = "test.fixedbin", colnames = "v")
    rng <- range(vals$v, na.rm = TRUE)

    bins <- gbins.summary("test.fixedbin", c(rng[1] - 100001, rng[2]), expr = expr, intervals = scope, iterator = "test.fixedbin")
    expect_equal(as.numeric(bins[1, c("Min", "Max")]), rng, tolerance = 1e-14)
    expect_equal(as.numeric(bins[1, "Sum"]), sum(vals$v, na.rm = TRUE), tolerance = 1e-12)

    in_memory <- gintervals.summary(expr, scope, iterator = "test.fixedbin")
    expect_equal(in_memory$Min, as.numeric(tapply(vals$v, vals$intervalID, min, na.rm = TRUE)), tolerance = 1e-14)

    nm <- paste0("test.isum_", sample(1:1e9, 1))
    on.exit(gintervals.rm(nm, force = TRUE), add = TRUE)
    withr::with_options(list(gbig.intervals.size = 1), {
        gintervals.summary(expr, scope, iterator = "test.fixedbin", intervals.set.out = nm)
    })
    expect_true(.gintervals.is_bigset(nm))
    saved <- gintervals.load(nm)
    cols <- c("Min", "Max", "Sum", "Mean", "Std dev")
    expect_equal(as.matrix(saved[, cols]), as.matrix(in_memory[, cols]), tolerance = 1e-14)
})

test_that("liftover orientation combines the target and source strands", {
    q <- data.frame(chrom = "chrS", start = 900, end = 910)
    mk <- function(strand, strandsrc) {
        gintervals.as_chain(data.frame(
            chrom = "chr1", start = 5000, end = 5100, strand = strand,
            chromsrc = "chrS", startsrc = 900, endsrc = 1000, strandsrc = strandsrc, chain_id = 1L, score = 1
        ))
    }
    lifted <- function(chain) {
        r <- gintervals.liftover(q, chain)
        c(r$start, r$end)
    }
    expect_equal(lifted(mk(1, 1)), c(5000, 5010))
    expect_equal(lifted(mk(-1, -1)), c(5000, 5010))
    expect_equal(lifted(mk(1, -1)), c(5090, 5100))
    expect_equal(lifted(mk(-1, 1)), c(5090, 5100))

    # chain files: the same alignment with the reversed strand on either side
    chr1_size <- gintervals.all()$end[as.character(gintervals.all()$chrom) == "chr1"]
    chain_file <- function(tstrand, tstart, qstrand, qstart) {
        f <- tempfile(fileext = ".chain")
        writeLines(c(sprintf(
            "chain 1 chrS 1000 %s %d %d chr1 %d %s %d %d 1",
            tstrand, tstart, tstart + 100, chr1_size, qstrand, qstart, qstart + 100
        ), "100", ""), f)
        f
    }
    ref_minus <- chain_file("-", 0, "+", 5000) # source chrS:900-1000 reversed against chr1:5000-5100
    query_minus <- chain_file("+", 900, "-", chr1_size - 5100) # the same alignment, written the usual way
    both_minus <- chain_file("-", 0, "-", chr1_size - 5100) # same direction
    on.exit(unlink(c(ref_minus, query_minus, both_minus)), add = TRUE)
    expect_equal(lifted(query_minus), c(5090, 5100))
    expect_equal(lifted(ref_minus), c(5090, 5100))
    expect_equal(lifted(both_minus), c(5000, 5010))

    loaded <- gintervals.load_chain(ref_minus)
    expect_equal(loaded$strand, -1)
    expect_equal(loaded$strandsrc, 1)
    expect_equal(
        loaded[, c("start", "end", "startsrc", "endsrc", "strand", "strandsrc")],
        gintervals.load_chain(query_minus)[, c("start", "end", "startsrc", "endsrc", "strand", "strandsrc")]
    )
})

test_that("gtrack.liftover follows the combined strand orientation", {
    # source genome: one 1000 bp contig chrS with a sparse track at 900-910
    src_root <- tempfile("src_db")
    on.exit(unlink(src_root, recursive = TRUE), add = TRUE)
    create_test_db(src_root, data.frame(chrom = "chrS", size = 1000))
    tgt_root <- get("GROOT", envir = .misha)
    withr::defer(gsetroot(tgt_root))
    gsetroot(src_root)
    gtrack.create_sparse("src_sp", "", gintervals("chrS", 900, 910), 7)
    gsetroot(tgt_root)

    chr1_size <- gintervals.all()$end[as.character(gintervals.all()$chrom) == "chr1"]
    chain_file <- function(tstrand, tstart, qstrand, qstart) {
        f <- tempfile(fileext = ".chain")
        writeLines(c(sprintf(
            "chain 1 chrS 1000 %s %d %d chr1 %d %s %d %d 1",
            tstrand, tstart, tstart + 100, chr1_size, qstrand, qstart, qstart + 100
        ), "100", ""), f)
        f
    }
    for (case in list(
        list(f = chain_file("-", 0, "+", 5000), start = 5090),
        list(f = chain_file("+", 900, "-", chr1_size - 5100), start = 5090),
        list(f = chain_file("-", 0, "-", chr1_size - 5100), start = 5000)
    )) {
        nm <- paste0("test.lifted_", sample(1:1e9, 1))
        gtrack.liftover(nm, "", file.path(src_root, "tracks", "src_sp.track"), case$f)
        r <- gextract(nm, gintervals("chr1", 4000, 6000))
        gtrack.rm(nm, force = TRUE)
        unlink(case$f)
        expect_equal(c(r$start, r$end, r[[nm]]), c(case$start, case$start + 10, 7))
    }
})

test_that("array slice stddev is accurate at a large offset", {
    # offsets representable in float at 1e8 (spacing 8)
    vals <- 1e8 + c(0, 8, 16, 24, 40)
    f <- tempfile()
    on.exit(unlink(f), add = TRUE)
    write.table(data.frame(chrom = "chr1", start = 100, end = 200, a = vals[1], b = vals[2], c = vals[3], d = vals[4], e = vals[5]),
        f,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    nm <- paste0("test.arr_sd_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(nm, force = TRUE))
    suppressMessages(gtrack.array.import(nm, "", f))

    withr::defer(remove_all_vtracks())
    gvtrack.create("arr_sd_all", nm)
    gvtrack.array.slice("arr_sd_all", NULL, "stddev")
    gvtrack.create("arr_sd_slice", nm)
    gvtrack.array.slice("arr_sd_slice", c("a", "b", "c", "d", "e"), "stddev")
    r <- gextract(c("arr_sd_all", "arr_sd_slice"), gintervals("chr1", 100, 200))
    expect_equal(r$arr_sd_all, sd(vals), tolerance = 1e-6)
    expect_equal(r$arr_sd_slice, sd(vals), tolerance = 1e-6)
})

test_that("gtrack.liftover applies min_score of a pre-loaded chain", {
    src_root <- tempfile("src_db")
    on.exit(unlink(src_root, recursive = TRUE), add = TRUE)
    create_test_db(src_root, data.frame(chrom = "chrS", size = 1000))
    tgt_root <- get("GROOT", envir = .misha)
    withr::defer(gsetroot(tgt_root))
    gsetroot(src_root)
    gtrack.create_sparse("src_sp", "", gintervals("chrS", 900, 910), 7)
    gsetroot(tgt_root)
    src <- file.path(src_root, "tracks", "src_sp.track")

    two <- data.frame(
        chrom = "chr1", start = c(1000, 3000), end = c(1100, 3100), strand = 1,
        chromsrc = "chrS", startsrc = 900, endsrc = 1000, strandsrc = 1, chain_id = c(1L, 2L), score = c(10, 1000)
    )
    lifted_starts <- function(...) {
        nm <- paste0("test.lifted_", sample(1:1e9, 1))
        gtrack.liftover(nm, "", src, ...)
        r <- gextract(nm, gintervals("chr1", 0, 5000))
        gtrack.rm(nm, force = TRUE)
        r$start
    }
    ch <- gintervals.as_chain(two, src_overlap_policy = "keep", tgt_overlap_policy = "keep", min_score = 100)
    expect_equal(lifted_starts(ch), 3000)
    # a chain without attributes takes min_score from the call
    expect_equal(lifted_starts(two, src_overlap_policy = "keep", tgt_overlap_policy = "keep", min_score = 100), 3000)

    # a chain row with an NA score is dropped by min_score
    na_two <- two
    na_two$score <- c(NA, 1000)
    res <- gintervals.liftover(data.frame(chrom = "chrS", start = 900, end = 910), na_two,
        src_overlap_policy = "keep", tgt_overlap_policy = "keep", min_score = 100
    )
    expect_equal(res$start, 3000)
})
