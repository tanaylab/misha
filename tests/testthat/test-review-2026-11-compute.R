# Regression tests for the compute and virtual-track findings of the 2026-11 review of 5.14.4.

create_isolated_test_db()

test_that("gcor is accurate at a large offset, single process and multitask", {
    iv <- gintervals(c(1, 2), 0, 2e5)
    x <- gextract("test.fixedbin", iv, iterator = "test.fixedbin")$test.fixedbin
    x <- x[!is.na(x)] + 1e8
    y <- 2 * x - 1e8
    for (mt in c(FALSE, TRUE)) {
        withr::local_options(list(gmultitasking = mt, gmin.scope4process = 1000, gmax.processes = 8))
        r <- gcor("test.fixedbin + 1e8", "test.fixedbin * 2 + 1e8", intervals = iv, iterator = "test.fixedbin", details = TRUE)
        expect_equal(r[1, "sd1"], sd(x), tolerance = 1e-6)
        expect_equal(r[1, "sd2"], sd(y), tolerance = 1e-6)
        expect_equal(r[1, "cov"], cov(x, y), tolerance = 1e-6)
        expect_equal(r[1, "cor"], 1, tolerance = 1e-9)
    }
})

test_that("summaries give a NaN standard deviation when a value is infinite, as sd() does", {
    expr <- "ifelse(test.fixedbin > 0.15, Inf, test.fixedbin)"
    iv <- gintervals(1, 0, 1000)
    x <- gextract(expr, iv, iterator = 100, colnames = "v")$v
    expect_true(any(is.infinite(x)))
    expect_true(is.nan(sd(x)))

    expect_true(is.nan(gsummary(expr, iv, iterator = 100)[["Std dev"]]))
    expect_true(is.nan(gintervals.summary(expr, iv, iterator = 100)[["Std dev"]]))
    expect_true(is.nan(gbins.summary("test.fixedbin", breaks = c(0, 0.5), expr = expr, intervals = iv, iterator = 100)[1, "Std dev"]))
    # a finite set still has a finite one
    expect_equal(gsummary("test.fixedbin", iv, iterator = 100)[["Std dev"]], sd(gextract("test.fixedbin", iv, iterator = 100)$test.fixedbin))
})

test_that("2D gintervals.annotate with tie_method and maxneighbors > 1 gives a clear error", {
    q <- gintervals.2d(1, 100, 200, 1, 300, 400)
    a <- gintervals.2d(1, c(500, 1000), c(600, 1100), 1, c(700, 1200), c(800, 1300))
    a$lab <- c("x", "y")
    expect_error(
        gintervals.annotate(q, a, annotation_columns = "lab", maxneighbors = 2, tie_method = "min.start"),
        "1D intervals only"
    )
    expect_equal(gintervals.annotate(q, a, annotation_columns = "lab", maxneighbors = 2)$lab, c("x", "y"))
})

test_that("first, last, sample, lse and their positions work on array tracks", {
    iv <- gintervals("chr1", c(100, 200, 300, 500, 700), c(150, 260, 400, 600, 750))
    a <- c(1, 2, 3, 4, 5)
    b <- c(11, 12, 13, 14, 15)
    vals <- (a + b) / 2
    f <- tempfile(fileext = ".tsv")
    on.exit(unlink(f), add = TRUE)
    write.table(cbind(iv, a = a, b = b), f, sep = "\t", quote = FALSE, row.names = FALSE)
    nm <- paste0("test.arr_fl_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(nm, force = TRUE))
    suppressMessages(gtrack.array.import(nm, "", f))
    withr::defer(remove_all_vtracks())

    q <- gintervals("chr1", c(0, 120, 230, 410), c(1000, 330, 240, 490))
    ext <- function(func) {
        gvtrack.create("arr_v", nm, func)
        gextract("arr_v", q, iterator = q, colnames = "v")$v
    }
    expect_equal(ext("first"), c(vals[1], vals[1], vals[2], NaN))
    expect_equal(ext("first.pos.abs"), c(100, 120, 230, NaN))
    expect_equal(ext("last"), c(vals[5], vals[3], vals[2], NaN))
    expect_equal(ext("last.pos.abs"), c(700, 300, 230, NaN))
    expect_equal(ext("last.pos.relative"), c(700, 180, 0, NaN))
    expect_equal(ext("lse"), c(log(sum(exp(vals))), log(sum(exp(vals[1:3]))), vals[2], NaN), tolerance = 1e-6)
    s <- ext("sample")
    expect_true(s[1] %in% vals && s[2] %in% vals[1:3])
    expect_equal(s[3:4], c(vals[2], NaN))
    sp <- ext("sample.pos.abs")
    expect_true(sp[1] %in% iv$start)
    expect_true(sp[2] %in% c(120, 200, 300))
    expect_equal(sp[3:4], c(230, NaN))
})

test_that("the minimal rectangle of a band intersection is exact", {
    sq <- gintervals.2d("chr1", 0, 100, "chr1", 0, 100)
    for (band in list(c(-1000, -50), c(-1000, 0), c(10, 1000), c(-20, 31))) {
        px <- expand.grid(x = 0:99, y = 0:99)
        px <- px[px$x - px$y >= band[1] & px$x - px$y < band[2], ]
        r <- gintervals.2d.band_intersect(sq, band = band)
        expect_equal(c(r$start1, r$end1, r$start2, r$end2), c(min(px$x), max(px$x) + 1, min(px$y), max(px$y) + 1), info = paste(band, collapse = ","))
    }
})

test_that("bands wider than 2^31 are accepted", {
    sq <- gintervals.2d("chr1", 0, 100, "chr1", 0, 100)
    r <- gintervals.2d.band_intersect(sq, band = c(10, 3e9))
    expect_equal(c(r$start1, r$end1, r$start2, r$end2), c(10, 100, 0, 90))
    r <- gintervals.2d.band_intersect(sq, band = c(-1e10, 1e10))
    expect_equal(c(r$start1, r$end1, r$start2, r$end2), c(0, 100, 0, 100))
})

test_that("gintervals.liftover returns an empty result when no chromosome is in the chain", {
    ch <- gintervals.as_chain(data.frame(
        chrom = "chr1", start = 1000, end = 2000, strand = 1, chromsrc = "chrA", startsrc = 0, endsrc = 1000,
        strandsrc = 1, chain_id = 1, score = 10
    ))
    expect_null(gintervals.liftover(data.frame(chrom = "chrC", start = 5, end = 60), ch))
    expect_null(gintervals.liftover(data.frame(chrom = "chrC", start = 5, end = 60, v = 1), ch, value_col = "v"))
    # as when the chromosome is in the chain but nothing overlaps it
    expect_null(gintervals.liftover(data.frame(chrom = "chrA", start = 5000, end = 6000), ch))
    r <- gintervals.liftover(data.frame(chrom = c("chrC", "chrA"), start = c(5, 500), end = c(60, 600)), ch)
    expect_equal(c(r$start, r$end, r$intervalID), c(1500, 1600, 2))
})

test_that("gintervals.2d.intersect gives every pairwise intersection over many rectangles", {
    set.seed(11)
    mk <- function(n) {
        s1 <- sample(0:9000, n, TRUE)
        s2 <- sample(0:9000, n, TRUE)
        gintervals.2d("chr1", s1, s1 + sample(1:300, n, TRUE), "chr1", s2, s2 + sample(1:300, n, TRUE))
    }
    a <- mk(1500)
    b <- mk(1000)
    # pairs with both axes overlapping, the a row varying fastest
    hit <- outer(a$start1, b$end1, "<") & outer(a$end1, b$start1, ">") & outer(a$start2, b$end2, "<") & outer(a$end2, b$start2, ">")
    idx <- which(hit, arr.ind = TRUE)
    exp_df <- data.frame(
        start1 = pmax(a$start1[idx[, 1]], b$start1[idx[, 2]]), end1 = pmin(a$end1[idx[, 1]], b$end1[idx[, 2]]),
        start2 = pmax(a$start2[idx[, 1]], b$start2[idx[, 2]]), end2 = pmin(a$end2[idx[, 1]], b$end2[idx[, 2]])
    )
    r <- gintervals.2d.intersect(a, b)
    key <- function(d) sort(paste(d$start1, d$end1, d$start2, d$end2))
    expect_equal(nrow(r), nrow(exp_df))
    expect_equal(key(r), key(exp_df))
})

test_that("sequence and value-based vtracks under a 2D iterator need a dimension projection", {
    withr::defer(remove_all_vtracks())
    pssm <- matrix(c(0.7, 0.1, 0.1, 0.1), 4, 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
    it2d <- gintervals.2d("chr1", 100, 200, "chr2", 10, 60)
    gvtrack.create("pwm2d", NULL, "pwm.max", pssm = pssm)
    expect_error(gextract("pwm2d", it2d, iterator = it2d), "without explicit dimension projection")
    gvtrack.create("kmer2d", NULL, "kmer.count", kmer = "ACG")
    expect_error(gextract("kmer2d", it2d, iterator = it2d), "without explicit dimension projection")
    gvtrack.create("val2d", src = data.frame(chrom = "chr1", start = 0, end = 300, score = 5), func = "avg")
    expect_error(gextract("val2d", it2d, iterator = it2d), "without explicit dimension projection")

    gvtrack.iterator("val2d", dim = 1)
    expect_equal(gextract("val2d", it2d, iterator = it2d)$val2d, 5)
})

test_that("gseq.kmer.dist(mask =) leaves out k-mers that reach into the mask", {
    iv <- gintervals("chr1", 0, 2000)
    mask <- gintervals("chr1", c(150, 1000), c(160, 1001))
    k <- 3
    s <- strsplit(toupper(gseq.extract(iv)), "")[[1]]
    masked <- rep(FALSE, length(s))
    masked[151:160] <- TRUE
    masked[1001] <- TRUE
    starts <- seq_len(length(s) - k + 1)
    keep <- vapply(starts, function(p) !any(masked[p:(p + k - 1)]) && all(s[p:(p + k - 1)] %in% c("A", "C", "G", "T")), logical(1))
    kmers <- vapply(starts[keep], function(p) paste(s[p:(p + k - 1)], collapse = ""), character(1))
    expected <- table(kmers)

    r <- gseq.kmer.dist(iv, k = k, mask = mask)
    expect_equal(sum(r$count), length(kmers))
    expect_equal(r$count[match(names(expected), r$kmer)], as.numeric(expected))
})

test_that("pwm and pwm.max are NA at a chromosome end that holds no anchor", {
    withr::defer(remove_all_vtracks())
    pssm <- matrix(c(0.7, 0.1, 0.1, 0.1), 6, 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
    sz <- gintervals.all()
    sz <- sz$end[sz$chrom == "chr1"]
    iv <- gintervals("chr1", sz - 4, sz)
    for (fn in c("pwm", "pwm.max", "pwm.max.pos")) {
        gvtrack.create("pwm_end", NULL, fn, pssm = pssm)
        expect_true(is.na(gextract("pwm_end", iv, iterator = iv)$pwm_end), info = fn)
    }
    # documented: pwm.count reports 0 there
    gvtrack.create("pwm_end", NULL, "pwm.count", pssm = pssm, score.thresh = -100)
    expect_equal(gextract("pwm_end", iv, iterator = iv)$pwm_end, 0)
})

test_that("pwm.edit_distance handles motifs whose prefilter block would be too large", {
    local_db_state()
    withr::defer(remove_all_vtracks())
    set.seed(4)
    bases <- c("A", "C", "G", "T")
    motif <- sample(bases, 32, TRUE)
    seqv <- sample(bases, 1000, TRUE)
    # the motif with one substitution at 100 and with two at 500 (0-based starts)
    seqv[101:132] <- motif
    seqv[110] <- setdiff(bases, motif[10])[1]
    seqv[501:532] <- motif
    seqv[c(505, 520)] <- c(setdiff(bases, motif[5])[1], setdiff(bases, motif[20])[1])
    db <- make_seq_db(paste(seqv, collapse = ""))
    suppressMessages(gsetroot(db))

    pssm <- t(vapply(motif, function(b) ifelse(bases == b, 0.97, 0.01), numeric(4)))
    colnames(pssm) <- bases
    thr <- sum(log(apply(pssm, 1, max))) - 0.1
    # 32 columns and max_edits = 1 make two prefilter blocks of 16 columns
    gvtrack.create("ed32", NULL, "pwm.edit_distance", pssm = pssm, score.thresh = thr, max_edits = 1, bidirect = FALSE, prior = 0)
    q <- gintervals("chr1", c(100, 500), c(101, 501))
    expect_equal(gextract("ed32", q, iterator = q)$ed32, c(1, NaN))
})

test_that("neighbors among tied intervals do not depend on the input order of the ties", {
    set.seed(5)
    starts <- sample(seq(1000, 100000, by = 1000), 200, TRUE)
    i2 <- data.frame(chrom = "chr1", start = starts, end = starts + 10, id = seq_along(starts))
    i1 <- data.frame(chrom = "chr1", start = seq(1500, 100500, by = 1000), end = seq(1500, 100500, by = 1000) + 1)
    unsorted <- gintervals.neighbors(i1, i2, maxneighbors = 1, mindist = 0)
    presorted <- gintervals.neighbors(i1, i2[order(i2$start, i2$id), ], maxneighbors = 1, mindist = 0)
    expect_equal(unsorted$id, presorted$id)
})

test_that("gseq.pwm_edits reports tied edits in motif order", {
    pssm <- matrix(c(0.7, 0.1, 0.1, 0.1), 20, 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
    thr <- 20 * log(0.1) + 5 * (log(0.7) - log(0.1)) - 0.01
    r <- gseq.pwm_edits(strrep("C", 20), pssm, score.thresh = thr, bidirect = FALSE, prior = 0)
    expect_equal(r$motif_col, 1:5)
    r <- gseq.pwm_edits(strrep("C", 20), pssm, score.thresh = thr, bidirect = FALSE, prior = 0, max_indels = 1)
    expect_equal(r$motif_col, 1:5)
})
