# Window skipping in the thresholded and maximum PWM scans (src/PwmFastScan.h).
#
# The scanner only decides which windows are scored; the survivors are scored by
# the unchanged code. So every result must be IDENTICAL - not just equal to a
# tolerance - to scoring every window, which MISHA_DISABLE_PWM_SKIP=1 forces.
# The PSSMs below are built to stress the bound: zero and near-zero
# probabilities, widths 1 to 30, N runs, lowercase, thresholds at and around an
# exact window score, spatial weights and both strands.

with_pwm_skip <- function(skip, expr) {
    old <- Sys.getenv("MISHA_DISABLE_PWM_SKIP", unset = NA)
    Sys.setenv(MISHA_DISABLE_PWM_SKIP = if (skip) "0" else "1")
    on.exit(if (is.na(old)) Sys.unsetenv("MISHA_DISABLE_PWM_SKIP") else Sys.setenv(MISHA_DISABLE_PWM_SKIP = old))
    force(expr)
}

both_ways <- function(expr_fn) {
    list(skip = with_pwm_skip(TRUE, expr_fn()), noskip = with_pwm_skip(FALSE, expr_fn()))
}

fast_scan_pssm <- function(w, kind) {
    m <- switch(kind,
        sparse = matrix(stats::rgamma(4 * w, 0.2), ncol = 4),
        flat = matrix(stats::runif(4 * w, 0.2, 1), ncol = 4),
        tiny = {
            x <- matrix(stats::runif(4 * w), ncol = 4)
            x[cbind(seq_len(w), sample(4, w, TRUE))] <- 1e-12
            x
        },
        zero = {
            x <- matrix(stats::runif(4 * w), ncol = 4)
            x[cbind(seq_len(w), sample(4, w, TRUE))] <- 0
            x
        },
        onehot = {
            x <- matrix(0.001, w, 4)
            x[cbind(seq_len(w), sample(4, w, TRUE))] <- 1
            x
        }
    )
    m <- m + 1e-300 * (rowSums(m) == 0)
    colnames(m) <- c("A", "C", "G", "T")
    m
}

fast_scan_seqs <- function(n, len) {
    vapply(seq_len(n), function(i) {
        s <- sample(c("A", "C", "G", "T", "a", "c", "g", "t"), len, TRUE, prob = c(rep(3, 4), rep(1, 4)))
        if (i %% 3 == 0) {
            at <- sample(len, 1)
            s[at:min(len, at + sample(0:20, 1))] <- "N"
        }
        if (i %% 5 == 0) s[sample(len, 3)] <- "N"
        paste(s, collapse = "")
    }, "")
}

test_that("gseq.pwm count/max/pos are identical with and without window skipping", {
    withr::local_options(gmultitasking = FALSE)
    set.seed(20261124)
    kinds <- c("sparse", "flat", "tiny", "zero", "onehot")
    seqs <- fast_scan_seqs(12, 400)
    for (it in 1:60) {
        w <- if (it <= 30) it else sample(1:30, 1)
        kind <- kinds[(it %% length(kinds)) + 1]
        pssm <- fast_scan_pssm(w, kind)
        prior <- if (kind == "zero") 0 else sample(c(0, 0.01, 0.1), 1)
        bidirect <- it %% 4 != 0
        strand <- if (it %% 8 == 0) -1L else 1L
        spat <- if (it %% 7 == 0) list(spat.factor = stats::runif(5, 0.1, 3), spat.bin = 37L) else list()
        args <- c(list(seqs = seqs, pssm = pssm, bidirect = bidirect, strand = strand, prior = prior), spat)

        mx <- do.call(gseq.pwm, c(args, mode = "max"))
        fin <- mx[is.finite(mx)]
        # thresholds at an exact window score (the maximum), just below and above it, and loose
        threshs <- if (length(fin)) c(max(fin), min(fin), max(fin) - 1e-6, max(fin) + 1e-6, stats::median(fin) - 3) else -1

        for (th in threshs) {
            r <- both_ways(function() do.call(gseq.pwm, c(args, mode = "count", score.thresh = th)))
            expect_identical(r$skip, r$noskip, info = paste("count", it, w, kind, th))
        }
        r <- both_ways(function() do.call(gseq.pwm, c(args, mode = "max")))
        expect_identical(r$skip, r$noskip, info = paste("max", it, w, kind))
        r <- both_ways(function() do.call(gseq.pwm, c(args, mode = "pos", return_strand = TRUE)))
        expect_identical(r$skip, r$noskip, info = paste("pos", it, w, kind))
    }
})

test_that("gseq.pwm skipping is exact with gaps and region bounds", {
    withr::local_options(gmultitasking = FALSE)
    set.seed(7)
    seqs <- fast_scan_seqs(6, 600)
    seqs <- vapply(seqs, function(s) {
        x <- strsplit(s, "")[[1]]
        x[sample(length(x), 20)] <- "-"
        paste(x, collapse = "")
    }, "", USE.NAMES = FALSE)
    for (w in c(3, 8, 15, 24)) {
        pssm <- fast_scan_pssm(w, "sparse")
        mx <- gseq.pwm(seqs, pssm, mode = "max")
        th <- max(mx[is.finite(mx)]) - 2
        for (skip_gaps in c(TRUE, FALSE)) {
            r <- both_ways(function() {
                list(
                    gseq.pwm(seqs, pssm, mode = "count", score.thresh = th, skip_gaps = skip_gaps, start_pos = 50, end_pos = 500, extend = TRUE),
                    gseq.pwm(seqs, pssm, mode = "max", skip_gaps = skip_gaps, start_pos = 30),
                    gseq.pwm(seqs, pssm, mode = "pos", skip_gaps = skip_gaps, end_pos = 400)
                )
            })
            expect_identical(r$skip, r$noskip)
        }
    }
})

test_that("gseq.pwm_hits are the windows gseq.pwm counts", {
    withr::local_options(gmultitasking = FALSE)
    set.seed(11)
    seqs <- fast_scan_seqs(10, 300)
    for (w in c(1, 4, 9, 17, 30)) {
        pssm <- fast_scan_pssm(w, "sparse")
        mx <- gseq.pwm(seqs, pssm, mode = "max")
        th <- stats::median(mx[is.finite(mx)])
        for (bidirect in c(TRUE, FALSE)) {
            for (strand in if (bidirect) 1L else c(1L, -1L)) {
                cnt <- gseq.pwm(seqs, pssm, mode = "count", score.thresh = th, bidirect = bidirect, strand = strand)
                r <- both_ways(function() gseq.pwm_hits(seqs, pssm, score.thresh = th, bidirect = bidirect, strand = strand))
                expect_identical(r$skip, r$noskip)
                h <- r$skip
                expect_equal(as.numeric(tabulate(h$seq_idx, length(seqs))), cnt)
                expect_true(all(h$score >= th))
                expect_equal(h$end - h$start + 1L, rep(w, nrow(h)))
                # each hit's window, scored on its own, is a window of the same score
                if (nrow(h)) {
                    k <- seq_len(min(nrow(h), 20))
                    sub <- substr(seqs[h$seq_idx[k]], h$start[k], h$end[k])
                    expect_equal(gseq.pwm(sub, pssm, mode = "max", bidirect = bidirect, strand = strand), h$score[k])
                    if (bidirect) {
                        fwd <- gseq.pwm(sub, pssm, mode = "max", bidirect = FALSE, strand = 1)
                        rev <- gseq.pwm(sub, pssm, mode = "max", bidirect = FALSE, strand = -1)
                        expect_equal(h$strand[k], ifelse(rev > fwd, -1L, 1L))
                    } else {
                        expect_true(all(h$strand == strand))
                    }
                }
            }
        }
    }
    expect_equal(nrow(gseq.pwm_hits("ACGT", fast_scan_pssm(3, "flat"), score.thresh = 100)), 0L)
    expect_error(gseq.pwm_hits("ACGT", fast_scan_pssm(3, "flat")), "score.thresh")
})

test_that("pwm.count, pwm.max and pwm.max.pos vtracks are identical with and without window skipping", {
    withr::local_options(gmultitasking = FALSE)
    old_root <- get("GROOT", envir = .misha)
    withr::defer(gsetroot(old_root))
    gdb.init_examples()
    withr::defer(gvtrack.clear())
    set.seed(3)

    tiles <- giterator.intervals(intervals = gintervals.all(), iterator = 5000)
    # non-overlapping intervals of very different sizes, one per 25 kb slot
    ragged <- do.call(rbind, lapply(seq(0, 475000, by = 25000), function(s0) {
        s <- s0 + sample(0:1000, 1)
        data.frame(chrom = "chr1", start = s, end = s + sample(c(10, 100, 2000, 23000), 1))
    }))

    for (it in 1:12) {
        w <- sample(c(1, 5, 9, 14, 20, 30), 1)
        pssm <- fast_scan_pssm(w, c("sparse", "tiny", "onehot", "flat")[(it %% 4) + 1])
        bidirect <- it %% 3 != 0
        strand <- if (it %% 6 == 0) -1 else 1
        spat <- if (it %% 5 == 0) list(spat_factor = stats::runif(4, 0.2, 2), spat_bin = 300L) else list()
        base <- c(list(pssm = pssm, bidirect = bidirect, strand = strand, prior = 0.01), spat)

        do.call(gvtrack.create, c(list("fs_max", NULL, "pwm.max"), base))
        mx <- gextract("fs_max", gintervals.all(), iterator = 50000)$fs_max
        th <- max(mx[is.finite(mx)]) - c(0.5, 3, 8)[(it %% 3) + 1]
        do.call(gvtrack.create, c(list("fs_count", NULL, "pwm.count"), base, score.thresh = th))
        do.call(gvtrack.create, c(list("fs_pos", NULL, "pwm.max.pos"), base))
        do.call(gvtrack.create, c(list("fs_count_sh", NULL, "pwm.count"), base, score.thresh = th))
        gvtrack.iterator("fs_count_sh", sshift = -700, eshift = 700)
        do.call(gvtrack.create, c(list("fs_max_sh", NULL, "pwm.max"), base))
        gvtrack.iterator("fs_max_sh", sshift = -700, eshift = 700)
        exprs <- c("fs_count", "fs_max", "fs_pos", "fs_count_sh", "fs_max_sh")

        for (iter in list(tiles, ragged, 1000)) {
            r <- both_ways(function() gextract(exprs, gintervals.all(), iterator = iter))
            expect_identical(r$skip, r$noskip, info = paste(it, w, bidirect, strand))
        }
    }
})

test_that("gintervals.pwm_hits returns exactly the windows pwm.count counts", {
    withr::local_options(gmultitasking = FALSE)
    old_root <- get("GROOT", envir = .misha)
    withr::defer(gsetroot(old_root))
    gdb.init_examples()
    withr::defer(gvtrack.clear())
    set.seed(5)
    intervs <- rbind(
        gintervals(1, c(0, 20000, 499000), c(15000, 21000, 5e5)),
        gintervals(2, 1000, 120000),
        gintervals("chrX", 199990, 2e5)
    )
    for (w in c(1, 6, 12, 25)) {
        pssm <- fast_scan_pssm(w, "sparse")
        for (bidirect in c(TRUE, FALSE)) {
            for (strand in if (bidirect) 1 else c(1, -1)) {
                for (extend in c(TRUE, FALSE)) {
                    gvtrack.create("fs_mx", NULL, "pwm.max", pssm = pssm, bidirect = bidirect, strand = strand, extend = extend)
                    mx <- gextract("fs_mx", intervs, iterator = intervs)$fs_mx
                    th <- max(mx[is.finite(mx)]) - 4
                    gvtrack.create("fs_cnt", NULL, "pwm.count", pssm = pssm, bidirect = bidirect, strand = strand, extend = extend, score.thresh = th)
                    cnt <- gextract("fs_cnt", intervs, iterator = intervs)$fs_cnt
                    r <- both_ways(function() gintervals.pwm_hits(intervs, pssm, score.thresh = th, bidirect = bidirect, strand = strand, extend = extend))
                    expect_identical(r$skip, r$noskip)
                    h <- r$skip
                    n_hits <- if (is.null(h)) rep(0, nrow(intervs)) else tabulate(h$intervalID, nrow(intervs))
                    expect_equal(n_hits, ifelse(is.na(cnt), 0, cnt), info = paste(w, bidirect, strand, extend))
                    if (!is.null(h)) {
                        expect_true(all(h$score >= th))
                        expect_true(all(h$end - h$start == w))
                        expect_true(all(h$start >= intervs$start[h$intervalID] & h$start < intervs$end[h$intervalID]))
                        # a hit is a window whose own pwm.max is its score, at its own start
                        k <- seq_len(min(nrow(h), 15))
                        gvtrack.create("fs_one", NULL, "pwm.max", pssm = pssm, bidirect = bidirect, strand = strand, extend = TRUE)
                        one <- gextract("fs_one", h[k, c("chrom", "start")] |> transform(end = start + 1), iterator = h[k, c("chrom", "start")] |> transform(end = start + 1))
                        expect_equal(one$fs_one, h$score[k], tolerance = 1e-6)
                        if (!bidirect) expect_true(all(h$strand == strand))
                    }
                }
            }
        }
    }
    expect_null(gintervals.pwm_hits(intervs, fast_scan_pssm(5, "flat"), score.thresh = 100))
    expect_error(gintervals.pwm_hits(intervs, fast_scan_pssm(5, "flat")), "score.thresh")
})

test_that("gseq.pwm_hits matches the rows pymisha's gseq_pwm_hits pins", {
    # tests/test_pwm_fast_scan.py in pymisha holds the same rows: both run the
    # same C++ (pwm_fast_scan::seq_hits), so they agree to the last bit
    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1, 0.1, 0.7, 0.1, 0.1,
        0.1, 0.6, 0.2, 0.1, 0.1, 0.1, 0.7, 0.1
    ), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
    seqs <- c("TTACCGTT", "AACGGTAA", "GGGGGGGG", "accgNACCGG")
    h <- gseq.pwm_hits(seqs, pssm, score.thresh = -4)
    expect_identical(h$seq_idx, c(1L, 1L, 1L, 2L, 2L, 2L, 4L, 4L, 4L))
    expect_identical(h$start, c(3L, 4L, 5L, 1L, 2L, 3L, 1L, 6L, 7L))
    expect_identical(h$end, h$start + 3L)
    expect_identical(h$strand, c(1L, -1L, -1L, 1L, 1L, -1L, 1L, 1L, 1L))
    # compared as 17 significant digits: R's parser does not always round a
    # decimal literal to the nearest double, so the literals cannot be compared
    # bit for bit, but these strings identify each double uniquely
    expect_identical(sprintf("%.17g", h$score), c(
        "-1.6773707465831633", "-2.7212820216154241", "-3.5391157588914117",
        "-3.5391157588914117", "-2.7212820216154241", "-1.6773707465831633",
        "-1.6773707465831633", "-1.6773707465831633", "-3.9166391510936496"
    ))
})
