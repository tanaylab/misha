# Regression tests for virtual track defects found in the November 2026 review of 5.14.4.
# Expected values are computed by hand or with base R, not with another misha code path.

vt_fix_db <- function() {
    set.seed(1)
    s1 <- paste(sample(c("A", "C", "G", "T"), 5000, TRUE), collapse = "")
    s2 <- paste(sample(c("A", "C", "G", "T"), 300, TRUE), collapse = "")
    setup_db(list(sprintf(">chr1\n%s\n", s1), sprintf(">chr2\n%s\n", s2)), return_db = TRUE)
}

vt_fix_pssm <- function() {
    set.seed(2)
    pssm <- matrix(runif(24), 6, 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")
    pssm
}

test_that("a NaN bin of a dense track gives no value to exists, size and the position functions", {
    local_db_state()
    vt_fix_db()
    gtrack.create_dense("dn", "x", gintervals("chr1", c(0, 50, 100), c(50, 100, 150)), c(1, NaN, 3), 50, NaN)
    withr::defer(gtrack.rm("dn", force = TRUE))
    gvtrack.create("ex", "dn", "exists")
    gvtrack.create("sz", "dn", "size")
    gvtrack.create("mp", "dn", "max.pos.abs")
    gvtrack.create("fp", "dn", "first.pos.relative")
    withr::defer(gvtrack.clear())

    # iterator equal to the bin size, all functions in one call
    r <- gextract(c("dn", "ex", "sz", "mp", "fp"), gintervals("chr1", 0, 150), iterator = 50)
    expect_equal(r$ex, c(1, 0, 1))
    expect_equal(r$sz, c(1, 0, 1))
    expect_equal(r$mp, c(0, NaN, 100))
    expect_equal(r$fp, c(0, NaN, 0))

    # a position function alone, and an iterator smaller than the bin
    r <- gextract("mp", gintervals("chr1", 50, 100), iterator = 50)
    expect_true(is.na(r$mp))
    r <- gextract(c("ex", "sz", "mp"), gintervals("chr1", 40, 110), iterator = 10)
    in_nan_bin <- r$start >= 50 & r$start < 100
    expect_true(all(r$ex[in_nan_bin] == 0))
    expect_true(all(r$sz[in_nan_bin] == 0))
    expect_true(all(is.na(r$mp[in_nan_bin])))
    expect_equal(r$mp[!in_nan_bin], r$start[!in_nan_bin])
})

test_that("value-based vtracks read only the query's chromosome", {
    local_db_state()
    vt_fix_db()
    vals <- data.frame(chrom = c("chr1", "chr2"), start = c(1000, 10), end = c(1100, 20), score = c(5, 7))
    gvtrack.create("v_avg", vals, "avg")
    gvtrack.create("v_near", vals, "nearest")
    gvtrack.create("v_ex", vals, "exists")
    withr::defer(gvtrack.clear())
    q <- gintervals(c("chr1", "chr1", "chr2"), c(1000, 3000, 10), c(1100, 3100, 20))
    r <- gextract(c("v_avg", "v_near", "v_ex"), q, iterator = q)
    r <- r[order(r$chrom, r$start), ]
    expect_equal(r$v_avg, c(5, NaN, 7))
    expect_equal(r$v_near, c(5, 5, 7))
    expect_equal(r$v_ex, c(1, 0, 1))
})

test_that("a filter that splits a bin or an interval does not count its value twice", {
    local_db_state()
    vt_fix_db()
    gtrack.create_dense("dn", "x", gintervals("chr1", c(0, 50, 100), c(50, 100, 150)), c(1, NaN, 3), 50, NaN)
    gtrack.create_sparse("sp", "x", gintervals("chr1", c(313, 400), c(332, 420)), c(-1.5, 2))
    withr::defer({
        gtrack.rm("dn", force = TRUE)
        gtrack.rm("sp", force = TRUE)
    })
    mask <- gintervals("chr1", c(20, 320), c(30, 322))
    for (f in c("sum", "size", "lse")) {
        gvtrack.create(paste0("d_", f), "dn", f, filter = mask)
        gvtrack.create(paste0("s_", f), "sp", f, filter = mask)
    }
    gvtrack.create("d_q", "dn", "quantile", 0.5, filter = mask)
    withr::defer(gvtrack.clear())

    # bin [0,50) is split into [0,20) and [30,50) by the mask; its value 1 counts once
    iv <- gintervals("chr1", 0, 150)
    r <- gextract(c("d_sum", "d_size", "d_lse", "d_q"), iv, iterator = iv)
    expect_equal(r$d_sum, 4)
    expect_equal(r$d_size, 2)
    expect_equal(r$d_lse, log(exp(1) + exp(3)), tolerance = 1e-6)
    expect_equal(r$d_q, median(c(1, 3)))

    # sparse interval [313,332) is split by [320,322); [400,420) is the second value
    iv <- gintervals("chr1", 300, 450)
    r <- gextract(c("s_sum", "s_size", "s_lse"), iv, iterator = iv)
    expect_equal(r$s_sum, 0.5)
    expect_equal(r$s_size, 2)
    expect_equal(r$s_lse, log(exp(-1.5) + exp(2)), tolerance = 1e-6)

    # the same for a value-based vtrack
    vv <- data.frame(chrom = "chr1", start = c(313, 400), end = c(332, 420), score = c(-1.5, 2))
    gvtrack.create("v_sum", vv, "sum", filter = mask)
    gvtrack.create("v_size", vv, "size", filter = mask)
    gvtrack.create("v_avg", vv, "avg", filter = mask)
    withr::defer(gvtrack.clear())
    r <- gextract(c("v_sum", "v_size", "v_avg"), iv, iterator = iv)
    expect_equal(r$v_sum, 0.5)
    expect_equal(r$v_size, 2)
    expect_equal(r$v_avg, 0.25)
})

test_that("a value-based quantile with a filter uses every unmasked part", {
    local_db_state()
    vt_fix_db()
    vv <- data.frame(chrom = "chr1", start = seq(100, 600, 100), end = seq(150, 650, 100), score = c(1, 9, 2, 8, 3, 7))
    gvtrack.create("v_q", vv, "quantile", 0.5, filter = gintervals("chr1", 360, 380))
    withr::defer(gvtrack.clear())
    iv <- gintervals("chr1", 0, 1000)
    expect_equal(gextract("v_q", iv, iterator = iv)$v_q, median(vv$score))
})

test_that("edit distance positions are offset when a filter leaves one part that does not start the interval", {
    local_db_state()
    vt_fix_db()
    pssm <- vt_fix_pssm()
    gvtrack.create("pm", NULL, "pwm.max", pssm = pssm)
    gvtrack.create("p", NULL, "pwm", pssm = pssm)
    withr::defer(gvtrack.clear())
    part <- gintervals("chr1", 110, 140)
    whole <- gintervals("chr1", 100, 140)
    scores <- gextract(c("pm", "p"), part, iterator = part)
    # thresholds half a unit above the current score, so that at least one edit is needed
    thresh <- c(pwm.edit_distance.pos = scores$pm + 0.5, pwm.edit_distance.lse.pos = scores$p + 0.5)

    for (f in names(thresh)) {
        gvtrack.create("ed", NULL, f, pssm = pssm, score.thresh = thresh[[f]])
        unfiltered <- gextract("ed", part, iterator = part)$ed
        expect_false(is.na(unfiltered), info = f)
        gvtrack.filter("ed", gintervals("chr1", 100, 110))
        filtered <- gextract("ed", whole, iterator = whole)$ed
        expect_equal(abs(filtered), abs(unfiltered) + 10, info = f)
        expect_equal(sign(filtered), sign(unfiltered), info = f)
        gvtrack.rm("ed")
    }
})

test_that("gtrack.array.extract reads columns given in any order", {
    asc <- gtrack.array.extract("test.array", c("col1", "col5"), gintervals(c(1, 2)))
    desc <- gtrack.array.extract("test.array", c("col5", "col1"), gintervals(c(1, 2)))
    expect_equal(desc$col1, asc$col1)
    expect_equal(desc$col5, asc$col5)
    expect_true(any(!is.na(desc$col1)))
})

test_that("min.pos on an array track reports the position of the minimum", {
    local_db_state()
    vt_fix_db()
    f <- withr::local_tempfile(fileext = ".tsv")
    writeLines(c(
        "chrom\tstart\tend\ta\tb",
        "chr1\t0\t10\t5\t5",
        "chr1\t20\t30\t1\t1",
        "chr1\t40\t50\t3\t3"
    ), f)
    suppressMessages(gtrack.array.import("arr", "x", f))
    withr::defer(gtrack.rm("arr", force = TRUE))
    gvtrack.create("a_minp", "arr", "min.pos.abs")
    gvtrack.create("a_minr", "arr", "min.pos.relative")
    gvtrack.create("a_maxp", "arr", "max.pos.abs")
    withr::defer(gvtrack.clear())
    iv <- gintervals("chr1", 0, 100)
    r <- gextract(c("a_minp", "a_minr", "a_maxp"), iv, iterator = iv)
    expect_equal(r$a_minp, 20)
    expect_equal(r$a_minr, 20)
    expect_equal(r$a_maxp, 0)
})

test_that("stddev with a filter is the sample sd of the unmasked values, each counted once", {
    local_db_state()
    vt_fix_db()
    starts <- seq(0, 950, 50)
    gtrack.create_dense("dn2", "x", gintervals("chr1", starts, starts + 50), as.numeric(1:20), 50, NaN)
    gtrack.create_sparse("sp2", "x", gintervals("chr1", c(10, 100, 200, 300), c(60, 120, 220, 320)), c(4, 1, 7, 2))
    withr::defer({
        gtrack.rm("dn2", force = TRUE)
        gtrack.rm("sp2", force = TRUE)
    })
    withr::defer(gvtrack.clear())
    # [20,30) splits bin 1 (and the first sparse interval); [500,600) removes bins 11 and 12; [190,230) removes 7
    mask <- gintervals("chr1", c(20, 190, 500), c(30, 230, 600))
    gvtrack.create("d_sd", "dn2", "stddev", filter = mask)
    gvtrack.create("d_avg", "dn2", "avg", filter = mask)
    gvtrack.create("s_sd", "sp2", "stddev", filter = mask)
    vv <- data.frame(chrom = "chr1", start = c(10, 100, 200, 300), end = c(60, 120, 220, 320), score = c(4, 1, 7, 2))
    gvtrack.create("v_sd", vv, "stddev", filter = mask)
    iv <- gintervals("chr1", 0, 1000)

    expected_dense <- sd(c(1:10, 13:20))
    expect_equal(gextract("d_sd", iv, iterator = iv)$d_sd, expected_dense, tolerance = 1e-5)
    expect_equal(gextract(c("d_sd", "d_avg"), iv, iterator = iv)$d_sd, expected_dense, tolerance = 1e-5)
    expect_equal(gextract("s_sd", iv, iterator = iv)$s_sd, sd(c(4, 1, 2)), tolerance = 1e-5)
    expect_equal(gextract("v_sd", iv, iterator = iv)$v_sd, sd(c(4, 1, 2)), tolerance = 1e-5)
})

test_that("positions of sparse, value-based and array values are clipped to the iterator start", {
    local_db_state()
    vt_fix_db()
    gtrack.create_sparse("sp3", "x", gintervals("chr1", c(313, 400), c(332, 420)), c(5, 1))
    f <- withr::local_tempfile(fileext = ".tsv")
    writeLines(c("chrom\tstart\tend\ta", "chr1\t313\t332\t5", "chr1\t400\t420\t1"), f)
    suppressMessages(gtrack.array.import("arr3", "x", f))
    withr::defer({
        gtrack.rm("sp3", force = TRUE)
        gtrack.rm("arr3", force = TRUE)
    })
    withr::defer(gvtrack.clear())
    vv <- data.frame(chrom = "chr1", start = c(313, 400), end = c(332, 420), score = c(5, 1))
    iv <- gintervals("chr1", 318, 450)
    one <- gintervals("chr1", 318, 330)

    for (src in c("sp3", "vv")) {
        source_obj <- if (src == "vv") vv else src
        for (fn in c("first", "last", "min", "max", "sample")) {
            gvtrack.create("pa", source_obj, paste0(fn, ".pos.abs"))
            gvtrack.create("pr", source_obj, paste0(fn, ".pos.relative"))
            r <- gextract(c("pa", "pr"), iv, iterator = iv)
            expected <- c(first = 318, last = 400, min = 400, max = 318)
            if (fn != "sample") {
                expect_equal(r$pa, expected[[fn]], info = paste(src, fn))
                expect_equal(r$pr, expected[[fn]] - 318, info = paste(src, fn))
            }
            r <- gextract(c("pa", "pr"), one, iterator = one)
            expect_equal(r$pa, 318, info = paste(src, fn, "one value"))
            expect_equal(r$pr, 0, info = paste(src, fn, "one value"))
        }
    }
    gvtrack.create("a_max", "arr3", "max.pos.relative")
    gvtrack.create("a_min", "arr3", "min.pos.relative")
    r <- gextract(c("a_max", "a_min"), iv, iterator = iv)
    expect_equal(r$a_max, 0)
    expect_equal(r$a_min, 400 - 318)
})

# Scores of the windows starting at 0-based offsets `anchors` of `bases`, read on the forward (strand 1) or the
# reverse (strand -1) strand. A window holding anything but A, C, G or T is not scored (-Inf).
vt_ref_window_scores <- function(bases, lp, anchors, strand) {
    L <- nrow(lp)
    comp <- c(A = "T", C = "G", G = "C", T = "A")
    vapply(anchors, function(a) {
        w <- bases[a + seq_len(L)]
        if (any(!w %in% c("A", "C", "G", "T"))) {
            return(-Inf)
        }
        if (strand == -1) w <- rev(comp[w])
        sum(lp[cbind(seq_len(L), match(w, c("A", "C", "G", "T")))])
    }, numeric(1))
}

vt_lse <- function(x) {
    m <- max(x)
    if (!is.finite(m)) {
        return(m)
    }
    m + log(sum(exp(x - m)))
}

test_that("pwm vtracks with bidirect = FALSE and strand = -1 count spatial bins and spat_min/spat_max from the 5' start", {
    local_db_state()
    vt_fix_db()
    withr::defer(gvtrack.clear())
    pssm <- vt_fix_pssm()
    lp <- log(pssm)
    L <- nrow(pssm)
    bases <- strsplit(toupper(gseq.extract(gintervals("chr1", 3000, 3300))), "")[[1]]

    # windows of 60 bp stepping by 3, so the sliding paths run too
    starts <- 3100 + seq(0, 57, 3)
    ivs <- gintervals("chr1", starts, starts + 60)
    spat <- c(1, 50, 3)
    gvtrack.create("m_pwm", NULL, "pwm", pssm = pssm, prior = 0, bidirect = FALSE, strand = -1, spat_factor = spat, spat_bin = 25L)
    gvtrack.create("m_max", NULL, "pwm.max", pssm = pssm, prior = 0, bidirect = FALSE, strand = -1, spat_factor = spat, spat_bin = 25L)
    gvtrack.create("m_pos", NULL, "pwm.max.pos", pssm = pssm, prior = 0, bidirect = FALSE, strand = -1, spat_factor = spat, spat_bin = 25L)
    gvtrack.create("m_rpos", NULL, "pwm.max.pos", pssm = pssm, prior = 0, bidirect = FALSE, strand = -1, spat_min = 1, spat_max = 20)
    gvtrack.create("m_rmax", NULL, "pwm.max", pssm = pssm, prior = 0, bidirect = FALSE, strand = -1, spat_min = 1, spat_max = 20)
    r <- gextract(c("m_pwm", "m_max", "m_pos", "m_rpos", "m_rmax"), ivs, iterator = ivs)
    r <- r[order(r$intervalID), ]

    for (k in seq_along(starts)) {
        anchors <- 0:59
        sc <- vt_ref_window_scores(bases, lp, starts[k] - 3000 + anchors, -1)
        w <- sc + log(spat[pmin(anchors %/% 25 + 1, length(spat))])
        expect_equal(r$m_pwm[k], vt_lse(w), tolerance = 1e-4)
        expect_equal(r$m_max[k], max(w), tolerance = 1e-4)
        expect_equal(r$m_pos[k], which.max(w))
        in_range <- anchors <= 20 - L
        expect_equal(r$m_rmax[k], max(sc[in_range]), tolerance = 1e-4)
        expect_equal(r$m_rpos[k], which.max(sc[in_range]))
    }
})

test_that("pwm vtracks never score a window holding N or another non-ACGT base, on any path", {
    local_db_state()
    set.seed(5)
    s <- sample(c("A", "C", "G", "T"), 400, TRUE)
    s[101:106] <- "N"
    s[201] <- "R"
    s[250:255] <- "Y"
    setup_db(list(sprintf(">chr1\n%s\n", paste(s, collapse = ""))), return_db = TRUE)
    withr::defer(gvtrack.clear())
    pssm <- vt_fix_pssm()
    lp <- log(pssm)
    L <- nrow(pssm)

    gvtrack.create("n_pwm", NULL, "pwm", pssm = pssm, prior = 0)
    gvtrack.create("n_max", NULL, "pwm.max", pssm = pssm, prior = 0)
    gvtrack.create("n_pos", NULL, "pwm.max.pos", pssm = pssm, prior = 0, bidirect = FALSE)
    gvtrack.create("n_fmax", NULL, "pwm.max", pssm = pssm, prior = 0, bidirect = FALSE)
    gvtrack.create("n_spw", NULL, "pwm.max", pssm = pssm, prior = 0, spat_factor = c(1, 2), spat_bin = 5L)
    gvtrack.create("n_sp", NULL, "pwm", pssm = pssm, prior = 0, spat_factor = c(1, 2), spat_bin = 5L)

    # intervals next to and over the N run, the R and the Y run; anchors whose windows reach past an
    # interval's end use the bases after it (extend = TRUE)
    starts <- c(90, 95, 100, 190, 196, 240, 245, 300)
    for (st in starts) {
        iv <- gintervals("chr1", st, st + 10)
        anchors <- st + 0:9
        fwd <- vt_ref_window_scores(s, lp, anchors, 1)
        rev <- vt_ref_window_scores(s, lp, anchors, -1)
        both <- mapply(function(a, b) vt_lse(c(a, b)), fwd, rev)
        spat_w <- log(c(1, 2)[pmin(0:9 %/% 5 + 1, 2)])
        r <- gextract(c("n_pwm", "n_max", "n_pos", "n_fmax", "n_spw", "n_sp"), iv, iterator = iv)
        expect_equal(r$n_pwm, vt_lse(both), tolerance = 1e-4, info = paste(st))
        expect_equal(r$n_max, max(both), tolerance = 1e-4, info = paste(st))
        expect_equal(r$n_fmax, max(fwd), tolerance = 1e-4, info = paste(st))
        if (is.finite(max(fwd))) {
            expect_equal(r$n_pos, which.max(fwd), info = paste(st))
        } else {
            expect_true(is.na(r$n_pos), info = paste(st))
        }
        expect_equal(r$n_spw, max(both + spat_w), tolerance = 1e-4, info = paste(st))
        expect_equal(r$n_sp, vt_lse(both + spat_w), tolerance = 1e-4, info = paste(st))

        # the spatial value of an interval does not depend on the interval scored before it
        two <- gintervals("chr1", c(10, st), c(20, st + 10))
        expect_equal(gextract("n_spw", two, iterator = two)$n_spw[2], r$n_spw, info = paste(st))
    }

    # the same over overlapping windows in one call, where the sliding caches run
    ivs <- gintervals("chr1", 85:120, 95:130)
    r <- gextract(c("n_pwm", "n_max", "n_spw", "n_sp"), ivs, iterator = ivs)
    r <- r[order(r$intervalID), ]
    spat_w <- log(c(1, 2)[pmin(0:9 %/% 5 + 1, 2)])
    for (k in seq_len(nrow(ivs))) {
        anchors <- ivs$start[k] + 0:9
        both <- mapply(function(a, b) vt_lse(c(a, b)), vt_ref_window_scores(s, lp, anchors, 1), vt_ref_window_scores(s, lp, anchors, -1))
        expect_equal(r$n_pwm[k], vt_lse(both), tolerance = 1e-4, info = paste(ivs$start[k]))
        expect_equal(r$n_max[k], max(both), tolerance = 1e-4, info = paste(ivs$start[k]))
        expect_equal(r$n_spw[k], max(both + spat_w), tolerance = 1e-4, info = paste(ivs$start[k]))
        expect_equal(r$n_sp[k], vt_lse(both + spat_w), tolerance = 1e-4, info = paste(ivs$start[k]))
    }
})

test_that("pwm.edit_distance.lse with bidirect = TRUE answers for the strand-combined LSE (brute force)", {
    local_db_state()
    vt_fix_db()
    withr::defer(gvtrack.clear())
    pssm <- vt_fix_pssm()
    lp <- log(pssm)
    L <- nrow(pssm)
    len <- 8
    iv <- gintervals("chr1", 3100, 3100 + len)
    bases <- strsplit(toupper(gseq.extract(gintervals("chr1", 3100, 3100 + len + L - 1))), "")[[1]]
    R <- length(bases)
    combined <- function(b) {
        vt_lse(c(vt_ref_window_scores(b, lp, 0:(len - 1), 1), vt_ref_window_scores(b, lp, 0:(len - 1), -1)))
    }
    gvtrack.create("e_pwm", NULL, "pwm", pssm = pssm, prior = 0)
    F0 <- combined(bases)
    expect_equal(gextract("e_pwm", iv, iterator = iv)$e_pwm, F0, tolerance = 1e-4)

    # every sequence one and two substitutions away, with the edited positions
    alt <- function(b) setdiff(c("A", "C", "G", "T"), b)
    one <- do.call(rbind, lapply(seq_len(R), function(q) {
        data.frame(q = q, F = vapply(alt(bases[q]), function(y) {
            b <- bases
            b[q] <- y
            combined(b)
        }, numeric(1)))
    }))
    pairs <- t(combn(R, 2))
    two <- unlist(lapply(seq_len(nrow(pairs)), function(k) {
        q1 <- pairs[k, 1]
        q2 <- pairs[k, 2]
        g <- expand.grid(y1 = alt(bases[q1]), y2 = alt(bases[q2]), stringsAsFactors = FALSE)
        mapply(function(y1, y2) {
            b <- bases
            b[q1] <- y1
            b[q2] <- y2
            combined(b)
        }, g$y1, g$y2)
    }))
    reach <- function(thresh, above) {
        hit <- function(x) if (above) any(x >= thresh) else any(x <= thresh)
        if (hit(F0)) 0 else if (hit(one$F)) 1 else if (hit(two)) 2 else NA
    }
    # thresholds in the gaps between reachable scores, away from any tie
    vals <- sort(unique(round(c(F0, one$F, two), 6)))
    gaps <- which(diff(vals) > 1e-3)
    pick <- gaps[unique(round(seq(1, length(gaps), length.out = 12)))]
    thresholds <- (vals[pick] + vals[pick + 1]) / 2

    for (thresh in thresholds) {
        for (dirn in c("above", "below")) {
            gvtrack.create("e_ed", NULL, "pwm.edit_distance.lse", pssm = pssm, prior = 0, score.thresh = thresh, max_edits = 2, direction = dirn)
            gvtrack.create("e_pos", NULL, "pwm.edit_distance.lse.pos", pssm = pssm, prior = 0, score.thresh = thresh, max_edits = 2, direction = dirn)
            r <- gextract(c("e_ed", "e_pos"), iv, iterator = iv)
            expected <- reach(thresh, dirn == "above")
            if (is.na(expected)) {
                expect_true(is.na(r$e_ed), info = paste(thresh, dirn))
            } else {
                expect_equal(r$e_ed, expected, info = paste(thresh, dirn))
            }
            if (identical(expected, 1)) {
                # the reported position is one where a single substitution reaches the threshold
                at <- one$F[one$q == abs(r$e_pos)]
                expect_true(if (dirn == "above") any(at >= thresh) else any(at <= thresh), info = paste(thresh, dirn))
            }
        }
    }
})

test_that("avg with a filter is the mean of the unmasked values, each counted once", {
    local_db_state()
    vt_fix_db()
    starts <- seq(0, 950, 50)
    gtrack.create_dense("dn4", "x", gintervals("chr1", starts, starts + 50), as.numeric(1:20), 50, NaN)
    gtrack.create_sparse("sp4", "x", gintervals("chr1", c(10, 100, 200, 300), c(60, 120, 220, 320)), c(4, 1, 7, 2))
    f <- withr::local_tempfile(fileext = ".tsv")
    writeLines(c("chrom\tstart\tend\ta", "chr1\t10\t60\t4", "chr1\t100\t120\t1", "chr1\t200\t220\t7", "chr1\t300\t320\t2"), f)
    suppressMessages(gtrack.array.import("arr4", "x", f))
    withr::defer({
        gtrack.rm("dn4", force = TRUE)
        gtrack.rm("sp4", force = TRUE)
        gtrack.rm("arr4", force = TRUE)
    })
    withr::defer(gvtrack.clear())
    # [20,30) splits bin 1 and the first interval; [190,230) removes the value 7; [500,600) removes bins 11, 12
    mask <- gintervals("chr1", c(20, 190, 500), c(30, 230, 600))
    gvtrack.create("d_avg", "dn4", "avg", filter = mask)
    gvtrack.create("s_avg", "sp4", "avg", filter = mask)
    gvtrack.create("a_avg", "arr4", "avg", filter = mask)
    iv <- gintervals("chr1", 0, 1000)
    expect_equal(gextract("d_avg", iv, iterator = iv)$d_avg, mean(c(1:10, 13:20)), tolerance = 1e-6)
    expect_equal(gextract("s_avg", iv, iterator = iv)$s_avg, mean(c(4, 1, 2)), tolerance = 1e-6)
    expect_equal(gextract("a_avg", iv, iterator = iv)$a_avg, mean(c(4, 1, 2)), tolerance = 1e-6)
})

vt_n_db <- function() {
    set.seed(5)
    s <- sample(c("A", "C", "G", "T"), 400, TRUE)
    s[101:106] <- "N"
    s[201] <- "R"
    s[300:330] <- "N"
    setup_db(list(sprintf(">chr1\n%s\n", paste(s, collapse = ""))), return_db = TRUE)
    s
}

test_that("gseq.pwm leaves out windows holding N by default, as the pwm vtracks do", {
    local_db_state()
    vt_n_db()
    withr::defer(gvtrack.clear())
    pssm <- vt_fix_pssm()
    for (f in c("pwm", "pwm.max", "pwm.max.pos")) gvtrack.create(gsub(".", "_", f, fixed = TRUE), NULL, f, pssm = pssm)
    gvtrack.create("pwm_count", NULL, "pwm.count", pssm = pssm, score.thresh = -7)
    for (st in c(95, 98, 195, 296)) {
        iv <- gintervals("chr1", st, st + 10)
        v <- gextract(c("pwm", "pwm_max", "pwm_max_pos", "pwm_count"), iv, iterator = iv)
        s <- gseq.extract(gintervals("chr1", st, st + 15))
        g <- function(m, ...) gseq.pwm(s, pssm, mode = m, end_pos = 10, extend = TRUE, ...)
        expect_equal(g("lse"), v$pwm, tolerance = 1e-5, info = paste(st))
        expect_equal(g("max"), v$pwm_max, tolerance = 1e-5, info = paste(st))
        if (is.finite(v$pwm_max)) {
            expect_equal(abs(g("pos")), abs(v$pwm_max_pos), info = paste(st))
        } else {
            expect_true(is.na(g("pos")), info = paste(st))
            expect_true(is.na(v$pwm_max_pos), info = paste(st))
        }
        expect_equal(g("count", score.thresh = -7), v$pwm_count, info = paste(st))
    }
    # nothing scorable
    expect_equal(gseq.pwm("NNNNNNNNNN", pssm, mode = "lse"), -Inf)
    expect_equal(gseq.pwm("NNNNNNNNNN", pssm, mode = "max"), -Inf)
    expect_true(is.na(gseq.pwm("NNNNNNNNNN", pssm, mode = "pos")))
    expect_true(is.na(gseq.pwm("RRRRRRRRRR", pssm, mode = "pos")))
    expect_equal(gseq.pwm("NNNNNNNNNN", pssm, mode = "count", score.thresh = -100), 0)
    # scoring N on request still charges the column mean
    lp <- log((pssm + 0.01) / rowSums(pssm + 0.01))
    expect_equal(
        gseq.pwm("ACGNNN", pssm, mode = "max", bidirect = FALSE, neutral_chars = "N"),
        sum(lp[cbind(1:3, 1:3)]) + sum(rowMeans(lp)[4:6]),
        tolerance = 1e-5
    )
})

# Fewest substitutions that bring one window's score to the threshold (NA when none do), from the
# per-column gains: the largest gains first, which is exact for a position weight matrix.
vt_ref_window_edits <- function(w, lp, thresh, above) {
    cols <- match(w, c("A", "C", "G", "T"))
    cur <- lp[cbind(seq_along(cols), cols)]
    score <- sum(cur)
    gains <- sort(if (above) apply(lp, 1, max) - cur else cur - apply(lp, 1, min), decreasing = TRUE)
    need <- if (above) thresh - score else score - thresh
    if (need <= 0) {
        return(0)
    }
    k <- which(cumsum(gains) >= need - 1e-9)
    if (length(k) == 0) NA else min(k)
}

test_that("edit distances never use a window holding N or another non-ACGT base", {
    local_db_state()
    s <- vt_n_db()
    withr::defer(gvtrack.clear())
    pssm <- vt_fix_pssm()
    lp <- log(pssm)
    L <- nrow(pssm)
    windows_of <- function(st) lapply(st + 0:9, function(a) s[a + seq_len(L)])
    acgt <- function(w) all(w %in% c("A", "C", "G", "T"))

    for (st in c(92, 97, 195)) {
        iv <- gintervals("chr1", st, st + 10)
        ws <- windows_of(st)
        ok <- vapply(ws, acgt, logical(1))
        scores <- vapply(ws, function(w) if (acgt(w)) sum(lp[cbind(1:L, match(w, c("A", "C", "G", "T")))]) else -Inf, numeric(1))
        for (dirn in c("above", "below")) {
            thresh <- if (dirn == "above") max(scores) + 1.5 else max(scores[ok]) - 1.5
            ref <- vapply(seq_along(ws), function(k) if (ok[k]) vt_ref_window_edits(ws[[k]], lp, thresh, dirn == "above") else NA_real_, numeric(1))
            for (f in c("pwm.edit_distance", "pwm.edit_distance.pos", "pwm.max.edit_distance")) {
                gvtrack.create("ed", NULL, f, pssm = pssm, prior = 0, bidirect = FALSE, score.thresh = thresh, direction = dirn)
                got <- gextract("ed", iv, iterator = iv)$ed
                expected <- switch(f,
                    pwm.edit_distance = if (all(is.na(ref))) NA else min(ref, na.rm = TRUE),
                    pwm.edit_distance.pos = if (all(is.na(ref))) NA else which.min(ref),
                    pwm.max.edit_distance = ref[which.max(scores)]
                )
                if (is.na(expected)) {
                    expect_true(is.na(got), info = paste(st, dirn, f))
                } else {
                    expect_equal(got, expected, info = paste(st, dirn, f))
                }
            }
            # the string version agrees
            if (!all(is.na(ref)) && dirn == "above") {
                e <- gseq.pwm_edits(gseq.extract(gintervals("chr1", st, st + 10 + L - 1)), pssm,
                    score.thresh = thresh, bidirect = FALSE, prior = 0, extend = FALSE
                )
                expect_equal(e$n_edits[1], min(ref, na.rm = TRUE), info = paste(st, "gseq.pwm_edits"))
            }
        }
    }

    # an interval whose every window holds an N: no candidate window, NA for the whole family
    iv <- gintervals("chr1", 302, 312)
    fns <- c(
        "pwm.edit_distance", "pwm.edit_distance.pos", "pwm.max.edit_distance", "pwm.n_mutations",
        "pwm.edit_distance.lse", "pwm.edit_distance.lse.pos"
    )
    for (f in fns) {
        for (dirn in c("above", "below")) {
            gvtrack.create("ed", NULL, f, pssm = pssm, score.thresh = -5, direction = dirn)
            expect_true(is.na(gextract("ed", iv, iterator = iv)$ed), info = paste(f, dirn))
        }
    }
})

test_that("pwm.edit_distance.lse leaves out windows holding N (brute force)", {
    local_db_state()
    s <- vt_n_db()
    withr::defer(gvtrack.clear())
    pssm <- vt_fix_pssm()
    lp <- log(pssm)
    L <- nrow(pssm)
    len <- 8
    st <- 89 # windows 89..96: the N run at 100..105 (0-based) cuts the last two
    iv <- gintervals("chr1", st, st + len)
    bases <- s[st + seq_len(len + L - 1)]
    editable <- which(bases %in% c("A", "C", "G", "T"))
    combined <- function(b) {
        vt_lse(c(vt_ref_window_scores(b, lp, 0:(len - 1), 1), vt_ref_window_scores(b, lp, 0:(len - 1), -1)))
    }
    gvtrack.create("e_pwm", NULL, "pwm", pssm = pssm, prior = 0)
    F0 <- combined(bases)
    expect_true(is.finite(F0))
    expect_equal(gextract("e_pwm", iv, iterator = iv)$e_pwm, F0, tolerance = 1e-4)
    alt <- function(b) setdiff(c("A", "C", "G", "T"), b)
    one <- unlist(lapply(editable, function(q) {
        vapply(alt(bases[q]), function(y) {
            b <- bases
            b[q] <- y
            combined(b)
        }, numeric(1))
    }))
    pairs <- t(combn(editable, 2))
    two <- unlist(lapply(seq_len(nrow(pairs)), function(k) {
        g <- expand.grid(y1 = alt(bases[pairs[k, 1]]), y2 = alt(bases[pairs[k, 2]]), stringsAsFactors = FALSE)
        mapply(function(y1, y2) {
            b <- bases
            b[pairs[k, 1]] <- y1
            b[pairs[k, 2]] <- y2
            combined(b)
        }, g$y1, g$y2)
    }))
    vals <- sort(unique(round(c(F0, one, two), 6)))
    gaps <- which(diff(vals) > 1e-3)
    pick <- gaps[unique(round(seq(1, length(gaps), length.out = 10)))]
    for (thresh in (vals[pick] + vals[pick + 1]) / 2) {
        for (dirn in c("above", "below")) {
            hit <- function(x) if (dirn == "above") any(x >= thresh) else any(x <= thresh)
            expected <- if (hit(F0)) 0 else if (hit(one)) 1 else if (hit(two)) 2 else NA
            gvtrack.create("e_ed", NULL, "pwm.edit_distance.lse", pssm = pssm, prior = 0, score.thresh = thresh, max_edits = 2, direction = dirn)
            got <- gextract("e_ed", iv, iterator = iv)$e_ed
            if (is.na(expected)) {
                expect_true(is.na(got), info = paste(thresh, dirn))
            } else {
                expect_equal(got, expected, info = paste(thresh, dirn))
            }
        }
    }
})

test_that("gtrack.create_pwm_energy leaves out windows holding N; an interval with none scorable is NaN", {
    local_db_state()
    s <- vt_n_db()
    pssm <- vt_fix_pssm()
    lp <- log(pssm)
    L <- nrow(pssm)
    dir.create(file.path(.misha$GROOT, "pssms"), showWarnings = FALSE)
    writeLines("0 x 0", file.path(.misha$GROOT, "pssms", "tset.key"))
    writeLines(
        sprintf("0\t%d\t%.10f\t%.10f\t%.10f\t%.10f", 0:(L - 1), pssm[, 1], pssm[, 2], pssm[, 3], pssm[, 4]),
        file.path(.misha$GROOT, "pssms", "tset.data")
    )
    gtrack.create_pwm_energy("en", "x", "tset", 0, 0, iterator = 10)
    withr::defer(gtrack.rm("en", force = TRUE))
    r <- gextract("en", gintervals("chr1", 60, 340))
    for (k in seq_len(nrow(r))) {
        # windows anchored at start..end of the bin (the sequence read is [start, end + L))
        anchors <- r$start[k]:r$end[k]
        anchors <- anchors[anchors + L <= 400]
        expected <- vt_lse(vt_ref_window_scores(s, lp, anchors, 1))
        # a track stores -Inf as NaN
        if (expected == -Inf) {
            expect_true(is.na(r$en[k]), info = paste(r$start[k]))
        } else {
            expect_equal(r$en[k], expected, tolerance = 1e-4, info = paste(r$start[k]))
        }
    }
    expect_true(any(is.na(r$en)))
})
