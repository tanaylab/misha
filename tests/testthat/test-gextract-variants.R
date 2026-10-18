create_isolated_test_db()

test_that("reach: a bin [s, e) of each sequence vtrack reads [s + lo, e + hi)", {
    define_variant_vtracks()
    pssm <- variant_test_pssm()
    gvtrack.create("r_noext", NULL, "pwm.max", params = list(pssm = pssm, bidirect = TRUE, extend = FALSE, prior = 0.01))
    gvtrack.create("r_kmer_fwd", NULL, "kmer.count", params = list(kmer = "ACG", extend = TRUE, strand = 1))
    gvtrack.create("r_kmer_noext", NULL, "kmer.count", params = list(kmer = "ACG", extend = FALSE, strand = 0))
    reach <- function(expr) misha:::.gvariants_reach(expr)
    expect_equal(reach("v_pwm"), cbind(lo = -5, hi = 3 + 3))
    expect_equal(reach("r_noext"), cbind(lo = 0, hi = 0))
    expect_equal(reach("v_potts"), cbind(lo = 4, hi = 4 + 5))
    expect_equal(reach("v_kmer"), cbind(lo = -2 - 2, hi = 2 + 2))
    expect_equal(reach("r_kmer_fwd"), cbind(lo = 0, hi = 2))
    expect_equal(reach("r_kmer_noext"), cbind(lo = 0, hi = 0))
    expect_equal(reach("v_masked"), cbind(lo = -3, hi = 3))
    expect_equal(reach("v_ed"), cbind(lo = 0, hi = 4 + 1 - 1))
    expect_equal(nrow(reach("v_pwm + 2 * v_potts + test.fixedbin")), 2L)
})

test_that("reach: an expression with no sequence vtrack is an error", {
    expect_error(misha:::.gvariants_reach("test.fixedbin * 2"), "no sequence virtual track")
})

test_that("engine: a substituted row equals the same bin on a really edited genome", {
    local_db_state()
    seq <- variant_test_seq()
    db <- make_seq_db(seq)
    subs <- data.frame(pos = c(100, 101, 150, 151, 152, 200), base = c("A", "C", "G", "T", "A", "C"))
    subs <- subs[toupper(substr(seq, subs$pos + 1, subs$pos + 1)) != subs$base, ]
    # A for loop, not lapply: make_seq_db() deletes its database when the
    # calling frame exits, and an lapply FUN frame exits at once.
    edited <- list()
    for (s in seq_len(nrow(subs))) {
        edited[[s]] <- make_seq_db(substitute_bases(seq, subs$pos[s], subs$base[s]))
    }
    bins <- 90:210
    # Bin-major: consecutive rows share a bin and differ in the substitution,
    # which is where a cache that ignores the substitution would leak.
    rows <- expand.grid(sub = seq_len(nrow(subs)), bin = bins)
    for (vt in VARIANT_VTRACKS) {
        suppressMessages(gsetroot(db))
        define_variant_vtracks()
        got <- misha:::.gvariants_eval(
            vt, rep("chr1", nrow(rows)), rows$bin, rows$bin + 1,
            subs$pos[rows$sub], match(subs$base[rows$sub], VARIANT_BASES) - 1L
        )
        for (s in seq_len(nrow(subs))) {
            expect_equal(got[rows$sub == s], values_on_db(edited[[s]], vt, bins),
                info = sprintf("%s, substitution %d", vt, s)
            )
        }
    }
})

test_that("engine: a row with no substitution after substituted rows reads the reference", {
    local_db_state()
    seq <- variant_test_seq(seed = 8L)
    db <- make_seq_db(seq)
    b <- 120
    alt <- setdiff(VARIANT_BASES, toupper(substr(seq, b + 1, b + 1)))[1]
    for (vt in VARIANT_VTRACKS) {
        want <- values_on_db(db, vt, b)
        got <- misha:::.gvariants_eval(
            vt, rep("chr1", 3), rep(b, 3), rep(b + 1, 3), rep(b, 3),
            c(match(alt, VARIANT_BASES) - 1L, -1L, -1L)
        )
        expect_equal(got[2:3], rep(want, 2), info = vt)
    }
})

test_that("engine: unsorted rows are an error", {
    local_db_state()
    db <- make_seq_db(variant_test_seq())
    suppressMessages(gsetroot(db))
    define_variant_vtracks()
    expect_error(misha:::.gvariants_eval("v_pwm", c("chr1", "chr1"), c(50, 40), c(51, 41), c(50, 40), c(0L, 0L)), "sorted")
})

test_that("saturation: every non-N base in reach, three alts each", {
    define_variant_vtracks()
    res <- gextract.variants("v_pwm", gintervals(1, 1000, 1001), iterator = 1, colnames = "s")
    expect_equal(colnames(res), c("chrom", "start", "end", "pos", "ref", "alt", "s.ref", "s.alt", "intervalID"))
    genome <- toupper(strsplit(gseq.extract(gintervals(1, 995, 1007)), "")[[1]])
    acgt <- (995:1006)[genome %in% VARIANT_BASES]
    expect_setequal(unique(res$pos), acgt)
    expect_equal(nrow(res), 3L * length(acgt))
    expect_equal(res$ref, genome[res$pos - 995 + 1])
    expect_true(all(res$alt != res$ref))
    expect_equal(unique(res$s.ref), gextract("v_pwm", gintervals(1, 1000, 1001), iterator = 1)$v_pwm)
})

test_that("values: .ref is the reference and .alt the really edited genome", {
    local_db_state()
    seq <- variant_test_seq(seed = 9L)
    db <- make_seq_db(seq)
    suppressMessages(gsetroot(db))
    define_variant_vtracks()
    res <- gextract.variants("v_pwm + v_potts", gintervals("chr1", 100, 104), iterator = 1, colnames = "s")
    set.seed(1)
    for (r in sample(nrow(res), 4)) {
        b <- res$start[r]
        expect_equal(res$s.ref[r], values_on_db(db, "v_pwm", b) + values_on_db(db, "v_potts", b))
        edb <- make_seq_db(substitute_bases(seq, res$pos[r], res$alt[r]))
        expect_equal(res$s.alt[r], values_on_db(edb, "v_pwm", b) + values_on_db(edb, "v_potts", b), info = r)
    }
})

test_that("reach is complete: changing every base outside it leaves the value unchanged", {
    local_db_state()
    seq <- variant_test_seq(seed = 10L)
    db <- make_seq_db(seq)
    b <- 200
    for (vt in VARIANT_VTRACKS) {
        suppressMessages(gsetroot(db))
        define_variant_vtracks()
        res <- gextract.variants(vt, gintervals("chr1", b, b + 1), iterator = 1)
        outside <- setdiff(170:240, res$pos)
        new_base <- c(A = "C", C = "G", G = "T", T = "A")[toupper(strsplit(seq, "")[[1]][outside + 1])]
        edb <- make_seq_db(substitute_bases(seq, outside, new_base))
        expect_equal(values_on_db(edb, vt, b), values_on_db(db, vt, b), info = vt)
    }
    # Single-vtrack expressions never take the batched sequence path (it runs
    # for >= 4 kmer/pwm/potts/masked vtracks), so check one such expression too.
    expr <- "v_kmer + 0 * v_pwm + 0 * v_potts + 0 * v_masked"
    value_on <- function(root) {
        suppressMessages(gsetroot(root))
        define_variant_vtracks()
        gextract(expr, gintervals("chr1", b, b + 1), iterator = 1, colnames = "x")$x
    }
    suppressMessages(gsetroot(db))
    define_variant_vtracks()
    res <- gextract.variants(expr, gintervals("chr1", b, b + 1), iterator = 1)
    outside <- setdiff(170:240, res$pos)
    new_base <- c(A = "C", C = "G", G = "T", T = "A")[toupper(strsplit(seq, "")[[1]][outside + 1])]
    edb <- make_seq_db(substitute_bases(seq, outside, new_base))
    expect_equal(value_on(edb), value_on(db), info = "batched expression")
})

test_that("chromosome ends clip the reach", {
    local_db_state()
    seq <- variant_test_seq(seed = 11L)
    db <- make_seq_db(seq)
    suppressMessages(gsetroot(db))
    define_variant_vtracks()
    res <- gextract.variants("v_pwm", gintervals("chr1", c(0, 398), c(2, 400)), iterator = 1)
    expect_true(all(res$pos >= 0 & res$pos < 400))
})

test_that("N bases are not substituted, and a listed variant on an N is a ref mismatch", {
    local_db_state()
    seq <- variant_test_seq(seed = 12L)
    substr(seq, 101, 103) <- "NNN"
    db <- make_seq_db(seq)
    suppressMessages(gsetroot(db))
    define_variant_vtracks()
    res <- gextract.variants("v_pwm", gintervals("chr1", 100, 101), iterator = 1)
    expect_false(any(res$pos %in% 100:102))
    v <- data.frame(chrom = "chr1", start = 100, end = 101, ref = "A", alt = "C")
    expect_warning(res2 <- gextract.variants("v_pwm", gintervals("chr1", 100, 101), variants = v, iterator = 1), "1 of 1 variants")
    expect_true(all(is.na(res2[["v_pwm.alt"]])))
})

test_that("a listed variant pairs with every bin whose reach contains it", {
    define_variant_vtracks()
    g <- toupper(strsplit(gseq.extract(gintervals(1, 1000, 1030)), "")[[1]])
    alt_at <- function(i) setdiff(VARIANT_BASES, g[i - 1000 + 1])[1]
    v <- data.frame(
        chrom = "chr1", start = c(1005, 1020, 1029), end = c(1006, 1021, 1030),
        ref = c(g[6], g[21], ifelse(g[30] == "A", "C", "A")), alt = c(alt_at(1005), alt_at(1020), "G")
    )
    v <- v[v$ref != v$alt, ]
    expect_warning(res <- gextract.variants("v_pwm", gintervals(1, 1000, 1010), variants = v, iterator = 1), "1 of 3 variants")
    # v_pwm reads [s - 5, s + 7): 1005 is read by bins 999..1010 (clipped to 1000..1009), 1020 by none
    # the third variant has a wrong ref: counted in the warning although no bin reads it
    expect_setequal(res$start[res$pos == 1005], 1000:1009)
    expect_false(any(res$pos %in% c(1020, 1029)))
})

test_that("variants are validated", {
    define_variant_vtracks()
    iv <- gintervals(1, 1000, 1001)
    bad <- function(...) data.frame(chrom = "chr1", start = 1000, end = 1001, ...)
    expect_error(gextract.variants("v_pwm", iv, variants = bad(ref = "A", alt = "AC")), "single-base")
    expect_error(gextract.variants("v_pwm", iv, variants = bad(ref = "A", alt = "C,G")), "single-base")
    expect_error(gextract.variants("v_pwm", iv, variants = bad(ref = "A", alt = "A")), "same")
    expect_error(gextract.variants("v_pwm", iv, variants = data.frame(chrom = "chr1", start = 1000)), "columns")
})

test_that("score.thresh keeps exactly the rows that cross it", {
    define_variant_vtracks()
    iv <- gintervals(1, 2000, 2050)
    all_rows <- gextract.variants("v_pwm", iv, iterator = 1, colnames = "s")
    # halfway across the largest change, so at least one row crosses
    i <- which.max(abs(all_rows$s.alt - all_rows$s.ref))
    t <- (all_rows$s.ref[i] + all_rows$s.alt[i]) / 2
    crossed <- all_rows[(all_rows$s.ref >= t) != (all_rows$s.alt >= t), ]
    got <- gextract.variants("v_pwm", iv, iterator = 1, colnames = "s", score.thresh = t)
    rownames(crossed) <- NULL
    rownames(got) <- NULL
    expect_equal(got, crossed)
})

test_that("a non-sequence track contributes the same value to .ref and .alt", {
    define_variant_vtracks()
    iv <- gintervals(1, 3000, 3010)
    a <- gextract.variants("v_pwm", iv, iterator = 1, colnames = "s")
    b <- gextract.variants("v_pwm + test.fixedbin", iv, iterator = 1, colnames = "s")
    keep <- !is.na(b$s.ref)
    expect_gt(sum(keep), 0)
    expect_equal((b$s.alt - b$s.ref)[keep], (a$s.alt - a$s.ref)[keep], tolerance = 1e-6)
})

test_that("intervalID matches gextract, and chunking does not change the result", {
    define_variant_vtracks()
    iv <- gintervals(1, c(4000, 4003, 4100), c(4006, 4010, 4103))
    ref <- gextract("v_pwm", iv, iterator = 1)
    res <- misha:::.gextract_variants_impl("v_pwm", iv, NULL, 1, NULL, "s", chunk = 3L)
    one <- misha:::.gextract_variants_impl("v_pwm", iv, NULL, 1, NULL, "s", chunk = 1e6)
    expect_equal(res, one)
    key <- paste(ref$start, ref$intervalID)
    expect_true(all(paste(res$start, res$intervalID) %in% key))
})

test_that("no rows gives NULL, and an expression with no sequence vtrack is an error", {
    define_variant_vtracks()
    v <- data.frame(chrom = "chr1", start = 90000, end = 90001, ref = "A", alt = "C")
    expect_null(suppressWarnings(gextract.variants("v_pwm", gintervals(1, 1000, 1001), variants = v, iterator = 1)))
    expect_error(gextract.variants("test.fixedbin", gintervals(1, 1000, 1001), iterator = 1), "no sequence virtual track")
})

test_that("intervals on two chromosomes sort the way the engine needs", {
    define_variant_vtracks()
    res <- gextract.variants("v_pwm", gintervals(c(1, 2), 1000, 1003), iterator = 1)
    expect_setequal(as.character(unique(res$chrom)), c("chr1", "chr2"))
})
