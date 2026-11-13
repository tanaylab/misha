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
    expect_equal(reach("v_kmer"), cbind(lo = -2, hi = 2 + 2))
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
    # K-mers take the batched path only with 32+ k-mer vtracks on bins of
    # 100 bp or more, so check one such expression too.
    pads <- paste0("v_pad", 1:31)
    define_vtracks <- function() {
        define_variant_vtracks()
        for (p in pads) gvtrack.create(p, NULL, "kmer.count", kmer = "CGT", extend = TRUE, strand = 0)
    }
    expr <- paste(c("v_kmer", "0 * v_pwm", "0 * v_potts", "0 * v_masked", paste0("0 * ", pads)), collapse = " + ")
    bin <- gintervals("chr1", b, b + 100)
    value_on <- function(root) {
        suppressMessages(gsetroot(root))
        define_vtracks()
        gextract(expr, bin, iterator = 100, colnames = "x")$x
    }
    suppressMessages(gsetroot(db))
    define_vtracks()
    res <- gextract.variants(expr, bin, iterator = 100)
    outside <- setdiff(150:350, res$pos)
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
    expect_setequal(unique(res$pos), c(0:7, 393:399))
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
    expect_equal(nrow(res2), 1L)
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

test_that("listed variants: each (bin, variant) with the variant in the bin's reach, duplicates and other chromosomes too", {
    define_variant_vtracks()
    iv <- gintervals(c(1, 1, 2), c(1000, 1003, 1000), c(1010, 1006, 1012))
    g1 <- toupper(strsplit(gseq.extract(gintervals(1, 990, 1030)), "")[[1]])
    g2 <- toupper(strsplit(gseq.extract(gintervals(2, 990, 1030)), "")[[1]])
    ref_at <- function(g, pos) g[pos - 990 + 1]
    alt_of <- function(r) setdiff(VARIANT_BASES, r)[1]
    pos1 <- c(995, 1000, 1004, 1004, 1016, 1020)
    pos2 <- c(1004, 1008, 1029)
    v <- data.frame(
        chrom = c(rep("chr1", length(pos1)), rep("chr2", length(pos2))),
        start = c(pos1, pos2),
        ref = c(ref_at(g1, pos1), ref_at(g2, pos2))
    )
    v <- v[v$ref %in% VARIANT_BASES, ]
    v$alt <- vapply(v$ref, alt_of, "")
    v$end <- v$start + 1
    # the same variant twice, and one with a wrong ref at the same position
    v <- rbind(v, v[2, ], transform(v[2, ], ref = alt_of(v$ref[2]), alt = v$ref[2]))
    res <- suppressWarnings(gextract.variants("v_pwm", iv, variants = v, iterator = 1))

    # v_pwm reads [s - 5, s + 7) from a 1 bp bin [s, s + 1)
    bins <- gextract("v_pwm", iv, iterator = 1)
    want <- do.call(rbind, lapply(seq_len(nrow(bins)), function(i) {
        hit <- v[as.character(v$chrom) == as.character(bins$chrom[i]) &
            v$start >= bins$start[i] - 5 & v$start < bins$end[i] + 6, ]
        if (!nrow(hit)) {
            return(NULL)
        }
        data.frame(chrom = as.character(bins$chrom[i]), start = bins$start[i], pos = hit$start, ref = hit$ref, alt = hit$alt)
    }))
    got <- data.frame(chrom = as.character(res$chrom), start = res$start, pos = res$pos, ref = res$ref, alt = res$alt)
    key <- function(d) d[do.call(order, unname(as.list(d))), ]
    expect_equal(nrow(got), nrow(want))
    expect_equal(key(got), key(want), ignore_attr = TRUE)
})

test_that("listed variants join only the bin positions that hold one", {
    define_variant_vtracks()
    # The join used to take every (bin, position) a bin reads: ~10^7 rows for 10^5 variants.
    rows_joined <- integer(0)
    record <- function(n) rows_joined <<- c(rows_joined, n)
    suppressMessages(trace("merge.data.frame", tracer = bquote(.(record)(nrow(x))), where = baseenv(), print = FALSE))
    withr::defer(suppressMessages(untrace("merge.data.frame", where = baseenv())))
    g <- toupper(strsplit(gseq.extract(gintervals(1, 5000, 5001)), "")[[1]])
    v <- data.frame(chrom = "chr1", start = 5000, end = 5001, ref = g, alt = setdiff(VARIANT_BASES, g)[1])
    res <- gextract.variants("v_pwm", gintervals(1, 4000, 6000), variants = v, iterator = 1)
    # bins 4994..5005 read 5000; the 2000 bins read ~24000 positions
    expect_equal(nrow(res), 12L)
    expect_gt(length(rows_joined), 0)
    expect_lte(max(rows_joined), nrow(res))
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
    threshold <- (all_rows$s.ref[i] + all_rows$s.alt[i]) / 2
    crossed <- all_rows[(all_rows$s.ref >= threshold) != (all_rows$s.alt >= threshold), ]
    got <- gextract.variants("v_pwm", iv, iterator = 1, colnames = "s", score.thresh = threshold)
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

test_that("a reach with zero span reads nothing and gives NULL", {
    define_variant_vtracks()
    gvtrack.create("v_empty", NULL, "pwm.max", params = list(pssm = variant_test_pssm(), extend = FALSE))
    gvtrack.iterator("v_empty", sshift = 0, eshift = -1)
    expect_null(gextract.variants("v_empty", gintervals(1, 1000, 1003), iterator = 1))
})

test_that("a reference byte >= 128 reads as no base, as N does", {
    local_db_state()
    seq <- charToRaw(variant_test_seq(seed = 11L))
    seq[101] <- charToRaw("N")
    db_n <- make_seq_db(rawToChar(seq))
    db_b <- make_seq_db(rawToChar(seq))
    # gdb.create() refuses such a byte in a FASTA, so it goes into the sequence file directly
    seq_file <- file.path(db_b, "seq", "genome.seq")
    expect_equal(file.size(seq_file), length(seq))
    seq[101] <- as.raw(0xe9)
    writeBin(seq, seq_file)
    run <- function(db) {
        suppressMessages(gsetroot(db))
        define_variant_vtracks()
        gextract.variants("v_kmer", gintervals("chr1", 96, 106), iterator = 1)
    }
    res_b <- run(db_b)
    expect_equal(charToRaw(gseq.extract(gintervals("chr1", 100, 101))), as.raw(0xe9))
    expect_false(100 %in% res_b$pos)
    expect_equal(res_b, run(db_n))
})
