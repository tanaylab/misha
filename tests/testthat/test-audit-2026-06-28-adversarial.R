# Adversarial regression tests for the 2026-06-28 audit fixes: each one attacks the
# NORMAL / boundary path a fix touched, cross-checking against an independent ground
# truth (serial vs parallel, file vs in-memory, bigset vs in-memory, consensus match),
# to catch a fix that over-corrected and broke the common case.

create_isolated_test_db()

test_that("ADV H7: track-parallel gextract == serial output (columns, order, values)", {
    exprs <- c(
        "test.fixedbin", "test.fixedbin+1", "test.fixedbin*2", "test.fixedbin-0.5",
        "test.fixedbin+10", "test.fixedbin+20", "test.fixedbin+30", "test.fixedbin/2"
    )
    scope <- gintervals(1, 0, 100000)
    for (nt in c(1L, 3L, 5L, 8L)) {
        sub <- exprs[seq_len(nt)]
        serial <- withr::with_options(
            list(gmultitasking = FALSE),
            gextract(sub, scope)
        )
        par <- withr::with_options(
            list(gmultitasking = TRUE, gmultitasking.strategy = "tracks", gmax.processes = 3),
            gextract(sub, scope)
        )
        rownames(serial) <- NULL
        rownames(par) <- NULL
        expect_identical(names(par), names(serial), info = sprintf("nt=%d names", nt))
        expect_equal(par, serial, ignore_attr = TRUE, info = sprintf("nt=%d values", nt))
    }
})

test_that("ADV H3: gextract(file=) round-trips to the in-memory result (incl. NaN)", {
    scope <- gintervals(1, 0, 50000)
    inmem <- gextract("test.fixedbin", "test.sparse", scope, iterator = 50)
    f <- tempfile()
    withr::defer(unlink(f))
    gextract("test.fixedbin", "test.sparse", scope, iterator = 50, file = f)
    ff <- read.table(f,
        header = TRUE, sep = "\t", na.strings = "nan",
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_equal(nrow(ff), nrow(inmem))
    expect_equal(ff$chrom, as.character(inmem$chrom))
    expect_equal(ff$start, inmem$start)
    expect_equal(ff$end, inmem$end)
    expect_equal(ff[["test.fixedbin"]], inmem[["test.fixedbin"]], tolerance = 1e-10)
    # NaN written as "nan" -> NA; positions must align with the in-memory NaN
    expect_equal(is.na(ff[["test.sparse"]]), is.nan(inmem[["test.sparse"]]))
    fin <- !is.nan(inmem[["test.sparse"]])
    expect_equal(ff[["test.sparse"]][fin], inmem[["test.sparse"]][fin], tolerance = 1e-10)
})

test_that("ADV H2: quantiles bigset == in-memory on the streaming path (1 and many percentiles)", {
    withr::local_options(list(gmultitasking = FALSE, gbig.intervals.size = 10))
    scope <- gintervals(c(1, 2), 0, 50000)
    for (pct in list(0.5, c(0.1, 0.5, 0.9), c(0.25, 0.5, 0.75, 0.95))) {
        tmp <- paste0("test.aq_", sample(1:1e9, 1))
        gintervals.rm(tmp, force = TRUE)
        gintervals.quantiles("test.fixedbin", pct, scope, intervals.set.out = tmp)
        from_disk <- gintervals.load(tmp)
        in_mem <- gintervals.quantiles("test.fixedbin", pct, scope)
        expect_equal(
            from_disk[order(from_disk$chrom, from_disk$start), ],
            in_mem[order(in_mem$chrom, in_mem$start), ],
            ignore_attr = TRUE, info = paste(pct, collapse = ",")
        )
        gintervals.rm(tmp, force = TRUE)
    }
})

test_that("ADV H4: pwm.edit_distance scores a perfect consensus as 0 at L=64/65/70 (heap-table boundary)", {
    remove_all_vtracks()
    consensus_pssm <- function(seqstr) {
        L <- nchar(seqstr)
        m <- matrix(0.02, L, 4)
        colnames(m) <- c("A", "C", "G", "T")
        for (i in seq_len(L)) m[i, substr(seqstr, i, i)] <- 0.94
        m / rowSums(m)
    }
    base <- 1000000 # clean (no-N) region
    for (L in c(64L, 65L, 70L)) {
        sq <- toupper(gseq.extract(gintervals(1, base, base + L)))
        skip_if(grepl("N", sq), "region has N bases")
        vt <- paste0("ed", L)
        gvtrack.create(vt, NULL, "pwm.edit_distance",
            pssm = consensus_pssm(sq), score.thresh = -40, score.min = -300,
            score.max = -5, max_edits = 5L, max_indels = 0L, bidirect = TRUE
        )
        # window starting at `base` is the exact consensus -> 0 edits needed
        perfect <- gextract(vt, gintervals(1, base, base + 1), iterator = gintervals(1, base, base + 1))[[vt]]
        # a far region does not reach the threshold within 5 edits -> NaN
        poor <- gextract(vt, gintervals(1, base + 5000, base + 5001), iterator = gintervals(1, base + 5000, base + 5001))[[vt]]
        expect_equal(perfect, 0, info = sprintf("L=%d perfect", L))
        expect_true(is.na(poor) || poor > 0, info = sprintf("L=%d poor", L))
    }
})

test_that("ADV H6: bedGraph export round-trips coordinates and values", {
    track <- paste0("test.exp_", sample(1:1e9, 1))
    gtrack.rm(track, force = TRUE)
    withr::defer(gtrack.rm(track, force = TRUE))
    ivs <- gintervals(c(1, 1, 1, 2), c(1e8, 2e8, 12345, 5000), c(1e8 + 100, 2e8 + 100, 12445, 5100))
    vals <- c(1.5, 2.25, 3.125, 4.0625)
    gtrack.create_sparse(track, "exp", ivs, vals)
    out <- tempfile(fileext = ".bedgraph")
    withr::defer(unlink(out))
    gtrack.export_bedgraph(track, out)
    bg <- read.table(out, skip = 1, sep = "\t", col.names = c("chrom", "start", "end", "value"), stringsAsFactors = FALSE)
    inmem <- gextract(track, gintervals.all())
    inmem <- inmem[!is.nan(inmem[[track]]), ]
    ord <- order(match(inmem$chrom, gintervals.all()$chrom), inmem$start)
    inmem <- inmem[ord, ]
    expect_false(any(grepl("e", c(bg$start, bg$end), ignore.case = TRUE))) # no scientific notation
    expect_equal(bg$start, inmem$start)
    expect_equal(bg$end, inmem$end)
    expect_equal(bg$value, inmem[[track]], tolerance = 1e-12)
})

test_that("ADV H8: a track create at the database root marks the cache clean (not dirty)", {
    misha:::.gdb.cache_mark_dirty(get("GROOT", envir = .misha))
    expect_true(misha:::.gdb.cache_is_dirty(get("GROOT", envir = .misha)))
    gtrack.create_sparse("rootclean", "r", gintervals(1, 0, 10000), 1)
    withr::defer(gtrack.rm("rootclean", force = TRUE))
    # at the root, the fix must still WRITE the cache and clear dirty (only subdirs mark dirty)
    expect_false(misha:::.gdb.cache_is_dirty(get("GROOT", envir = .misha)))
})

test_that("ADV M9: gintervals.annotate keeps real annotations for found neighbors", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- gintervals(1, c(1050, 5400), c(1150, 5500))
    ann$remark <- c("hit", "far")
    r <- gintervals.annotate(intervs, ann,
        annotation_columns = "remark",
        mindist = -200, maxdist = 200, na_value = "none"
    )
    # 1000-1100 has neighbor "hit" within range -> must keep "hit", not "none"
    expect_equal(r$remark[r$start == 1000], "hit")
    # 5000-5050 has no neighbor within 200 -> "none"
    expect_equal(r$remark[r$start == 5000], "none")
})

test_that("ADV M7: gtrack.ls(db=) returns all tracks for GROOT and NULL for a bogus db", {
    expect_setequal(gtrack.ls(db = get("GROOT", envir = .misha)), gtrack.ls())
    expect_null(gtrack.ls(db = tempfile(pattern = "no_such_db_")))
})

test_that("ADV M2: gintervals.canonic 1D mapping and 2D identity are unchanged", {
    # 1D (untouched by the fix): documented example "2 1 2 2"
    d1 <- data.frame(
        chrom = "chr1", start = c(11000, 100, 10000, 10500),
        end = c(12000, 200, 13000, 10600), stringsAsFactors = FALSE
    )
    res1 <- gintervals.canonic(d1)
    m1 <- as.numeric(attr(res1, "mapping"))
    expect_equal(m1, c(2, 1, 2, 2))

    # 2D already sorted -> identity mapping
    d2 <- gintervals.2d(c(1, 1), c(0, 2000), c(1000, 3000), c(1, 1), c(0, 0), c(1000, 1000))
    res2 <- gintervals.canonic(d2)
    expect_equal(as.numeric(attr(res2, "mapping")), seq_len(nrow(d2)))
})

test_that("ADV L5: multiple ASCII and non-ASCII track attributes all round-trip", {
    gtrack.create_sparse("attrtrk", "t", gintervals(1, 0, 10000), 1)
    withr::defer(gtrack.rm("attrtrk", force = TRUE))
    vals <- list(plain = "simple ascii", num = "12345", accent = "café über", mix = "a\tb-c 50%")
    for (nm in names(vals)) gtrack.attr.set("attrtrk", nm, vals[[nm]])
    for (nm in names(vals)) expect_identical(gtrack.attr.get("attrtrk", nm), vals[[nm]])
})
