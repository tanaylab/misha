# Regression tests for the medium-severity 2026-06-28 audit fixes.
# See dev/notes/2026-06-28_full-audit.md (M1-M11).
# M1/M3/M4 (C++ OOB/defensive) are exercised by the liftover/autocorr/gextract
# suites. M5 (gtrack.convert registration) has no runnable trigger here - it
# needs a genuinely old-format track, and every fixture reports "does not require
# conversion"; it is verified by code review.

create_isolated_test_db()

test_that("M2: gintervals.canonic 2D mapping uses the original->result convention", {
    # 1D writes mapping[original]=result (documented); 2D wrote the inverse.
    df <- data.frame(
        chrom1 = c("chr3", "chr1", "chr2"), start1 = 0, end1 = 1000,
        chrom2 = c("chr3", "chr1", "chr2"), start2 = 0, end2 = 1000,
        stringsAsFactors = FALSE
    )
    res <- gintervals.canonic(df)
    m <- attr(res, "mapping")
    # sorted order chr1<chr2<chr3 => a 3-cycle distinguishes the correct mapping
    # (3,1,2) from the inverted one (2,3,1).
    expect_equal(as.numeric(m), c(3, 1, 2))
    for (i in seq_len(nrow(df))) {
        expect_equal(as.character(df$chrom1[i]), as.character(res$chrom1[m[i]]))
    }
})

test_that("M7: gtrack.ls(db = GROOT) returns tracks with a single database", {
    groot <- get("GROOT", envir = .misha)
    all_tracks <- gtrack.ls()
    by_db <- gtrack.ls(db = groot)
    expect_false(is.null(by_db))
    expect_setequal(by_db, all_tracks)
})

test_that("M8: gtrack.var.ls accepts a call-expression track argument", {
    gtrack.create_sparse("m8trk", "t", gintervals(1, 0, 10000), 1)
    withr::defer(gtrack.rm("m8trk", force = TRUE))
    gtrack.var.set("m8trk", "m8var", 42)
    tracks <- c("m8trk")

    # a call expression (tracks[1]) used to fail length(substitute(track)) != 1
    expect_true("m8var" %in% gtrack.var.ls(tracks[1]))
})

test_that("M9: gintervals.annotate fills na_value for per-row no-neighbor cases", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")
    r <- gintervals.annotate(intervs, ann,
        annotation_columns = "remark",
        mindist = -80, maxdist = 80, na_value = "none"
    )
    # the 5000-5050 query has no neighbor within +/-80 -> must be "none", not NA
    expect_true(all(!is.na(r$remark)))
    expect_equal(r$remark[r$start == 5000], "none")
})

test_that("M10: max-data-size error formats the double option (no sprintf '%d' failure)", {
    # gmax.data.size is a (multi-GB) double; a fractional value reproduces the
    # sprintf("%d", ...) failure without needing a billion-row result.
    withr::local_options(list(gmax.data.size = 100.5))
    expect_error(
        misha:::.gverify_max_data_size(200),
        "size exceeded the maximal allowed"
    )
})

test_that("M11: gsetroot with a missing chrom_sizes.txt leaves the current DB loaded", {
    cur_groot <- get("GROOT", envir = .misha)
    bad <- tempfile(pattern = "baddb_")
    dir.create(file.path(bad, "tracks"), recursive = TRUE)
    dir.create(file.path(bad, "seq"))
    withr::defer(unlink(bad, recursive = TRUE))

    expect_error(gsetroot(bad), "chrom_sizes")
    # current session must be intact (baseline wiped it before the check)
    expect_identical(get("GROOT", envir = .misha), cur_groot)
    expect_false(is.null(get("ALLGENOME", envir = .misha)))
})
