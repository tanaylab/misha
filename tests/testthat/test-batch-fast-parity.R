# Parity tests for the batched fast path (gscreen/gsummary/gquantiles
# fast=TRUE) against the legacy slow path, on controlled tracks that exercise
# the cases the example DB doesn't: negative values, NaN gaps, and chromosome
# sizes that aren't a multiple of the iterator. These guard the fixes for:
#   #1 SUM aggregate bounds (window_max/min * over-estimated count was not a
#      valid bound -> missed/spurious gscreen intervals, wrong gquantiles)
#   #2 all-NaN MIN/MAX certain-pass (an all-NaN window certain-passed instead
#      of breaking the run -> fused intervals / phantom hits)
#   #5 gscreen run end past chrom size when chrom size % iterator != 0
#   #3 bare sparse track is not fast-path eligible without an explicit iterator
#   #4 only grid-aligned scopes are fast-path eligible

# Build a controlled genome with a signed track containing NaN gaps and a
# fully-covered positive track.
.setup_parity_db <- function(chrom_size = 600L, binsize = 10L) {
    root <- tempfile("batchparity_")
    dir.create(root)
    dir.create(file.path(root, "tracks"))
    dir.create(file.path(root, "seq"))
    writeLines(sprintf("chr1\t%d", chrom_size), file.path(root, "chrom_sizes.txt"))
    file.create(file.path(root, "seq", "chr1.seq"))
    gsetroot(root)

    nbins <- chrom_size %/% binsize
    # Signed track "t": data in bins 1:20 and 31:50, NaN gaps elsewhere.
    vals <- rep(NA_real_, nbins)
    vals[1:20] <- rep(c(2, -2, 1, -1), length.out = 20)
    vals[31:50] <- rep(c(3, -3), length.out = 20)
    idx <- which(!is.na(vals))
    gtrack.create_dense(
        "t", "signed w/ NaN gaps",
        intervals = data.frame(chrom = "chr1", start = (idx - 1) * binsize, end = idx * binsize),
        values = vals[idx], binsize = binsize, defval = NaN, func = "weighted.mean"
    )

    # Fully-covered positive track "tfull": every bin = 0.5.
    gtrack.create_dense(
        "tfull", "all 0.5",
        intervals = data.frame(chrom = "chr1", start = (seq_len(nbins) - 1) * binsize, end = seq_len(nbins) * binsize),
        values = rep(0.5, nbins), binsize = binsize, defval = NaN, func = "weighted.mean"
    )
    invisible(root)
}

# Normalize a gscreen result for comparison (chrom -> character, sort).
.norm_iv <- function(d) {
    if (is.null(d) || nrow(d) == 0) {
        return(data.frame(chrom = character(0), start = numeric(0), end = numeric(0)))
    }
    out <- data.frame(
        chrom = as.character(d$chrom),
        start = as.numeric(d$start),
        end = as.numeric(d$end),
        stringsAsFactors = FALSE
    )
    out <- out[order(out$chrom, out$start, out$end), ]
    rownames(out) <- NULL
    out
}

# Compare a single screen expression: legacy (single-expr) vs fast (multi-expr,
# filtered back to that expression). The filler expression reuses the SAME lhs
# (track/vtrack) so the multi-expr set shares one (func, sshift, eshift) tuple
# and stays fast-path eligible.
.screen_parity <- function(expr, iterator) {
    lhs <- sub("\\s*(<=|>=|==|<|>).*$", "", expr)
    extra <- paste0(lhs, " < 1e12")
    leg <- gscreen(expr, iterator = iterator)
    fast <- gscreen(c(expr, extra), iterator = iterator, fast = TRUE)
    f1 <- fast[fast$track == expr, , drop = FALSE]
    list(leg = .norm_iv(leg), fast = .norm_iv(f1))
}

test_that("SUM screen fast path matches legacy on a signed track (#1)", {
    old <- options(gmultitasking = FALSE, misha.quiet_dispatch = TRUE)
    on.exit(options(old), add = TRUE)
    .setup_parity_db()
    remove_all_vtracks()
    gvtrack.create("vsum", "t", "sum")
    gvtrack.iterator("vsum", sshift = 0, eshift = 0)

    # The old SUM bounds (window_min*count) over-pruned LT and mis-handled GT
    # on negative data. Sweep both sides and several thresholds.
    for (thr in c(-5, -3, -1, 0, 1, 3)) {
        for (op in c("<", ">")) {
            expr <- sprintf("vsum %s %g", op, thr)
            r <- .screen_parity(expr, iterator = 10)
            expect_equal(r$fast, r$leg, info = expr)
        }
    }
})

test_that("MIN/MAX screen fast path does not fuse runs across NaN gaps (#2)", {
    old <- options(gmultitasking = FALSE, misha.quiet_dispatch = TRUE)
    on.exit(options(old), add = TRUE)
    .setup_parity_db()
    remove_all_vtracks()
    gvtrack.create("vmin", "t", "min")
    gvtrack.iterator("vmin", sshift = 0, eshift = 0)
    gvtrack.create("vmax", "t", "max")
    gvtrack.iterator("vmax", sshift = 0, eshift = 0)

    # vmin > t : the non-informative MIN upper bound (+inf) used to certain-pass
    # an all-NaN window, fusing the two data runs into one. Must stay split.
    r <- .screen_parity("vmin > -10", iterator = 10)
    expect_equal(r$fast, r$leg)
    expect_equal(nrow(r$fast), 2L) # two data regions, NaN gap between them

    # Symmetric MAX < t (MAX lower bound is -inf).
    r2 <- .screen_parity("vmax < 10", iterator = 10)
    expect_equal(r2$fast, r2$leg)
    expect_equal(nrow(r2$fast), 2L)
})

test_that("gscreen fast path clamps run end to chrom size (#5)", {
    old <- options(gmultitasking = FALSE, misha.quiet_dispatch = TRUE)
    on.exit(options(old), add = TRUE)
    .setup_parity_db(chrom_size = 600L) # 600 is not a multiple of iterator 70
    remove_all_vtracks()

    r <- .screen_parity("tfull > 0", iterator = 70L)
    expect_equal(r$fast, r$leg)
    # No emitted interval may exceed the chromosome size.
    expect_true(all(r$fast$end <= 600))
})

test_that("fast-path eligibility declines sparse-without-iterator and non-aligned scopes (#3,#4)", {
    old <- options(gmultitasking = FALSE, misha.quiet_dispatch = TRUE)
    on.exit(options(old), add = TRUE)
    root <- .setup_parity_db()
    # A sparse track.
    gtrack.create_sparse(
        "sp", "sparse",
        intervals = data.frame(chrom = "chr1", start = c(50L, 150L, 350L), end = c(60L, 160L, 360L)),
        values = c(1, 2, 3)
    )

    # #3: bare sparse track with iterator = NULL is NOT eligible (its implicit
    # iterator is the irregular sparse iterator, not a grid).
    expect_null(misha:::.detect_fast_path("sp", NULL, NULL, NULL))
    # ... but eligible with an explicit iterator.
    expect_false(is.null(misha:::.detect_fast_path("sp", 10L, NULL, NULL)))

    # Dense bare track with iterator = NULL stays eligible (bin grid).
    expect_false(is.null(misha:::.detect_fast_path("t", NULL, NULL, NULL)))

    # #4: a non-grid-aligned explicit scope is NOT eligible.
    aligned <- data.frame(chrom = "chr1", start = 100L, end = 200L) # multiples of 10
    nonaligned <- data.frame(chrom = "chr1", start = 25L, end = 125L) # not multiples of 10
    expect_false(is.null(misha:::.detect_fast_path("t", 10L, aligned, NULL)))
    expect_null(misha:::.detect_fast_path("t", 10L, nonaligned, NULL))

    # A scope ending exactly at chrom size is allowed even if not a multiple.
    full_chrom <- data.frame(chrom = "chr1", start = 0L, end = 600L)
    expect_false(is.null(misha:::.detect_fast_path("t", 70L, full_chrom, NULL)))
})
