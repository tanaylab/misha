# Tests for gscreen fast path (Phase 4) on gdb.init_examples().

test_that("gscreen single-expr legacy behavior preserved", {
    gdb.init_examples()
    r <- gscreen("dense_track > 0.15")
    expect_s3_class(r, "data.frame")
    expect_true(all(c("chrom", "start", "end") %in% colnames(r)))
})

test_that("gscreen with character vector returns long data.frame", {
    gdb.init_examples()
    res <- gscreen(c("dense_track > 0.15", "dense_track < 0.05"),
        fast = TRUE
    )
    expect_s3_class(res, "data.frame")
    expect_equal(colnames(res), c("chrom", "start", "end", "track"))
    expect_true(all(res$track %in% c(
        "dense_track > 0.15",
        "dense_track < 0.05"
    )))
})

test_that("gscreen multi-track fast path matches per-track legacy calls", {
    gdb.init_examples()
    exprs <- c("dense_track > 0.15", "subdir.dense_track2 > 0.3")
    fast <- gscreen(exprs, fast = TRUE)

    ref <- do.call(rbind, lapply(exprs, function(e) {
        iv <- gscreen(e)
        if (nrow(iv) == 0) {
            data.frame(
                chrom = character(0), start = numeric(0),
                end = numeric(0), track = character(0),
                stringsAsFactors = FALSE
            )
        } else {
            iv$track <- e
            iv[, c("chrom", "start", "end", "track")]
        }
    }))

    # Sort both the same way (fast path returns per-track, then per-chrom;
    # legacy returns per-chrom).
    fast$chrom <- as.character(fast$chrom)
    ref$chrom <- as.character(ref$chrom)
    ord <- function(df) df[order(df$track, df$chrom, df$start, df$end), ]
    fast_s <- ord(fast)
    rownames(fast_s) <- NULL
    ref_s <- ord(ref)
    rownames(ref_s) <- NULL
    expect_equal(nrow(fast_s), nrow(ref_s))
    expect_equal(fast_s$chrom, ref_s$chrom)
    expect_equal(fast_s$start, ref_s$start, tolerance = 0)
    expect_equal(fast_s$end, ref_s$end, tolerance = 0)
    expect_equal(fast_s$track, ref_s$track)
})

test_that("gscreen supports all five comparison operators", {
    gdb.init_examples()
    for (op in c("<", "<=", "==", ">=", ">")) {
        # Use a threshold that's likely to pass or not trivially, and a
        # second expression so the fast path is exercised.
        exprs <- c(
            paste("dense_track", op, "0.1"),
            paste("subdir.dense_track2", op, "0.2")
        )
        df <- tryCatch(gscreen(exprs, fast = TRUE),
            error = function(e) NULL
        )
        expect_false(is.null(df), info = op)
        expect_s3_class(df, "data.frame")
        expect_equal(colnames(df), c("chrom", "start", "end", "track"))
    }
})

test_that("gscreen flushes at interval-mask boundary (no spurious fusion)", {
    gdb.init_examples()
    # Two non-adjacent intervals on chrom 1. The boundary() hook must
    # prevent a passing run on the first from fusing with a passing run
    # on the second.
    ivs <- data.frame(
        chrom = c("1", "1"),
        start = c(0, 200000),
        end = c(100000, 300000),
        stringsAsFactors = FALSE
    )
    res <- gscreen(c("dense_track > 0.15", "subdir.dense_track2 > 0.3"),
        intervals = ivs, fast = TRUE
    )
    # No returned interval may straddle [100000, 200000).
    bad <- with(res, as.character(chrom) == "1" &
        start < 100000 & end > 200000)
    expect_false(any(bad))
})

test_that("gscreen handles multiple predicates on the SAME underlying track", {
    # Regression for Phase 4 review C1: when two or more expressions share
    # the same underlying source track, the R-side track-column rewrite
    # previously collapsed the per-row mapping because it relied on
    # the underlying-track name string and the rows were indistinguishable.
    gdb.init_examples()
    exprs <- c("dense_track > 0.05", "dense_track > 0.15")
    fast <- gscreen(exprs, fast = TRUE)
    # Both input expressions must appear in the track column.
    expect_true(all(c("dense_track > 0.05", "dense_track > 0.15") %in%
        unique(as.character(fast$track))))
    # The stricter threshold produces a subset of the looser one; every
    # row tagged as the stricter predicate must also satisfy the looser.
    ref_looser <- gscreen("dense_track > 0.05")
    ref_stricter <- gscreen("dense_track > 0.15")
    expect_equal(sum(fast$track == "dense_track > 0.05"), nrow(ref_looser))
    expect_equal(sum(fast$track == "dense_track > 0.15"), nrow(ref_stricter))
})

test_that("gscreen multi-expr with fast=FALSE errors cleanly (Phase 5 scope)", {
    gdb.init_examples()
    expect_error(
        gscreen(c("dense_track > 0.15", "dense_track < 0.05"), fast = FALSE),
        regexp = "multi-expression slow path not yet implemented"
    )
})

test_that("gscreen with intervals.set.out + multi-expr rejects fast path", {
    gdb.init_examples()
    # intervals.set.out is incompatible with multi-expr fast path; must
    # fall through (and then error from the slow-path stub).
    expect_error(
        gscreen(c("dense_track > 0.15", "dense_track < 0.05"),
            intervals.set.out = "foo", fast = TRUE
        ),
        regexp = "multi-expression slow path not yet implemented"
    )
})
