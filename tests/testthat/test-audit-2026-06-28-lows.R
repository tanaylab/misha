# Regression tests for the low-severity 2026-06-28 audit fixes.
# See dev/notes/2026-06-28_full-audit.md (L1-L12). L1/L3/L9/L10 have no runnable
# trigger here (benign PROTECT, corrupt-file hardening, dead code, numeric-consistency
# deferral); L4 (indexed 2D finalizer) is exercised by test-gintervals-2d-indexed.R;
# L6 (directional dist==0) and L11 (already fixed with H6) are not retested here.

create_isolated_test_db()

test_that("L5: non-ASCII track attribute values round-trip (signed-char getc)", {
    gtrack.create_sparse("lowtrk", "t", gintervals(1, 0, 10000), 1)
    withr::defer(gtrack.rm("lowtrk", force = TRUE))
    value <- "café über" # bytes >= 0x80
    gtrack.attr.set("lowtrk", "note", value)
    expect_identical(gtrack.attr.get("lowtrk", "note"), value)
})

test_that("L12: gtrack.array.extract rejects NA slice indices", {
    arr <- grep("array", gtrack.ls(), value = TRUE)
    skip_if(length(arr) == 0, "no array track fixture")
    expect_error(
        gtrack.array.extract(arr[1], c(1, NA), gintervals(1, 0, 10000)),
        "must not be NA"
    )
})

test_that("L7: gtrack.2d.create writes a sensible created.by (not a deparsed value)", {
    ivs2d <- gintervals.2d(1, 0, 1000, 1, 0, 1000)
    gtrack.2d.create("low2d", "t", ivs2d, 5)
    withr::defer(gtrack.rm("low2d", force = TRUE))
    cb <- gtrack.attr.get("low2d", "created.by")
    expect_true(startsWith(cb, "gtrack.2d.create("))
    expect_false(grepl("structure(", cb, fixed = TRUE)) # the deparse-after-reassign bug
})
