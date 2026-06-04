# Regression for the N-count skip false-negative in PWMEditDistanceScorer
# (introduced by commit cb1f287a, "perf: sliding-window N-count skip for edit
# distance").
#
# evaluate_windows skips any window whose N-count exceeds max_edits, on the
# rationale that "each N forces a mandatory edit". That holds only for
# direction="above" with substitutions only. For direction="below" an N is NOT
# a mandatory edit (precompute_tables sets m_mandatory_table[i][N]=false and
# scores N at col_max), so an N-heavy window can already be at/below threshold
# (0 edits). The unconditional skip wrongly returned NA / missed the true
# minimum for below queries with max_edits set. (For indels the fixed-width
# N window is also not a sound bound.) The fix gates the skip to
# !is_below() && max_indels == 0.
#
# This needs sequence containing real N bases, so we build a tiny genome whose
# .seq is all N (misha .seq files are raw bytes, position == coordinate).

make_all_n_genome <- function(root, chrom_size = 120L) {
    dir.create(root)
    dir.create(file.path(root, "tracks"))
    dir.create(file.path(root, "seq"))
    writeLines(sprintf("chr1\t%d", chrom_size), file.path(root, "chrom_sizes.txt"))
    # Raw-byte .seq: every base is N.
    writeBin(
        charToRaw(paste(rep("N", chrom_size), collapse = "")),
        file.path(root, "seq", "chr1.seq")
    )
    gsetroot(root)
    invisible(root)
}

test_that("pwm.edit_distance below: N-heavy window already below threshold returns 0, not NA (max_edits set)", {
    root <- tempfile("nskip_below_")
    make_all_n_genome(root, chrom_size = 120L)
    tryCatch(remove_all_vtracks(), error = function(e) NULL)

    # 4-position motif: more positions than max_edits=1, so an all-N window has
    # n_count = 4 > max_edits = 1 and the buggy skip fires.
    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(1, 10, 50)

    # Threshold extremely high: any window already scores below it -> 0 edits.
    # True regardless of how N is scored (col_max <= 0 < 100).
    threshold <- 100.0

    gvtrack.create("edist_n_maxedits", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf, direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0,
        max_edits = 1
    )
    res_maxedits <- gextract("edist_n_maxedits", test_intervals, iterator = test_intervals)

    # The N window is already below the threshold -> 0 edits. Pre-fix the window
    # was skipped (n_count=4 > max_edits=1) and the result was NA.
    expect_false(is.na(res_maxedits$edist_n_maxedits[1]))
    expect_equal(res_maxedits$edist_n_maxedits[1], 0, tolerance = 1e-6)

    # Cross-check: with no max_edits the skip never fires (m_max_edits = -1), so
    # this path was always correct. The max_edits path must now agree with it.
    gvtrack.create("edist_n_nomax", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf, direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    res_nomax <- gextract("edist_n_nomax", test_intervals, iterator = test_intervals)
    expect_equal(res_maxedits$edist_n_maxedits[1], res_nomax$edist_n_nomax[1], tolerance = 1e-6)
})

test_that("pwm.max.edit_distance below over an N region is not poisoned by the N-count skip", {
    root <- tempfile("nskip_below2_")
    make_all_n_genome(root, chrom_size = 120L)
    tryCatch(remove_all_vtracks(), error = function(e) NULL)

    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(1, 10, 60)
    threshold <- 100.0

    gvtrack.create("medist_max", NULL,
        func = "pwm.max.edit_distance",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf, direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0,
        max_edits = 1
    )
    gvtrack.create("medist_nomax", NULL,
        func = "pwm.max.edit_distance",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf, direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    res <- gextract("medist_max", "medist_nomax", test_intervals, iterator = test_intervals)

    # With max_edits set, the result must match the (always-correct) no-max path
    # rather than collapsing to NA because every N window got skipped.
    expect_equal(res$medist_max[1], res$medist_nomax[1], tolerance = 1e-6)
})
