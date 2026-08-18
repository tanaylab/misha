# gvtrack.iterator() sshift/eshift can push the iterator interval off the
# chromosome. Partially: the window is clamped to the chromosome. Completely
# (or into an empty window): every virtual track function returns NaN.

shift_pssm <- function() {
    data.frame(
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
}

# Create every vtrack flavor that has its own out-of-range path, and return the
# extracted values for the given shift.
shifted_values <- function(sshift, eshift, scope) {
    remove_all_vtracks()
    src <- gintervals("chr1", c(3000100, 3000500), c(3000200, 3000600))

    gvtrack.create("s_pwm", NULL, func = "pwm", pssm = shift_pssm(), prior = 0.01)
    gvtrack.create("s_pwm_max", NULL, func = "pwm.max", pssm = shift_pssm(), prior = 0.01)
    gvtrack.create("s_kmer", NULL, func = "kmer.count", kmer = "ACGT", strand = 1)
    gvtrack.create("s_masked", NULL, func = "masked.count")
    gvtrack.create("s_coverage", src, func = "coverage")
    gvtrack.create("s_neighbor", src, func = "neighbor.count")
    gvtrack.create("s_distance", src, func = "distance")
    gvtrack.create("s_track", "test.fixedbin", "avg")
    vtracks <- c(
        "s_pwm", "s_pwm_max", "s_kmer", "s_masked",
        "s_coverage", "s_neighbor", "s_distance", "s_track"
    )

    for (v in vtracks) {
        gvtrack.iterator(v, sshift = sshift, eshift = eshift)
    }

    # 4+ sequence vtracks in one call also exercises the batched sequence path
    r <- gextract(vtracks, scope, iterator = scope)
    unlist(r[1, vtracks])
}

test_that("a shift that leaves the chromosome yields NaN for every vtrack function", {
    withr::defer(remove_all_vtracks())
    scope <- gintervals("chr1", 3000000, 3001000)

    baseline <- shifted_values(0, 0, scope)
    expect_false(any(is.na(baseline)))

    # window pushed 5Mb to the left of the chromosome start
    expect_true(all(is.na(shifted_values(-5e6, -5e6, scope))))

    # window pushed past the chromosome end
    chr1_end <- gintervals.all()$end[gintervals.all()$chrom == "chr1"]
    tail_scope <- gintervals("chr1", chr1_end - 1000, chr1_end)
    expect_true(all(is.na(shifted_values(10000, 10000, tail_scope))))

    # sshift/eshift that collapse the window to nothing
    expect_true(all(is.na(shifted_values(500, -500, scope))))
})

test_that("a shift that only partly leaves the chromosome clamps the window", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())
    scope <- gintervals("chr1", c(0, 1000), c(1000, 2000))

    gvtrack.create("clamped", "test.fixedbin", "avg")
    unshifted <- gextract("clamped", scope, iterator = scope)$clamped

    gvtrack.iterator("clamped", sshift = -1000, eshift = 0)
    shifted <- gextract("clamped", scope, iterator = scope)$clamped

    # the first bin has nothing to its left, so it keeps its unshifted value;
    # the second one widens to [0, 2000) and changes
    expect_equal(shifted[1], unshifted[1])
    expect_false(isTRUE(all.equal(shifted[2], unshifted[2])))
    expect_false(any(is.na(shifted)))
})
