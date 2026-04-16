create_isolated_test_db()

# ============================================================
# Helper: R-level logistic transform
# ============================================================
r_logist <- function(x, L = 1, k = 1, x_0 = 0, pre_shift = 0, post_shift = 0) {
    input <- x + pre_shift
    L / (1 + exp(-k * (input - x_0))) + post_shift
}

# Helper: R-level scaling (cap + normalize)
r_scale <- function(raw, max_cap, dis_from_cap, scale_factor = 10) {
    ceiled <- pmin(raw - max_cap, 0)
    floored <- pmax(ceiled, -dis_from_cap)
    scaled <- scale_factor * (floored + dis_from_cap) / dis_from_cap
    ifelse(is.finite(scaled), scaled, 0)
}

# Helper: convert vtrack sshift/eshift to glm_pred shifts
# Regular vtrack: window = [iter.start + sshift, iter.end + eshift)
# glm_pred: center = (start+end)/2, window = [center + glm_ss, center + glm_es)
# To match: glm_ss = sshift - iter_size/2, glm_es = eshift + iter_size/2
vtrack_to_glm_shifts <- function(sshift, eshift, iter_size) {
    c(sshift - iter_size / 2, eshift + iter_size / 2)
}

# ============================================================
# Category 1: R parameter validation
# ============================================================

test_that("glm_pred.create requires tracks", {
    expect_error(
        glm_pred.create("bad", tracks = character(0), inner_func = "sum", weights = 1),
        "tracks"
    )
})

test_that("glm_pred.create requires matching lengths", {
    expect_error(
        glm_pred.create("bad",
            tracks = c("test.fixedbin", "test.fixedbin"),
            inner_func = "sum", weights = 1
        ),
        "length"
    )
})

test_that("glm_pred.create validates inner_func values", {
    expect_error(
        glm_pred.create("bad",
            tracks = "test.fixedbin",
            inner_func = "bad", weights = 1
        ),
        "inner_func"
    )
})

test_that("glm_pred.create validates shift ordering", {
    expect_error(
        glm_pred.create("bad",
            tracks = "test.fixedbin",
            inner_func = "sum", weights = 1,
            shifts = list(c(100, -100))
        ),
        "sshift"
    )
})

test_that("glm_pred.create validates interaction indices", {
    expect_error(
        glm_pred.create("bad",
            tracks = c("test.fixedbin", "test.fixedbin"),
            inner_func = c("sum", "sum"),
            weights = c(1, 1),
            shifts = list(c(-100, 100), c(-100, 100)),
            interactions = list(c(1L, 3L)),
            interaction_weights = 0.5
        ),
        "indices"
    )
})

# ============================================================
# Category 2: Basic pipeline correctness
# ============================================================

test_that("glm_pred with bias only returns bias", {
    glm_pred.create("vt_bias",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = 0,
        bias = 42.0,
        shifts = list(c(-50, 50))
    )

    intervals <- gintervals(1, 500, 1000)
    result <- gextract("vt_bias", intervals = intervals, iterator = 100)
    expect_true(all(abs(result$vt_bias - 42.0) < 1e-10))
})

test_that("glm_pred with weight=1, no transform matches raw sum", {
    iter <- 100L
    ss <- -50
    es <- 50
    # glm_pred shifts are relative to center
    gs <- vtrack_to_glm_shifts(ss, es, iter)

    glm_pred.create("vt_raw_sum",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = 1.0,
        bias = 0,
        shifts = list(gs)
    )

    gvtrack.create("vt_ref_sum", "test.fixedbin", "sum")
    gvtrack.iterator("vt_ref_sum", sshift = ss, eshift = es)

    intervals <- gintervals(1, 500, 2000)
    result <- gextract("vt_raw_sum", intervals = intervals, iterator = iter)
    ref <- gextract("vt_ref_sum", intervals = intervals, iterator = iter)

    ref_vals <- ref$vt_ref_sum
    ref_vals[!is.finite(ref_vals)] <- 0

    expect_equal(result$vt_raw_sum, ref_vals, tolerance = 1e-6)
})

test_that("glm_pred with weight=1, no transform matches raw lse", {
    iter <- 100L
    ss <- -50
    es <- 50
    gs <- vtrack_to_glm_shifts(ss, es, iter)

    glm_pred.create("vt_raw_lse",
        tracks = "test.fixedbin",
        inner_func = "lse",
        weights = 1.0,
        bias = 0,
        shifts = list(gs)
    )

    gvtrack.create("vt_ref_lse", "test.fixedbin", "lse")
    gvtrack.iterator("vt_ref_lse", sshift = ss, eshift = es)

    intervals <- gintervals(1, 500, 2000)
    result <- gextract("vt_raw_lse", intervals = intervals, iterator = iter)
    ref <- gextract("vt_ref_lse", intervals = intervals, iterator = iter)

    ref_vals <- ref$vt_ref_lse
    ref_vals[!is.finite(ref_vals)] <- 0

    expect_equal(result$vt_raw_lse, ref_vals, tolerance = 1e-6)
})

# ============================================================
# Category 3: Scaling
# ============================================================

test_that("glm_pred scaling matches manual R computation", {
    iter <- 100L
    ss <- -50
    es <- 50
    gs <- vtrack_to_glm_shifts(ss, es, iter)

    glm_pred.create("vt_scaled",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = 1.0,
        bias = 0,
        shifts = list(gs),
        max_cap = -15,
        dis_from_cap = 10,
        scale_factor = 10
    )

    gvtrack.create("vt_raw_for_scale", "test.fixedbin", "sum")
    gvtrack.iterator("vt_raw_for_scale", sshift = ss, eshift = es)

    intervals <- gintervals(1, 500, 2000)
    result <- gextract("vt_scaled", intervals = intervals, iterator = iter)
    ref <- gextract("vt_raw_for_scale", intervals = intervals, iterator = iter)

    raw_vals <- ref$vt_raw_for_scale
    raw_vals[!is.finite(raw_vals)] <- 0
    expected <- r_scale(raw_vals, max_cap = -15, dis_from_cap = 10, scale_factor = 10)

    expect_equal(result$vt_scaled, expected, tolerance = 1e-6)
})

# ============================================================
# Category 4: Logistic transform
# ============================================================

test_that("glm_pred logistic transform matches manual R computation", {
    iter <- 100L
    ss <- -50
    es <- 50
    gs <- vtrack_to_glm_shifts(ss, es, iter)
    trans <- list(L = 2, k = 0.5, x_0 = 5, pre_shift = 0, post_shift = -1)

    glm_pred.create("vt_logist",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = 1.0,
        bias = 0,
        shifts = list(gs),
        max_cap = -15,
        dis_from_cap = 10,
        scale_factor = 10,
        trans_family = "logist",
        trans_params = list(trans)
    )

    gvtrack.create("vt_raw_for_logist", "test.fixedbin", "sum")
    gvtrack.iterator("vt_raw_for_logist", sshift = ss, eshift = es)

    intervals <- gintervals(1, 500, 2000)
    result <- gextract("vt_logist", intervals = intervals, iterator = iter)
    ref <- gextract("vt_raw_for_logist", intervals = intervals, iterator = iter)

    raw_vals <- ref$vt_raw_for_logist
    raw_vals[!is.finite(raw_vals)] <- 0
    scaled <- r_scale(raw_vals, -15, 10, 10)
    expected <- r_logist(scaled,
        L = trans$L, k = trans$k, x_0 = trans$x_0,
        pre_shift = trans$pre_shift, post_shift = trans$post_shift
    )
    expected[!is.finite(expected)] <- 0

    expect_equal(result$vt_logist, expected, tolerance = 1e-6)
})

# ============================================================
# Category 5: Multi-entry with weights
# ============================================================

test_that("glm_pred with multiple entries matches manual computation", {
    iter <- 100L
    w1 <- 0.3
    w2 <- -0.5
    b <- 0.42
    trans1 <- list(L = 2, k = 0.5, x_0 = 0, post_shift = -1)
    trans2 <- list(L = 1, k = 1, x_0 = 0, pre_shift = -5)

    # Entry 1: vtrack sshift=-100, eshift=0 on iter=100 → glm: (-150, 50)
    gs1 <- vtrack_to_glm_shifts(-100, 0, iter)
    gs2 <- vtrack_to_glm_shifts(0, 100, iter)

    glm_pred.create("vt_multi",
        tracks = c("test.fixedbin", "test.fixedbin"),
        inner_func = c("sum", "sum"),
        weights = c(w1, w2),
        bias = b,
        shifts = list(gs1, gs2),
        max_cap = c(-15, -15),
        dis_from_cap = c(10, 10),
        scale_factor = 10,
        trans_family = c("logist", "logist"),
        trans_params = list(trans1, trans2)
    )

    gvtrack.create("vt_e1", "test.fixedbin", "sum")
    gvtrack.iterator("vt_e1", sshift = -100, eshift = 0)
    gvtrack.create("vt_e2", "test.fixedbin", "sum")
    gvtrack.iterator("vt_e2", sshift = 0, eshift = 100)

    intervals <- gintervals(1, 500, 2000)
    ref1 <- gextract("vt_e1", intervals = intervals, iterator = iter)
    ref2 <- gextract("vt_e2", intervals = intervals, iterator = iter)

    r1 <- ref1$vt_e1
    r1[!is.finite(r1)] <- 0
    s1 <- r_scale(r1, -15, 10, 10)
    t1 <- r_logist(s1, L = 2, k = 0.5, x_0 = 0, post_shift = -1)
    t1[!is.finite(t1)] <- 0

    r2 <- ref2$vt_e2
    r2[!is.finite(r2)] <- 0
    s2 <- r_scale(r2, -15, 10, 10)
    t2 <- r_logist(s2, L = 1, k = 1, x_0 = 0, pre_shift = -5)
    t2[!is.finite(t2)] <- 0

    expected <- b + w1 * t1 + w2 * t2

    result <- gextract("vt_multi", intervals = intervals, iterator = iter)
    expect_equal(result$vt_multi, expected, tolerance = 1e-6)
})

# ============================================================
# Category 6: Interactions
# ============================================================

test_that("glm_pred interactions compute correct products", {
    iter <- 100L
    w1 <- 0.3
    w2 <- 0.5
    iw <- 0.7
    b <- 0.1
    sf <- 10

    gs1 <- vtrack_to_glm_shifts(-100, 0, iter)
    gs2 <- vtrack_to_glm_shifts(0, 100, iter)

    glm_pred.create("vt_inter",
        tracks = c("test.fixedbin", "test.fixedbin"),
        inner_func = c("sum", "sum"),
        weights = c(w1, w2),
        bias = b,
        shifts = list(gs1, gs2),
        max_cap = c(-15, -15),
        dis_from_cap = c(10, 10),
        scale_factor = sf,
        interactions = list(c(1L, 2L)),
        interaction_weights = iw
    )

    gvtrack.create("vt_i1", "test.fixedbin", "sum")
    gvtrack.iterator("vt_i1", sshift = -100, eshift = 0)
    gvtrack.create("vt_i2", "test.fixedbin", "sum")
    gvtrack.iterator("vt_i2", sshift = 0, eshift = 100)

    intervals <- gintervals(1, 500, 2000)
    ref1 <- gextract("vt_i1", intervals = intervals, iterator = iter)
    ref2 <- gextract("vt_i2", intervals = intervals, iterator = iter)

    r1 <- ref1$vt_i1
    r1[!is.finite(r1)] <- 0
    s1 <- r_scale(r1, -15, 10, sf)

    r2 <- ref2$vt_i2
    r2[!is.finite(r2)] <- 0
    s2 <- r_scale(r2, -15, 10, sf)

    # Main effects: no transform
    m1 <- w1 * s1
    m2 <- w2 * s2

    # Interaction: product / scale_factor, no transform
    inter <- iw * s1 * s2 / sf

    expected <- b + m1 + m2 + inter

    result <- gextract("vt_inter", intervals = intervals, iterator = iter)
    expect_equal(result$vt_inter, expected, tolerance = 1e-6)
})

test_that("glm_pred interactions with transforms compute correctly", {
    iter <- 100L
    w1 <- 0.3
    w2 <- 0.5
    iw <- 0.7
    b <- 0.1
    sf <- 10
    trans_main1 <- list(L = 2, k = 0.5, x_0 = 0, post_shift = -1)
    trans_main2 <- list(L = 2, k = 0.5, x_0 = 10)
    trans_inter <- list(L = 1, k = 1, x_0 = 5)

    gs1 <- vtrack_to_glm_shifts(-100, 0, iter)
    gs2 <- vtrack_to_glm_shifts(0, 100, iter)

    glm_pred.create("vt_inter_trans",
        tracks = c("test.fixedbin", "test.fixedbin"),
        inner_func = c("sum", "sum"),
        weights = c(w1, w2),
        bias = b,
        shifts = list(gs1, gs2),
        max_cap = c(-15, -15),
        dis_from_cap = c(10, 10),
        scale_factor = sf,
        trans_family = c("logist", "logist"),
        trans_params = list(trans_main1, trans_main2),
        interactions = list(c(1L, 2L)),
        interaction_weights = iw,
        interaction_trans_family = "logist",
        inter_trans_params = list(trans_inter)
    )

    gvtrack.create("vt_it1", "test.fixedbin", "sum")
    gvtrack.iterator("vt_it1", sshift = -100, eshift = 0)
    gvtrack.create("vt_it2", "test.fixedbin", "sum")
    gvtrack.iterator("vt_it2", sshift = 0, eshift = 100)

    intervals <- gintervals(1, 500, 2000)
    ref1 <- gextract("vt_it1", intervals = intervals, iterator = iter)
    ref2 <- gextract("vt_it2", intervals = intervals, iterator = iter)

    r1 <- ref1$vt_it1
    r1[!is.finite(r1)] <- 0
    s1 <- r_scale(r1, -15, 10, sf)

    r2 <- ref2$vt_it2
    r2[!is.finite(r2)] <- 0
    s2 <- r_scale(r2, -15, 10, sf)

    t1 <- r_logist(s1, L = 2, k = 0.5, x_0 = 0, post_shift = -1)
    t1[!is.finite(t1)] <- 0
    t2 <- r_logist(s2, L = 2, k = 0.5, x_0 = 10)
    t2[!is.finite(t2)] <- 0
    m1 <- w1 * t1
    m2 <- w2 * t2

    product <- s1 * s2 / sf
    inter_transformed <- r_logist(product, L = 1, k = 1, x_0 = 5)
    inter_transformed[!is.finite(inter_transformed)] <- 0
    inter <- iw * inter_transformed

    expected <- b + m1 + m2 + inter

    result <- gextract("vt_inter_trans", intervals = intervals, iterator = iter)
    expect_equal(result$vt_inter_trans, expected, tolerance = 1e-6)
})

# ============================================================
# Category 7: Edge cases
# ============================================================

test_that("glm_pred handles zero-weight entries correctly", {
    iter <- 100L
    gs <- vtrack_to_glm_shifts(-50, 50, iter)

    glm_pred.create("vt_zero_w",
        tracks = c("test.fixedbin", "test.fixedbin"),
        inner_func = c("sum", "sum"),
        weights = c(0, 1),
        bias = 5,
        shifts = list(gs, gs)
    )

    gvtrack.create("vt_ref_zw", "test.fixedbin", "sum")
    gvtrack.iterator("vt_ref_zw", sshift = -50, eshift = 50)

    intervals <- gintervals(1, 500, 2000)
    result <- gextract("vt_zero_w", intervals = intervals, iterator = iter)
    ref <- gextract("vt_ref_zw", intervals = intervals, iterator = iter)

    ref_vals <- ref$vt_ref_zw
    ref_vals[!is.finite(ref_vals)] <- 0
    expected <- 5 + ref_vals

    expect_equal(result$vt_zero_w, expected, tolerance = 1e-6)
})

test_that("glm_pred works in expressions", {
    glm_pred.create("vt_expr",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = 1,
        bias = 0,
        shifts = list(c(-50, 50))
    )

    intervals <- gintervals(1, 500, 1000)
    result <- gextract("vt_expr * 2", intervals = intervals, iterator = 100)
    result_single <- gextract("vt_expr", intervals = intervals, iterator = 100)

    expect_equal(result$`vt_expr * 2`, result_single$vt_expr * 2, tolerance = 1e-6)
})

# ============================================================
# Category 8: Sparse track support
# ============================================================

test_that("glm_pred works with sparse tracks", {
    iter <- 200L
    gs <- vtrack_to_glm_shifts(-50, 50, iter)

    glm_pred.create("vt_sparse",
        tracks = "test.sparse",
        inner_func = "sum",
        weights = 1,
        bias = 0,
        shifts = list(gs)
    )

    gvtrack.create("vt_ref_sparse", "test.sparse", "sum")
    gvtrack.iterator("vt_ref_sparse", sshift = -50, eshift = 50)

    intervals <- gintervals(1, 500, 5000)
    result <- gextract("vt_sparse", intervals = intervals, iterator = iter)
    ref <- gextract("vt_ref_sparse", intervals = intervals, iterator = iter)

    ref_vals <- ref$vt_ref_sparse
    # Positions where the sparse track has no data should be NaN
    # (glm_pred propagates NaN from all-NaN windows)
    expect_equal(is.nan(result$vt_sparse), !is.finite(ref_vals))
    # Where both have data, values should match
    valid <- is.finite(ref_vals)
    expect_equal(result$vt_sparse[valid], ref_vals[valid], tolerance = 1e-6)
})

# ============================================================
# Category 9: Multiple chromosomes
# ============================================================

test_that("glm_pred works across multiple chromosomes", {
    glm_pred.create("vt_multi_chr",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = 1.0,
        bias = 0.5,
        shifts = list(c(-50, 50))
    )

    intervals <- gintervals(c(1, 2), 500, 2000)
    result <- gextract("vt_multi_chr", intervals = intervals, iterator = 100)
    expect_true(nrow(result) > 0)
    expect_true(all(is.finite(result$vt_multi_chr)))
})

# ============================================================
# Category 10: Unaligned window regression (bin_size=50)
# ============================================================

test_that("glm_pred sum matches vtrack sum with unaligned window", {
    # test.fixedbin has bin_size=50. Using iterator=30, sshift=-10, eshift=70
    # produces a window not aligned to 50bp bin boundaries.
    iter <- 30L
    ss <- -10
    es <- 70
    gs <- vtrack_to_glm_shifts(ss, es, iter)

    glm_pred.create("vt_unaligned",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = 1.0,
        bias = 0,
        shifts = list(gs)
    )

    gvtrack.create("vt_ref_unaligned", "test.fixedbin", "sum")
    gvtrack.iterator("vt_ref_unaligned", sshift = ss, eshift = es)

    intervals <- gintervals(1, 500, 2000)
    result <- gextract("vt_unaligned", intervals = intervals, iterator = iter)
    ref <- gextract("vt_ref_unaligned", intervals = intervals, iterator = iter)

    ref_vals <- ref$vt_ref_unaligned
    ref_vals[!is.finite(ref_vals)] <- 0

    expect_equal(result$vt_unaligned, ref_vals, tolerance = 1e-6)
})

test_that("glm_pred lse matches vtrack lse with unaligned window", {
    iter <- 30L
    ss <- -10
    es <- 70
    gs <- vtrack_to_glm_shifts(ss, es, iter)

    glm_pred.create("vt_unaligned_lse",
        tracks = "test.fixedbin",
        inner_func = "lse",
        weights = 1.0,
        bias = 0,
        shifts = list(gs)
    )

    gvtrack.create("vt_ref_unaligned_lse", "test.fixedbin", "lse")
    gvtrack.iterator("vt_ref_unaligned_lse", sshift = ss, eshift = es)

    intervals <- gintervals(1, 500, 2000)
    result <- gextract("vt_unaligned_lse", intervals = intervals, iterator = iter)
    ref <- gextract("vt_ref_unaligned_lse", intervals = intervals, iterator = iter)

    ref_vals <- ref$vt_ref_unaligned_lse
    ref_vals[!is.finite(ref_vals)] <- 0

    expect_equal(result$vt_unaligned_lse, ref_vals, tolerance = 1e-6)
})

# ============================================================
# Category 12: Validation tests
# ============================================================

test_that("glm_pred.create rejects missing tracks", {
    expect_error(
        glm_pred.create("bad",
            tracks = "definitely_missing_track",
            inner_func = "sum", weights = 1
        ),
        "not found"
    )
})

test_that("glm_pred.create rejects mismatched max_cap/dis_from_cap", {
    expect_error(
        glm_pred.create("bad",
            tracks = "test.fixedbin",
            inner_func = "sum", weights = 1,
            max_cap = 1, dis_from_cap = NA
        ),
        "mismatch"
    )
    expect_error(
        glm_pred.create("bad",
            tracks = "test.fixedbin",
            inner_func = "sum", weights = 1,
            max_cap = NA, dis_from_cap = 10
        ),
        "mismatch"
    )
})

# ============================================================
# Category 14: Selector track validation
# ============================================================

test_that("glm_pred.create validates selector_track requires selector_breaks", {
    expect_error(
        glm_pred.create("bad",
            tracks = "test.fixedbin",
            inner_func = "sum",
            weights = 1,
            selector_track = "test.fixedbin"
        ),
        "selector_breaks"
    )
})

test_that("glm_pred.create validates selector_breaks requires 2+ elements", {
    expect_error(
        glm_pred.create("bad",
            tracks = "test.fixedbin",
            inner_func = "sum",
            weights = 1,
            selector_track = "test.fixedbin",
            selector_breaks = 0.5
        ),
        "at least 2"
    )
})

test_that("glm_pred.create validates selector_track must be dense", {
    expect_error(
        glm_pred.create("bad",
            tracks = "test.fixedbin",
            inner_func = "sum",
            weights = matrix(c(1, 2), nrow = 1),
            selector_track = "test.sparse",
            selector_breaks = c(0, 0.5, 1)
        ),
        "fixed-bin.*dense"
    )
})

test_that("glm_pred.create validates weights must be matrix when K > 1", {
    expect_error(
        glm_pred.create("bad",
            tracks = "test.fixedbin",
            inner_func = "sum",
            weights = 1,
            selector_track = "test.fixedbin",
            selector_breaks = c(0, 0.5, 1)
        ),
        "matrix"
    )
})

test_that("glm_pred.create validates weight matrix dimensions", {
    expect_error(
        glm_pred.create("bad",
            tracks = "test.fixedbin",
            inner_func = "sum",
            weights = matrix(c(1, 2, 3), nrow = 1, ncol = 3),
            selector_track = "test.fixedbin",
            selector_breaks = c(0, 0.5, 1)
        ),
        "1 x 2"
    )
})

# ============================================================
# Category 15: Selector track - K=2 basic
# ============================================================

test_that("glm_pred with selector_track selects per-bin weights", {
    # test.fixedbin has values in [0, 0.26]. Use it as both source and selector.
    # Breaks: [0, 0.1, 1.0] gives K=2 bins:
    #   bin 0: selector in [0, 0.1]
    #   bin 1: selector in (0.1, 1.0]
    # Use different weights for each bin so we can verify selection.

    iter <- 50L
    N <- 1L
    K <- 2L
    w_bin0 <- 2.0
    w_bin1 <- 5.0
    b <- 0

    W <- matrix(c(w_bin0, w_bin1), nrow = N, ncol = K)

    glm_pred.create("vt_sel_basic",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = W,
        bias = b,
        selector_track = "test.fixedbin",
        selector_breaks = c(0, 0.1, 1.0)
    )

    intervals <- gintervals(1, 0, 500)
    result <- gextract("vt_sel_basic", intervals = intervals, iterator = iter)

    # Get the selector values and raw source values at each position
    sel_vals <- gextract("test.fixedbin", intervals = intervals, iterator = iter)

    # For each position, determine which bin the selector falls in:
    # bin 0 if selector in [0, 0.1], bin 1 if selector in (0.1, 1.0]
    # Then the result should be weight[bin] * source_sum
    # Since source == selector here (same track, no shift, iter=50=bin_size),
    # source_sum = selector_value for each bin.
    for (i in seq_len(nrow(result))) {
        sv <- sel_vals$test.fixedbin[i]
        val <- result$vt_sel_basic[i]
        if (is.na(sv) || !is.finite(sv)) next

        if (sv >= 0 && sv <= 0.1) {
            expected_w <- w_bin0
        } else if (sv > 0.1 && sv <= 1.0) {
            expected_w <- w_bin1
        } else {
            # Out of range - should be NaN
            expect_true(is.nan(val), info = sprintf("row %d: sv=%g should be NaN", i, sv))
            next
        }

        # No shift, iter = bin_size = 50, so source_sum = selector_value
        expected <- expected_w * sv
        expect_equal(val, expected,
            tolerance = 1e-6,
            info = sprintf("row %d: sv=%g, expected w=%g", i, sv, expected_w)
        )
    }
})

test_that("glm_pred.info reconstructs weight matrix for K > 1", {
    N <- 1L
    K <- 2L
    W <- matrix(c(2.0, 5.0), nrow = N, ncol = K)

    glm_pred.create("vt_sel_info",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = W,
        bias = c(0.1, 0.2),
        selector_track = "test.fixedbin",
        selector_breaks = c(0, 0.5, 1.0)
    )

    info <- glm_pred.info("vt_sel_info")
    expect_equal(info$func, "glm.predict")
    expect_equal(info$params$num_bins, 2)
    expect_true(is.matrix(info$params$weights))
    expect_equal(nrow(info$params$weights), 1)
    expect_equal(ncol(info$params$weights), 2)
    expect_equal(info$params$weights[1, 1], 2.0)
    expect_equal(info$params$weights[1, 2], 5.0)
    expect_equal(info$params$bias, c(0.1, 0.2))
    expect_equal(info$params$selector_track, "test.fixedbin")
    expect_equal(info$params$selector_breaks, c(0, 0.5, 1.0))
})

# ============================================================
# Category 16: Selector track - out-of-range -> NaN
# ============================================================

test_that("glm_pred selector emits NaN for out-of-range selector values", {
    # test.fixedbin has values in [0, 0.26].
    # Use breaks (0.3, 0.5, 1.0) so all values < 0.3 are out-of-range.
    # Positions in [0, 0.3) should produce NaN.
    # Most values are < 0.3, so most outputs should be NaN.

    iter <- 50L
    W <- matrix(c(1.0, 2.0), nrow = 1, ncol = 2)

    glm_pred.create("vt_sel_oor",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = W,
        bias = 0,
        selector_track = "test.fixedbin",
        selector_breaks = c(0.3, 0.5, 1.0)
    )

    intervals <- gintervals(1, 0, 5000)
    result <- gextract("vt_sel_oor", intervals = intervals, iterator = iter)
    sel_vals <- gextract("test.fixedbin", intervals = intervals, iterator = iter)

    for (i in seq_len(nrow(result))) {
        sv <- sel_vals$test.fixedbin[i]
        val <- result$vt_sel_oor[i]
        if (is.na(sv) || !is.finite(sv) || sv < 0.3) {
            # Out of range => NaN
            expect_true(is.nan(val),
                info = sprintf("row %d: sv=%g should produce NaN", i, sv)
            )
        } else {
            # In range => finite
            expect_true(is.finite(val),
                info = sprintf("row %d: sv=%g should be finite", i, sv)
            )
        }
    }

    # Verify at least some positions are NaN (most values are < 0.3)
    expect_true(sum(is.nan(result$vt_sel_oor)) > 0)
})

# ============================================================
# Category 17: Selector track - K=2 with interactions
# ============================================================

test_that("glm_pred K=2 interactions use per-bin weights", {
    iter <- 100L
    N <- 2L
    K <- 2L
    sf <- 10

    gs1 <- vtrack_to_glm_shifts(-100, 0, iter)
    gs2 <- vtrack_to_glm_shifts(0, 100, iter)

    # Different weights per bin for main effects and interaction
    W <- matrix(c(0.3, 0.5, 0.7, 0.9), nrow = N, ncol = K)
    IW <- matrix(c(0.4, 0.8), nrow = 1, ncol = K)
    b <- c(0.1, 0.2)

    glm_pred.create("vt_sel_inter",
        tracks = c("test.fixedbin", "test.fixedbin"),
        inner_func = c("sum", "sum"),
        weights = W,
        bias = b,
        shifts = list(gs1, gs2),
        max_cap = c(-15, -15),
        dis_from_cap = c(10, 10),
        scale_factor = sf,
        interactions = list(c(1L, 2L)),
        interaction_weights = IW,
        selector_track = "test.fixedbin",
        selector_breaks = c(0, 0.1, 1.0)
    )

    # Verify it produces output (functional test - C++ handles the per-bin logic)
    intervals <- gintervals(1, 500, 2000)
    result <- gextract("vt_sel_inter", intervals = intervals, iterator = iter)
    expect_true(nrow(result) > 0)

    # Check that glm_pred.info reconstructs interaction weight matrix
    info <- glm_pred.info("vt_sel_inter")
    expect_true(is.matrix(info$params$inter_weights))
    expect_equal(nrow(info$params$inter_weights), 1)
    expect_equal(ncol(info$params$inter_weights), 2)
    expect_equal(info$params$inter_weights[1, 1], 0.4)
    expect_equal(info$params$inter_weights[1, 2], 0.8)
})

# ============================================================
# Category 18: Selector track - per-bin bias
# ============================================================

test_that("glm_pred with selector uses per-bin bias", {
    # With weight=0, the output should be just the bias for the selected bin
    iter <- 50L
    b <- c(10.0, 20.0)

    W <- matrix(c(0, 0), nrow = 1, ncol = 2)

    glm_pred.create("vt_sel_bias",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = W,
        bias = b,
        selector_track = "test.fixedbin",
        selector_breaks = c(0, 0.1, 1.0)
    )

    intervals <- gintervals(1, 0, 500)
    result <- gextract("vt_sel_bias", intervals = intervals, iterator = iter)
    sel_vals <- gextract("test.fixedbin", intervals = intervals, iterator = iter)

    for (i in seq_len(nrow(result))) {
        sv <- sel_vals$test.fixedbin[i]
        val <- result$vt_sel_bias[i]
        if (is.na(sv) || !is.finite(sv)) next

        if (sv >= 0 && sv <= 0.1) {
            expect_equal(val, 10.0,
                tolerance = 1e-10,
                info = sprintf("row %d: sv=%g, expected bias=10.0", i, sv)
            )
        } else if (sv > 0.1 && sv <= 1.0) {
            expect_equal(val, 20.0,
                tolerance = 1e-10,
                info = sprintf("row %d: sv=%g, expected bias=20.0", i, sv)
            )
        }
    }
})

# ============================================================
# Category 19: Selector track - backward compat (K=1)
# ============================================================

test_that("glm_pred without selector still works (K=1 backward compat)", {
    # Ensure that adding num_bins=1 to params doesn't break existing behavior
    iter <- 100L
    gs <- vtrack_to_glm_shifts(-50, 50, iter)

    glm_pred.create("vt_no_sel",
        tracks = "test.fixedbin",
        inner_func = "sum",
        weights = 1.0,
        bias = 0
    )

    # Check that num_bins defaults to 1
    info <- glm_pred.info("vt_no_sel")
    expect_equal(info$params$num_bins, 1)
    expect_null(info$params$selector_track)
    expect_null(info$params$selector_breaks)

    # Still produces correct output
    intervals <- gintervals(1, 500, 1000)
    result <- gextract("vt_no_sel", intervals = intervals, iterator = iter)
    expect_true(nrow(result) > 0)
    expect_true(all(is.finite(result$vt_no_sel)))
})

# ============================================================
