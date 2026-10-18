create_isolated_test_db()

test_that("reach: a bin [s, e) of each sequence vtrack reads [s + lo, e + hi)", {
    define_variant_vtracks()
    pssm <- variant_test_pssm()
    gvtrack.create("r_noext", NULL, "pwm.max", params = list(pssm = pssm, bidirect = TRUE, extend = FALSE, prior = 0.01))
    reach <- function(expr) misha:::.gvariants_reach(expr)
    expect_equal(reach("v_pwm"), cbind(lo = -5, hi = 3 + 3))
    expect_equal(reach("r_noext"), cbind(lo = 0, hi = 0))
    expect_equal(reach("v_potts"), cbind(lo = 4, hi = 4 + 5))
    expect_equal(reach("v_kmer"), cbind(lo = -2, hi = 2 + 2))
    expect_equal(reach("v_masked"), cbind(lo = -3, hi = 3))
    expect_equal(reach("v_ed"), cbind(lo = 0, hi = 4 + 1 - 1))
    expect_equal(nrow(reach("v_pwm + 2 * v_potts + test.fixedbin")), 2L)
})

test_that("reach: an expression with no sequence vtrack is an error", {
    expect_error(misha:::.gvariants_reach("test.fixedbin * 2"), "no sequence virtual track")
})
