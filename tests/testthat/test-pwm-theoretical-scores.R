test_that("gseq.pwm_score_theoretical returns correct data frame structure", {
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.05, 0.8, 0.1, 0.05,
        0.1, 0.1, 0.7, 0.1
    ), nrow = 3, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))

    res <- gseq.pwm_score_theoretical(pssm, sub_rate = 0.05, n = 100, bidirect = FALSE)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 100)
    expect_true(all(c("score", "n_subs", "n_indels") %in% names(res)))
    expect_true(all(is.numeric(res$score)))
    expect_true(all(res$n_indels == 0L))
})

test_that("analytical mode attaches PMF attribute", {
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.05, 0.8, 0.1, 0.05
    ), nrow = 2, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))

    res <- gseq.pwm_score_theoretical(pssm,
        sub_rate = 0.05, n = 100,
        start = "consensus", bidirect = FALSE
    )
    pmf <- attr(res, "pmf")
    expect_s3_class(pmf, "data.frame")
    expect_true(all(c("score", "prob") %in% names(pmf)))
    expect_equal(sum(pmf$prob), 1, tolerance = 1e-6)
    expect_true(all(pmf$prob >= 0))
})

test_that("sub_rate = 0 gives consensus score", {
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.05, 0.8, 0.1, 0.05,
        0.1, 0.1, 0.7, 0.1
    ), nrow = 3, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))

    # Compute expected consensus score
    prior <- 0.01
    pssm_prob <- (pssm + prior) / rowSums(pssm + prior)
    pssm_log <- log(pssm_prob)
    consensus_score <- sum(apply(pssm_log, 1, max))

    # Analytical mode
    res <- gseq.pwm_score_theoretical(pssm, sub_rate = 0, n = 50, bidirect = FALSE)
    expect_true(all(abs(res$score - consensus_score) < 0.002))
    expect_true(all(res$n_subs == 0))

    # Simulation mode (bidirect)
    res2 <- gseq.pwm_score_theoretical(pssm, sub_rate = 0, n = 50, bidirect = TRUE)
    # Bidirect score >= forward-only consensus score
    expect_true(all(res2$score >= consensus_score - 1e-6))
})

test_that("higher sub_rate decreases mean score", {
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.05, 0.8, 0.1, 0.05,
        0.1, 0.1, 0.7, 0.1,
        0.02, 0.02, 0.02, 0.94
    ), nrow = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))

    set.seed(123)
    res_low <- gseq.pwm_score_theoretical(pssm, sub_rate = 0.01, n = 5000, bidirect = FALSE)
    res_high <- gseq.pwm_score_theoretical(pssm, sub_rate = 0.2, n = 5000, bidirect = FALSE)
    expect_gt(mean(res_low$score), mean(res_high$score))
})

test_that("consensus vs pssm start produce different distributions", {
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.05, 0.8, 0.1, 0.05,
        0.1, 0.1, 0.7, 0.1
    ), nrow = 3, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))

    set.seed(42)
    res_cons <- gseq.pwm_score_theoretical(pssm,
        sub_rate = 0.05, n = 5000,
        start = "consensus", bidirect = FALSE
    )
    set.seed(42)
    res_pssm <- gseq.pwm_score_theoretical(pssm,
        sub_rate = 0.05, n = 5000,
        start = "pssm", bidirect = FALSE
    )

    # Consensus start should have higher mean score (starting from best)
    expect_gt(mean(res_cons$score), mean(res_pssm$score))
})

test_that("simulation mode with indels works", {
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.05, 0.8, 0.1, 0.05,
        0.1, 0.1, 0.7, 0.1
    ), nrow = 3, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))

    set.seed(42)
    res <- gseq.pwm_score_theoretical(pssm,
        sub_rate = 0.05, indel_rate = 0.05,
        n = 200
    )
    expect_equal(nrow(res), 200)
    expect_true(any(res$n_indels > 0))
    # PMF should NOT be attached in simulation mode
    expect_null(attr(res, "pmf"))
})

test_that("rate_matrix parameter works", {
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.05, 0.8, 0.1, 0.05
    ), nrow = 2, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))

    # Transition-biased rate matrix (A<->G and C<->T more likely)
    rm <- matrix(c(
        0, 0.1, 0.8, 0.1,
        0.1, 0, 0.1, 0.8,
        0.8, 0.1, 0, 0.1,
        0.1, 0.8, 0.1, 0
    ), nrow = 4, byrow = TRUE)

    # Should run without error
    set.seed(42)
    res <- gseq.pwm_score_theoretical(pssm,
        sub_rate = 0.1, n = 500,
        rate_matrix = rm, bidirect = FALSE
    )
    expect_equal(nrow(res), 500)
})

test_that("bidirect mode runs and produces valid output", {
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.05, 0.8, 0.1, 0.05,
        0.1, 0.1, 0.7, 0.1
    ), nrow = 3, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))

    set.seed(42)
    res_bi <- gseq.pwm_score_theoretical(pssm,
        sub_rate = 0.1, n = 500,
        bidirect = TRUE
    )
    expect_equal(nrow(res_bi), 500)
    expect_true(all(is.finite(res_bi$score)))
    # No PMF in simulation mode (bidirect forces simulation)
    expect_null(attr(res_bi, "pmf"))
})

test_that("very short PSSM (L=1) works", {
    pssm <- matrix(c(0.9, 0.03, 0.03, 0.04),
        nrow = 1,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    res <- gseq.pwm_score_theoretical(pssm, sub_rate = 0.1, n = 100, bidirect = FALSE)
    expect_equal(nrow(res), 100)
    expect_true(all(is.finite(res$score)))
})

test_that("sub_rate = 1 mutates all positions", {
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.05, 0.8, 0.1, 0.05
    ), nrow = 2, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))

    set.seed(42)
    res <- gseq.pwm_score_theoretical(pssm,
        sub_rate = 1, n = 500,
        bidirect = FALSE
    )
    # All positions should be mutated
    expect_true(all(res$n_subs == 2))
})

test_that("scores are within theoretical range", {
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.05, 0.8, 0.1, 0.05,
        0.1, 0.1, 0.7, 0.1
    ), nrow = 3, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))

    prior <- 0.01
    pssm_prob <- (pssm + prior) / rowSums(pssm + prior)
    pssm_log <- log(pssm_prob)
    theoretical_max <- sum(apply(pssm_log, 1, max))
    theoretical_min <- sum(apply(pssm_log, 1, min))

    res <- gseq.pwm_score_theoretical(pssm, sub_rate = 0.1, n = 1000, bidirect = FALSE)
    # All scores should be within [min, max] (with small tolerance for discretization)
    expect_true(all(res$score >= theoretical_min - 0.01))
    expect_true(all(res$score <= theoretical_max + 0.01))
})

test_that("input validation works", {
    pssm <- matrix(c(0.9, 0.03, 0.03, 0.04),
        nrow = 1,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    expect_error(gseq.pwm_score_theoretical(pssm, sub_rate = -0.1))
    expect_error(gseq.pwm_score_theoretical(pssm, sub_rate = 1.5))
    expect_error(gseq.pwm_score_theoretical(pssm, sub_rate = 0.1, indel_rate = -1))
    expect_error(gseq.pwm_score_theoretical(pssm, sub_rate = 0.1, n = 0))
    expect_error(gseq.pwm_score_theoretical(pssm, sub_rate = 0.1, start = "invalid"))
    expect_error(gseq.pwm_score_theoretical(pssm, sub_rate = 0.1, rate_matrix = matrix(1:9, 3, 3)))
})
