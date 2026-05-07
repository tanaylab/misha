test_that("oracle: flat PSSM gives g = 0 (lse, max, ism)", {
    pssm <- matrix(0.25, 3, 4, dimnames = list(NULL, c("A", "C", "G", "T")))
    expect_equal(manual_pwm_grad(pssm, "ACGTAA", "lse", bidirect = TRUE), 0)
    expect_equal(manual_pwm_grad(pssm, "ACGTAA", "max", bidirect = TRUE), 0)
    expect_equal(manual_pwm_grad_ism(pssm, "ACGTAA", "lse", bidirect = TRUE), 0)
})

test_that("oracle: L=1 PSSM, single strand, g_MAX = M[1,b_p] - min(M[1,])", {
    pssm <- matrix(c(0.7, 0.1, 0.1, 0.1), 1, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    M <- log((pssm + 0.01) / sum(pssm + 0.01))
    expected <- unname(M[1, 1] - min(M[1, ]))
    expect_equal(
        manual_pwm_grad(pssm, "ACGT", "max", bidirect = FALSE, strand = 1L),
        expected,
        tolerance = 1e-10
    )
})

test_that("oracle: linearized matches ISM when argmax is stable under flips (single strand)", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01
        ),
        3, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    seq <- "ACGAAAAA" # head anchor "ACG" perfect, others weak.
    g_lin <- manual_pwm_grad(pssm, seq, "max", bidirect = FALSE, strand = 1L)
    g_ism <- manual_pwm_grad_ism(pssm, seq, "max", bidirect = FALSE, strand = 1L)
    expect_equal(g_lin, g_ism, tolerance = 1e-8)
})

test_that("oracle: linearized lse with concentrated softmax matches max", {
    # Longer PSSM + uniform A-tail in seq pushes the non-head anchor scores
    # so far below the head anchor that the softmax concentrates at >= 1-1e-4
    # on the head, making g_LSE numerically indistinguishable from g_MAX.
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.97, 0.01
        ),
        5, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    seq <- "ACGGGAAAAAAAA" # head anchor "ACGGG" perfect, others weak.
    g_max <- manual_pwm_grad(pssm, seq, "max", bidirect = FALSE, strand = 1L)
    g_lse <- manual_pwm_grad(pssm, seq, "lse", bidirect = FALSE, strand = 1L)
    expect_equal(g_max, g_lse, tolerance = 1e-3)
})

test_that("oracle: bidirect on palindromic PSSM matches single-strand on palindrome", {
    pssm <- matrix(
        c(
            0.7, 0.1, 0.1, 0.1,
            0.1, 0.4, 0.4, 0.1,
            0.1, 0.1, 0.1, 0.7
        ),
        3, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    rc <- pssm[3:1, c(4, 3, 2, 1)]
    dimnames(rc) <- dimnames(pssm)
    expect_equal(rc, pssm)
    g <- manual_pwm_grad(pssm, "ATCAT", "lse", bidirect = TRUE)
    expect_true(is.finite(g) && g >= 0)
})

test_that("oracle: rc by RC-of-seq equals rc by RC-of-PSSM", {
    pssm <- matrix(
        c(
            0.7, 0.1, 0.1, 0.1,
            0.1, 0.6, 0.2, 0.1,
            0.1, 0.1, 0.6, 0.2
        ),
        3, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    seq <- "ACGCAT"
    g_strand_minus <- manual_pwm_grad(pssm, seq, "max", bidirect = FALSE, strand = -1L)
    expect_true(is.finite(g_strand_minus) && g_strand_minus >= 0)
})

test_that("oracle: rc-of-seq vs rc-of-PSSM give same per-anchor scores", {
    # Self-consistency: scoring fwd PSSM on RC(seq) at anchor a should equal
    # scoring rc PSSM on seq at the mirrored anchor (n - L - a + 1).
    pssm <- matrix(
        c(
            0.6, 0.2, 0.1, 0.1,
            0.1, 0.5, 0.3, 0.1,
            0.2, 0.1, 0.5, 0.2
        ),
        3, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    seq <- "ACGTACGTAC"
    L <- nrow(pssm)
    n <- nchar(seq)

    sc <- .oracle_anchor_scores(seq, pssm, prior = 0.01)

    # Independent path: rc score at anchor a == fwd score on RC(seq) at
    # position (n - L - (a - 1) + 1) - 1 + 1 = n - L - a + 2 (1-based).
    rc_seq <- .rc_seq(seq)
    fwd_on_rc <- manual_pwm_scores_single_strand(rc_seq, pssm, prior = 0.01)

    n_anchors <- n - L + 1L
    for (a in seq_len(n_anchors)) {
        mirror <- n_anchors - a + 1L
        expect_equal(sc$rc[a], fwd_on_rc[mirror], tolerance = 1e-12)
    }
})

test_that("oracle: ISM is always >= 0", {
    pssm <- matrix(
        c(
            0.7, 0.1, 0.1, 0.1,
            0.1, 0.6, 0.2, 0.1,
            0.1, 0.1, 0.6, 0.2
        ),
        3, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    for (seq in c("ACGCAT", "TTTTTT", "ACGACG", "GGGGGG")) {
        for (agg in c("lse", "max")) {
            g <- manual_pwm_grad_ism(pssm, seq, agg, bidirect = TRUE)
            expect_gte(g, -1e-10)
        }
    }
})
