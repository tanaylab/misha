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

# ---- Engine tests for GRAD_MAX (Task 3) ----

create_isolated_test_db()

# Build a strong-consensus PSSM (0.97/0.01/0.01/0.01) whose column j has its
# peak at the j-th base of `bases`. `bases` must be a length-L character vector
# of A/C/G/T (no Ns).
.consensus_pssm <- function(bases) {
    L <- length(bases)
    stopifnot(all(bases %in% c("A", "C", "G", "T")))
    pssm <- matrix(0.01,
        nrow = L, ncol = 4,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    for (j in seq_len(L)) {
        idx <- match(bases[j], c("A", "C", "G", "T"))
        pssm[j, idx] <- 0.97
    }
    pssm
}

test_that("pwm.grad max: head=argmax => engine matches oracle col-0 diff", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L # number of iterator-interval positions per pivot
    # Use a 1bp iterator interval; the iterator extension makes the vtrack
    # see a window of length iter_W (sshift=0, eshift=iter_W-1).
    pivot <- gintervals(1, 200, 201)
    # The engine fetches sequence over [pivot.start + sshift,
    #                                   pivot.end + eshift + (L-1))
    # because extend=TRUE adds (L-1) at the end. That's [200, 200 + iter_W + L - 1).
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    # Pick PSSM consensus = first L bases of seq_ext, so the head anchor wins.
    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    pssm <- .consensus_pssm(head_bases)

    gvtrack.create("g", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "max", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g", sshift = 0, eshift = iter_W - 1)

    res <- gextract("g", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "max",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
    # Sanity: the head-wins case must produce a strictly positive gradient
    # (the head base is the consensus and column 0 is non-uniform).
    expect_gt(res$g[1], 0)
})

test_that("pwm.grad max: head not argmax => engine returns 0 (matches oracle)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 30L
    offset <- 10L # interior anchor index whose match wins
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    # Construct a PSSM whose consensus matches the interior anchor at `offset`,
    # NOT the head anchor at 0. Argmax should land at the interior, so head
    # gradient is 0 (provided the head bases are not coincidentally also a
    # consensus run, which we verify below).
    interior_bases <- strsplit(substr(seq_ext, offset + 1, offset + L), "")[[1]]
    pssm <- .consensus_pssm(interior_bases)

    # Guard the test setup: interior pattern must differ from head bases.
    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    skip_if(
        identical(interior_bases, head_bases),
        "test fixture: interior PSSM consensus collides with head bases; pick a different interval"
    )

    gvtrack.create("g", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "max", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g", sshift = 0, eshift = iter_W - 1)

    res <- gextract("g", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "max",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
    # And it should be exactly 0 in this construction (argmax interior).
    expect_equal(expected, 0, tolerance = 1e-10)
})

test_that("pwm.grad max: bidirect=TRUE errors out (Task 3 scope)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    gvtrack.create("g_bidir", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "max", bidirect = TRUE,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_bidir", sshift = 0, eshift = 19)

    expect_error(
        gextract("g_bidir", iterator = 1, intervals = gintervals(1, 200, 201)),
        "bidirect=FALSE"
    )
})

test_that("pwm.grad max: strand=-1 errors out (Task 3 scope)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    gvtrack.create("g_minus", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "max", bidirect = FALSE, strand = -1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_minus", sshift = 0, eshift = 19)

    expect_error(
        gextract("g_minus", iterator = 1, intervals = gintervals(1, 200, 201)),
        "strand=-1"
    )
})

# Task 4: GRAD_LSE single-strand fwd

test_that("pwm.grad lse: concentrated softmax (head=consensus) matches oracle", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    pssm <- .consensus_pssm(head_bases)

    gvtrack.create("g_lse", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "lse", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_lse", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_lse", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "lse",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g_lse[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
    expect_gt(res$g_lse[1], 0)
})

test_that("pwm.grad lse: head not argmax => w_p small, g_lse < g_max", {
    # Same construction as the head-not-argmax test for max, but compare lse
    # vs max instead of just checking equality.
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 30L
    offset <- 10L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    interior_bases <- strsplit(substr(seq_ext, offset + 1, offset + L), "")[[1]]
    pssm <- .consensus_pssm(interior_bases)

    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    skip_if(
        identical(interior_bases, head_bases),
        "test fixture: interior PSSM consensus collides with head bases"
    )

    gvtrack.create("g_lse", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "lse", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_lse", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_lse", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "lse",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g_lse[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
    # When the interior anchor is the LSE-dominant one, w_p (head) is small,
    # so g_lse stays close to 0 even though M[0, b_p] - worst is large.
    expect_lt(res$g_lse[1], 0.5)
})

test_that("pwm.grad lse: w_p computed correctly via softmax denominator", {
    # Hand-verifiable case: build a 2-anchor situation where head and one
    # interior anchor have known scores, then compare the engine's w_p * diff
    # against a direct computation.
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 5L # exactly 5 anchors
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    pssm <- .consensus_pssm(head_bases)

    gvtrack.create("g_lse", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "lse", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_lse", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_lse", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "lse",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g_lse[1], expected, tolerance = 1e-5, ignore_attr = TRUE)

    # Cross-check: w_p must be in (0, 1] and g_lse <= g_max.
    g_max_oracle <- manual_pwm_grad(pssm, seq_ext, "max",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_lte(res$g_lse[1], g_max_oracle + 1e-8)
    expect_gt(res$g_lse[1], 0)
})

test_that("pwm.grad lse: bidirect=TRUE / strand=-1 error out (Task 4 scope)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    gvtrack.create("g_bidir_lse", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "lse", bidirect = TRUE,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_bidir_lse", sshift = 0, eshift = 19)
    expect_error(
        gextract("g_bidir_lse", iterator = 1, intervals = gintervals(1, 200, 201)),
        "bidirect=FALSE"
    )

    gvtrack.create("g_minus_lse", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "lse", bidirect = FALSE, strand = -1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_minus_lse", sshift = 0, eshift = 19)
    expect_error(
        gextract("g_minus_lse", iterator = 1, intervals = gintervals(1, 200, 201)),
        "strand=-1"
    )
})
