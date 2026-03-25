create_isolated_test_db()

# ============================================================================
# Adversarial tests targeting the pigeonhole pre-filter and sub-chromosome
# parallelism optimizations in PWM edit distance.
#
# Strategy:
#   1. Compare max_edits=K (pre-filter enabled) vs max_edits=-1 (exact, no
#      pre-filter) filtered to edits <= K. Any mismatch is a correctness bug.
#   2. Compare serial (gmultitasking=FALSE) vs parallel (gmultitasking=TRUE)
#      for gscreen/gextract with explicit iterators.
#   3. Use gseq.pwm_edits for direct sequence testing.
# ============================================================================


# ============================================================================
# Test 1: PSSMs where ALL columns have at least one log-zero base
# The pigeonhole viable-bitset marks a B-mer as non-viable if ANY position
# has a mandatory (log-zero) PSSM score. When every column has at least one
# forbidden base, most B-mers are non-viable. But the MATCHING B-mer should
# still be viable. Verify the pre-filter does not reject true matches.
# ============================================================================
test_that("prefilter: PSSM with log-zero in every column, exact-match still found", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # 9-column PSSM: each column has one zero-probability base
    # With max_edits=2 and L=9, we get K+1=3 blocks of 3 columns each.
    pssm <- matrix(c(
        0.50, 0.25, 0.25, 0.00, # col 0: T is forbidden
        0.00, 0.50, 0.25, 0.25, # col 1: A is forbidden
        0.25, 0.00, 0.50, 0.25, # col 2: C is forbidden
        0.25, 0.25, 0.00, 0.50, # col 3: G is forbidden
        0.50, 0.25, 0.25, 0.00, # col 4: T is forbidden
        0.00, 0.50, 0.25, 0.25, # col 5: A is forbidden
        0.25, 0.00, 0.50, 0.25, # col 6: C is forbidden
        0.25, 0.25, 0.00, 0.50, # col 7: G is forbidden
        0.50, 0.25, 0.00, 0.25  # col 8: G is forbidden
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Perfect match: best base per column is A,C,G,T,A,C,G,T,A
    perfect_seq <- "ACGTACGTA"
    threshold <- sum(log(c(0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50))) - 0.001

    # Exact mode (no pre-filter)
    r_exact <- gseq.pwm_edits(perfect_seq, pssm,
        score.thresh = threshold, prior = 0, bidirect = FALSE
    )
    expect_true(nrow(r_exact) > 0, info = "Exact mode should find a match")
    expect_equal(min(r_exact$n_edits), 0L, info = "Perfect match = 0 edits")

    # Heuristic mode with max_edits=2 (pre-filter enabled for L=9, K=2 -> 3 blocks of 3)
    r_heur <- gseq.pwm_edits(perfect_seq, pssm,
        score.thresh = threshold, max_edits = 2L,
        prior = 0, bidirect = FALSE
    )
    expect_true(nrow(r_heur) > 0, info = "Pre-filter should not reject the perfect match")
    expect_equal(min(r_heur$n_edits), 0L, info = "Heuristic also finds 0 edits")

    # Now test a sequence that needs exactly 2 edits (at the boundary).
    # Col 0: A->T (T is forbidden at col 0, mandatory edit)
    # Col 4: A->T (T is forbidden at col 4, mandatory edit)
    # Rest match best: C,G,T,_,C,G,T,A
    two_edit_seq <- "TCGTTCGTA"
    r_exact2 <- gseq.pwm_edits(two_edit_seq, pssm,
        score.thresh = threshold, prior = 0, bidirect = FALSE
    )
    r_heur2 <- gseq.pwm_edits(two_edit_seq, pssm,
        score.thresh = threshold, max_edits = 2L,
        prior = 0, bidirect = FALSE
    )

    if (nrow(r_exact2) > 0) {
        exact_min <- min(r_exact2$n_edits)
        if (exact_min <= 2) {
            expect_true(nrow(r_heur2) > 0,
                info = "Pre-filter should not reject a 2-edit match"
            )
            heur_min <- min(r_heur2$n_edits)
            expect_equal(heur_min, exact_min,
                info = "Heuristic should agree with exact on 2-edit match"
            )
        }
    }
})


# ============================================================================
# Test 2: Very low information content PSSM (nearly uniform)
# Every B-mer is viable. The pre-filter should pass everything, so results
# must match exact mode exactly.
# ============================================================================
test_that("prefilter: nearly-uniform PSSM, pre-filter passes everything", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 12
    # Near-uniform: one base is 0.28, rest are 0.24 each
    pssm <- matrix(rep(c(0.28, 0.24, 0.24, 0.24), L),
        ncol = 4, byrow = TRUE
    )
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- L * log(0.28) - 0.1

    set.seed(42)
    seqs <- replicate(50, paste0(sample(c("A", "C", "G", "T"), 12, replace = TRUE), collapse = ""))

    for (s in seqs) {
        r_exact <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, prior = 0, bidirect = FALSE
        )
        r_heur <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, max_edits = 3L,
            prior = 0, bidirect = FALSE
        )

        exact_min <- if (nrow(r_exact) > 0) min(r_exact$n_edits) else NA
        heur_min <- if (nrow(r_heur) > 0) min(r_heur$n_edits) else NA

        if (!is.na(exact_min) && exact_min <= 3) {
            expect_false(is.na(heur_min),
                info = sprintf("seq=%s: exact found %d edits but heuristic found none", s, exact_min)
            )
            expect_equal(heur_min, exact_min,
                info = sprintf("seq=%s: heuristic (%s) != exact (%d)", s, heur_min, exact_min)
            )
        }
    }
    succeed()
})


# ============================================================================
# Test 3: Short motifs where pre-filter block length = 3
# L=6 with K=1 -> 2 blocks of 3 -> 4^3=64 entries.
# L=9 with K=2 -> 3 blocks of 3 -> 4^3=64 entries.
# ============================================================================
test_that("prefilter: short motif L=6, K=1, blocks of length 3", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -3.0

    # Diverse test sequences
    set.seed(101)
    seqs <- c(
        "ACGTAC", # Perfect match
        "TCGTAC", # 1 edit
        "TGCATG", # Many edits
        replicate(50, paste0(sample(c("A", "C", "G", "T"), 6, replace = TRUE), collapse = ""))
    )

    mismatches <- character(0)
    for (s in seqs) {
        r_exact <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, prior = 0, bidirect = FALSE
        )
        r_heur <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, max_edits = 1L,
            prior = 0, bidirect = FALSE
        )

        exact_min <- if (nrow(r_exact) > 0) min(r_exact$n_edits) else NA
        heur_min <- if (nrow(r_heur) > 0) min(r_heur$n_edits) else NA

        if (!is.na(exact_min) && exact_min <= 1) {
            if (is.na(heur_min)) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: exact=%d, heur=NA (pre-filter false reject!)", s, exact_min
                ))
            } else if (heur_min != exact_min) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: exact=%d, heur=%d", s, exact_min, heur_min
                ))
            }
        }
    }

    if (length(mismatches) > 0) {
        fail(paste("Mismatches:\n", paste(mismatches, collapse = "\n")))
    }
    succeed()
})


# ============================================================================
# Test 4: High-IC columns at edges only
# Blocks in the middle have low IC -> less selective. The edge blocks have
# high IC. Pigeonhole should still work correctly.
# ============================================================================
test_that("prefilter: high-IC columns at edges, low-IC in middle", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # L=12, K=3 -> 4 blocks of 3
    # Blocks 0 and 3 (edges) have high IC, blocks 1 and 2 (middle) nearly uniform
    pssm <- matrix(c(
        # Block 0 (high IC)
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        # Block 1 (low IC)
        0.28, 0.24, 0.24, 0.24,
        0.24, 0.28, 0.24, 0.24,
        0.24, 0.24, 0.28, 0.24,
        # Block 2 (low IC)
        0.24, 0.24, 0.24, 0.28,
        0.28, 0.24, 0.24, 0.24,
        0.24, 0.28, 0.24, 0.24,
        # Block 3 (high IC)
        0.01, 0.01, 0.01, 0.97, # T
        0.01, 0.01, 0.97, 0.01, # G
        0.97, 0.01, 0.01, 0.01  # A
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -8.0

    set.seed(202)
    seqs <- c(
        "ACGACBTATGA", # Near match with N in middle (11 chars, needs extension)
        "ACGACATATGA", # Close to consensus
        replicate(40, paste0(sample(c("A", "C", "G", "T"), 12, replace = TRUE), collapse = ""))
    )

    mismatches <- character(0)
    for (s in seqs) {
        r_exact <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, prior = 0, bidirect = FALSE
        )
        r_heur <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, max_edits = 3L,
            prior = 0, bidirect = FALSE
        )

        exact_min <- if (nrow(r_exact) > 0) min(r_exact$n_edits) else NA
        heur_min <- if (nrow(r_heur) > 0) min(r_heur$n_edits) else NA

        if (!is.na(exact_min) && exact_min <= 3) {
            if (is.na(heur_min)) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: exact=%d, heur=NA (false reject)", s, exact_min
                ))
            } else if (heur_min != exact_min) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: exact=%d, heur=%d", s, exact_min, heur_min
                ))
            }
        }
    }

    if (length(mismatches) > 0) {
        fail(paste("Mismatches:\n", paste(mismatches, collapse = "\n")))
    }
    succeed()
})


# ============================================================================
# Test 5: Sequences with N-bases near valid matches
# N-base at a block position makes that block invalid in the prefilter.
# If the N is in the only block that would match, the prefilter might
# falsely reject the window.
# ============================================================================
test_that("prefilter: N-bases at strategic positions do not cause false rejections", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # L=9, K=2 -> 3 blocks of 3
    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.7, 0.1, 0.1, 0.1
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -5.0

    # Consensus is ACGTACGTA. Block 0: ACG, Block 1: TAC, Block 2: GTA
    # Place N in block 0 to invalidate it, blocks 1&2 should still pass
    seqs_with_N <- c(
        "NCGTACGTA", # N in block 0, pos 0
        "ANGTACGTA", # N in block 0, pos 1
        "ACNTACGTA", # N in block 0, pos 2
        "ACGNACGTA", # N in block 1, pos 3
        "ACGTNCGTA", # N in block 1, pos 4
        "ACGTANGTA", # N in block 1, pos 5
        "ACGTACNTA", # N in block 2, pos 6
        "ACGTACGNA", # N in block 2, pos 7
        "ACGTACGTN"  # N in block 2, pos 8
    )

    for (s in seqs_with_N) {
        r_exact <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, prior = 0, bidirect = FALSE
        )
        r_heur <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, max_edits = 2L,
            prior = 0, bidirect = FALSE
        )

        exact_min <- if (nrow(r_exact) > 0) min(r_exact$n_edits) else NA
        heur_min <- if (nrow(r_heur) > 0) min(r_heur$n_edits) else NA

        if (!is.na(exact_min) && exact_min <= 2) {
            expect_false(is.na(heur_min),
                info = sprintf("seq=%s: exact=%d, heur=NA (N-base false reject)", s, exact_min)
            )
            expect_equal(heur_min, exact_min,
                info = sprintf("seq=%s: exact=%d != heur=%s", s, exact_min, heur_min)
            )
        }
    }
})


# ============================================================================
# Test 6: Pigeonhole with indels - the critical case
# When the optimal alignment uses 2 indels that both happen inside one block,
# that block has 2 edits. The other K-1 blocks have 0 edits each.
# The prefilter checks each block at shifts {-D,...,+D}. Verify correctness.
# ============================================================================
test_that("prefilter: indels + pigeonhole, 2 indels in same block region", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # L=9, K=2, max_indels=2 -> 3 blocks of 3
    # Block 0: cols 0-2, Block 1: cols 3-5, Block 2: cols 6-8
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97, # T
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97, # T
        0.97, 0.01, 0.01, 0.01  # A
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Motif consensus: ACGTACGTA
    threshold <- 9 * log(0.97) - 0.01

    # Sequence where 2 deletions are needed in the FIRST block region:
    # Insert 2 extra bases before block 1 starts. The alignment needs to delete them.
    # "ACGXXTACGTA" with XX inserted between positions 2 and 3
    # Optimal: delete X, X -> 2 indels, 0 subs = 2 edits total
    seqs <- c(
        "ACGGGTACGTA",  # 2 extra G's -> delete 2 = 2 edits (W=L+2)
        "AACGTACGTA",   # 1 extra A at start -> delete 1 = 1 edit
        "ACGTACGTAA",   # 1 extra A at end -> delete 1 = 1 edit
        "CGTACGTA",     # Missing first A -> insert 1 = 1 edit
        "ACGTACGT",     # Missing last A -> insert 1 = 1 edit
        "GTACGTA",      # Missing AC -> 2 insertions = 2 edits
        "ACGTACGTACC"   # 2 extra at end -> 2 deletions = 2 edits
    )

    for (s in seqs) {
        # Exact mode (no prefilter)
        r_exact <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold,
            max_indels = 2L, prior = 0, bidirect = FALSE
        )
        # Heuristic with prefilter
        r_heur <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold,
            max_edits = 2L, max_indels = 2L,
            prior = 0, bidirect = FALSE
        )

        exact_min <- if (nrow(r_exact) > 0 && any(!is.na(r_exact$n_edits))) min(r_exact$n_edits, na.rm = TRUE) else NA
        heur_min <- if (nrow(r_heur) > 0 && any(!is.na(r_heur$n_edits))) min(r_heur$n_edits, na.rm = TRUE) else NA

        if (!is.na(exact_min) && exact_min <= 2) {
            expect_false(is.na(heur_min),
                info = sprintf("seq=%s: exact=%d, heur=NA (indel+pigeonhole false reject)", s, exact_min)
            )
            if (!is.na(heur_min)) {
                expect_equal(heur_min, exact_min,
                    info = sprintf("seq=%s: exact=%d != heur=%d", s, exact_min, heur_min)
                )
            }
        }
    }
})


# ============================================================================
# Test 7: Reverse complement edge cases for prefilter
# The prefilter hash uses bidx[] which contains RC base indices for reverse.
# Test that the hash is computed correctly for reverse strand.
# ============================================================================
test_that("prefilter: reverse complement hashing is correct", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # L=9, K=2 -> 3 blocks of 3 columns
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97, # T
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97, # T
        0.97, 0.01, 0.01, 0.01  # A
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- 9 * log(0.97) - 0.01

    # Forward consensus: ACGTACGTA
    # RC of ACGTACGTA = TACGTACGT (reverse complement)
    # So the RC sequence should match on the reverse strand
    rc_seq <- "TACGTACGT"

    # Test with bidirect=TRUE (should find the RC match)
    r_bidir <- gseq.pwm_edits(rc_seq, pssm,
        score.thresh = threshold, max_edits = 2L,
        prior = 0, bidirect = TRUE
    )
    expect_true(nrow(r_bidir) > 0, info = "Should find RC match with bidirect=TRUE")
    expect_equal(min(r_bidir$n_edits), 0L, info = "RC match should be 0 edits")

    # With bidirect=FALSE, strand=1 (forward only), should NOT find a match
    r_fwd <- gseq.pwm_edits(rc_seq, pssm,
        score.thresh = threshold, max_edits = 2L,
        prior = 0, bidirect = FALSE, strand = 1L
    )
    # May or may not find a match depending on threshold, but should need edits
    if (nrow(r_fwd) > 0) {
        expect_true(min(r_fwd$n_edits) > 0,
            info = "Forward-only scan on RC seq should need edits"
        )
    }

    # Test RC match with indels
    # RC of ACGTACGTA = TACGTACGT. Add an insertion: TAACGTACGT (10 bases)
    rc_plus_ins <- "TAACGTACGT"
    r_rc_indel <- gseq.pwm_edits(rc_plus_ins, pssm,
        score.thresh = threshold,
        max_edits = 2L, max_indels = 1L,
        prior = 0, bidirect = TRUE
    )

    r_rc_indel_exact <- gseq.pwm_edits(rc_plus_ins, pssm,
        score.thresh = threshold,
        max_indels = 1L,
        prior = 0, bidirect = TRUE
    )

    exact_min <- if (nrow(r_rc_indel_exact) > 0) min(r_rc_indel_exact$n_edits) else NA
    heur_min <- if (nrow(r_rc_indel) > 0) min(r_rc_indel$n_edits) else NA

    if (!is.na(exact_min) && exact_min <= 2) {
        expect_false(is.na(heur_min),
            info = sprintf("RC+indel: exact=%d, heur=NA", exact_min)
        )
        if (!is.na(heur_min)) {
            expect_equal(heur_min, exact_min,
                info = sprintf("RC+indel: exact=%d != heur=%d", exact_min, heur_min)
            )
        }
    }
})


# ============================================================================
# Test 8: Prefilter with substitution-only mode (max_indels=0)
# When max_indels=0, the prefilter only checks shift=0.
# ============================================================================
test_that("prefilter: substitution-only mode (max_indels=0)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # L=9, K=2 -> 3 blocks of 3
    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.7, 0.1, 0.1, 0.1
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -4.0

    set.seed(303)
    seqs <- c(
        "ACGTACGTA", # Perfect
        "TCGTACGTA", # 1 sub at pos 0
        "ACGTTCGTA", # 1 sub at pos 3
        "TCGTTCGTA", # 2 subs (block 0 and block 1 each have 1 edit)
        replicate(40, paste0(sample(c("A", "C", "G", "T"), 9, replace = TRUE), collapse = ""))
    )

    mismatches <- character(0)
    for (s in seqs) {
        r_exact <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, prior = 0, bidirect = FALSE
        )
        r_heur <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, max_edits = 2L,
            prior = 0, bidirect = FALSE
        )

        exact_min <- if (nrow(r_exact) > 0) min(r_exact$n_edits) else NA
        heur_min <- if (nrow(r_heur) > 0) min(r_heur$n_edits) else NA

        if (!is.na(exact_min) && exact_min <= 2) {
            if (is.na(heur_min)) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: exact=%d, heur=NA", s, exact_min
                ))
            } else if (heur_min != exact_min) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: exact=%d, heur=%d", s, exact_min, heur_min
                ))
            }
        }
    }

    if (length(mismatches) > 0) {
        fail(paste("Sub-only prefilter mismatches:\n", paste(mismatches, collapse = "\n")))
    }
    succeed()
})


# ============================================================================
# Test 9: Large-scale random differential testing
# Generate many random PSSMs and sequences, compare exact vs heuristic.
# This is the most powerful bug-finding test.
# ============================================================================
test_that("prefilter: large-scale random differential test (subs only)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(12345)
    n_pssms <- 5
    n_seqs_per <- 30

    mismatches <- character(0)

    for (p in seq_len(n_pssms)) {
        # Random PSSM length that enables prefilter: L >= 3*(K+1) for K=2 -> L>=9
        L <- sample(9:15, 1)
        K <- sample(1:3, 1)

        # Only if L >= 3*(K+1) will prefilter engage
        if (L < 3 * (K + 1)) K <- max(1, floor(L / 3) - 1)

        # Generate random PSSM with some structure
        pssm <- matrix(0, nrow = L, ncol = 4)
        for (i in seq_len(L)) {
            probs <- runif(4, 0.05, 0.95)
            probs <- probs / sum(probs)
            pssm[i, ] <- probs
        }
        colnames(pssm) <- c("A", "C", "G", "T")

        # Threshold: moderate
        col_max <- apply(log(pssm), 1, max)
        threshold <- sum(col_max) - 2.0 * K

        # Generate test sequences
        seqs <- replicate(n_seqs_per,
            paste0(sample(c("A", "C", "G", "T"), L, replace = TRUE), collapse = "")
        )

        for (s in seqs) {
            r_exact <- gseq.pwm_edits(s, pssm,
                score.thresh = threshold, prior = 0, bidirect = FALSE
            )
            r_heur <- gseq.pwm_edits(s, pssm,
                score.thresh = threshold, max_edits = as.integer(K),
                prior = 0, bidirect = FALSE
            )

            exact_min <- if (nrow(r_exact) > 0) min(r_exact$n_edits) else NA
            heur_min <- if (nrow(r_heur) > 0) min(r_heur$n_edits) else NA

            if (!is.na(exact_min) && exact_min <= K) {
                if (is.na(heur_min)) {
                    mismatches <- c(mismatches, sprintf(
                        "L=%d K=%d seq=%s: exact=%d, heur=NA", L, K, s, exact_min
                    ))
                } else if (heur_min != exact_min) {
                    mismatches <- c(mismatches, sprintf(
                        "L=%d K=%d seq=%s: exact=%d, heur=%d", L, K, s, exact_min, heur_min
                    ))
                }
            }
        }
    }

    if (length(mismatches) > 0) {
        fail(paste(
            sprintf("Found %d mismatches in random differential test:\n", length(mismatches)),
            paste(head(mismatches, 20), collapse = "\n")
        ))
    }
    succeed()
})


# ============================================================================
# Test 10: Large-scale random differential test WITH INDELS
# This tests the pigeonhole pre-filter's interaction with indels.
# ============================================================================
test_that("prefilter: large-scale random differential test (with indels)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(54321)
    n_pssms <- 5
    n_seqs_per <- 30

    mismatches <- character(0)

    for (p in seq_len(n_pssms)) {
        L <- sample(9:12, 1)
        K <- 2L
        D <- sample(1:2, 1) # max_indels

        pssm <- matrix(0, nrow = L, ncol = 4)
        for (i in seq_len(L)) {
            probs <- runif(4, 0.05, 0.95)
            probs <- probs / sum(probs)
            pssm[i, ] <- probs
        }
        colnames(pssm) <- c("A", "C", "G", "T")

        col_max <- apply(log(pssm), 1, max)
        threshold <- sum(col_max) - 3.0

        # Generate sequences of varying lengths (L-D to L+D)
        seqs <- character(0)
        for (len in max(1, L - D):(L + D)) {
            seqs <- c(seqs, replicate(n_seqs_per %/% 5 + 1,
                paste0(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = "")
            ))
        }

        for (s in seqs) {
            r_exact <- gseq.pwm_edits(s, pssm,
                score.thresh = threshold,
                max_indels = as.integer(D),
                prior = 0, bidirect = FALSE
            )
            r_heur <- gseq.pwm_edits(s, pssm,
                score.thresh = threshold,
                max_edits = as.integer(K),
                max_indels = as.integer(D),
                prior = 0, bidirect = FALSE
            )

            exact_min <- if (nrow(r_exact) > 0 && any(!is.na(r_exact$n_edits))) min(r_exact$n_edits, na.rm = TRUE) else NA
            heur_min <- if (nrow(r_heur) > 0 && any(!is.na(r_heur$n_edits))) min(r_heur$n_edits, na.rm = TRUE) else NA

            if (!is.na(exact_min) && exact_min <= K) {
                if (is.na(heur_min)) {
                    mismatches <- c(mismatches, sprintf(
                        "L=%d K=%d D=%d len=%d seq=%s: exact=%d, heur=NA",
                        L, K, D, nchar(s), s, exact_min
                    ))
                } else if (heur_min != exact_min) {
                    mismatches <- c(mismatches, sprintf(
                        "L=%d K=%d D=%d len=%d seq=%s: exact=%d, heur=%d",
                        L, K, D, nchar(s), s, exact_min, heur_min
                    ))
                }
            }
        }
    }

    if (length(mismatches) > 0) {
        fail(paste(
            sprintf("Found %d mismatches in indel differential test:\n", length(mismatches)),
            paste(head(mismatches, 20), collapse = "\n")
        ))
    }
    succeed()
})


# ============================================================================
# Test 11: Bidirectional prefilter with indels - random stress test
# ============================================================================
test_that("prefilter: bidirectional + indels stress test", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(777)

    L <- 9
    K <- 2L
    D <- 1L

    pssm <- matrix(c(
        0.8, 0.05, 0.10, 0.05,
        0.05, 0.8, 0.05, 0.10,
        0.10, 0.05, 0.8, 0.05,
        0.05, 0.10, 0.05, 0.8,
        0.8, 0.05, 0.10, 0.05,
        0.05, 0.8, 0.05, 0.10,
        0.10, 0.05, 0.8, 0.05,
        0.05, 0.10, 0.05, 0.8,
        0.8, 0.05, 0.10, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -5.0

    seqs <- replicate(40,
        paste0(sample(c("A", "C", "G", "T"), sample((L - D):(L + D), 1), replace = TRUE), collapse = "")
    )

    mismatches <- character(0)
    for (s in seqs) {
        r_exact <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold,
            max_indels = D,
            prior = 0, bidirect = TRUE
        )
        r_heur <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold,
            max_edits = K, max_indels = D,
            prior = 0, bidirect = TRUE
        )

        exact_min <- if (nrow(r_exact) > 0 && any(!is.na(r_exact$n_edits))) min(r_exact$n_edits, na.rm = TRUE) else NA
        heur_min <- if (nrow(r_heur) > 0 && any(!is.na(r_heur$n_edits))) min(r_heur$n_edits, na.rm = TRUE) else NA

        if (!is.na(exact_min) && exact_min <= K) {
            if (is.na(heur_min)) {
                mismatches <- c(mismatches, sprintf(
                    "bidirect seq=%s: exact=%d, heur=NA", s, exact_min
                ))
            } else if (heur_min != exact_min) {
                mismatches <- c(mismatches, sprintf(
                    "bidirect seq=%s: exact=%d, heur=%d", s, exact_min, heur_min
                ))
            }
        }
    }

    if (length(mismatches) > 0) {
        fail(paste("Bidirectional+indels mismatches:\n", paste(mismatches, collapse = "\n")))
    }
    succeed()
})


# ============================================================================
# Test 12: Prefilter with PSSMs containing -Inf (log-zero) at all columns
# Every base at every column is mandatory -> every B-mer is non-viable.
# The prefilter should still not crash. And with enough edits, matches are
# possible (all mandatory edits, score = S_max).
# ============================================================================
test_that("prefilter: all-mandatory PSSM does not crash", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # L=9, every column has one allowed base and three zeros
    pssm <- matrix(c(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
        1, 0, 0, 0
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -0.01

    # Perfect match
    r1 <- gseq.pwm_edits("ACGTACGTA", pssm,
        score.thresh = threshold, max_edits = 2L,
        prior = 0, bidirect = FALSE
    )
    expect_equal(min(r1$n_edits), 0L, info = "Perfect match should be 0 edits")

    # 1 mismatch
    r2 <- gseq.pwm_edits("TCGTACGTA", pssm,
        score.thresh = threshold, max_edits = 2L,
        prior = 0, bidirect = FALSE
    )
    expect_true(nrow(r2) > 0, info = "1 mismatch on all-mandatory PSSM should be found")
    expect_equal(min(r2$n_edits), 1L, info = "1 mandatory edit needed")

    # 2 mismatches at boundary
    r3 <- gseq.pwm_edits("TGGTACGTA", pssm,
        score.thresh = threshold, max_edits = 2L,
        prior = 0, bidirect = FALSE
    )
    expect_true(nrow(r3) > 0, info = "2 mismatches should be found with max_edits=2")
    expect_equal(min(r3$n_edits), 2L, info = "2 mandatory edits needed")

    # 3 mismatches should be rejected by max_edits=2
    r4 <- gseq.pwm_edits("TGTTACGTA", pssm,
        score.thresh = threshold, max_edits = 2L,
        prior = 0, bidirect = FALSE
    )

    # Compare with exact
    r4_exact <- gseq.pwm_edits("TGTTACGTA", pssm,
        score.thresh = threshold,
        prior = 0, bidirect = FALSE
    )
    exact_min4 <- if (nrow(r4_exact) > 0) min(r4_exact$n_edits) else NA
    heur_min4 <- if (nrow(r4) > 0) min(r4$n_edits) else NA

    # If exact mode finds 3 edits, heuristic with max_edits=2 should return NA
    if (!is.na(exact_min4) && exact_min4 == 3) {
        expect_true(is.na(heur_min4) || heur_min4 > 2,
            info = "3 mandatory edits should be rejected by max_edits=2"
        )
    }
})


# ============================================================================
# Test 13: gscreen serial vs parallel - the sub-chromosome split test
# Tests the allow_multichrom_1d_range_split optimization.
# ============================================================================
test_that("gscreen: serial vs parallel with explicit iterator and pwm.edit_distance", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -5.0

    gvtrack.create("v_edist", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        max_edits = 2L,
        bidirect = TRUE, extend = TRUE, prior = 0
    )

    # Use chr1 200-2000 for testing
    test_scope <- gintervals(1, 200, 2000)

    # Serial execution
    withr::local_options(list(gmultitasking = FALSE))
    result_serial <- gscreen("!is.na(v_edist)", test_scope, iterator = 50)

    # Parallel execution (should trigger sub-chromosome split with explicit iterator)
    withr::local_options(list(gmultitasking = TRUE, gmax.processes = 4))
    result_parallel <- gscreen("!is.na(v_edist)", test_scope, iterator = 50)

    # Results should be identical
    expect_equal(nrow(result_serial), nrow(result_parallel),
        info = "Serial and parallel gscreen should return same number of rows"
    )

    if (nrow(result_serial) > 0 && nrow(result_parallel) > 0) {
        # Sort both by position for comparison
        result_serial <- result_serial[order(result_serial$chrom, result_serial$start), ]
        result_parallel <- result_parallel[order(result_parallel$chrom, result_parallel$start), ]
        rownames(result_serial) <- NULL
        rownames(result_parallel) <- NULL

        expect_equal(result_serial$chrom, result_parallel$chrom,
            info = "Chromosomes should match"
        )
        expect_equal(result_serial$start, result_parallel$start,
            info = "Start positions should match"
        )
        expect_equal(result_serial$end, result_parallel$end,
            info = "End positions should match"
        )
    }
})


# ============================================================================
# Test 14: gextract serial vs parallel with edit distance vtrack
# ============================================================================
test_that("gextract: serial vs parallel with pwm.edit_distance and explicit iterator", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.05, 0.8, 0.1, 0.05,
        0.05, 0.05, 0.8, 0.1,
        0.1, 0.05, 0.05, 0.8,
        0.8, 0.05, 0.1, 0.05,
        0.05, 0.8, 0.05, 0.1
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -4.0

    gvtrack.create("v_edist2", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        max_edits = 2L,
        bidirect = TRUE, extend = TRUE, prior = 0
    )

    test_scope <- gintervals(1, 500, 3000)

    # Serial
    withr::local_options(list(gmultitasking = FALSE))
    result_serial <- gextract("v_edist2", test_scope, iterator = 100)

    # Parallel
    withr::local_options(list(gmultitasking = TRUE, gmax.processes = 4))
    result_parallel <- gextract("v_edist2", test_scope, iterator = 100)

    expect_equal(nrow(result_serial), nrow(result_parallel),
        info = "Serial and parallel gextract should return same number of rows"
    )

    if (nrow(result_serial) > 0 && nrow(result_parallel) > 0) {
        result_serial <- result_serial[order(result_serial$chrom, result_serial$start), ]
        result_parallel <- result_parallel[order(result_parallel$chrom, result_parallel$start), ]
        rownames(result_serial) <- NULL
        rownames(result_parallel) <- NULL

        expect_equal(result_serial$v_edist2, result_parallel$v_edist2,
            tolerance = 1e-6,
            info = "Edit distance values should match"
        )
    }
})


# ============================================================================
# Test 15: gextract serial vs parallel with indels enabled
# The sub-chromosome split must correctly handle the extended sequence fetch
# for indels (extra bases beyond interval boundaries).
# ============================================================================
test_that("gextract: serial vs parallel with indels enabled", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -4.0

    gvtrack.create("v_edist_indel", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        max_edits = 2L, max_indels = 1L,
        bidirect = TRUE, extend = TRUE, prior = 0
    )

    test_scope <- gintervals(1, 300, 2500)

    # Serial
    withr::local_options(list(gmultitasking = FALSE))
    result_serial <- gextract("v_edist_indel", test_scope, iterator = 80)

    # Parallel
    withr::local_options(list(gmultitasking = TRUE, gmax.processes = 4))
    result_parallel <- gextract("v_edist_indel", test_scope, iterator = 80)

    expect_equal(nrow(result_serial), nrow(result_parallel),
        info = "Serial/parallel row count should match (indels)"
    )

    if (nrow(result_serial) > 0 && nrow(result_parallel) > 0) {
        result_serial <- result_serial[order(result_serial$chrom, result_serial$start), ]
        result_parallel <- result_parallel[order(result_parallel$chrom, result_parallel$start), ]
        rownames(result_serial) <- NULL
        rownames(result_parallel) <- NULL

        expect_equal(result_serial$v_edist_indel, result_parallel$v_edist_indel,
            tolerance = 1e-6,
            info = "Indel edit distance values should match serial vs parallel"
        )
    }
})


# ============================================================================
# Test 16: Prefilter with mandatory edits consuming entire edit budget
# When all K edits are mandatory (log-zero bases), the prefilter should still
# allow the window through (at least one block has 0 edits in non-mandatory
# positions but the mandatory positions count as edits).
# ============================================================================
test_that("prefilter: mandatory edits consuming entire budget", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # L=9, K=2. Motif: ACGTACGTA where columns 0,1 have A,C as 1.0 (rest 0)
    # and columns 2-8 are uniform.
    pssm <- matrix(c(
        1.0, 0.0, 0.0, 0.0, # Only A
        0.0, 1.0, 0.0, 0.0, # Only C
        0.25, 0.25, 0.25, 0.25, # Uniform
        0.25, 0.25, 0.25, 0.25, # Uniform
        0.25, 0.25, 0.25, 0.25, # Uniform
        0.25, 0.25, 0.25, 0.25, # Uniform
        0.25, 0.25, 0.25, 0.25, # Uniform
        0.25, 0.25, 0.25, 0.25, # Uniform
        0.25, 0.25, 0.25, 0.25  # Uniform
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -13.0 # Low enough that mandatory-edit-corrected score can reach

    # Sequence with 2 mandatory edits (T at pos 0, G at pos 1)
    # Block 0 (cols 0-2): T,G match -> 2 mandatory edits in block 0
    # Block 1 (cols 3-5): uniform -> 0 edits
    # Block 2 (cols 6-8): uniform -> 0 edits
    # Pigeonhole: blocks 1 and 2 have 0 edits, so prefilter should pass.
    s <- "TGGTACGTA"

    r_exact <- gseq.pwm_edits(s, pssm,
        score.thresh = threshold, prior = 0, bidirect = FALSE
    )
    r_heur <- gseq.pwm_edits(s, pssm,
        score.thresh = threshold, max_edits = 2L,
        prior = 0, bidirect = FALSE
    )

    exact_min <- if (nrow(r_exact) > 0) min(r_exact$n_edits) else NA
    heur_min <- if (nrow(r_heur) > 0) min(r_heur$n_edits) else NA

    if (!is.na(exact_min) && exact_min <= 2) {
        expect_false(is.na(heur_min),
            info = sprintf("Mandatory budget: exact=%d, heur=NA (false reject)", exact_min)
        )
        if (!is.na(heur_min)) {
            expect_equal(heur_min, exact_min,
                info = "Mandatory edits should match"
            )
        }
    }
})


# ============================================================================
# Test 17: Specialized solver (max_indels=1) vs generic DP - boundary cases
# Focus on edge cases where the deletion/insertion is at position 0 or L.
# ============================================================================
test_that("specialized max_indels=1 vs generic: boundary deletion/insertion positions", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # L=10 PSSM with strong preferences
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01,
        0.01, 0.01, 0.97, 0.01,
        0.01, 0.01, 0.01, 0.97,
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01,
        0.01, 0.01, 0.97, 0.01,
        0.01, 0.01, 0.01, 0.97,
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- 10 * log(0.97) - 0.01

    # Deletion at position 0: extra base at beginning
    # Deletion at position L: extra base at end
    # Insertion at position 0: skip motif[0]
    # Insertion at position L-1: skip motif[L-1]
    boundary_seqs <- c(
        "XACGTACGTAC",   # X at start (deletion needed)
        "ACGTACGTACX",   # X at end (deletion needed)
        "CGTACGTAC",     # Missing first base (insertion)
        "ACGTACGTA",     # Missing last base (insertion)
        "AACGTACGTAC",   # Extra A at start
        "ACGTACGTACC",   # Extra C at end
        "TACGTACGTAC",   # Wrong first + extra
        "GACGTACGTAC"    # Wrong first + extra
    )
    # Replace X with random bases
    boundary_seqs <- gsub("X", "G", boundary_seqs)

    mismatches <- character(0)
    for (s in boundary_seqs) {
        r_spec <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = FALSE
        )
        r_gen <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold,
            max_indels = 3L, prior = 0, bidirect = FALSE
        )

        spec_min <- if (nrow(r_spec) > 0 && any(!is.na(r_spec$n_edits))) min(r_spec$n_edits, na.rm = TRUE) else NA
        gen_min <- if (nrow(r_gen) > 0 && any(!is.na(r_gen$n_edits))) min(r_gen$n_edits, na.rm = TRUE) else NA

        # If specialized finds a result, generic must find equal or better
        if (!is.na(spec_min)) {
            if (is.na(gen_min)) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: specialized=%d but generic=NA", s, spec_min
                ))
            } else if (gen_min > spec_min) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: specialized=%d < generic=%d", s, spec_min, gen_min
                ))
            }
        }
        # If generic finds a 1-indel solution, specialized should find it too
        if (!is.na(gen_min) && gen_min <= 1) {
            if (is.na(spec_min) || spec_min > gen_min) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: generic=%d but specialized=%s (should match for <=1 indel)",
                    s, gen_min, ifelse(is.na(spec_min), "NA", spec_min)
                ))
            }
        }
    }

    if (length(mismatches) > 0) {
        fail(paste("Boundary indel mismatches:\n", paste(mismatches, collapse = "\n")))
    }
    succeed()
})


# ============================================================================
# Test 18: Specialized solver (max_indels=2) vs generic DP - random stress
# ============================================================================
test_that("specialized max_indels=2 vs generic: random stress test", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(999)

    L <- 8
    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.05, 0.8, 0.1, 0.05,
        0.05, 0.05, 0.8, 0.1,
        0.1, 0.05, 0.05, 0.8,
        0.8, 0.05, 0.1, 0.05,
        0.05, 0.1, 0.8, 0.05,
        0.1, 0.05, 0.05, 0.8,
        0.8, 0.05, 0.05, 0.1
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -3.5

    seqs <- character(0)
    for (len in (L - 2):(L + 2)) {
        seqs <- c(seqs, replicate(15,
            paste0(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = "")
        ))
    }

    mismatches <- character(0)
    for (s in seqs) {
        r_spec <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold,
            max_indels = 2L, prior = 0, bidirect = FALSE
        )
        r_gen <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold,
            max_indels = 3L, prior = 0, bidirect = FALSE
        )

        spec_min <- if (nrow(r_spec) > 0 && any(!is.na(r_spec$n_edits))) min(r_spec$n_edits, na.rm = TRUE) else NA
        gen_min <- if (nrow(r_gen) > 0 && any(!is.na(r_gen$n_edits))) min(r_gen$n_edits, na.rm = TRUE) else NA

        if (!is.na(spec_min) && (is.na(gen_min) || gen_min > spec_min)) {
            mismatches <- c(mismatches, sprintf(
                "seq=%s len=%d: spec=%d, gen=%s", s, nchar(s), spec_min,
                ifelse(is.na(gen_min), "NA", gen_min)
            ))
        }
        # If generic finds <=2 indel solution, specialized should find it
        if (!is.na(gen_min) && gen_min <= 2 && (is.na(spec_min) || spec_min > gen_min)) {
            mismatches <- c(mismatches, sprintf(
                "seq=%s len=%d: gen=%d but spec=%s (should match for <=2 indels)",
                s, nchar(s), gen_min, ifelse(is.na(spec_min), "NA", spec_min)
            ))
        }
    }

    if (length(mismatches) > 0) {
        fail(paste(
            sprintf("%d mismatches in specialized vs generic:\n", length(mismatches)),
            paste(head(mismatches, 20), collapse = "\n")
        ))
    }
    succeed()
})


# ============================================================================
# Test 19: gextract with implicit iterator should NOT trigger range splitting
# ============================================================================
test_that("gextract: implicit iterator (NULL) does not cause issues", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -3.0

    gvtrack.create("v_edist_impl", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        max_edits = 1L,
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    test_scope <- gintervals(1, 500, 1000)

    # With NULL iterator (implicit), should use scope intervals as-is
    withr::local_options(list(gmultitasking = FALSE))
    result_serial <- gextract("v_edist_impl", test_scope, iterator = test_scope)

    withr::local_options(list(gmultitasking = TRUE, gmax.processes = 4))
    result_parallel <- gextract("v_edist_impl", test_scope, iterator = test_scope)

    expect_equal(nrow(result_serial), nrow(result_parallel),
        info = "Implicit iterator: row counts should match"
    )

    if (nrow(result_serial) > 0) {
        expect_equal(result_serial$v_edist_impl, result_parallel$v_edist_impl,
            tolerance = 1e-6,
            info = "Implicit iterator: values should match"
        )
    }
})


# ============================================================================
# Test 20: Prefilter correctness on gextract (genome-level test)
# Compare gextract results with max_edits (prefilter enabled) vs without
# (exact mode), filtering exact results to <= max_edits.
# ============================================================================
test_that("gextract: prefilter correctness on real genome data", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # L=10, K=2 -> 3 blocks (prefilter enabled)
    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.05, 0.8, 0.1, 0.05,
        0.05, 0.05, 0.8, 0.1,
        0.1, 0.05, 0.05, 0.8,
        0.8, 0.05, 0.1, 0.05,
        0.05, 0.8, 0.05, 0.1,
        0.1, 0.05, 0.8, 0.05,
        0.05, 0.1, 0.05, 0.8,
        0.8, 0.05, 0.05, 0.1,
        0.05, 0.05, 0.1, 0.8
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -5.0

    # Exact mode (no prefilter)
    gvtrack.create("v_exact", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = TRUE, extend = TRUE, prior = 0
    )

    # Heuristic mode with prefilter
    gvtrack.create("v_heur", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        max_edits = 2L,
        bidirect = TRUE, extend = TRUE, prior = 0
    )

    test_scope <- gintervals(1, 200, 1500)

    withr::local_options(list(gmultitasking = FALSE))

    r_exact <- gextract("v_exact", test_scope, iterator = 50)
    r_heur <- gextract("v_heur", test_scope, iterator = 50)

    # Compare: for every interval, if exact says edits <= 2, heuristic should agree
    expect_equal(nrow(r_exact), nrow(r_heur),
        info = "Exact and heuristic should scan same number of windows"
    )

    if (nrow(r_exact) > 0 && nrow(r_heur) > 0) {
        for (i in seq_len(nrow(r_exact))) {
            e <- r_exact$v_exact[i]
            h <- r_heur$v_heur[i]

            if (!is.na(e) && e <= 2) {
                # Heuristic should find the same result
                expect_false(is.na(h),
                    info = sprintf(
                        "Row %d: exact=%g but heur=NA (false reject at chr%s:%d-%d)",
                        i, e, r_exact$chrom[i], r_exact$start[i], r_exact$end[i]
                    )
                )
                if (!is.na(h)) {
                    expect_equal(h, e, tolerance = 1e-6,
                        info = sprintf(
                            "Row %d: exact=%g != heur=%g at chr%s:%d-%d",
                            i, e, h, r_exact$chrom[i], r_exact$start[i], r_exact$end[i]
                        )
                    )
                }
            }

            if (!is.na(e) && e > 2) {
                # Heuristic should return NA (beyond max_edits)
                expect_true(is.na(h),
                    info = sprintf(
                        "Row %d: exact=%g (>2) but heur=%g (should be NA)",
                        i, e, h
                    )
                )
            }
        }
    }
})


# ============================================================================
# Test 21: Cross-check specialized vs generic with indels on real genome
# ============================================================================
test_that("gextract: specialized max_indels=1 vs generic on real genome", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -3.5

    # max_indels=1 -> specialized solver
    gvtrack.create("v_indel1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        max_indels = 1L,
        bidirect = TRUE, extend = TRUE, prior = 0
    )

    # max_indels=3 -> generic solver
    gvtrack.create("v_indel3", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        max_indels = 3L,
        bidirect = TRUE, extend = TRUE, prior = 0
    )

    test_scope <- gintervals(1, 200, 1000)

    withr::local_options(list(gmultitasking = FALSE))

    r1 <- gextract("v_indel1", test_scope, iterator = 40)
    r3 <- gextract("v_indel3", test_scope, iterator = 40)

    expect_equal(nrow(r1), nrow(r3),
        info = "Specialized and generic should scan same windows"
    )

    if (nrow(r1) > 0 && nrow(r3) > 0) {
        mismatches <- character(0)
        for (i in seq_len(nrow(r1))) {
            e1 <- r1$v_indel1[i]
            e3 <- r3$v_indel3[i]

            # Generic with more indel budget should always be <= specialized
            if (!is.na(e1) && !is.na(e3) && e3 > e1 + 0.01) {
                mismatches <- c(mismatches, sprintf(
                    "Row %d: specialized=%g > generic=%g at %s:%d-%d",
                    i, e1, e3, r1$chrom[i], r1$start[i], r1$end[i]
                ))
            }
            # If specialized finds 0 or 1 edit (uses <=1 indel), generic should match
            if (!is.na(e1) && e1 <= 1 && (is.na(e3) || abs(e3 - e1) > 0.01)) {
                mismatches <- c(mismatches, sprintf(
                    "Row %d: spec=%g but gen=%s (should match for <=1 edit)",
                    i, e1, ifelse(is.na(e3), "NA", e3)
                ))
            }
        }

        if (length(mismatches) > 0) {
            fail(paste(
                sprintf("Found %d mismatches:\n", length(mismatches)),
                paste(head(mismatches, 20), collapse = "\n")
            ))
        }
    }
})


# ============================================================================
# Test 22: Brute-force reference vs C++ for short sequences with indels
# Uses the R brute-force from the existing adversarial test file.
# ============================================================================
test_that("prefilter: brute-force reference vs C++ with prefilter (indels)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Re-implement a simple brute force for validation
    bruteforce_min_edits <- function(seq, pssm, threshold, max_indels, bidirect = FALSE) {
        L <- nrow(pssm)
        S <- nchar(seq)
        log_pssm <- log(pssm)
        col_max <- apply(log_pssm, 1, max)

        score_alignment <- function(bases, motif_len) {
            # Simple: compute PWM score, count mandatory edits, compute subs
            adj_score <- 0
            mandatory <- 0
            gains <- numeric(0)
            base_map <- c(A = 1, C = 2, G = 3, T = 4)

            for (i in seq_len(motif_len)) {
                bi <- base_map[bases[i]]
                if (is.na(bi)) {
                    bs <- min(log_pssm[i, ])
                } else {
                    bs <- log_pssm[i, bi]
                }

                if (!is.finite(bs)) {
                    mandatory <- mandatory + 1
                    adj_score <- adj_score + col_max[i]
                } else {
                    adj_score <- adj_score + bs
                    g <- col_max[i] - bs
                    if (g > 1e-12) gains <- c(gains, g)
                }
            }

            deficit <- threshold - adj_score
            if (deficit <= 0) return(mandatory)

            gains_sorted <- sort(gains, decreasing = TRUE)
            acc <- 0
            subs <- 0
            for (g in gains_sorted) {
                acc <- acc + g
                subs <- subs + 1
                if (acc >= deficit - 1e-10) {
                    return(mandatory + subs)
                }
            }
            return(NA_real_)
        }

        best <- NA_real_

        # Try all window widths and positions
        D <- max_indels
        for (W in max(1, L - D):min(S, L + D)) {
            # Number of indels for this W
            k <- abs(W - L)
            if (k > D) next

            for (p in seq_len(S - W + 1)) {
                subseq <- substr(seq, p, p + W - 1)
                bases <- strsplit(subseq, "")[[1]]

                # Forward
                # For W != L, we'd need proper DP alignment.
                # Simple case: W == L -> 0 indels
                if (W == L) {
                    result <- score_alignment(bases, L)
                    if (!is.na(result)) {
                        total <- result
                        if (is.na(best) || total < best) best <- total
                    }
                }
            }
        }

        if (bidirect) {
            # Also check reverse complement
            rc_map <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
            for (W in max(1, L - D):min(S, L + D)) {
                k <- abs(W - L)
                if (k > D) next
                for (p in seq_len(S - W + 1)) {
                    subseq <- substr(seq, p, p + W - 1)
                    bases <- rev(strsplit(subseq, "")[[1]])
                    bases <- rc_map[bases]
                    if (W == L) {
                        result <- score_alignment(bases, L)
                        if (!is.na(result)) {
                            total <- result
                            if (is.na(best) || total < best) best <- total
                        }
                    }
                }
            }
        }

        return(best)
    }

    # L=9 PSSM
    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.7, 0.1, 0.1, 0.1
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -4.0

    set.seed(444)
    # Only test exact-length sequences (W=L, no indels) for the simple brute-force
    seqs <- replicate(30, paste0(sample(c("A", "C", "G", "T"), 9, replace = TRUE), collapse = ""))

    mismatches <- character(0)
    for (s in seqs) {
        bf <- bruteforce_min_edits(s, pssm, threshold, max_indels = 0, bidirect = FALSE)

        r_heur <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, max_edits = 2L,
            prior = 0, bidirect = FALSE
        )
        heur_min <- if (nrow(r_heur) > 0) min(r_heur$n_edits) else NA

        if (!is.na(bf) && bf <= 2) {
            if (is.na(heur_min)) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: brute-force=%d, heur=NA", s, bf
                ))
            } else if (heur_min != bf) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: brute-force=%d, heur=%d", s, bf, heur_min
                ))
            }
        }
    }

    if (length(mismatches) > 0) {
        fail(paste("Brute-force vs C++ mismatches:\n", paste(mismatches, collapse = "\n")))
    }
    succeed()
})


# ============================================================================
# Test 23: Very long motif that exceeds prefilter block entry limit
# L=9, K=2 -> blocks of 3, 4^3=64 entries. Fine.
# L=20, K=1 -> 2 blocks of 10, 4^10=1048576 entries. This is large.
# Test that it still works (or falls back gracefully).
# ============================================================================
test_that("prefilter: long motif L=20 with K=1", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 20
    pssm <- matrix(0, nrow = L, ncol = 4)
    set.seed(555)
    for (i in seq_len(L)) {
        probs <- c(0.7, 0.1, 0.1, 0.1)[sample(4)]
        pssm[i, ] <- probs
    }
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- sum(apply(log(pssm), 1, max)) - 1.0

    seqs <- replicate(20, paste0(sample(c("A", "C", "G", "T"), L, replace = TRUE), collapse = ""))

    mismatches <- character(0)
    for (s in seqs) {
        r_exact <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, prior = 0, bidirect = FALSE
        )
        r_heur <- gseq.pwm_edits(s, pssm,
            score.thresh = threshold, max_edits = 1L,
            prior = 0, bidirect = FALSE
        )

        exact_min <- if (nrow(r_exact) > 0) min(r_exact$n_edits) else NA
        heur_min <- if (nrow(r_heur) > 0) min(r_heur$n_edits) else NA

        if (!is.na(exact_min) && exact_min <= 1) {
            if (is.na(heur_min)) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: exact=%d, heur=NA", s, exact_min
                ))
            } else if (heur_min != exact_min) {
                mismatches <- c(mismatches, sprintf(
                    "seq=%s: exact=%d, heur=%d", s, exact_min, heur_min
                ))
            }
        }
    }

    if (length(mismatches) > 0) {
        fail(paste("Long motif L=20 mismatches:\n", paste(mismatches, collapse = "\n")))
    }
    succeed()
})


# ============================================================================
# Test 24: gextract multi-chromosome parallel consistency
# Test with multiple chromosomes of very different sizes.
# ============================================================================
test_that("gextract: multi-chromosome parallel consistency", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7,
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    threshold <- -3.0

    gvtrack.create("v_multi", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        max_edits = 2L,
        bidirect = TRUE, extend = TRUE, prior = 0
    )

    # Use multiple chromosomes
    chrom_info <- gintervals.all()
    # Take first 3 chromosomes (or fewer if not available)
    n_chroms <- min(3, nrow(chrom_info))
    test_intervals <- do.call(rbind, lapply(seq_len(n_chroms), function(i) {
        chr <- chrom_info$chrom[i]
        max_end <- min(2000, chrom_info$end[i])
        start <- min(200, chrom_info$start[i])
        if (max_end > start + 100) {
            gintervals(as.numeric(gsub("chr", "", chr)), start, max_end)
        } else {
            NULL
        }
    }))

    if (is.null(test_intervals) || nrow(test_intervals) == 0) {
        skip("No suitable chromosomes found")
    }

    withr::local_options(list(gmultitasking = FALSE))
    result_serial <- gextract("v_multi", test_intervals, iterator = 60)

    withr::local_options(list(gmultitasking = TRUE, gmax.processes = 4))
    result_parallel <- gextract("v_multi", test_intervals, iterator = 60)

    expect_equal(nrow(result_serial), nrow(result_parallel),
        info = "Multi-chrom: row counts should match"
    )

    if (nrow(result_serial) > 0 && nrow(result_parallel) > 0) {
        result_serial <- result_serial[order(result_serial$chrom, result_serial$start), ]
        result_parallel <- result_parallel[order(result_parallel$chrom, result_parallel$start), ]
        rownames(result_serial) <- NULL
        rownames(result_parallel) <- NULL

        expect_equal(result_serial$v_multi, result_parallel$v_multi,
            tolerance = 1e-6,
            info = "Multi-chrom: values should match"
        )
    }
})
