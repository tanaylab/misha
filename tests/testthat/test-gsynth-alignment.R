# Regression tests for misha#94: gsynth.sample must honor model$iterator when
# intervals are not aligned to the iterator bin boundary, instead of silently
# inferring iter_size from the first same-chrom diff in iter_starts.

# --- shared helpers -----------------------------------------------------------

.gsa_train_small_model <- function(intervals, iterator = 200) {
    # Caller is responsible for removing "test_vt" (the vtrack is referenced
    # again at sample time via model$dim_specs).
    if ("test_vt" %in% gvtrack.ls()) gvtrack.rm("test_vt")
    gvtrack.create("test_vt", "dense_track", "avg")
    rng <- gsummary("dense_track", intervals = intervals)
    gsynth.train(
        list(
            expr = "test_vt",
            breaks = seq(rng["Min"], rng["Max"], length.out = 6)
        ),
        intervals = intervals,
        iterator = iterator
    )
}

.gsa_force_emit_A <- function(model) {
    # Make every CDF row emit base 0 (A) deterministically. The sampler picks
    # base b as the smallest b with unif_rand() < cdf[b]; setting every cell to
    # 1.0 means base 0 always wins.
    for (b in seq_along(model$model_data$cdf)) {
        model$model_data$cdf[[b]][] <- 1.0
    }
    model
}

# --- tests --------------------------------------------------------------------

test_that("gsynth.sample honors model$iterator on unaligned intervals (#94)", {
    gdb.init_examples()

    train_intervals <- gintervals(1, 0, 50000)
    model <- .gsa_train_small_model(train_intervals, iterator = 200)
    model <- .gsa_force_emit_A(model)

    # Unaligned to iterator=200: starts at 64, not a multiple of 200. Before
    # the fix, iter_size was inferred as (first same-chrom iter_start diff),
    # which for this interval is 136 — positions 336..399, 536..599 would
    # fall through to uniform-random sampling and emit non-A bases.
    unaligned <- gintervals(1, 64, 664)
    seqs <- gsynth.sample(
        model,
        output_format = "vector",
        intervals = unaligned,
        seed = 60427
    )

    expect_length(seqs, 1L)
    s <- seqs[[1]]
    expect_equal(nchar(s), 600L)

    # Drop the first k=5 seeded positions (those are uniform random).
    post_seed <- substring(s, 6)
    # With the bug active, expect ~128/595 non-A bases (~16%). With the fix,
    # zero non-A bases.
    non_A <- sum(strsplit(post_seed, "")[[1]] != "A")
    expect_equal(non_A, 0L)

    gvtrack.rm("test_vt")
})

test_that("aligned and unaligned intervals produce matching forbidden-kmer stats", {
    gdb.init_examples()

    train_intervals <- gintervals(1, 0, 50000)
    model <- .gsa_train_small_model(train_intervals, iterator = 200)

    # Forbid "CG": zero every cdf cell where (state ends in C) AND (next=G).
    # State row r (0..1023) encodes a 5-mer; its last base is r %% 4.
    # Base C=1, G=2. After zeroing, renormalize the state's row.
    state_ends_in_C <- (seq_len(1024) - 1L) %% 4L == 1L
    for (b in seq_along(model$model_data$cdf)) {
        cdf <- model$model_data$cdf[[b]]
        # Recover per-row probabilities from the cumulative.
        probs <- cbind(cdf[, 1], cdf[, 2] - cdf[, 1], cdf[, 3] - cdf[, 2], 1 - cdf[, 3])
        probs[state_ends_in_C, 3L] <- 0  # C -> G blocked
        rs <- rowSums(probs)
        nz <- rs > 0
        probs[nz, ] <- probs[nz, ] / rs[nz]
        new_cdf <- t(apply(probs, 1L, cumsum))
        new_cdf[, 4L] <- 1
        model$model_data$cdf[[b]] <- new_cdf
    }

    aligned <- gintervals(1, 0, 2000)
    unaligned <- gintervals(1, 64, 2064)

    seq_a <- gsynth.sample(model, output_format = "vector",
                           intervals = aligned, seed = 60427)[[1]]
    seq_u <- gsynth.sample(model, output_format = "vector",
                           intervals = unaligned, seed = 60427)[[1]]

    # Ignore the first k=5 seeded bases in each (they are uniform random so
    # can form "CG" across the seeding boundary).
    tail_a <- substring(seq_a, 6)
    tail_u <- substring(seq_u, 6)

    cg_count <- function(s) {
        m <- gregexpr("CG", s, fixed = TRUE)[[1]]
        if (length(m) == 1L && m == -1L) 0L else length(m)
    }

    # Both must be exactly zero beyond seeding. Before the fix, the unaligned
    # interval picks up CG bigrams from the uniform-fallback regions.
    expect_equal(cg_count(tail_a), 0L)
    expect_equal(cg_count(tail_u), 0L)

    gvtrack.rm("test_vt")
})

test_that("gsynth.sample errors clearly on non-positive model$iterator", {
    gdb.init_examples()

    train_intervals <- gintervals(1, 0, 50000)
    model <- .gsa_train_small_model(train_intervals, iterator = 200)

    # Zero or negative iterator is a malformed model. The error can come from
    # gextract (iterator must be positive) or from the C validator — either is
    # fine, both signal a bad model loud and early.
    model$iterator <- 0L
    expect_error(
        gsynth.sample(model, output_format = "vector",
                      intervals = gintervals(1, 0, 200), seed = 1)
    )

    gvtrack.rm("test_vt")
})
