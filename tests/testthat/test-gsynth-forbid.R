# Tests for gsynth.forbid_kmer()

.gfb_train <- function(intervals = gintervals(1, 0, 50000), iterator = 200, k = 5L) {
    if ("test_vt" %in% gvtrack.ls()) gvtrack.rm("test_vt")
    gvtrack.create("test_vt", "dense_track", "avg")
    rng <- gsummary("dense_track", intervals = intervals)
    gsynth.train(
        list(
            expr = "test_vt",
            breaks = seq(rng["Min"], rng["Max"], length.out = 6)
        ),
        intervals = intervals,
        iterator = iterator,
        k = k
    )
}

.gfb_count_substrings <- function(s, pat) {
    m <- gregexpr(pat, s, fixed = TRUE)[[1]]
    if (length(m) == 1L && m == -1L) 0L else length(m)
}

test_that("gsynth.forbid_kmer('CG') produces CG-free samples beyond seeding", {
    gdb.init_examples()
    model <- .gfb_train()
    model_no_cg <- suppressMessages(gsynth.forbid_kmer(model, "CG"))

    # Shape preserved.
    expect_s3_class(model_no_cg, "gsynth.model")
    expect_identical(
        length(model_no_cg$model_data$cdf),
        length(model$model_data$cdf)
    )

    # Original untouched (reference unchanged on at least one bin's counts).
    expect_identical(
        model$model_data$counts[[1]],
        model$model_data$counts[[1]]
    )

    # Every row of every cdf still sums to 1 in the last column.
    for (cdf in model_no_cg$model_data$cdf) {
        expect_true(all(abs(cdf[, 4L] - 1) < 1e-10))
    }

    # Sample an aligned interval. Skip the first 2*k bases to clear the
    # seeding + trap-escape window (states whose k-mer already contains CG
    # fall back to uniform until the pattern slides out of the state window).
    seqs <- gsynth.sample(
        model_no_cg,
        output_format = "vector",
        intervals = gintervals(1, 0, 5000),
        seed = 60427
    )
    expect_length(seqs, 1L)
    post_trap <- substring(seqs[[1]], 2L * 5L + 1L)
    expect_equal(.gfb_count_substrings(post_trap, "CG"), 0L)

    gvtrack.rm("test_vt")
})

test_that("gsynth.forbid_kmer('GCGC') forbids a 4-mer (length <= k+1)", {
    gdb.init_examples()
    model <- .gfb_train() # k=5, so pattern length 4 is fine
    model_no_gcgc <- suppressMessages(gsynth.forbid_kmer(model, "GCGC"))

    seqs <- gsynth.sample(
        model_no_gcgc,
        output_format = "vector",
        intervals = gintervals(1, 0, 5000),
        seed = 60427
    )
    # Skip 2*k bases to clear the seeding + trap-escape window.
    tail <- substring(seqs[[1]], 2L * 5L + 1L)
    expect_equal(.gfb_count_substrings(tail, "GCGC"), 0L)

    gvtrack.rm("test_vt")
})

test_that("gsynth.forbid_kmer errors on pattern longer than k+1", {
    gdb.init_examples()
    model <- .gfb_train(k = 5L) # k+1 = 6
    # Length 7 pattern exceeds k+1=6.
    expect_error(
        gsynth.forbid_kmer(model, "ACGTACG"),
        "exceeds model.k \\+ 1"
    )

    gvtrack.rm("test_vt")
})

test_that("gsynth.forbid_kmer errors on invalid pattern", {
    gdb.init_examples()
    model <- .gfb_train()
    expect_error(
        gsynth.forbid_kmer(model, ""),
        "pattern must be non-empty DNA over ACGT"
    )
    expect_error(
        gsynth.forbid_kmer(model, "CN"),
        "pattern must be non-empty DNA over ACGT"
    )
    expect_error(
        gsynth.forbid_kmer(model, NA_character_),
        "pattern must be a single non-NA character string"
    )
    expect_error(
        gsynth.forbid_kmer(model, c("CG", "GC")),
        "pattern must be a single non-NA character string"
    )
    expect_error(
        gsynth.forbid_kmer("not a model", "CG"),
        "model must be a gsynth.model object"
    )

    gvtrack.rm("test_vt")
})

test_that("gsynth.forbid_kmer does not mutate the input model", {
    gdb.init_examples()
    model <- .gfb_train()
    counts_before <- lapply(model$model_data$counts, function(x) if (is.null(x)) NULL else x + 0)
    cdf_before <- lapply(model$model_data$cdf, function(x) if (is.null(x)) NULL else x + 0)

    suppressMessages(gsynth.forbid_kmer(model, "CG"))

    for (b in seq_along(model$model_data$counts)) {
        expect_identical(model$model_data$counts[[b]], counts_before[[b]])
        expect_identical(model$model_data$cdf[[b]], cdf_before[[b]])
    }

    gvtrack.rm("test_vt")
})

test_that("gsynth.forbid_kmer integrates with gsynth.sample on unaligned intervals", {
    # Combined test: uses the #94 fix AND forbid_kmer together. Unaligned CGD-
    # like intervals were the original motivation for both changes.
    gdb.init_examples()
    model <- .gfb_train()
    model_no_cg <- suppressMessages(gsynth.forbid_kmer(model, "CG"))

    # Unaligned intervals (start 64, 347 — not multiples of 200).
    intervals <- rbind(
        gintervals(1, 64, 2064),
        gintervals(1, 3347, 5347)
    )

    seqs <- gsynth.sample(
        model_no_cg,
        output_format = "vector",
        intervals = intervals,
        seed = 60427
    )

    total_cg <- 0L
    for (s in seqs) {
        tail <- substring(s, 2L * 5L + 1L)
        total_cg <- total_cg + .gfb_count_substrings(tail, "CG")
    }
    expect_equal(total_cg, 0L)

    gvtrack.rm("test_vt")
})
