test_that("gsynth.train resolves chromids via full chromkey, not input subset", {
    # Regression test for positional chromid bug in R/synth.R.
    #
    # BUG: chrom_ids/iter_chroms were computed via
    #   match(as.character(intervals$chrom), chrom_sizes$chrom) - 1L
    # where `chrom_sizes$chrom` (as coerced to character by match()) only
    # contains chromosomes present in the input subset. When the input is
    # missing a chromosome that sorts earlier in the chromkey, every later
    # chromosome's ID shifts down, and the C++ code reads sequences from the
    # WRONG chromosome (silently).
    #
    # FIX: match against `levels(chrom_sizes$chrom)`, which the C++ side
    # populates with the FULL chromkey in ID order (see
    # GenomeIntervalUtils.cpp::gintervals_chrom_sizes).
    #
    # Test scenario: the test DB has chr1, chr2, chrX. Passing intervals on
    # chrX only excludes chr1 and chr2 from the input subset. Under the bug,
    # chrX-input chromid resolves to 0 (chr1's chromkey ID), so gsynth.train
    # reads chr1 sequence instead of chrX, and the resulting model matches
    # a chr1-only-input model. After the fix, the two models differ.

    gdb.init_examples()

    chr1_intervs <- gintervals("chr1", 0, 100000)
    chrX_intervs <- gintervals("chrX", 0, 100000)

    # Sanity: chr1 and chrX carry different sequences in the test DB
    expect_false(
        gseq.extract(gintervals("chr1", 0, 1000)) ==
            gseq.extract(gintervals("chrX", 0, 1000))
    )

    # Zero-dim model: no dim_specs, so gsynth.train just counts k-mers
    # across the input intervals. This isolates the chromid resolution
    # from track-extraction logic.
    model_chr1 <- gsynth.train(
        intervals = chr1_intervs,
        iterator = 100
    )
    model_chrX <- gsynth.train(
        intervals = chrX_intervs,
        iterator = 100
    )

    # Both models should be valid
    expect_s3_class(model_chr1, "gsynth.model")
    expect_s3_class(model_chrX, "gsynth.model")

    # Under the bug: chrX-input loads chr1's sequence → the two CDFs match.
    # After the fix: each loads its own chromosome → the CDFs differ.
    cdf_chr1 <- model_chr1$model_data$cdf[[1]]
    cdf_chrX <- model_chrX$model_data$cdf[[1]]
    expect_false(isTRUE(all.equal(cdf_chr1, cdf_chrX, tolerance = 1e-6)))
})
