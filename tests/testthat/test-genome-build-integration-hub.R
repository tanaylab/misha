test_that("install_intervals from ucsc-hub installs rmsk against a real assembly", {
    skip_on_cran()
    skip_if_offline()
    skip_if_not(
        Sys.getenv("MISHA_RUN_SLOW_TESTS") == "1",
        "Skipping slow ucsc-hub install (downloads ~152MB; set MISHA_RUN_SLOW_TESTS=1 to enable)."
    )

    # Hits a real UCSC hub: downloads chromAlias (~190KB) and the assembly's
    # .repeatMasker.out.gz (~152MB for rabbit GCF_000003625.3) before failing
    # chromAlias detection. Gated because the download is heavy for CI.
    groot <- testthat::test_path("fixtures", "tiny-hub-groot")
    # The fixture's chrom names (SYNTHETIC_CHR_*) won't appear in any chromAlias
    # column, so we expect a hard error — that is a positive end-to-end signal
    # that the detection path is wired up correctly.
    expect_error(
        gdb.install_intervals(
            groot = groot,
            source = list(source = "ucsc-hub", accession = "GCF_000003625.3"),
            sets = "rmsk",
            verbose = FALSE
        ),
        "100%"
    )
})

test_that("install_intervals from ucsc-hub succeeds when chromAlias matches", {
    skip_on_cran()
    skip_if_offline()
    skip_if_not(
        Sys.getenv("MISHA_TEST_HUB_GROOT_PATH") != "",
        "MISHA_TEST_HUB_GROOT_PATH not set; skipping ucsc-hub e2e test."
    )

    # The user provides a path to a real groot built from a hub assembly's
    # FASTA. We then install rmsk onto it.
    groot <- Sys.getenv("MISHA_TEST_HUB_GROOT_PATH")
    accession <- Sys.getenv("MISHA_TEST_HUB_ACCESSION", "GCF_000003625.3")
    gdb.install_intervals(
        groot = groot,
        source = list(source = "ucsc-hub", accession = accession),
        sets = "rmsk",
        prefix = "test_install.",
        verbose = FALSE
    )
    expect_true(gintervals.exists("test_install.rmsk"))
})
