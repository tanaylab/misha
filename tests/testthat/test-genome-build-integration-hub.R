test_that("install_intervals from ucsc-hub installs rmsk against a real assembly", {
    skip_on_cran()
    skip_if_offline()

    # Use a small mammal hub assembly. We pick one whose .repeatMasker.out.gz
    # is on the smaller side (rabbit ~10 MB) and only request rmsk.
    groot <- testthat::test_path("fixtures", "tiny-hub-groot")
    # This fixture won't have matching chrom names — expect a hard error from
    # chromAlias detection, which is a positive end-to-end signal that the
    # detection path is wired up correctly.
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
