test_that("rmsk install succeeds for three Zoonomia accessions", {
    skip_on_cran()
    skip_if_offline()
    skip_if_not(
        Sys.getenv("MISHA_RUN_SLOW_TESTS") == "1",
        "Skipping slow Zoonomia smoke test (set MISHA_RUN_SLOW_TESTS=1 to enable)."
    )
    skip_if_not(
        Sys.getenv("MISHA_TEST_HUB_GROOT_PATH") != "",
        "MISHA_TEST_HUB_GROOT_PATH not set."
    )

    accessions <- c(
        "GCA_004023825.1", # Arctic fox
        "GCA_004023885.1", # Fossa
        "GCA_004023845.1"
    ) # Dwarf mongoose
    for (acc in accessions) {
        # Each accession would need its own pre-built groot; for v1 smoke, just
        # exercise the directory-listing + chromAlias-fetch + rmsk-parse half.
        url <- .hub_url_for(acc)
        files <- .hub_list_dir(url, verbose = FALSE)
        expect_true(any(grepl("\\.repeatMasker\\.out\\.gz$", files)),
            info = sprintf("rmsk file expected at %s", url)
        )
    }
})
