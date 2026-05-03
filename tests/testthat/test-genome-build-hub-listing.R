test_that(".hub_list_dir parses Apache index from a captured fixture", {
    f <- testthat::test_path("fixtures", "hub-index.html")
    body <- readLines(f, warn = FALSE)
    body <- paste(body, collapse = "\n")
    files <- .hub_list_dir_parse(body)
    expect_true(any(grepl("\\.chromAlias\\.txt", files)))
    expect_true(any(grepl("\\.repeatMasker\\.out\\.gz", files)))
    expect_true(any(grepl("\\.fa\\.gz$|\\.2bit$", files)))
    # Does not include sort/parent-dir noise.
    expect_false(any(startsWith(files, "?")))
    expect_false(any(files == "../"))
    expect_false(any(files == "/"))
})

test_that(".hub_list_dir returns NULL on non-existent URL (404)", {
    skip_on_cran()
    skip_if_offline()
    expect_null(.hub_list_dir(
        "https://hgdownload.soe.ucsc.edu/hubs/GCA/999/999/999/GCA_999999999.1/",
        verbose = FALSE
    ))
})
