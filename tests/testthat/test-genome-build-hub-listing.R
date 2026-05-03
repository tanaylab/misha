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

test_that(".hub_url_for builds correct path from accession", {
    expect_equal(
        .hub_url_for("GCA_004023825.1"),
        "https://hgdownload.soe.ucsc.edu/hubs/GCA/004/023/825/GCA_004023825.1/"
    )
    expect_equal(
        .hub_url_for("GCF_000003625.3"),
        "https://hgdownload.soe.ucsc.edu/hubs/GCF/000/003/625/GCF_000003625.3/"
    )
})

test_that(".pick_gtf returns first matching GTF in priority order", {
    files <- c(
        "genes/GCA_X.augustus.gtf.gz",
        "genes/GCA_X.ncbiRefSeq.gtf.gz",
        "genes/GCA_X.ensGene.gtf.gz"
    )
    res <- .pick_gtf(files,
        priority = c("ncbiRefSeq", "bestRefSeq", "ensGene", "augustus", "xenoRefGene")
    )
    expect_equal(res$source, "ncbiRefSeq")
    expect_match(res$file, "ncbiRefSeq")
})

test_that(".pick_gtf returns NULL when nothing matches", {
    files <- c("genes/GCA_X.unknown.gtf.gz")
    res <- .pick_gtf(files,
        priority = c("ncbiRefSeq", "ensGene")
    )
    expect_null(res)
})
