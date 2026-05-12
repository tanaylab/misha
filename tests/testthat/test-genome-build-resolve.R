test_that("pattern fallback for GC[FA]_* synthesizes source: ucsc-hub", {
    res <- .resolve_genome("GCA_004023825.1", registry = NULL)
    expect_equal(res$recipe$source, "ucsc-hub")
    expect_equal(res$recipe$accession, "GCA_004023825.1")
    expect_match(res$resolved_from, "pattern fallback")
})

test_that("registry entries with source: ncbi still resolve to ncbi", {
    yaml_path <- tempfile(fileext = ".yaml")
    on.exit(unlink(yaml_path))
    writeLines(c(
        "genome:",
        "  my_acc:",
        "    source: ncbi",
        "    accession: GCA_004023825.1"
    ), yaml_path)
    res <- .resolve_genome("my_acc", registry = yaml_path)
    expect_equal(res$recipe$source, "ncbi")
})
