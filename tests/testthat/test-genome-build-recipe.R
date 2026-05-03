test_that(".validate_recipe accepts ucsc-hub with accession", {
    r <- list(source = "ucsc-hub", accession = "GCA_004023825.1")
    expect_silent(.validate_recipe(r, "x"))
})

test_that(".validate_recipe rejects ucsc-hub without accession", {
    r <- list(source = "ucsc-hub")
    expect_error(.validate_recipe(r, "x"), "accession")
})

test_that(".validate_recipe rejects malformed accession", {
    r <- list(source = "ucsc-hub", accession = "GCAfoo")
    expect_error(.validate_recipe(r, "x"), "accession")
})

test_that(".validate_recipe accepts optional gtf_priority", {
    r <- list(
        source = "ucsc-hub", accession = "GCA_004023825.1",
        gtf_priority = c("ncbiRefSeq", "ensGene")
    )
    expect_silent(.validate_recipe(r, "x"))
})
