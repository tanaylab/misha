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

test_that(".validate_recipe accepts arbitrary chrom_naming for ucsc-hub (validated at build time)", {
    # Pre-fix this rejected anything outside ucsc/accession/sequence_name.
    # Now we accept any non-empty string; mismatches against alias columns
    # surface at build time with an actionable error.
    r <- list(
        source = "ucsc-hub", accession = "GCA_004023825.1",
        chrom_naming = "genbank"
    )
    expect_silent(.validate_recipe(r, "x"))
})

test_that(".validate_recipe rejects empty / non-string chrom_naming for ucsc-hub", {
    r <- list(
        source = "ucsc-hub", accession = "GCA_004023825.1",
        chrom_naming = ""
    )
    expect_error(.validate_recipe(r, "x"), "chrom_naming")
    r2 <- list(
        source = "ucsc-hub", accession = "GCA_004023825.1",
        chrom_naming = c("ucsc", "genbank")
    )
    expect_error(.validate_recipe(r2, "x"), "chrom_naming")
})
