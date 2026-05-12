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

test_that("gdb.build_genome accepts target_chroms and validates type", {
    d <- tempfile()
    dir.create(d)
    on.exit(unlink(d, recursive = TRUE))

    expect_error(
        gdb.build_genome("hg38", path = d, target_chroms = 1:5),
        "target_chroms"
    )
    expect_error(
        gdb.build_genome("hg38", path = d, target_chroms = c("chr1", NA)),
        "target_chroms"
    )
    # Valid: char vector. Existing-path error wins.
    expect_error(
        gdb.build_genome("hg38", path = d, target_chroms = c("chr1", "chr2")),
        "already exists"
    )
})

test_that("gdb.build_genome accepts min_coverage and validates range", {
    # Signature/argument plumbing only -- no resolution / network call.
    # `path` already exists -> early error before any network is touched, so
    # we exercise the param-validation path cleanly.
    d <- tempfile()
    dir.create(d)
    on.exit(unlink(d, recursive = TRUE))

    # Out-of-range values rejected upfront.
    expect_error(
        gdb.build_genome("hg38", path = d, min_coverage = 0),
        "min_coverage"
    )
    expect_error(
        gdb.build_genome("hg38", path = d, min_coverage = 1.5),
        "min_coverage"
    )
    # Valid in-range value: existing-path error wins (proves min_coverage was
    # accepted by the signature without complaint).
    expect_error(
        gdb.build_genome("hg38", path = d, min_coverage = 0.99),
        "already exists"
    )
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
