# Unit tests for gdb.build_genome() and friends. All offline — no network,
# no UCSC binary required. End-to-end build tests are deferred (would need
# skip_on_cran + skip_if_offline + a real upstream).

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

write_temp_yaml <- function(text) {
    f <- tempfile(fileext = ".yaml")
    writeLines(text, f)
    f
}

# ---------------------------------------------------------------------------
# Registry parsing / normalization / validation
# ---------------------------------------------------------------------------

test_that(".normalize_recipe accepts string (legacy local form)", {
    r <- misha:::.normalize_recipe("/tmp/groot", "hg19")
    expect_equal(r$source, "local")
    expect_equal(r$path, "/tmp/groot")
})

test_that(".normalize_recipe preserves a dict and requires 'source'", {
    r <- misha:::.normalize_recipe(list(source = "ucsc", assembly = "hg38"), "hg38")
    expect_equal(r$source, "ucsc")
    expect_equal(r$assembly, "hg38")

    expect_error(
        misha:::.normalize_recipe(list(assembly = "hg38"), "hg38"),
        "missing required 'source:'"
    )
    expect_error(
        misha:::.normalize_recipe(list(source = "made_up"), "x"),
        "unknown source"
    )
    expect_error(
        misha:::.normalize_recipe(42, "x"),
        "invalid registry entry"
    )
})

test_that(".validate_recipe enforces per-source required fields", {
    expect_error(misha:::.validate_recipe(list(source = "ucsc"), "x"), "assembly")
    expect_error(misha:::.validate_recipe(list(source = "ncbi"), "x"), "accession")
    expect_error(
        misha:::.validate_recipe(list(source = "ncbi", accession = "garbage"), "x"),
        "does not match"
    )
    expect_error(misha:::.validate_recipe(list(source = "local"), "x"), "path")
    expect_error(misha:::.validate_recipe(list(source = "manual"), "x"), "fasta")

    # OK paths
    expect_silent(misha:::.validate_recipe(list(source = "ucsc", assembly = "hg38"), "x"))
    expect_silent(misha:::.validate_recipe(list(source = "ncbi", accession = "GCF_000001405.40"), "x"))
    expect_silent(misha:::.validate_recipe(list(source = "manual", fasta = "url"), "x"))
})

test_that(".parse_genome_registry handles valid + invalid YAML", {
    f <- write_temp_yaml(c(
        "genome:",
        "  hg38:",
        "    source: ucsc",
        "    assembly: hg38",
        "  custom: /tmp/local"
    ))
    on.exit(unlink(f))
    entries <- misha:::.parse_genome_registry(f)
    expect_named(entries, c("hg38", "custom"), ignore.order = TRUE)
    expect_equal(entries$custom, "/tmp/local")
    expect_equal(entries$hg38$source, "ucsc")

    # No genome: key -> empty list, not error
    f2 <- write_temp_yaml(c("other: stuff"))
    on.exit(unlink(f2), add = TRUE)
    expect_equal(length(misha:::.parse_genome_registry(f2)), 0)

    # Missing file
    expect_error(
        misha:::.parse_genome_registry("/nonexistent/registry.yaml"),
        "does not exist"
    )
})

# ---------------------------------------------------------------------------
# Resolution chain
# ---------------------------------------------------------------------------

test_that(".resolve_genome explicit registry overrides built-in", {
    f <- write_temp_yaml(c(
        "genome:",
        "  hg38:",
        "    source: manual",
        "    fasta: file:///tmp/custom.fa"
    ))
    on.exit(unlink(f))
    res <- misha:::.resolve_genome("hg38", registry = f)
    expect_equal(res$recipe$source, "manual")
    expect_match(res$resolved_from, "registry arg")
})

test_that(".resolve_genome falls back to built-in registry", {
    res <- misha:::.resolve_genome("mm10")
    expect_equal(res$recipe$source, "ucsc")
    expect_equal(res$recipe$assembly, "mm10")
    expect_match(res$resolved_from, "built-in")
})

test_that(".resolve_genome pattern fallback for NCBI accessions", {
    res <- misha:::.resolve_genome("GCF_009806435.1")
    expect_equal(res$recipe$source, "ncbi")
    expect_equal(res$recipe$accession, "GCF_009806435.1")
    expect_match(res$resolved_from, "pattern fallback")

    res2 <- misha:::.resolve_genome("GCA_000001405.29")
    expect_equal(res2$recipe$source, "ncbi")
})

test_that(".resolve_genome errors clearly for unknown name", {
    expect_error(misha:::.resolve_genome("totally_made_up_genome"), "not found")
})

# ---------------------------------------------------------------------------
# UCSC URL builder
# ---------------------------------------------------------------------------

test_that(".ucsc_urls constructs the right URLs for hg38", {
    urls <- misha:::.ucsc_urls("hg38",
        annotations = c("genes", "rmsk", "cpgIsland", "cytoband")
    )
    expect_match(urls$fasta, "/goldenPath/hg38/bigZips/hg38\\.fa\\.gz$")
    expect_match(urls$genes, "/database/ncbiRefSeq\\.txt\\.gz$")
    expect_match(urls$annots, "/database/ncbiRefSeqLink\\.txt\\.gz$")
    expect_match(urls$rmsk, "/database/rmsk\\.txt\\.gz$")
    expect_match(urls$cpgIsland, "/database/cpgIslandExt\\.txt\\.gz$")
    expect_match(urls$cytoband, "/database/cytoBandIdeo\\.txt\\.gz$")
})

test_that(".ucsc_urls omits annotations not requested", {
    urls <- misha:::.ucsc_urls("mm10", annotations = c("genes"))
    expect_true(!is.null(urls$genes))
    expect_null(urls$rmsk)
    expect_null(urls$cpgIsland)
    expect_null(urls$cytoband)
})

# ---------------------------------------------------------------------------
# Architecture detection
# ---------------------------------------------------------------------------

test_that(".detect_arch returns one of the known platforms on this host", {
    arch <- misha:::.detect_arch()
    expect_true(arch %in% c("linux_x86_64", "macos_x86_64", "macos_arm64"))
})

# ---------------------------------------------------------------------------
# Binary path resolver
# ---------------------------------------------------------------------------

test_that(".gff3_to_genepred_path honors MISHA_GFF3_TO_GENEPRED", {
    skip_if_not_installed("withr")
    stub <- tempfile(fileext = ".sh")
    writeLines("#!/bin/sh\necho 'stub'", stub)
    Sys.chmod(stub, "0755")
    on.exit(unlink(stub))

    withr::with_envvar(
        c(MISHA_GFF3_TO_GENEPRED = stub),
        {
            p <- misha:::.gff3_to_genepred_path()
            expect_equal(normalizePath(p), normalizePath(stub))
        }
    )
})

test_that(".gff3_to_genepred_path errors if env var points at nonexistent file", {
    skip_if_not_installed("withr")
    withr::with_envvar(
        c(MISHA_GFF3_TO_GENEPRED = "/nonexistent/binary"),
        {
            expect_error(misha:::.gff3_to_genepred_path(), "non-existent")
        }
    )
})

test_that(".gff3_to_genepred_path returns NULL when nothing is cached and no env var", {
    skip_if_not_installed("withr")
    # Guard: make sure no env var is set.
    withr::with_envvar(
        c(MISHA_GFF3_TO_GENEPRED = ""),
        {
            # Move any existing cache out of the way for the duration.
            cached <- misha:::.gff3_to_genepred_cache_path()
            if (file.exists(cached)) {
                bak <- paste0(cached, ".bak.", Sys.getpid())
                file.rename(cached, bak)
                on.exit(file.rename(bak, cached), add = TRUE)
            }
            expect_null(misha:::.gff3_to_genepred_path())
        }
    )
})

# ---------------------------------------------------------------------------
# UCSC table parsers
# ---------------------------------------------------------------------------

test_that(".parse_ucsc_rmsk handles a minimal record", {
    f <- tempfile(fileext = ".txt")
    on.exit(unlink(f))
    # 17 cols, two rows
    writeLines(c(
        paste0(
            "0\t1000\t100\t10\t10\t", # bin/swScore/milliDiv/milliDel/milliIns
            "chr1\t100\t200\t-1000\t+\t", # genoName/genoStart/genoEnd/genoLeft/strand
            "L1HS\tLINE\tL1\t", # repName/repClass/repFamily
            "1\t100\t-50\t1"
        ), # repStart/repEnd/repLeft/id
        paste0(
            "0\t800\t100\t10\t10\t",
            "chr2\t500\t800\t-500\t-\t",
            "AluSx\tSINE\tAlu\t",
            "1\t300\t-100\t2"
        )
    ), f)
    df <- misha:::.parse_ucsc_rmsk(f)
    expect_equal(nrow(df), 2)
    expect_equal(df$chrom, c("chr1", "chr2"))
    expect_equal(df$strand, c(1L, -1L))
    expect_equal(df$class, c("LINE", "SINE"))
    expect_equal(df$family, c("L1", "Alu"))
})

test_that(".parse_ucsc_cpg_island handles a minimal record", {
    f <- tempfile(fileext = ".txt")
    on.exit(unlink(f))
    # 11 cols
    writeLines("0\tchr1\t100\t300\tCpG: 25\t200\t25\t150\t25.0\t75.0\t0.85", f)
    df <- misha:::.parse_ucsc_cpg_island(f)
    expect_equal(nrow(df), 1)
    expect_equal(df$chrom, "chr1")
    expect_equal(df$start, 100)
    expect_equal(df$end, 300)
})

test_that(".parse_ucsc_cytoband handles a minimal record", {
    f <- tempfile(fileext = ".txt")
    on.exit(unlink(f))
    # 5 cols
    writeLines(c(
        "chr1\t0\t2300000\tp36.33\tgneg",
        "chr1\t2300000\t5400000\tp36.32\tgpos25"
    ), f)
    df <- misha:::.parse_ucsc_cytoband(f)
    expect_equal(nrow(df), 2)
    expect_equal(df$stain, c("gneg", "gpos25"))
})

# ---------------------------------------------------------------------------
# gdb.list_genomes / gdb.genome_info
# ---------------------------------------------------------------------------

test_that("gdb.list_genomes includes built-in entries", {
    df <- gdb.list_genomes()
    expect_true("hg38" %in% df$name)
    expect_true("mm10" %in% df$name)
    expect_equal(unique(df$source[df$name %in% c("hg38", "mm10")]), "ucsc")
})

test_that("gdb.list_genomes accepts an explicit registry", {
    f <- write_temp_yaml(c(
        "genome:",
        "  fake: {source: manual, fasta: file:///tmp/x.fa}"
    ))
    on.exit(unlink(f))
    df <- gdb.list_genomes(registry = f)
    expect_true("fake" %in% df$name)
    expect_equal(df$source[df$name == "fake"], "manual")
})

test_that("gdb.genome_info returns a recipe", {
    info <- gdb.genome_info("hg38")
    expect_equal(info$recipe$source, "ucsc")
    expect_equal(info$recipe$assembly, "hg38")
})

test_that("gdb.genome_info on NCBI accession uses pattern fallback", {
    info <- gdb.genome_info("GCF_009806435.1")
    expect_equal(info$recipe$source, "ncbi")
    expect_equal(info$recipe$accession, "GCF_009806435.1")
})

# ---------------------------------------------------------------------------
# Stub-based dispatch sanity check
# ---------------------------------------------------------------------------

test_that("gdb.build_genome refuses to overwrite existing path", {
    d <- tempfile()
    dir.create(d)
    on.exit(unlink(d, recursive = TRUE))
    expect_error(gdb.build_genome("hg38", path = d), "already exists")
})
