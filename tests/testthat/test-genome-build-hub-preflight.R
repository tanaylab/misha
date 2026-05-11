# Unit tests for .hub_preflight_coverage. All network calls are mocked.

# Synthetic 3-row chromAlias: ucsc has full coverage, genbank misses chrM.
.preflight_alias <- function() {
    data.frame(
        ucsc = c("chr1", "chr2", "chrM"),
        genbank = c("CM00001.1", "CM00002.1", ""),
        refseq = c("NC_00001.1", "NC_00002.1", "NC_M.1"),
        stringsAsFactors = FALSE
    )
}

# Synthetic chrom.sizes (FASTA's source col is refseq for this fixture).
.preflight_sizes <- "NC_00001.1\t100000000\nNC_00002.1\t100000000\nNC_M.1\t16000\n"

with_preflight_mocks <- function(downloads_recorded, code) {
    testthat::local_mocked_bindings(
        .hub_url_for = function(accession) "https://hub.example/test/",
        .hub_list_dir = function(url, verbose = TRUE) {
            c(
                "GCF_TEST.chromAlias.txt",
                "GCF_TEST.chrom.sizes.txt",
                "GCF_TEST.fa.gz"
            )
        },
        .download_to = function(url, dest, verbose = TRUE) {
            downloads_recorded$urls <- c(downloads_recorded$urls, url)
            if (grepl("chromAlias", url)) {
                writeLines("# ucsc\tgenbank\trefseq", dest)
            } else if (grepl("chrom\\.sizes", url)) {
                cat(.preflight_sizes, file = dest)
            }
            invisible(dest)
        },
        .parse_ucsc_chromalias = function(path) .preflight_alias(),
        .package = "misha"
    )
    force(code)
}

test_that(".hub_preflight_coverage fails without fetching FASTA when min_coverage unmet", {
    rec <- new.env()
    rec$urls <- character()
    with_preflight_mocks(rec, {
        wd <- tempfile("preflight_")
        dir.create(wd)
        on.exit(unlink(wd, recursive = TRUE), add = TRUE)
        expect_error(
            .hub_preflight_coverage(
                accession = "GCF_TEST.1",
                target_chroms = NULL,
                chrom_naming = "genbank",
                min_coverage = 1.0,
                workdir = wd,
                verbose = FALSE
            ),
            "no column with 100"
        )
    })
    expect_false(any(grepl("\\.fa\\.gz$", rec$urls)),
        info = paste("URLs touched:", paste(rec$urls, collapse = ", "))
    )
    expect_true(any(grepl("chromAlias", rec$urls)))
    expect_true(any(grepl("chrom\\.sizes", rec$urls)))
})

test_that(".hub_preflight_coverage returns prefetched alias on success", {
    rec <- new.env()
    rec$urls <- character()
    with_preflight_mocks(rec, {
        wd <- tempfile("preflight_")
        dir.create(wd)
        on.exit(unlink(wd, recursive = TRUE), add = TRUE)
        out <- .hub_preflight_coverage(
            accession = "GCF_TEST.1",
            target_chroms = NULL,
            chrom_naming = "ucsc",
            min_coverage = 1.0,
            workdir = wd,
            verbose = FALSE
        )
        expect_named(out, c(
            "df", "row_lengths", "alias_file",
            "sizes_file", "report_file", "canonical_col"
        ),
        ignore.order = TRUE
        )
        expect_equal(out$canonical_col, "ucsc")
        expect_equal(nrow(out$df), 3L)
        expect_equal(length(out$row_lengths), 3L)
    })
    expect_false(any(grepl("\\.fa\\.gz$", rec$urls)))
})

test_that(".hub_preflight_coverage rescues a canonical with empty cells via target_lengths", {
    # target_chroms picks the genbank column (2/3 of names match) but that
    # column has an empty cell for chrM. Without rescue the strict gate
    # would fail (16 kb missing). With target_lengths supplied and
    # match_by_length=TRUE, the empty cell is matched to "chrM" by unique
    # length pairing -> preflight passes.
    rec <- new.env()
    rec$urls <- character()
    with_preflight_mocks(rec, {
        wd <- tempfile("preflight_")
        dir.create(wd)
        on.exit(unlink(wd, recursive = TRUE), add = TRUE)
        out <- .hub_preflight_coverage(
            accession       = "GCF_TEST.1",
            target_chroms   = c("CM00001.1", "CM00002.1", "chrM"),
            target_lengths  = c(1e8, 1e8, 16000),
            chrom_naming    = "genbank",
            min_coverage    = 1.0,
            match_by_length = TRUE,
            workdir         = wd,
            verbose         = FALSE
        )
        expect_equal(out$canonical_col, "genbank")
    })
    expect_false(any(grepl("\\.fa\\.gz$", rec$urls)))
})

test_that(".hub_preflight_coverage merges assembly_report.txt and rescues coverage", {
    # chromAlias alone covers 2 of 3 groot contigs by name (66.7%); the
    # mirrored assembly_report.txt fills in the missing 3rd row (chrUn_X)
    # via its UCSC-style-name column, so the strict gate passes at 1.0.
    rec <- new.env()
    rec$urls <- character()
    testthat::local_mocked_bindings(
        .hub_url_for = function(accession) "https://hub.example/test/",
        .hub_list_dir = function(url, verbose = TRUE) {
            c(
                "GCF_TEST.chromAlias.txt",
                "GCF_TEST.chrom.sizes.txt",
                "GCF_TEST_assembly_report.txt",
                "GCF_TEST.fa.gz"
            )
        },
        .download_to = function(url, dest, verbose = TRUE) {
            rec$urls <- c(rec$urls, url)
            if (grepl("chromAlias", url)) {
                # chromAlias has only 2 rows -- chrUn_X is missing.
                writeLines(c(
                    "# ucsc\tgenbank\trefseq",
                    "chr1\tCM00001.1\tNC_00001.1",
                    "chr2\tCM00002.1\tNC_00002.1"
                ), dest)
            } else if (grepl("assembly_report", url)) {
                # Report has all 3 rows including chrUn_X.
                writeLines(c(
                    "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name",
                    "1\tassembled-molecule\t1\tChromosome\tCM00001.1\t=\tNC_00001.1\tPrimary Assembly\t100000000\tchr1",
                    "2\tassembled-molecule\t2\tChromosome\tCM00002.1\t=\tNC_00002.1\tPrimary Assembly\t100000000\tchr2",
                    "unp1\tunplaced-scaffold\tna\tna\tJH99.1\t=\tNW_99.1\tPrimary Assembly\t500\tchrUn_X"
                ), dest)
            } else if (grepl("chrom\\.sizes", url)) {
                # Sizes keyed on refseq, including the unplaced row.
                cat("NC_00001.1\t100000000\nNC_00002.1\t100000000\nNW_99.1\t500\n",
                    file = dest
                )
            }
            invisible(dest)
        },
        .parse_ucsc_chromalias = function(path) {
            data.frame(
                ucsc = c("chr1", "chr2"),
                genbank = c("CM00001.1", "CM00002.1"),
                refseq = c("NC_00001.1", "NC_00002.1"),
                stringsAsFactors = FALSE
            )
        },
        .package = "misha"
    )
    wd <- tempfile("preflight_")
    dir.create(wd)
    on.exit(unlink(wd, recursive = TRUE), add = TRUE)
    # Without the assembly_report merge this would fail (chromAlias only
    # covers 2 of 3 contigs); with it, alias_df grows to 3 rows and the
    # ucsc column covers all of them.
    out <- .hub_preflight_coverage(
        accession = "GCF_TEST.1",
        target_chroms = NULL,
        chrom_naming = "ucsc",
        min_coverage = 1.0,
        workdir = wd,
        verbose = FALSE
    )
    expect_equal(nrow(out$df), 3L)
    expect_true("chrUn_X" %in% out$df$ucsc)
    expect_false(any(grepl("\\.fa\\.gz$", rec$urls)))
    expect_true(any(grepl("assembly_report", rec$urls)))
})

test_that(".hub_preflight_coverage errors before FASTA when a target_chrom is unplaceable", {
    # 'HALxyz' isn't in any chromAlias column and its length (42) doesn't
    # appear on the alias side -> can't be placed -> preflight errors before
    # the FASTA download.
    rec <- new.env()
    rec$urls <- character()
    with_preflight_mocks(rec, {
        wd <- tempfile("preflight_")
        dir.create(wd)
        on.exit(unlink(wd, recursive = TRUE), add = TRUE)
        expect_error(
            .hub_preflight_coverage(
                accession       = "GCF_TEST.1",
                target_chroms   = c("chr1", "chr2", "HALxyz"),
                target_lengths  = c(1e8, 1e8, 42),
                chrom_naming    = "ucsc",
                min_coverage    = 1.0,
                match_by_length = TRUE,
                workdir         = wd,
                verbose         = FALSE
            ),
            "Cannot align"
        )
    })
    expect_false(any(grepl("\\.fa\\.gz$", rec$urls)))
})

test_that(".hub_preflight_coverage still fails strict when match_by_length=FALSE", {
    # Same fixture as the rescue test but with match_by_length disabled - the
    # strict gate must fire (no rescue eligibility).
    rec <- new.env()
    rec$urls <- character()
    with_preflight_mocks(rec, {
        wd <- tempfile("preflight_")
        dir.create(wd)
        on.exit(unlink(wd, recursive = TRUE), add = TRUE)
        expect_error(
            .hub_preflight_coverage(
                accession       = "GCF_TEST.1",
                target_chroms   = c("CM00001.1", "CM00002.1", "chrM"),
                target_lengths  = c(1e8, 1e8, 16000),
                chrom_naming    = "genbank",
                min_coverage    = 1.0,
                match_by_length = FALSE,
                workdir         = wd,
                verbose         = FALSE
            ),
            "no column with 100"
        )
    })
})

test_that("gdb.build_genome ucsc-hub fails before fetching FASTA on coverage failure", {
    rec <- new.env()
    rec$urls <- character()
    fake_recipe <- list(
        source = "ucsc-hub", accession = "GCF_TEST.1",
        chrom_naming = "ucsc"
    )
    testthat::local_mocked_bindings(
        .resolve_genome = function(name, registry = NULL) {
            list(recipe = fake_recipe, resolved_from = "test-fixture")
        },
        .hub_url_for = function(accession) "https://hub.example/test/",
        .hub_list_dir = function(url, verbose = TRUE) {
            c(
                "GCF_TEST.chromAlias.txt",
                "GCF_TEST.chrom.sizes.txt",
                "GCF_TEST.fa.gz"
            )
        },
        .download_to = function(url, dest, verbose = TRUE) {
            rec$urls <- c(rec$urls, url)
            if (grepl("chromAlias", url)) {
                writeLines("# ucsc\tgenbank\trefseq", dest)
            } else if (grepl("chrom\\.sizes", url)) {
                cat(.preflight_sizes, file = dest)
            }
            invisible(dest)
        },
        .parse_ucsc_chromalias = function(path) {
            data.frame(
                ucsc = c("chr1", "chr2", ""),
                genbank = c("CM00001.1", "", ""),
                refseq = c("NC_00001.1", "NC_00002.1", "NC_M.1"),
                stringsAsFactors = FALSE
            )
        },
        .package = "misha"
    )
    out <- tempfile("misha_test_groot_")
    on.exit(unlink(out, recursive = TRUE), add = TRUE)
    expect_error(
        gdb.build_genome(
            name         = "ucsc-hub-test-fixture",
            path         = out,
            sets         = "rmsk",
            min_coverage = 1.0,
            verbose      = FALSE
        ),
        "no column with 100"
    )
    expect_false(any(grepl("\\.fa\\.gz$", rec$urls)),
        info = paste("Touched URLs:", paste(rec$urls, collapse = ", "))
    )
})
