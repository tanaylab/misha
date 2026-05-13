# Tests for gdb.install_intervals partial-failure handling.

make_tiny_groot_for_install <- function() {
    groot <- tempfile("misha_test_install_")
    dir.create(groot)
    dir.create(file.path(groot, "tracks"))
    dir.create(file.path(groot, "seq"))
    writeLines(c("chr1\t1000", "chr2\t1000"), file.path(groot, "chrom_sizes.txt"))
    file.create(file.path(groot, "seq", "chr1.seq"))
    file.create(file.path(groot, "seq", "chr2.seq"))
    groot
}

test_that("gdb.install_intervals errors when a requested set has no asset (default)", {
    groot <- make_tiny_groot_for_install()
    on.exit(unlink(groot, recursive = TRUE), add = TRUE)
    gdb.init(groot, rescan = TRUE)

    # `manual` source with no URLs: .manual_fetch_assets returns no entries
    # for any of the requested sets. Currently the function returns silently;
    # we want a hard error.
    expect_error(
        suppressWarnings(suppressMessages(gdb.install_intervals(
            groot = groot,
            source = list(source = "manual", fasta = "n/a"),
            sets = c("genes", "rmsk"),
            verbose = FALSE
        ))),
        "genes.*rmsk|rmsk.*genes"
    )
})

test_that("gdb.install_intervals with force=TRUE warns instead of erroring on missing sets", {
    groot <- make_tiny_groot_for_install()
    on.exit(unlink(groot, recursive = TRUE), add = TRUE)
    gdb.init(groot, rescan = TRUE)

    expect_warning(
        suppressMessages(gdb.install_intervals(
            groot = groot,
            source = list(source = "manual", fasta = "n/a"),
            sets = c("genes", "rmsk"),
            verbose = FALSE,
            force = TRUE
        )),
        "genes.*rmsk|rmsk.*genes"
    )
})

# Hybrid-naming case: when no single alias column covers the groot at
# `min_coverage` but the per-row name-override pass does, the post-rescue
# canonical should reach the threshold. Exercises the gate-after-rescue
# logic without going through a live fetcher.
test_that("post-rescue: name-override rescues hybrid HAL naming (UCSC + bare GenBank)", {
    # Phylo447 shape in miniature: assembled chrom matches `ucsc` column,
    # unplaced bare-accession matches `genbank` column. Neither column
    # alone reaches 100%, but per-row name-override lifts both rows to
    # canonical-in-groot values.
    alias_df <- data.frame(
        ucsc = c("chr1", "chrFOO"),
        genbank = c("CM000663.2", "GL000001.1"),
        stringsAsFactors = FALSE
    )
    groot_chroms <- c("chr1", "GL000001.1")
    groot_lengths <- c(1000, 500)

    # Pre-rescue single-column scores: ucsc=66.7%, genbank=33.3%. The new
    # semantics pick best column unconditionally (min_coverage = 0):
    best <- .detect_alias_column(alias_df, groot_chroms,
        min_coverage = 0, chrom_lengths = groot_lengths
    )
    expect_equal(as.character(best), "ucsc")

    # After name-match-override across the row's other columns, every
    # alias row gets canonicalized to a groot chrom.
    canonical <- .name_match_override(
        alias_df$ucsc, alias_df, "ucsc", groot_chroms
    )
    expect_setequal(canonical, c("chr1", "GL000001.1"))
    expect_equal(
        .canonical_coverage(canonical, groot_chroms, groot_lengths),
        1.0
    )
})

test_that("post-rescue: gate fires when canonical still doesn't cover groot", {
    # Single alias row can't rescue a missing groot contig. Post-rescue
    # canonical covers 1000 / 1500 = 66.7% bp; min_coverage = 1.0 must
    # reject this.
    alias_df <- data.frame(
        ucsc = "chr1",
        genbank = "CM000663.2",
        stringsAsFactors = FALSE
    )
    groot_chroms <- c("chr1", "GL000001.1")
    groot_lengths <- c(1000, 500)
    canonical <- .name_match_override(
        alias_df$ucsc, alias_df, "ucsc", groot_chroms
    )
    expect_lt(
        .canonical_coverage(canonical, groot_chroms, groot_lengths),
        1.0
    )
})
