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
