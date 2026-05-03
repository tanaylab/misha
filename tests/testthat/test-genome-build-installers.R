# Helper: build a tiny in-memory groot with two contigs.
make_test_groot <- function() {
    groot <- tempfile()
    dir.create(groot)
    dir.create(file.path(groot, "tracks"))
    dir.create(file.path(groot, "seq"))
    # Minimal chrom_sizes
    writeLines(c("chr1\t1000", "chr2\t1000"), file.path(groot, "chrom_sizes.txt"))
    # Empty .seq files (tests that don't touch sequence won't notice)
    file.create(file.path(groot, "seq", "chr1.seq"))
    file.create(file.path(groot, "seq", "chr2.seq"))
    groot
}

test_that(".save_intervals filters to ALLGENOME and writes the set", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    gdb.init(groot, rescan = TRUE)
    df <- data.frame(
        chrom = c("chr1", "chr1", "chrX"),
        start = c(0L, 100L, 0L),
        end = c(50L, 150L, 50L),
        stringsAsFactors = FALSE
    )
    .save_intervals("foo", df, overwrite = FALSE, verbose = FALSE)
    expect_true(gintervals.exists("foo"))
    expect_equal(nrow(gintervals.load("foo")), 2L)
})

test_that(".save_intervals errors on existing set when overwrite=FALSE", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    gdb.init(groot, rescan = TRUE)
    df <- data.frame(chrom = "chr1", start = 0L, end = 50L, stringsAsFactors = FALSE)
    .save_intervals("foo", df, overwrite = FALSE, verbose = FALSE)
    expect_error(
        .save_intervals("foo", df, overwrite = FALSE, verbose = FALSE),
        "already exists"
    )
})

test_that(".save_intervals overwrites when overwrite=TRUE", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    gdb.init(groot, rescan = TRUE)
    df1 <- data.frame(chrom = "chr1", start = 0L, end = 50L, stringsAsFactors = FALSE)
    df2 <- data.frame(
        chrom = "chr1",
        start = c(0L, 200L),
        end = c(50L, 300L),
        stringsAsFactors = FALSE
    )
    .save_intervals("foo", df1, overwrite = FALSE, verbose = FALSE)
    .save_intervals("foo", df2, overwrite = TRUE, verbose = FALSE)
    expect_equal(nrow(gintervals.load("foo")), 2L)
})

test_that(".save_intervals skips empty frame with message", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    gdb.init(groot, rescan = TRUE)
    df <- data.frame(
        chrom = character(0),
        start = integer(0),
        end = integer(0),
        stringsAsFactors = FALSE
    )
    expect_message(
        .save_intervals("empty", df, overwrite = FALSE, verbose = TRUE),
        "0 rows"
    )
    expect_false(gintervals.exists("empty"))
})

test_that(".install_rmsk_set produces combined + per-class subsets", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    gdb.init(groot, rescan = TRUE)
    df <- data.frame(
        chrom = c("chr1", "chr1", "chr1", "chr1", "chr1"),
        start = c(0L, 100L, 200L, 300L, 400L),
        end = c(50L, 150L, 250L, 350L, 450L),
        strand = c(1L, -1L, 1L, 1L, -1L),
        name = c("MIR", "L1", "Alu", "Charlie", "Repeat?"),
        class = c("SINE", "LINE", "SINE", "DNA", "DNA?"),
        family = c("MIR", "L1", "Alu", "hAT-Charlie", NA),
        stringsAsFactors = FALSE
    )
    .install_rmsk_set(df, prefix = "", overwrite = FALSE, verbose = FALSE)
    expect_true(gintervals.exists("rmsk"))
    expect_true(gintervals.exists("rmsk_sine"))
    expect_true(gintervals.exists("rmsk_line"))
    expect_true(gintervals.exists("rmsk_dna"))
    expect_true(gintervals.exists("rmsk_dna_qmark"))
    expect_equal(nrow(gintervals.load("rmsk")), 5L)
    expect_equal(nrow(gintervals.load("rmsk_sine")), 2L)
    expect_equal(nrow(gintervals.load("rmsk_dna_qmark")), 1L)
})

test_that(".install_rmsk_set respects prefix", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    # Dotted prefix = misha namespace hierarchy; pre-create the subdirectory.
    dir.create(file.path(groot, "tracks", "intervs"), showWarnings = FALSE)
    dir.create(file.path(groot, "tracks", "intervs", "global"), showWarnings = FALSE)
    gdb.init(groot, rescan = TRUE)
    df <- data.frame(
        chrom = "chr1", start = 0L, end = 50L, strand = 1L,
        name = "x", class = "SINE", family = "MIR",
        stringsAsFactors = FALSE
    )
    .install_rmsk_set(df, prefix = "intervs.global.", overwrite = FALSE, verbose = FALSE)
    expect_true(gintervals.exists("intervs.global.rmsk"))
    expect_true(gintervals.exists("intervs.global.rmsk_sine"))
})

test_that(".install_cgi_set saves cgi (not cpgIsland)", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    gdb.init(groot, rescan = TRUE)
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L,
        name = "cgi-1", length = 100L, cpgNum = 10L,
        perCpg = 0.5, perGc = 0.6, obsExp = 0.7,
        stringsAsFactors = FALSE
    )
    .install_cgi_set(df, prefix = "", overwrite = FALSE, verbose = FALSE)
    expect_true(gintervals.exists("cgi"))
    expect_false(gintervals.exists("cpgIsland"))
})

test_that(".install_cytoband_set saves cytoband", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    gdb.init(groot, rescan = TRUE)
    df <- data.frame(
        chrom = "chr1", start = 0L, end = 1000L,
        name = "p1", stain = "gneg",
        stringsAsFactors = FALSE
    )
    .install_cytoband_set(df, prefix = "", overwrite = FALSE, verbose = FALSE)
    expect_true(gintervals.exists("cytoband"))
})

test_that(".install_genes_set with default gene_sets produces tss/exons/utr3/utr5", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    # gintervals.import_genes needs <groot>/downloads.
    dir.create(file.path(groot, "downloads"), showWarnings = FALSE)
    gdb.init(groot, rescan = TRUE)

    # Minimal 12-col genePred: name chrom strand txStart txEnd cdsStart cdsEnd
    #                          exonCount exonStarts exonEnds score name2
    gp <- tempfile(fileext = ".genePred")
    writeLines(c(
        "tx1\tchr1\t+\t100\t500\t150\t450\t1\t100,\t500,\t0\tgene1"
    ), gp)

    .install_genes_set(
        list(file = gp, format = "genepred"),
        prefix = "",
        gene_sets = c(tss = "tss", exons = "exons", utr3 = "utr3", utr5 = "utr5"),
        overwrite = FALSE, verbose = FALSE
    )
    expect_true(gintervals.exists("tss"))
    expect_true(gintervals.exists("exons"))
    expect_true(gintervals.exists("utr3"))
    expect_true(gintervals.exists("utr5"))
})

test_that(".install_genes_set with NA in gene_sets skips that role", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    dir.create(file.path(groot, "downloads"), showWarnings = FALSE)
    gdb.init(groot, rescan = TRUE)

    gp <- tempfile(fileext = ".genePred")
    writeLines(c(
        "tx1\tchr1\t+\t100\t500\t150\t450\t1\t100,\t500,\t0\tgene1"
    ), gp)

    .install_genes_set(
        list(file = gp, format = "genepred"),
        prefix = "",
        gene_sets = c(tss = "tss", exons = "exon", utr3 = NA, utr5 = NA),
        overwrite = FALSE, verbose = FALSE
    )
    expect_true(gintervals.exists("tss"))
    expect_true(gintervals.exists("exon")) # renamed
    expect_false(gintervals.exists("exons")) # not present
    expect_false(gintervals.exists("utr3"))
    expect_false(gintervals.exists("utr5"))
})

test_that(".install_genes_set creates intermediate dirs for dotted prefix", {
    groot <- make_test_groot()
    on.exit(unlink(groot, recursive = TRUE))
    dir.create(file.path(groot, "downloads"), showWarnings = FALSE)
    gdb.init(groot, rescan = TRUE)

    gp <- tempfile(fileext = ".genePred")
    writeLines(c(
        "tx1\tchr1\t+\t100\t500\t150\t450\t1\t100,\t500,\t0\tgene1"
    ), gp)

    .install_genes_set(
        list(file = gp, format = "genepred"),
        prefix = "intervs.global.",
        gene_sets = c(tss = "tss", exons = "exon", utr3 = NA, utr5 = NA),
        overwrite = FALSE, verbose = FALSE
    )
    expect_true(file.exists(file.path(.misha$GROOT, "tracks/intervs/global/tss.interv")))
    expect_true(file.exists(file.path(.misha$GROOT, "tracks/intervs/global/exon.interv")))
})
