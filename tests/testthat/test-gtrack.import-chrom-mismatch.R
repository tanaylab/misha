# Chromosome-name mismatch on WIG / bedGraph / tab import.
#
# gtrack.import used to hand the parser ignore_unknown_chroms = TRUE and never
# look at the outcome, so a file whose chrom names match nothing in the database
# produced a complete, all-NaN track and reported success. The BED branch of the
# same function has always hard-errored on that input.
#
# Policy under test:
#   - zero chromosomes matched  -> error
#   - some matched, some not    -> warning, matched data imported
#   - all of the file's names matched (even if the database has more) -> silent
#   - names the alias chain resolves (chr1 <-> 1) count as matched

ensure_valid_groot()

# Database with UCSC-style names; "scaffold_*" resolves to nothing in it.
setup_mismatch_db <- function() {
    create_test_db(
        "mismatch_db",
        chrom_sizes = data.frame(
            chrom = c("chr1", "chr2"),
            size = c(1000L, 1000L)
        )
    )
    gdb.init("mismatch_db")
}

write_lines_to <- function(name, lines) {
    path <- file.path(getwd(), name)
    writeLines(lines, path)
    path
}

test_that("bedGraph matching zero chromosomes errors instead of importing an empty track", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        bg <- write_lines_to("nomatch.bedgraph", c(
            "scaffold_1\t0\t100\t1.0",
            "scaffold_1\t100\t200\t2.5",
            "scaffold_7\t10\t30\t-1.0"
        ))

        expect_error(
            gtrack.import("nomatch_bg", "no match", bg, binsize = 10),
            "scaffold_1"
        )
        # the message must show both sides of the mismatch and the way out
        err <- tryCatch(gtrack.import("nomatch_bg", "no match", bg, binsize = 10),
            error = function(e) conditionMessage(e)
        )
        expect_match(err, "scaffold_7")
        expect_match(err, "chr1")
        expect_match(err, "chrom_aliases.tsv")

        # and no track is left behind
        expect_false(gtrack.exists("nomatch_bg"))
    })
})

test_that("fixedStep WIG matching zero chromosomes errors", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        wig <- write_lines_to("nomatch.wig", c(
            "fixedStep chrom=scaffold_1 start=1 step=1",
            "1.5", "2.0", "1.8"
        ))

        expect_error(
            gtrack.import("nomatch_wig", "no match", wig, binsize = 1),
            "scaffold_1"
        )
        expect_false(gtrack.exists("nomatch_wig"))
    })
})

test_that("variableStep WIG matching zero chromosomes errors for a sparse track too", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        wig <- write_lines_to("nomatch_var.wig", c(
            "variableStep chrom=ctgA span=10",
            "1\t1.5",
            "21\t2.0"
        ))

        expect_error(
            gtrack.import("nomatch_var", "no match", wig, binsize = 0),
            "ctgA"
        )
        expect_false(gtrack.exists("nomatch_var"))
    })
})

test_that("tab file matching zero chromosomes errors", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        tab <- write_lines_to("nomatch.tab", c(
            "chrom\tstart\tend\tvalue",
            "scaffold_1\t0\t100\t1.0",
            "scaffold_1\t100\t200\t2.5"
        ))

        expect_error(
            gtrack.import("nomatch_tab", "no match", tab, binsize = 10),
            "scaffold_1"
        )
        expect_false(gtrack.exists("nomatch_tab"))
    })
})

test_that("alias-resolvable chromosome names still import", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        # Ensembl-style names against a UCSC-named database: the alias chain
        # maps 1 -> chr1, so this is a legitimate import, not a mismatch.
        bg <- write_lines_to("ensembl.bedgraph", c(
            "1\t0\t100\t1.0",
            "1\t100\t200\t2.5",
            "2\t10\t30\t-1.0"
        ))

        expect_no_warning(gtrack.import("ens_bg", "ensembl names", bg, binsize = 10))
        v <- gextract("ens_bg", gintervals("chr1", 0, 200), colnames = "v")
        expect_equal(unique(v$v[v$start < 100]), 1.0)
        expect_equal(unique(v$v[v$start >= 100]), 2.5)
    })
})

test_that("a file covering only part of the database imports silently", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        # chr1-only file into a chr1 + chr2 database: every name in the file
        # matched, so nothing was dropped and nothing is reported.
        bg <- write_lines_to("chr1only.bedgraph", c(
            "chr1\t0\t100\t1.0",
            "chr1\t100\t200\t2.5"
        ))

        expect_no_warning(gtrack.import("chr1_bg", "chr1 only", bg, binsize = 10))
        v <- gextract("chr1_bg", gintervals("chr1", 0, 200), colnames = "v")
        expect_equal(unique(v$v[v$start < 100]), 1.0)
        expect_true(all(is.na(gextract("chr1_bg", gintervals("chr2"), colnames = "v")$v)))
    })
})

test_that("partially matching file warns, names the skipped chromosomes and imports the rest", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        bg <- write_lines_to("partial.bedgraph", c(
            "chr1\t0\t100\t1.0",
            "chr1\t100\t200\t2.5",
            "scaffold_9\t0\t100\t5.0"
        ))

        expect_warning(
            gtrack.import("partial_bg", "partial", bg, binsize = 10),
            "scaffold_9"
        )
        expect_true(gtrack.exists("partial_bg"))
        v <- gextract("partial_bg", gintervals("chr1", 0, 200), colnames = "v")
        expect_equal(unique(v$v[v$start < 100]), 1.0)
        expect_equal(unique(v$v[v$start >= 100]), 2.5)
    })
})

test_that("partially matching tab file warns", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        tab <- write_lines_to("partial.tab", c(
            "chrom\tstart\tend\tvalue",
            "chr1\t0\t100\t1.0",
            "scaffold_9\t0\t100\t5.0"
        ))

        expect_warning(
            gtrack.import("partial_tab", "partial", tab, binsize = 10),
            "scaffold_9"
        )
        v <- gextract("partial_tab", gintervals("chr1", 0, 100), colnames = "v")
        expect_equal(unique(v$v), 1.0)
    })
})

test_that("gtrack.import_set reports a zero-match file as failed and keeps importing", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        dir.create("set_in")
        writeLines(c("chr1\t0\t100\t1.0"), file.path("set_in", "good.bedgraph"))
        writeLines(c("scaffold_1\t0\t100\t1.0"), file.path("set_in", "bad.bedgraph"))

        res <- suppressMessages(gtrack.import_set(
            "set import", file.path(getwd(), "set_in", "*.bedgraph"),
            binsize = 10
        ))
        expect_true("good.bedgraph" %in% res$files.imported)
        expect_true("bad.bedgraph" %in% res$files.failed)
        expect_true(gtrack.exists("good"))
        expect_false(gtrack.exists("bad"))
    })
})
