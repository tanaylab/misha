# Chromosome-name mismatch on WIG / bedGraph / tab import.
#
# gtrack.import used to hand the parser ignore_unknown_chroms = TRUE and never
# look at the outcome, so a file whose chrom names match nothing in the database
# produced a complete, all-NaN track and reported success. The BED branch of the
# same function has always hard-errored on that input.
#
# Policy under test:
#   - zero chromosomes matched  -> error
#   - some matched, some not    -> matched data imported, skipped names reported:
#     a message for contigs the database does not have, a warning when a skipped
#     name looks like a primary chromosome (chr7, X, MT)
#   - all of the file's names matched (even if the database has more) -> silent
#   - no chromosome records in the file at all -> silent (deliberate carve-out)
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

test_that("skipped scaffolds are reported as a message, not a warning, and the rest is imported", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        bg <- write_lines_to("partial.bedgraph", c(
            "chr1\t0\t100\t1.0",
            "chr1\t100\t200\t2.5",
            "scaffold_9\t0\t100\t5.0"
        ))

        # A contig the database does not have is the everyday case (whole-genome
        # bigWigs are full of them), so it must not consume the warning channel.
        expect_no_warning(
            expect_message(
                gtrack.import("partial_bg", "partial", bg, binsize = 10),
                "scaffold_9"
            )
        )
        expect_true(gtrack.exists("partial_bg"))
        v <- gextract("partial_bg", gintervals("chr1", 0, 200), colnames = "v")
        expect_equal(unique(v$v[v$start < 100]), 1.0)
        expect_equal(unique(v$v[v$start >= 100]), 2.5)
    })
})

test_that("skipped scaffolds in a tab file are reported as a message", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        tab <- write_lines_to("partial.tab", c(
            "chrom\tstart\tend\tvalue",
            "chr1\t0\t100\t1.0",
            "scaffold_9\t0\t100\t5.0"
        ))

        expect_no_warning(
            expect_message(
                gtrack.import("partial_tab", "partial", tab, binsize = 10),
                "scaffold_9"
            )
        )
        v <- gextract("partial_tab", gintervals("chr1", 0, 100), colnames = "v")
        expect_equal(unique(v$v), 1.0)
    })
})

test_that("a skipped primary chromosome warns and points at the alias mechanism", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        # chr3 is not a scaffold the database chose to omit - it looks like a
        # chromosome that should have been there, so this one is a warning.
        bg <- write_lines_to("primary.bedgraph", c(
            "chr1\t0\t100\t1.0",
            "chr3\t0\t100\t5.0"
        ))

        w <- capture_warnings(gtrack.import("primary_bg", "primary", bg, binsize = 10))
        expect_length(w, 1)
        expect_match(w, "chr3")
        expect_match(w, "chrom_aliases.tsv")
        v <- gextract("primary_bg", gintervals("chr1", 0, 100), colnames = "v")
        expect_equal(unique(v$v), 1.0)
    })
})

test_that("mixed scaffolds and a primary chromosome warn about the primary one", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        bg <- write_lines_to("mixed.bedgraph", c(
            "chr1\t0\t100\t1.0",
            "chr1_gl000191_random\t0\t100\t5.0",
            "chrUn_KI270302v1\t0\t100\t5.0",
            "chrX\t0\t100\t5.0"
        ))

        w <- capture_warnings(gtrack.import("mixed_bg", "mixed", bg, binsize = 10))
        expect_length(w, 1)
        expect_match(w, "chrX")
        # the scaffolds are counted, but they are not what the warning is about
        expect_match(w, "^3 chromosome name")
        expect_false(grepl("gl000191", w))
    })
})

test_that("a skipped mitochondrial name is a message, not a warning", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        # The alias chain generates every mito spelling for whatever mito contig a
        # database has, so an unmatched MT proves the database has no mito contig at
        # all - unfixable by an alias, and the everyday case for hg19 / mm10 as built
        # here. It belongs in the message channel with the other absent contigs.
        bg <- write_lines_to("mito.bedgraph", c(
            "chr1\t0\t100\t1.0",
            "MT\t0\t100\t5.0"
        ))

        expect_no_warning(
            expect_message(
                gtrack.import("mito_bg", "mito", bg, binsize = 10),
                "MT"
            )
        )
        v <- gextract("mito_bg", gintervals("chr1", 0, 100), colnames = "v")
        expect_equal(unique(v$v), 1.0)
    })
})

test_that("chrM in a whole-genome file does not push the report into the warning channel", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        bg <- write_lines_to("wg.bedgraph", c(
            "chr1\t0\t100\t1.0",
            "chr1_gl000191_random\t0\t100\t5.0",
            "chrM\t0\t100\t5.0"
        ))

        expect_no_warning(expect_message(
            gtrack.import("wg_bg", "whole genome", bg, binsize = 10),
            "chrM"
        ))
    })
})

test_that("a primary chromosome after the tracking cap is still noticed", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        # 120 junk names exhaust the distinct-name cap; chr7 comes after all of them
        # and must still route the report to the warning channel.
        junk <- sprintf("READ%05d\t0\t100\t1.0", 1:120)
        bg <- write_lines_to("late_primary.bedgraph", c(
            "chr1\t0\t100\t1.0", junk, "chr7\t0\t100\t5.0"
        ))

        w <- capture_warnings(gtrack.import("late_bg", "late", bg, binsize = 10))
        expect_length(w, 1)
        expect_match(w, "chr7")
        # the junk names blew past the cap, so the count is reported as approximate
        expect_match(w, "^100\\+ chromosome name")
    })
})

test_that("the skipped-name count is only marked approximate when names were dropped", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        # exactly at the cap, with one name repeated: nothing was dropped, so no "+"
        junk <- sprintf("READ%05d\t0\t100\t1.0", 1:100)
        bg <- write_lines_to("exact_cap.bedgraph", c(
            "chr1\t0\t100\t1.0", junk, "READ00001\t200\t300\t1.0"
        ))

        m <- capture_messages(gtrack.import("cap_bg", "cap", bg, binsize = 10))
        m <- paste(m, collapse = "")
        expect_match(m, "100 chromosome name")
        expect_false(grepl("100+ chromosome name", m, fixed = TRUE))
    })
})

test_that("a warning naming more primary chromosomes than it can print says how many", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        bg <- write_lines_to("many_primary.bedgraph", c(
            "chr1\t0\t100\t1.0",
            sprintf("chr%d\t0\t100\t5.0", 3:10)
        ))

        w <- capture_warnings(gtrack.import("many_bg", "many", bg, binsize = 10))
        expect_length(w, 1)
        # 8 primary names dropped, 5 printed
        expect_match(w, "(and 3 more)", fixed = TRUE)
    })
})

test_that("a file with no chromosome records at all still imports silently", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        # Deliberate carve-out: nothing was skipped and no name mismatched, so
        # there is nothing to report. An empty sample in a batch keeps working.
        wig <- write_lines_to("empty.wig", c("track type=wiggle_0 name=\"empty\""))

        expect_no_warning(expect_no_message(
            gtrack.import("empty_wig", "empty", wig, binsize = 10)
        ))
        expect_true(gtrack.exists("empty_wig"))
        expect_true(all(is.na(gextract("empty_wig", gintervals.all(), colnames = "v")$v)))
    })
})

test_that("diagnostics name the file the user asked for, not the temporary copy", {
    skip_if(Sys.which("gunzip") == "", "gunzip not on PATH")
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        plain <- write_lines_to("gzipped.bedgraph", c("scaffold_1\t0\t100\t1.0"))
        system(sprintf("gzip -f %s", shQuote(plain)))
        gz <- paste0(plain, ".gz")
        expect_true(file.exists(gz))

        # gtrack.import unzips into a temp file that is gone by the time the user
        # reads the error; the message must still name the .gz they passed.
        err <- tryCatch(suppressMessages(gtrack.import("gz_nomatch", "gz", gz, binsize = 10)),
            error = function(e) conditionMessage(e)
        )
        expect_match(err, "gzipped.bedgraph.gz", fixed = TRUE)
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

test_that("a skipped-contig message queued inside gintervals.mapply's FUN is delivered, and does not displace the outer call's warning", {
    local_db_state()
    withr::local_options(gmultitasking = FALSE)
    withr::with_tempdir({
        setup_mismatch_db()
        base <- write_lines_to("base.bedgraph", c(
            "chr1\t0\t1000\t1.0",
            "chr2\t0\t1000\t2.0"
        ))
        gtrack.import("base", "base", base, binsize = 100)

        # A partial match: chr1 lands, scaffold_9 is skipped and reported as a message.
        partial <- write_lines_to("partial.bedgraph", c(
            "chr1\t0\t100\t1.0",
            "scaffold_9\t0\t100\t5.0"
        ))

        # Overlapping iterator for the outer call, so it queues a warning of its own
        # while FUN - and the import nested inside it - runs.
        pk <- gintervals("chr1", c(0, 50), c(100, 300))
        notes <- character()
        # The inner call's diagnostics have to be caught inside FUN: gintervals.mapply
        # evaluates FUN through R_tryEval, which resets R's handler stack.
        expect_warning(
            gintervals.mapply(
                function(x) {
                    withCallingHandlers(
                        gtrack.import("partial", "partial", partial, binsize = 10),
                        message = function(m) {
                            notes <<- c(notes, conditionMessage(m))
                            invokeRestart("muffleMessage")
                        }
                    )
                    1
                },
                "base", gintervals("chr1", 0, 300),
                iterator = pk
            ),
            "were merged into"
        )
        expect_equal(sum(grepl("scaffold_9", notes)), 1)
        expect_true(gtrack.exists("partial"))
    })
})

test_that("under options(warn = 2) a partial import aborts and leaves no track behind", {
    local_db_state()
    withr::with_tempdir({
        setup_mismatch_db()
        # chrX is primary-shaped, so it goes to the warning channel; warn = 2 turns
        # that into an error, and .gtrack.create_atomic rolls the half-built track back.
        bg <- write_lines_to("warn2.bedgraph", c(
            "chr1\t0\t100\t1.0",
            "chrX\t0\t100\t5.0"
        ))
        withr::local_options(warn = 2)
        expect_error(gtrack.import("warn2_bg", "warn2", bg, binsize = 10), "chrX")
        expect_false(gtrack.exists("warn2_bg"))
    })
})
