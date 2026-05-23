# Tests for bedGraph import on databases whose contig names do not
# start with "chr" (Ensembl-style, non-mammalian, organelles, etc.).
# Regression: the Wig parser used to gate bedGraph detection on a
# "chr" prefix check, which silently mis-parsed any 4-field line
# from such genomes as an INVALID_FORMAT error.

test_that("gtrack.import parses bedGraph with non-chr-prefixed contig names", {
    local_db_state()
    withr::with_tempdir({
        # Build a tiny DB with Ensembl-style chrom names.
        create_test_db(
            "nonchr_db",
            chrom_sizes = data.frame(
                chrom = c("1", "Y", "MT"),
                size = c(1000L, 500L, 200L)
            )
        )
        gdb.init("nonchr_db")

        # Write a bedGraph using those bare contig names. The whole
        # tempdir (including the DB and any tracks created inside it)
        # is removed by withr::with_tempdir at end of block, so we
        # don't need an explicit gtrack.rm cleanup.
        bg <- file.path(getwd(), "nonchr.bedgraph")
        writeLines(
            c(
                "track type=bedGraph name=\"nonchr\"",
                "1\t0\t100\t1.0",
                "1\t100\t200\t2.5",
                "Y\t10\t30\t-1.0",
                "MT\t0\t50\t3.5"
            ),
            bg
        )

        # Import as a dense binsize=10 track. Before the fix, this
        # threw "Invalid format of WIG file" on the first data row
        # because the parser rejected the lack of a chr prefix.
        expect_no_error(
            gtrack.import("nonchr_bg", "non-chr bedGraph", bg, binsize = 10)
        )
        expect_true(gtrack.exists("nonchr_bg"))

        # Spot-check values via gextract on chrom "1".
        vals <- gextract("nonchr_bg",
            intervals = gintervals("1", 0, 200),
            iterator = "nonchr_bg", colnames = "v"
        )
        # First 100bp should be 1.0, next 100bp 2.5
        expect_equal(unique(vals$v[vals$start < 100]), 1.0)
        expect_equal(unique(vals$v[vals$start >= 100 & vals$start < 200]), 2.5)
    })
})
