# End-to-end tests for BAM auto-detect in gtrack.import_mappedseq.
# All tests skip cleanly when samtools is not on PATH.

create_isolated_test_db()

test_that("gtrack.import_mappedseq imports BAM as sparse track", {
    bam <- default_bam_path()

    tmptrack <- random_track_name("test")
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.rm(tmptrack, force = TRUE)

    stats <- gtrack.import_mappedseq(
        tmptrack, "BAM sparse", bam,
        cols.order = NULL, remove.dups = TRUE
    )
    # 2 mapped records on chrom "chr1", 1 unmapped record (FLAG 4 -> chrom "*").
    expect_equal(stats[[1]][["total.mapped"]], 2)
    expect_equal(stats[[1]][["total.unmapped"]], 1)
})

test_that("gtrack.import_mappedseq imports BAM as dense pileup", {
    bam <- default_bam_path()

    tmptrack <- random_track_name("test")
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.rm(tmptrack, force = TRUE)

    stats <- gtrack.import_mappedseq(
        tmptrack, "BAM dense", bam,
        pileup = 10, binsize = 10,
        cols.order = NULL, remove.dups = TRUE
    )
    expect_equal(stats[[1]][["total.mapped"]], 2)

    intervs <- gintervals("chr1", c(100, 200), c(110, 210))
    r <- gextract(tmptrack, intervs, iterator = 10, colnames = "v")
    # Each read covers exactly its bin -> v == 1.
    expect_equal(as.numeric(r$v), c(1, 1))
})

test_that("gtrack.import_mappedseq auto-switches default cols.order for BAM", {
    bam <- default_bam_path()

    tmptrack <- random_track_name("test")
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.rm(tmptrack, force = TRUE)

    # Caller does NOT pass cols.order. The legacy default c(9, 11, 13, 14)
    # is the tab-delimited layout; samtools view emits SAM (columns 10/3/4/2).
    # The wrapper must treat the default as "user didn't pass anything" and
    # switch to NULL; otherwise the import returns 0 mapped reads.
    stats <- gtrack.import_mappedseq(
        tmptrack, "BAM default cols.order", bam,
        remove.dups = TRUE
    )
    expect_equal(stats[[1]][["total.mapped"]], 2)
})

test_that("gtrack.import_mappedseq errors on explicit cols.order with BAM", {
    bam <- default_bam_path()

    tmptrack <- random_track_name("test")
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    expect_error(
        gtrack.import_mappedseq(
            tmptrack, "BAM explicit cols.order", bam,
            cols.order = c(9, 11, 13, 14), remove.dups = TRUE
        ),
        "BAM input forces SAM column layout"
    )
})

test_that("gtrack.import_mappedseq surfaces samtools-not-found for BAM", {
    bam <- default_bam_path()

    tmptrack <- random_track_name("test")
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Clear PATH so the child shell can't find samtools.
    err <- tryCatch(
        withr::with_envvar(
            c(PATH = "/nonexistent"),
            gtrack.import_mappedseq(
                tmptrack, "BAM no samtools", bam,
                cols.order = NULL, remove.dups = TRUE
            )
        ),
        error = function(e) conditionMessage(e)
    )
    expect_match(err, "samtools is not on PATH")
    expect_match(err, "conda install")
})

test_that("gtrack.import_mappedseq accepts gzipped SAM", {
    tmptrack <- random_track_name("test")
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.rm(tmptrack, force = TRUE)

    sam_gz <- tempfile(fileext = ".sam.gz")
    con <- gzfile(sam_gz, "w")
    writeLines(default_sam_text(), con)
    close(con)

    stats <- gtrack.import_mappedseq(
        tmptrack, "gzipped SAM", sam_gz,
        cols.order = NULL, remove.dups = TRUE
    )
    expect_equal(stats[[1]][["total.mapped"]], 2)
    expect_equal(stats[[1]][["total.unmapped"]], 1)
})
