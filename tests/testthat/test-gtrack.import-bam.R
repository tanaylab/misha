# End-to-end tests for BAM auto-detect in gtrack.import_mappedseq.
# All tests skip cleanly when samtools is not on PATH.

create_isolated_test_db()

test_that("gtrack.import_mappedseq imports BAM as sparse track", {
    skip_if_no_samtools()
    bam <- make_test_bam(default_sam_text())

    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.rm(tmptrack, force = TRUE)

    stats <- gtrack.import_mappedseq(
        tmptrack, "BAM sparse", bam,
        cols.order = NULL, remove.dups = TRUE
    )
    # 2 mapped records on chrom "1", 1 unmapped record (FLAG 4 -> chrom "*").
    expect_equal(stats[[1]][["total.mapped"]], 2)
    expect_equal(stats[[1]][["total.unmapped"]], 1)
})

test_that("gtrack.import_mappedseq imports BAM as dense pileup", {
    skip_if_no_samtools()
    bam <- make_test_bam(default_sam_text())

    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
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
    skip_if_no_samtools()
    bam <- make_test_bam(default_sam_text())

    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.rm(tmptrack, force = TRUE)

    # NB: caller does NOT pass cols.order. The legacy default
    # c(9, 11, 13, 14) is the tab-delimited layout; samtools view emits SAM
    # (columns 10/3/4/2). The wrapper must detect BAM and silently switch
    # to NULL (SAM mode); otherwise the import returns 0 mapped reads.
    stats <- gtrack.import_mappedseq(
        tmptrack, "BAM default cols.order", bam,
        remove.dups = TRUE
    )
    expect_equal(stats[[1]][["total.mapped"]], 2)
})

test_that("gtrack.import_mappedseq surfaces samtools-not-found for BAM", {
    skip_if_no_samtools()
    bam <- make_test_bam(default_sam_text())

    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Clear PATH so the child shell can't find samtools.
    withr::with_envvar(
        c(PATH = "/nonexistent"),
        expect_error(
            gtrack.import_mappedseq(
                tmptrack, "BAM no samtools", bam,
                cols.order = NULL, remove.dups = TRUE
            ),
            "samtools is not on PATH"
        )
    )
})

test_that("gtrack.import_mappedseq accepts gzipped SAM", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
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
