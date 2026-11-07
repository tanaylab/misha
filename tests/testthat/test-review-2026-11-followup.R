# Defects in 5.16.0 found by comparing R misha with pymisha

restore_groot_on_exit()

# A small database of its own (chr1, chr2), per-chromosome or indexed
small_db <- function(indexed) {
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    if (!indexed) {
        create_test_db(dir)
        suppressMessages(gsetroot(dir))
        return(normalizePath(dir))
    }
    fa <- file.path(dir, "g.fa")
    set.seed(3)
    writeLines(paste0(">", c("chr1", "chr2"), "\n", replicate(2, paste(sample(c("A", "C", "G", "T"), 10000, TRUE), collapse = ""))), fa)
    db <- file.path(dir, "db")
    suppressMessages(gdb.create(groot = db, fasta = fa, format = "indexed", verbose = FALSE))
    suppressMessages(gsetroot(db))
    normalizePath(db)
}

# Runs expr in a forked process, so that a crash fails the test instead of the test run
in_child <- function(expr) {
    job <- parallel::mcparallel(withr::with_options(list(gmultitasking = FALSE), tryCatch(expr, error = function(e) e)))
    res <- parallel::mccollect(job)[[1]]
    if (is.null(res)) stop("the child process died")
    if (inherits(res, "try-error")) stop("the child process failed: ", res)
    res
}

test_that("a sequence virtual track next to a sparse track picks the sparse track's iterator", {
    skip_on_cran()
    skip_on_os("windows")
    gvtrack.create("followup_kmer", NULL, "kmer.count", kmer = "ACG")
    withr::defer(gvtrack.rm("followup_kmer"))
    intervs <- gintervals(1, 0, 100000)
    expected <- gextract("test.sparse + followup_kmer", intervs, iterator = "test.sparse")
    for (expr in c("test.sparse + followup_kmer", "followup_kmer + test.sparse")) {
        res <- in_child(gextract(expr, intervs))
        expect_s3_class(res, "data.frame")
        expect_equal(res[[expr]], expected[["test.sparse + followup_kmer"]])
    }
})

test_that("a sequence virtual track next to a 2D track gives an error, not a crash", {
    skip_on_cran()
    skip_on_os("windows")
    local_db_state()
    db <- withr::local_tempdir()
    create_test_db(db)
    suppressMessages(gsetroot(db))
    gtrack.2d.create("rects", "", data.frame(chrom1 = "chr1", start1 = 100, end1 = 300, chrom2 = "chr1", start2 = 200, end2 = 400), 5)
    gvtrack.create("followup_kmer2d", NULL, "kmer.count", kmer = "ACG")
    withr::defer(gvtrack.rm("followup_kmer2d"))
    for (expr in c("rects + followup_kmer2d", "followup_kmer2d + rects")) {
        res <- in_child(gextract(expr, gintervals.2d.all()))
        expect_s3_class(res, "error")
        expect_match(conditionMessage(res), "2D iterator is used without explicit dimension projection")
    }
})

test_that("a numeric strand other than -1, 0 or 1 is an error", {
    intervs <- data.frame(chrom = "chr1", start = 200, end = 210)
    fwd <- gseq.extract(cbind(intervs, strand = 1))
    expect_equal(gseq.extract(cbind(intervs, strand = 1L)), fwd)
    expect_equal(gseq.extract(cbind(intervs, strand = 0)), fwd)
    for (s in c(255, 511, 2, 0.5, -2)) {
        expect_error(gseq.extract(cbind(intervs, strand = s)), "Invalid strand value")
    }
    expect_error(gseq.extract(cbind(intervs, strand = 255L)), "Invalid strand value")
})

test_that("gtrack.import_parquet reads empty shards written by pyarrow and pandas", {
    local_db_state()
    db <- withr::local_tempdir()
    create_test_db(db)
    suppressMessages(gsetroot(db))
    # (fixtures/parquet/make_fixtures.py: one row group of no rows, column chunks with
    # data_page_offset 0; pandas also types the chrom column null)
    fixtures <- test_path("fixtures", "parquet")
    gtrack.import_parquet("", file.path(fixtures, "shard_chr1.parquet"), 32, track.prefix = "one_")
    for (empty in c("empty_pyarrow.parquet", "empty_pandas.parquet")) {
        dir <- withr::local_tempdir()
        file.copy(file.path(fixtures, c("shard_chr1.parquet", empty)), dir)
        prefix <- paste0(sub("\\.parquet$", "", empty), "_")
        gtrack.import_parquet("", file.path(dir, "*.parquet"), 32, track.prefix = prefix)
        expect_equal(
            gextract(paste0(prefix, "a"), gintervals.all(), colnames = "v")$v,
            gextract("one_a", gintervals.all(), colnames = "v")$v
        )
    }
    expect_equal(gextract("one_a", gintervals("chr1", 0, 64))$one_a, c(0.25, 1.25))
})

test_that("gtrack.create and gtrack.lookup of a 2D expression that give no rectangles create no track", {
    local_db_state()
    for (indexed in c(FALSE, TRUE)) {
        db <- small_db(indexed)
        gtrack.2d.create("rA", "", data.frame(chrom1 = "chr1", start1 = 100, end1 = 300, chrom2 = "chr1", start2 = 200, end2 = 400), 5)
        gtrack.2d.create("rB", "", data.frame(chrom1 = "chr1", start1 = 200, end1 = 250, chrom2 = "chr1", start2 = 300, end2 = 350), 11)
        for (mt in c(FALSE, TRUE)) {
            withr::local_options(gmultitasking = mt)
            # no rectangle of rA is in this band
            expect_error(gtrack.create("e1", "", "rB", iterator = "rA", band = c(10000, 20000)), "gave no rectangles")
            expect_false(dir.exists(file.path(db, "tracks", "e1.track")))
            expect_false("e1" %in% gtrack.ls())
            expect_error(gtrack.lookup("l1", "", 1:3, "rB", c(0, 1, 2, 1e9), iterator = "rA", band = c(10000, 20000)), "gave no rectangles")
            expect_false(dir.exists(file.path(db, "tracks", "l1.track")))
            expect_false("l1" %in% gtrack.ls())
        }
        # (a rectangle in the scope still gives a track)
        gtrack.create("e2", "", "rB", iterator = "rA")
        expect_equal(nrow(gextract("e2", gintervals.2d.all())), 1)
    }
})
