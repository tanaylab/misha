# Regression tests for the database-state findings of the 2026-11 review of 5.14.4: interval sets,
# datasets, track copy/mv and the session state around them.

restore_groot_on_exit()

# A small indexed database of its own (chr1, chr2, chrX, chrM), built from a FASTA
review_db <- function(dir, name = "db") {
    fa <- file.path(dir, paste0(name, ".fa"))
    chroms <- c("chr1", "chr2", "chrX", "chrM")
    set.seed(17)
    seqs <- vapply(chroms, function(x) paste(sample(c("A", "C", "G", "T"), 5000, TRUE), collapse = ""), "")
    writeLines(paste0(">", chroms, "\n", seqs), fa)
    db <- file.path(dir, name)
    suppressMessages(gdb.create(groot = db, fasta = fa, format = "indexed", verbose = FALSE))
    suppressMessages(gsetroot(db))
    normalizePath(db)
}

test_that("gintervals.save stores alias chromosome names under their canonical names", {
    local_db_state()
    review_db(withr::local_tempdir())
    d <- data.frame(chrom = c("1", "X", "MT", "chr2"), start = c(0, 10, 20, 30), end = c(5, 15, 25, 35))

    gintervals.save("al_small", d)
    res <- gintervals.load("al_small")
    res <- res[order(res$start), ]
    expect_equal(as.character(res$chrom), c("chr1", "chrX", "chrM", "chr2"))
    expect_equal(res$end, d$end)

    # a big set (one file per chromosome) kept only the chr2 rows
    withr::local_options(gmax.data.size = 2)
    gintervals.save("al_big", rbind(d, d))
    withr::local_options(gmax.data.size = 1e8)
    res <- gintervals.load("al_big")
    expect_equal(nrow(res), 8)
    expect_equal(sort(as.character(res$chrom)), sort(rep(c("chr1", "chrX", "chrM", "chr2"), 2)))

    expect_error(gintervals.save("al_bad", data.frame(chrom = "chrNOPE", start = 0, end = 5)), "chrNOPE")
})

test_that("a big set saved in indexed format from a character chrom column reads per chromosome", {
    local_db_state()
    review_db(withr::local_tempdir())
    d <- data.frame(chrom = rep(c("chr1", "chr2", "chrX", "chrM"), each = 20), start = rep(seq(0, 1900, 100), 4), stringsAsFactors = FALSE)
    d$end <- d$start + 50
    withr::local_options(gmulticontig.indexed_format = TRUE, gmax.data.size = 30)
    gintervals.save("ixc", d)
    expect_true(file.exists(file.path(.misha$GROOT, "tracks", "ixc.interv", "intervals.idx")))

    res <- gintervals.load("ixc", chrom = "chrX")
    expect_equal(nrow(res), 20)
    expect_equal(res$start, seq(0, 1900, 100))
    withr::local_options(gmax.data.size = 1e8)
    expect_true(is.factor(gintervals.load("ixc")$chrom))
})

test_that("gintervals.update takes a data frame with a character chrom column, and aliases", {
    local_db_state()
    review_db(withr::local_tempdir())
    gintervals.save("up", gintervals(c("chr1", "chr2"), 0, 100))

    gintervals.update("up", data.frame(chrom = "chr2", start = 200, end = 300, stringsAsFactors = FALSE), chrom = "chr2")
    res <- gintervals.load("up")
    expect_equal(res$start[as.character(res$chrom) == "chr2"], 200)

    # by an alias of the chromosome, in a big set
    withr::local_options(gmax.data.size = 2)
    gintervals.save("upb", gintervals(c("chr1", "chrX", "chr2"), 0, 100))
    gintervals.update("upb", data.frame(chrom = "X", start = c(500, 700), end = c(600, 800)), chrom = "chrX")
    withr::local_options(gmax.data.size = 1e8)
    res <- gintervals.load("upb")
    expect_equal(res$start[as.character(res$chrom) == "chrX"], c(500, 700))
    expect_equal(nrow(res), 4)
})

test_that("gintervals.update refuses an interval set of a loaded dataset", {
    local_db_state()
    td <- withr::local_tempdir()
    ds <- review_db(td, "ds")
    gintervals.save("dsmall", gintervals("chrX", 300, 400))
    withr::local_options(gmax.data.size = 2)
    gintervals.save("dbig", gintervals(c("chr1", "chrX", "chr2"), 300, 400))
    withr::local_options(gmax.data.size = 1e8)
    work <- review_db(td, "work")
    gdataset.load(ds)

    expect_error(gintervals.update("dsmall", gintervals("chrX", 1000, 2000), chrom = "chrX"), "read-only dataset")
    expect_error(gintervals.update("dbig", gintervals("chrX", 1000, 2000), chrom = "chrX"), "read-only dataset")
    expect_false(file.exists(file.path(work, "tracks", "dsmall.interv")))
    expect_false(file.exists(file.path(work, "tracks", "dbig.interv")))
    expect_equal(gintervals.load("dsmall")$start, 300)
})

test_that("interval set attributes survive gintervals.update moving a set between small and big", {
    local_db_state()
    db <- review_db(withr::local_tempdir())
    withr::local_options(gmax.data.size = 30)
    gintervals.save("at1", gintervals("chr1", seq(0, 900, 100), seq(50, 950, 100)))
    expect_false(dir.exists(file.path(db, "tracks", "at1.interv")))
    gintervals.attr.set("at1", "source", "lab")

    # small to big
    gintervals.update("at1", gintervals("chrX", seq(0, 2900, 100), seq(50, 2950, 100)), chrom = "chrX")
    expect_true(dir.exists(file.path(db, "tracks", "at1.interv")))
    expect_equal(gintervals.attr.get("at1", "source"), "lab")
    expect_false(file.exists(file.path(db, "tracks", "at1.iattr")))

    # big to small
    gintervals.update("at1", NULL, chrom = "chrX")
    expect_false(dir.exists(file.path(db, "tracks", "at1.interv")))
    expect_equal(gintervals.attr.get("at1", "source"), "lab")
})

test_that("gtrack.mv and gtrack.copy refuse the name of an interval set", {
    local_db_state()
    db <- review_db(withr::local_tempdir())
    gtrack.create_sparse("t1", "", gintervals("chr1", 0, 10), 1)
    gintervals.save("iset", gintervals("chr1", 0, 100))

    expect_error(gtrack.mv("t1", "iset"), "iset already exists")
    expect_true(dir.exists(file.path(db, "tracks", "t1.track")))
    expect_false(dir.exists(file.path(db, "tracks", "iset.track")))
    expect_error(gtrack.copy("t1", "iset"), "iset already exists")
    expect_false(dir.exists(file.path(db, "tracks", "iset.track")))

    # the database still opens
    suppressMessages(gsetroot(db, rescan = TRUE))
    expect_true("t1" %in% gtrack.ls())
})

test_that("gdataset.load refuses a track that is an interval set elsewhere, leaving the session as it was", {
    local_db_state()
    td <- withr::local_tempdir()
    ds <- review_db(td, "ds")
    gtrack.create_sparse("x", "", gintervals("chr1", 0, 10), 1)
    work <- review_db(td, "work")
    gintervals.save("x", gintervals("chr1", 0, 100))
    gtrack.create_sparse("mine", "", gintervals("chr1", 0, 10), 2)
    datasets <- gdataset.ls()

    expect_error(gdataset.load(ds), "'x' would be both a track and an interval set")
    expect_equal(gtrack.ls(), "mine")
    expect_equal(gintervals.ls(), "x")
    expect_equal(gdataset.ls(), datasets)
    expect_equal(gextract("mine", gintervals("chr1", 0, 10))$mine, 2)
})

test_that("ggenome.implant from a donor database leaves the session as it was", {
    local_db_state()
    td <- withr::local_tempdir()
    donor <- review_db(td, "donor")
    ds <- review_db(td, "ds")
    gtrack.create_sparse("dtrack", "", gintervals("chr1", 0, 10), 5)
    work <- review_db(td, "work")
    gdataset.load(ds)
    datasets <- gdataset.ls()
    expect_true(ds %in% datasets)
    fa <- file.path(td, "ref.fa")
    gdb.export_fasta(fa)

    # failure
    expect_error(ggenome.implant(data.frame(chrom = "chrNOPE", start = 0, end = 10), donor, file.path(td, "o1.fa"), genome_fasta = fa, create_trackdb = FALSE))
    expect_equal(.misha$GROOT, work)
    expect_equal(gdataset.ls(), datasets)

    # success
    ggenome.implant(data.frame(chrom = "chr1", start = 0, end = 10), donor, file.path(td, "o2.fa"), genome_fasta = fa, create_trackdb = FALSE)
    expect_equal(.misha$GROOT, work)
    expect_equal(gdataset.ls(), datasets)
    expect_equal(gextract("dtrack", gintervals("chr1", 0, 10))$dtrack, 5)

    # gdb.export_fasta of another database, the same way
    gdb.export_fasta(file.path(td, "donor_export.fa"), groot = donor)
    expect_equal(.misha$GROOT, work)
    expect_equal(gdataset.ls(), datasets)
})

test_that("a failed creation over the name of a dataset track leaves that track readable", {
    local_db_state()
    td <- withr::local_tempdir()
    ds <- review_db(td, "ds")
    gtrack.create_sparse("dtrack", "", gintervals("chr1", 0, 10), 5)
    work <- review_db(td, "work")
    gdataset.load(ds)

    expect_error(gtrack.create("dtrack", "", "no_such_track + 1"))
    expect_equal(gtrack.dataset("dtrack"), ds)
    expect_equal(gextract("dtrack", gintervals("chr1", 0, 10))$dtrack, 5)
    expect_false(dir.exists(file.path(work, "tracks", "dtrack.track")))

    # a successful one writes into the working database, which then wins
    gtrack.create_sparse("dtrack", "", gintervals("chr1", 0, 10), 7)
    expect_equal(gtrack.dataset("dtrack"), work)
    expect_equal(gextract("dtrack", gintervals("chr1", 0, 10))$dtrack, 7)
    expect_equal(sum(gtrack.ls() == "dtrack"), 1)
})
