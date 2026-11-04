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
