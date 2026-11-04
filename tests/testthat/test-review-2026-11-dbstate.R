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

test_that("removing a working interval set shows the dataset's set of that name at once", {
    local_db_state()
    td <- withr::local_tempdir()
    d1 <- review_db(td, "d1")
    gintervals.save("shared", gintervals("chr1", 100, 200))
    d2 <- review_db(td, "d2")
    gintervals.save("shared", gintervals("chr1", 300, 400))
    work <- review_db(td, "work")
    gintervals.save("shared", gintervals("chr1", 500, 600))
    gdataset.load(d1, force = TRUE)
    gdataset.load(d2, force = TRUE)
    expect_equal(gintervals.load("shared")$start, 500)

    gintervals.rm("shared", force = TRUE)
    # the dataset loaded last wins, as gdb.reload() has it
    expect_true(gintervals.exists("shared"))
    expect_equal(gintervals.load("shared")$start, 300)
    suppressMessages(gdb.reload())
    expect_equal(gintervals.load("shared")$start, 300)
})

test_that("removing a working track shows the track of the dataset loaded last", {
    local_db_state()
    td <- withr::local_tempdir()
    d1 <- review_db(td, "d1")
    gtrack.create_sparse("tt", "", gintervals("chr1", 0, 10), 1)
    d2 <- review_db(td, "d2")
    gtrack.create_sparse("tt", "", gintervals("chr1", 0, 10), 2)
    work <- review_db(td, "work")
    gtrack.create_sparse("tt", "", gintervals("chr1", 0, 10), 3)
    gdataset.load(d1, force = TRUE)
    gdataset.load(d2, force = TRUE)

    gtrack.rm("tt", force = TRUE)
    expect_equal(gextract("tt", gintervals("chr1", 0, 10))$tt, 2)
    suppressMessages(gdb.reload())
    expect_equal(gextract("tt", gintervals("chr1", 0, 10))$tt, 2)
})

test_that("gdataset.save keeps the attributes of a small interval set", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- review_db(td, "src")
    gintervals.save("ii", gintervals("chr1", 0, 100))
    gintervals.attr.set("ii", "who", "DA")
    gdataset.save(file.path(td, "DS"), "x", intervals = "ii")
    expect_true(file.exists(file.path(td, "DS", "tracks", "ii.iattr")))

    review_db(td, "other")
    gdataset.load(file.path(td, "DS"))
    expect_equal(gintervals.attr.get("ii", "who"), "DA")
})

test_that("gdb.convert_to_indexed converts a 1D set over a contig whose name has a dash as 1D", {
    local_db_state()
    db <- create_test_db(file.path(withr::local_tempdir(), "db"), data.frame(chrom = c("chr1", "HLA-A", "scaf-2"), size = c(10000, 5000, 5000)))
    suppressMessages(gsetroot(db))
    d <- data.frame(chrom = rep(c("chr1", "HLA-A", "scaf-2"), each = 5), start = rep(seq(0, 400, 100), 3))
    d$end <- d$start + 50
    withr::with_options(list(gmax.data.size = 4), gintervals.save("dash", d))
    expect_true(file.exists(file.path(db, "tracks", "dash.interv", "HLA-A")))

    expect_no_warning(suppressMessages(gdb.convert_to_indexed(convert_intervals = TRUE, force = TRUE, remove_old_files = TRUE)))
    expect_true(file.exists(file.path(db, "tracks", "dash.interv", "intervals.idx")))
    suppressMessages(gsetroot(db))
    expect_equal(nrow(gintervals.load("dash", chrom = "HLA-A")), 5)
    expect_equal(nrow(gintervals.load("dash")), 15)
})

# gtrack.copy(overwrite = TRUE) SIGKILLed in a forked child, as an OOM kill or a job limit: no R
# handler runs. Whatever point it dies at, the real name holds the old or the new track complete, or
# the old one is in a directory no later cleanup removes.
copy_killed_at <- function(kill_at) {
    td <- withr::local_tempdir(.local_envir = parent.frame())
    db <- review_db(td)
    gtrack.create_sparse("orig", "", gintervals("chr1", 0, 10), 42)
    gtrack.create_sparse("src", "", gintervals("chr1", 100, 110), 7)
    tracks <- file.path(db, "tracks")
    orig_bytes <- readBin(file.path(tracks, "orig.track", "track.dat"), "raw", 1e6)

    job <- parallel::mcparallel({
        rename <- base::file.rename
        unlockBinding("file.rename", .BaseNamespaceEnv)
        assign("file.rename", function(from, to) {
            # the first rename out of the real name (the track set aside), or into it (the copy put in place)
            if ((kill_at == "aside" && basename(from[1]) == "orig.track") ||
                (kill_at == "swap" && basename(to[1]) == "orig.track")) {
                res <- rename(from, to)
                tools::pskill(Sys.getpid(), tools::SIGKILL)
            }
            rename(from, to)
        }, envir = .BaseNamespaceEnv)
        if (kill_at == "copy") {
            copy <- base::file.copy
            unlockBinding("file.copy", .BaseNamespaceEnv)
            assign("file.copy", function(...) {
                res <- copy(...)
                tools::pskill(Sys.getpid(), tools::SIGKILL)
                res
            }, envir = .BaseNamespaceEnv)
        }
        gtrack.copy("src", "orig", overwrite = TRUE)
    })
    expect_warning(parallel::mccollect(job), "did not deliver a result")

    # a later write in the directory sweeps what the dead process left
    suppressMessages(gsetroot(db))
    gtrack.create_sparse("another", "", gintervals("chr1", 0, 10), 1)
    leftovers <- setdiff(list.files(tracks, all.files = TRUE, no.. = TRUE), c("orig.track", "src.track", "another.track"))
    list(db = db, tracks = tracks, leftovers = leftovers, orig_bytes = orig_bytes)
}

test_that("gtrack.copy(overwrite = TRUE) killed during the copy leaves the old track in place", {
    skip_on_cran()
    skip_on_os("windows")
    local_db_state()
    st <- copy_killed_at("copy")
    expect_equal(gextract("orig", gintervals("chr1", 0, 200))$orig, 42)
    expect_length(st$leftovers, 0)
})

test_that("gtrack.copy(overwrite = TRUE) killed between its renames keeps the old track where no cleanup removes it", {
    skip_on_cran()
    skip_on_os("windows")
    local_db_state()
    st <- copy_killed_at("aside")
    expect_false(dir.exists(file.path(st$tracks, "orig.track")))
    expect_length(st$leftovers, 1)
    kept <- file.path(st$tracks, st$leftovers)
    expect_equal(readBin(file.path(kept, "track.dat"), "raw", 1e6), st$orig_bytes)
    # put back by hand, it is the track as it was
    file.rename(kept, file.path(st$tracks, "orig.track"))
    suppressMessages(gdb.reload())
    expect_equal(gextract("orig", gintervals("chr1", 0, 200))$orig, 42)
})

test_that("gtrack.copy(overwrite = TRUE) killed once the copy is in place leaves the new track", {
    skip_on_cran()
    skip_on_os("windows")
    local_db_state()
    st <- copy_killed_at("swap")
    suppressMessages(gdb.reload())
    res <- gextract("orig", gintervals("chr1", 0, 200))
    expect_equal(res$start, 100)
    expect_equal(res$orig, 7)
})

# Another session on the same database: a forked child, which starts with a copy of this session's
# state and writes to the same files
in_other_session <- function(expr) {
    job <- parallel::mcparallel(withr::with_options(list(gmultitasking = FALSE), {
        expr
        TRUE
    }))
    res <- parallel::mccollect(job)[[1]]
    if (!isTRUE(res)) stop("the other session failed: ", format(res))
    invisible()
}

read_cache <- function(db) {
    f <- file(file.path(db, ".db.cache"), "rb")
    on.exit(close(f))
    unserialize(f)
}

test_that(".db.cache keeps the tracks and interval sets another session created or removed", {
    skip_on_os("windows")
    local_db_state()
    db <- review_db(withr::local_tempdir())
    gtrack.create_sparse("t1", "", gintervals("chr1", 0, 10), 1)

    in_other_session({
        gtrack.create_sparse("t2", "", gintervals("chr1", 0, 10), 2)
        gintervals.save("iB", gintervals("chr1", 0, 10))
    })
    # this session's lists predate t2 and iB
    gtrack.create_sparse("t3", "", gintervals("chr1", 0, 10), 3)
    gintervals.save("iA", gintervals("chr1", 0, 10))
    expect_setequal(read_cache(db)[[1]], c("t1", "t2", "t3"))
    expect_setequal(read_cache(db)[[2]], c("iA", "iB"))
    expect_false(.gdb.cache_is_dirty(db))

    # a new session reads the cache
    suppressMessages(gsetroot(db))
    expect_setequal(gtrack.ls(), c("t1", "t2", "t3"))
    expect_setequal(gintervals.ls(), c("iA", "iB"))

    # the other session removes t2, then this one writes
    in_other_session(gtrack.rm("t2", force = TRUE))
    gtrack.create_sparse("t4", "", gintervals("chr1", 0, 10), 4)
    expect_setequal(read_cache(db)[[1]], c("t1", "t3", "t4"))
})

test_that("a write leaves a .db.cache that another session marked dirty dirty, and a rescan clears it", {
    local_db_state()
    db <- review_db(withr::local_tempdir())
    dirty <- file.path(db, ".db.cache.dirty")
    writeLines("set by another session", dirty)

    gtrack.create_sparse("t5", "", gintervals("chr1", 0, 10), 5)
    expect_equal(readLines(dirty), "set by another session")
    gintervals.rm("nothing_here", force = TRUE)
    expect_true(file.exists(dirty))

    suppressMessages(gdb.reload())
    expect_false(file.exists(dirty))
    expect_true("t5" %in% read_cache(db)[[1]])
})
