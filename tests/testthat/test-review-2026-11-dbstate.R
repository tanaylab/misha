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
    # a big set saved in per-chromosome files, to be converted
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
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

    flags <- function() list.files(db, pattern = "^\\.db\\.cache\\.dirty", all.files = TRUE)
    gtrack.create_sparse("t5", "", gintervals("chr1", 0, 10), 5)
    # a flag of this writer's own is added: whoever set the other cannot clear the cache while t5
    # is missing too
    expect_true(file.exists(dirty))
    expect_length(flags(), 2)
    expect_false(.gdb.cache_clear_dirty(db, "set by another session"))
    gintervals.rm("nothing_here", force = TRUE)
    expect_length(flags(), 3)

    # a rescan removes the flags that were there when it started
    suppressMessages(gdb.reload())
    expect_length(flags(), 0)
    expect_true("t5" %in% read_cache(db)[[1]])
})

test_that("a forced re-conversion of an indexed set refuses to drop what is only in its index", {
    local_db_state()
    db <- review_db(withr::local_tempdir())
    d <- data.frame(chrom = rep(c("chr1", "chr2", "chrX", "chrM"), each = 20), start = rep(seq(0, 1900, 100), 4))
    d$end <- d$start + 50
    withr::with_options(list(gmax.data.size = 50), gintervals.save("u3", d))
    gintervals.convert_to_indexed("u3", remove.old = TRUE)
    dir <- file.path(db, "tracks", "u3.interv")
    before <- tools::md5sum(file.path(dir, c("intervals.dat", "intervals.idx")))

    expect_error(gintervals.convert_to_indexed("u3", force = TRUE), "Nothing was changed")
    expect_equal(tools::md5sum(file.path(dir, c("intervals.dat", "intervals.idx"))), before)
    expect_equal(nrow(gintervals.load("u3")), 80)

    # one chromosome updated since (its own file), the others only in the index
    withr::with_options(list(gmax.data.size = 50), gintervals.update("u3", gintervals("chrX", c(0, 500), c(10, 510)), chrom = "chrX"))
    expect_error(gintervals.convert_to_indexed("u3", force = TRUE), "chr1")
    expect_equal(nrow(gintervals.load("u3")), 62)

    # 2D
    d2 <- data.frame(chrom1 = rep(c("chr1", "chr2"), each = 10), start1 = rep(seq(0, 900, 100), 2), chrom2 = rep(c("chr1", "chrX"), each = 10), start2 = 0)
    d2$end1 <- d2$start1 + 50
    d2$end2 <- d2$start2 + 50
    d2 <- d2[, c("chrom1", "start1", "end1", "chrom2", "start2", "end2")]
    withr::with_options(list(gmax.data.size = 15), gintervals.save("u2d", d2))
    gintervals.2d.convert_to_indexed("u2d", remove.old = TRUE)
    expect_error(gintervals.2d.convert_to_indexed("u2d", force = TRUE), "Nothing was changed")
    expect_equal(nrow(gintervals.load("u2d")), 20)
})

test_that("2D computations over the deferred whole-genome 2D scope stop and name gmulticontig.2d.mode", {
    local_db_state()
    withr::local_options(list(gmulticontig.2d.threshold = 2, gmulticontig.2d.mode = NULL))
    review_db(withr::local_tempdir())
    r <- data.frame(
        chrom1 = c("chr1", "chr1", "chr2", "chrX"), start1 = c(0, 0, 0, 0), end1 = c(100, 100, 100, 100),
        chrom2 = c("chr1", "chr2", "chrX", "chrM"), start2 = c(200, 0, 0, 0), end2 = c(300, 100, 100, 100)
    )
    gtrack.2d.create("mr", "", r, c(10, 20, 30, 40))
    gtrack.create_sparse("t1d", "", gintervals(c("chr1", "chr2", "chrX"), 0, 100), c(1, 2, 3))
    msg <- 'options(gmulticontig.2d.mode = "full")'

    for (mt in c(TRUE, FALSE)) {
        withr::local_options(list(gmultitasking = mt))
        # gquantiles gave the median of the cis value only (10)
        expect_error(gquantiles("mr", 0.5), msg, fixed = TRUE)
        expect_error(gextract("mr", .misha$ALLGENOME), msg, fixed = TRUE)
        expect_error(gsummary("mr"), msg, fixed = TRUE)
        expect_error(gscreen("mr > 0"), msg, fixed = TRUE)
        expect_error(gdist("mr", c(0, 100)), msg, fixed = TRUE)

        # 1D computations over the whole-genome scope are not affected (gsummary, gscreen and gdist
        # failed on the placeholder)
        expect_equal(gsummary("t1d")[["Sum"]], 6)
        expect_equal(nrow(gscreen("t1d > 0")), 3)
        expect_equal(sum(gdist("t1d", c(0, 2, 4))), 3)
        expect_equal(unname(gquantiles("t1d", 0.5)), 2)
        expect_equal(nrow(gextract("t1d", .misha$ALLGENOME)), 3)
    }

    # gintervals.2d.all() itself gives the cis pairs, with a warning; what is made from them works
    expect_warning(all2d <- gintervals.2d.all(), msg, fixed = TRUE)
    expect_equal(nrow(all2d), 4)
    expect_null(attr(all2d, "deferred"))
    expect_warning(expect_equal(gextract("mr", gintervals.2d.all())$mr, 10), msg, fixed = TRUE)
    some <- suppressWarnings(head(gintervals.2d.all(), 2))
    expect_equal(nrow(gintervals.intersect(some, some)), 2)
    gintervals.save("some2d", some)
    expect_equal(nrow(gintervals.load("some2d")), 2)
    gtrack.2d.create("some_track", "", some, c(1, 2))
    expect_equal(gextract("some_track", some)$some_track, c(1, 2))

    # track creation over the whole genome too
    expect_error(gtrack.create("c2d", "", "mr"), msg, fixed = TRUE)

    # in a mode that was set, the whole-genome scope is built in it, everywhere, and only for a 2D
    # computation
    built <- 0
    real_generate <- .generate_2d_on_demand
    local_mocked_bindings(.generate_2d_on_demand = function(...) {
        built <<- built + 1
        real_generate(...)
    })
    withr::with_options(list(gmulticontig.2d.mode = "full"), {
        expect_equal(gsummary("t1d")[["Sum"]], 6)
        expect_equal(nrow(gextract("t1d", .misha$ALLGENOME)), 3)
        expect_equal(unname(gquantiles("t1d", 0.5)), 2)
        expect_equal(nrow(gscreen("t1d > 0")), 3)
        gtrack.create("c1d", "", "t1d", iterator = 50)
        expect_equal(built, 0)
        gtrack.create("c2d", "", "mr")
        expect_equal(sort(gextract("c2d", r)$c2d), c(10, 20, 30, 40))
        expect_gt(built, 0)
    })
    withr::with_options(list(gmulticontig.2d.mode = "full"), {
        suppressWarnings({
            expect_equal(nrow(gextract("mr", gintervals.2d.all())), 4)
            expect_equal(unname(gquantiles("mr", 0.5)), median(c(10, 20, 30, 40)))
            expect_equal(gsummary("mr")[["Sum"]], 100)
            expect_equal(nrow(gscreen("mr > 0")), 4)
            expect_equal(sum(gdist("mr", c(0, 100))), 4)
            expect_equal(nrow(gextract("mr", .misha$ALLGENOME)), 4)
            expect_no_error(gcis_decay("mr", c(0, 1000), gintervals.all(), gintervals.all()))
        })
    })
    withr::with_options(list(gmulticontig.2d.mode = "diagonal"), {
        expect_equal(gextract("mr", gintervals.2d.all())$mr, 10)
        expect_equal(gsummary("mr")[["Sum"]], 10)
    })
})

test_that("1D track creation works on a database above gmulticontig.2d.threshold, in every 2D mode", {
    local_db_state()
    withr::local_options(list(gmulticontig.2d.threshold = 2))
    td <- withr::local_tempdir()
    review_db(td)
    bg <- file.path(td, "a.bedgraph")
    writeLines(c("chr1\t0\t100\t1", "chr2\t0\t100\t2"), bg)
    for (mode in list(NULL, "full", "diagonal")) {
        withr::local_options(list(gmulticontig.2d.mode = mode))
        tag <- if (is.null(mode)) "unset" else mode
        sp <- paste0("sp_", tag)
        gtrack.create_sparse(sp, "", gintervals(c("chr1", "chr2"), 0, 100), c(1, 2))
        gtrack.create(paste0("cr_", tag), "", sp, iterator = 50)
        expect_equal(gextract(paste0("cr_", tag), gintervals("chr2", 0, 50))[[paste0("cr_", tag)]], 2)
        gtrack.lookup(paste0("lk_", tag), "", c(10, 20), sp, c(0, 1.5, 3))
        expect_equal(gextract(paste0("lk_", tag), gintervals("chr2", 0, 50), iterator = 50)[[paste0("lk_", tag)]], 20)
        gtrack.smooth(paste0("sm_", tag), "", paste0("cr_", tag), 100, alg = "MEAN")
        expect_true(gtrack.exists(paste0("sm_", tag)))
        suppressMessages(gtrack.import(paste0("im_", tag), "", bg, 0))
        expect_equal(gextract(paste0("im_", tag), gintervals.all())[[paste0("im_", tag)]], c(1, 2))
    }
})

test_that("tracks and interval sets another session re-created with another layout are read right", {
    skip_on_os("windows")
    local_db_state()
    db <- review_db(withr::local_tempdir())
    gtrack.create_sparse("st", "", gintervals(c("chr1", "chr2"), 0, 10), c(1, 2))
    gtrack.create("dt", "", "st", iterator = 10)
    gtrack.2d.create("rt", "", data.frame(chrom1 = "chr1", start1 = 0, end1 = 10, chrom2 = "chr2", start2 = 0, end2 = 10), 1)
    withr::with_options(list(gmax.data.size = 1), gintervals.save("bs", gintervals(c("chr1", "chr2"), 0, 10)))
    gintervals.convert_to_indexed("bs", remove.old = TRUE)
    expect_true(file.exists(file.path(db, "tracks", "st.track", "track.idx")))
    expect_true(file.exists(file.path(db, "tracks", "bs.interv", "intervals.idx")))
    # this session caches their indexes
    gextract("st", gintervals.all())
    gextract("dt", gintervals.all())
    gextract("rt", gintervals.2d("chr1", 0, 100, "chr2", 0, 100))
    gextract("st", intervals = "bs")

    in_other_session({
        gtrack.rm("st", force = TRUE)
        gtrack.create_sparse("st", "", gintervals(c("chr1", "chr1", "chr1", "chr2"), c(0, 20, 40, 0), c(10, 30, 50, 10)), c(5, 6, 7, 8))
        gtrack.rm("dt", force = TRUE)
        gtrack.create("dt", "", "st", iterator = 7)
        gtrack.rm("rt", force = TRUE)
        gtrack.2d.create("rt", "", data.frame(
            chrom1 = c("chr1", "chr1", "chr1"), start1 = c(0, 20, 40), end1 = c(10, 30, 50),
            chrom2 = c("chr1", "chr2", "chr2"), start2 = c(0, 0, 20), end2 = c(10, 10, 30)
        ), c(3, 4, 5))
        gintervals.rm("bs", force = TRUE)
        withr::with_options(list(gmax.data.size = 1), gintervals.save("bs", gintervals(c("chr1", "chr1", "chr2"), c(0, 40, 0), c(10, 50, 10))))
        gintervals.convert_to_indexed("bs", remove.old = TRUE)
    })

    expect_equal(gextract("st", gintervals.all())$st, c(5, 6, 7, 8))
    expect_equal(gextract("dt", gintervals("chr1", 0, 14))$dt, c(5, 5))
    expect_equal(gextract("rt", gintervals.2d("chr1", 0, 100, "chr2", 0, 100))$rt, c(4, 5))
    expect_equal(gextract("st", intervals = "bs")$st, c(5, 7, 8))
})

test_that("a rescan does not clear the dirty flag over a track created during its scan", {
    skip_on_os("windows")
    local_db_state()
    db <- review_db(withr::local_tempdir())
    .gdb.cache_mark_dirty(db)
    real_gcall <- misha:::.gcall
    local_mocked_bindings(.gcall = function(...) {
        res <- real_gcall(...)
        # another session, loaded before (its lists of this database), creates a track right after
        # the scan (a forked child of this reload would start from its half-reloaded state)
        if (identical(..1, "gfind_tracks_n_intervals") && !dir.exists(file.path(db, "tracks", "late.track"))) {
            in_other_session({
                assign("GTRACKS", character(0), envir = .misha)
                assign("GINTERVS", character(0), envir = .misha)
                assign("GTRACKS_SRC", db, envir = .misha)
                gtrack.create_sparse("late", "", gintervals("chr1", 0, 10), 1)
            })
        }
        res
    })
    suppressMessages(gdb.reload())
    expect_true(dir.exists(file.path(db, "tracks", "late.track")))
    expect_true("late" %in% read_cache(db)[[1]] || .gdb.cache_is_dirty(db))
})

test_that("concurrent sessions leave every track in .db.cache or the cache dirty", {
    skip_on_os("windows")
    local_db_state()
    db <- review_db(withr::local_tempdir())
    jobs <- lapply(1:8, function(k) {
        parallel::mcparallel(withr::with_options(list(gmultitasking = FALSE), {
            for (i in 1:5) gtrack.create_sparse(sprintf("c%d_%d", k, i), "", gintervals("chr1", 0, 10), i)
            TRUE
        }))
    })
    expect_true(all(vapply(parallel::mccollect(jobs), isTRUE, logical(1))))
    tracks <- sprintf("c%d_%d", rep(1:8, each = 5), 1:5)
    expect_true(all(dir.exists(file.path(db, "tracks", paste0(tracks, ".track")))))
    expect_true(all(tracks %in% read_cache(db)[[1]]) || .gdb.cache_is_dirty(db))
    expect_false(file.exists(file.path(db, ".db.cache.lock")))
    suppressMessages(gsetroot(db))
    expect_true(all(tracks %in% gtrack.ls()))
})

test_that("a stale .db.cache lock is broken, and a held one makes a writer mark the cache dirty", {
    local_db_state()
    db <- review_db(withr::local_tempdir())
    lock <- file.path(db, ".db.cache.lock")

    # left by a killed writer a minute ago
    writeLines("n0 1 dead", lock)
    Sys.setFileTime(lock, Sys.time() - 60)
    gtrack.create_sparse("after_stale", "", gintervals("chr1", 0, 10), 1)
    expect_false(file.exists(lock))
    expect_true("after_stale" %in% read_cache(db)[[1]])
    expect_false(.gdb.cache_is_dirty(db))

    # held now by another writer
    writeLines("n0 2 alive", lock)
    gtrack.create_sparse("while_held", "", gintervals("chr1", 0, 10), 1)
    expect_equal(readLines(lock), "n0 2 alive")
    expect_true(.gdb.cache_is_dirty(db))
    unlink(lock)
})

# Replaces a base function for the duration of the calling test
mock_base <- function(name, fun, env = parent.frame()) {
    real <- get(name, envir = baseenv())
    unlockBinding(name, .BaseNamespaceEnv)
    assign(name, fun, envir = .BaseNamespaceEnv)
    lockBinding(name, .BaseNamespaceEnv)
    withr::defer(
        {
            unlockBinding(name, .BaseNamespaceEnv)
            assign(name, real, envir = .BaseNamespaceEnv)
            lockBinding(name, .BaseNamespaceEnv)
        },
        envir = env
    )
    real
}

test_that("a failed gtrack.copy(overwrite = TRUE) keeps what another session put at the name", {
    local_db_state()
    db <- review_db(withr::local_tempdir())
    tracks <- file.path(db, "tracks")
    gtrack.create_sparse("orig", "", gintervals("chr1", 0, 10), 42)
    gtrack.create_sparse("src", "", gintervals("chr1", 100, 110), 7)
    gtrack.create_sparse("other", "", gintervals("chr1", 200, 210), 99)
    other_bytes <- readBin(file.path(tracks, "other.track", "track.dat"), "raw", 1e6)

    # another session's copy lands at the name between this one's two renames
    rename <- mock_base("file.rename", function(from, to) {
        if (length(to) == 1 && basename(to) == "orig.track" && grepl("\\.tmp\\.", basename(from)) && !dir.exists(to)) {
            dir.create(to)
            file.copy(list.files(file.path(tracks, "other.track"), full.names = TRUE, all.files = TRUE, no.. = TRUE), to, recursive = TRUE)
        }
        rename(from, to)
    })
    orig_bytes <- readBin(file.path(tracks, "orig.track", "track.dat"), "raw", 1e6)
    expect_warning(expect_error(gtrack.copy("src", "orig", overwrite = TRUE), "Failed to rename"), "is kept at")
    expect_equal(readBin(file.path(tracks, "orig.track", "track.dat"), "raw", 1e6), other_bytes)
    # the track that was there is kept, under its hidden name: what took the name may be incomplete
    # (an older misha writes in place) or go again (a copy that fails after it was placed)
    kept <- list.files(tracks, pattern = "^\\.[0-9a-f]+$", all.files = TRUE, full.names = TRUE)
    expect_length(kept, 1)
    expect_equal(readBin(file.path(kept, "track.dat"), "raw", 1e6), orig_bytes)
})

test_that("gtrack.copy stops on a short copy instead of committing it", {
    local_db_state()
    db <- review_db(withr::local_tempdir())
    tracks <- file.path(db, "tracks")
    gtrack.create_sparse("orig", "", gintervals("chr1", 0, 10), 42)
    gtrack.create_sparse("src", "", gintervals(c("chr1", "chr2"), 100, 110), c(7, 8))
    before <- readBin(file.path(tracks, "orig.track", "track.dat"), "raw", 1e6)

    # file.copy() succeeds, but the data reached the disk short (as a close error it does not report)
    copy <- mock_base("file.copy", function(from, to, ...) {
        res <- copy(from, to, ...)
        if (basename(from[1]) == "track.dat") {
            dest <- file.path(to, "track.dat")
            writeBin(readBin(dest, "raw", 8), dest)
        }
        res
    })
    expect_error(gtrack.copy("src", "orig", overwrite = TRUE), "incomplete")
    expect_equal(readBin(file.path(tracks, "orig.track", "track.dat"), "raw", 1e6), before)
})

test_that("an array track another session re-created with another layout is read right", {
    skip_on_os("windows")
    local_db_state()
    td <- withr::local_tempdir()
    db <- review_db(td)
    write_arr <- function(f, iv) write.table(cbind(iv, a = seq_len(nrow(iv)), b = 10 * seq_len(nrow(iv))), f, sep = "\t", quote = FALSE, row.names = FALSE)
    f1 <- file.path(td, "a1.tsv")
    write_arr(f1, gintervals(c("chr1", "chr2"), 0, 10))
    # (gtrack.array.import writes per-chromosome files even in an indexed database)
    suppressMessages(gtrack.array.import("arr", "", f1))
    suppressMessages(gtrack.convert_to_indexed("arr"))
    expect_true(file.exists(file.path(db, "tracks", "arr.track", "track.idx")))
    gtrack.array.extract("arr", c("a", "b"), gintervals.all())

    f2 <- file.path(td, "a2.tsv")
    write_arr(f2, gintervals(c("chr1", "chr1", "chr1", "chr2"), c(0, 20, 40, 0), c(10, 30, 50, 10)))
    in_other_session({
        gtrack.rm("arr", force = TRUE)
        suppressMessages(gtrack.array.import("arr", "", f2))
        suppressMessages(gtrack.convert_to_indexed("arr"))
    })
    res <- gtrack.array.extract("arr", c("a", "b"), gintervals.all())
    expect_equal(res$a, 1:4)
    expect_equal(res$b, 10 * (1:4))
})

test_that("an older rescan does not overwrite a newer cache", {
    skip_on_os("windows")
    local_db_state()
    db <- review_db(withr::local_tempdir())
    real_gcall <- misha:::.gcall
    other <- FALSE
    local_mocked_bindings(.gcall = function(...) {
        res <- real_gcall(...)
        # while this rescan is slow, another session rescans, writes, and creates a track (set before
        # the fork: the other session's own rescan does not start one more)
        if (identical(..1, "gfind_tracks_n_intervals") && !other) {
            other <<- TRUE
            in_other_session({
                suppressMessages(gdb.reload())
                gtrack.create_sparse("newtrack", "", gintervals("chr1", 0, 10), 1)
            })
        }
        res
    })
    suppressMessages(gdb.reload())
    expect_true("newtrack" %in% read_cache(db)[[1]] || .gdb.cache_is_dirty(db))
})

test_that("a mark made by another session while a writer merges is not lost", {
    skip_on_os("windows")
    local_db_state()
    db <- review_db(withr::local_tempdir())
    expect_false(.gdb.cache_is_dirty(db))
    real_mark <- misha:::.gdb.cache_mark_dirty
    other <- FALSE
    local_mocked_bindings(.gdb.cache_mark_dirty = function(groot = NULL, ...) {
        # between this writer's look at the flag and its own mark, another session marks the cache
        if (!other) {
            other <<- TRUE
            in_other_session(real_mark(groot))
        }
        real_mark(groot, ...)
    })
    gtrack.create_sparse("t6", "", gintervals("chr1", 0, 10), 1)
    expect_true(other)
    expect_true(.gdb.cache_is_dirty(db))
})

test_that("an unbreakable stale lock or a root that cannot be written does not make a write hang", {
    skip_on_os("windows")
    skip_if(Sys.info()[["user"]] == "root")
    local_db_state()
    db <- review_db(withr::local_tempdir())
    lock <- file.path(db, ".db.cache.lock")
    writeLines("n0 1 dead", lock)
    Sys.setFileTime(lock, Sys.time() - 60)
    Sys.chmod(lock, "000")
    withr::defer(Sys.chmod(lock, "644"))
    # in a forked child, given 60 s (it hung, busy-waiting)
    job <- parallel::mcparallel(withr::with_options(list(gmultitasking = FALSE), {
        gtrack.create_sparse("t7", "", gintervals("chr1", 0, 10), 1)
        TRUE
    }))
    res <- parallel::mccollect(job, wait = FALSE, timeout = 60)
    if (is.null(res)) {
        tools::pskill(job$pid, tools::SIGKILL)
        parallel::mccollect(job)
    }
    expect_true(isTRUE(res[[1]]))
    expect_true(.gdb.cache_is_dirty(db))

    # a root that cannot be written (tracks/ can): no lock, no flag, no warning
    Sys.chmod(lock, "644")
    unlink(lock)
    unlink(list.files(db, pattern = "^\\.db\\.cache\\.dirty", all.files = TRUE, full.names = TRUE))
    Sys.chmod(db, "555")
    withr::defer(Sys.chmod(db, "755"))
    # (the cache there stays stale: said)
    expect_warning(gtrack.create_sparse("t8", "", gintervals("chr1", 0, 10), 1), "Failed to mark cache dirty")
    expect_true(dir.exists(file.path(db, "tracks", "t8.track")))
})

test_that("a rescan interrupted while it holds the lock does not leave it", {
    local_db_state()
    db <- review_db(withr::local_tempdir())
    real_mark <- misha:::.gdb.cache_mark_dirty
    local_mocked_bindings(.gdb.cache_mark_dirty = function(...) stop("interrupted here"))
    expect_error(suppressMessages(gdb.reload()), "interrupted here")
    expect_false(file.exists(file.path(db, ".db.cache.lock")))
})

test_that("a mark also makes the flag earlier versions read, which a rescan that saw it removes", {
    local_db_state()
    db <- review_db(withr::local_tempdir())
    legacy <- file.path(db, ".db.cache.dirty")
    flags <- function() list.files(db, pattern = "^\\.db\\.cache\\.dirty", all.files = TRUE)
    expect_length(flags(), 0)

    token <- .gdb.cache_mark_dirty(db)
    expect_setequal(flags(), c(".db.cache.dirty", paste0(".db.cache.dirty.", token)))

    # a rescan that saw both at its start removes both
    suppressMessages(gdb.reload())
    expect_length(flags(), 0)

    # a write that merges its name leaves neither (its own mark while it works is not a lasting one)
    gtrack.create_sparse("t9", "", gintervals("chr1", 0, 10), 1)
    expect_length(flags(), 0)
    expect_true("t9" %in% read_cache(db)[[1]])
})

test_that("interval operations and cartesian iterators over the deferred whole-genome 2D scope build it or stop", {
    local_db_state()
    withr::local_options(list(gmulticontig.2d.threshold = 2, gmulticontig.2d.mode = NULL))
    review_db(withr::local_tempdir())
    r <- data.frame(
        chrom1 = c("chr1", "chr1", "chr2", "chrX"), start1 = c(0, 0, 0, 0), end1 = c(100, 100, 100, 100),
        chrom2 = c("chr1", "chr2", "chrX", "chrM"), start2 = c(200, 0, 0, 0), end2 = c(300, 100, 100, 100)
    )
    gtrack.2d.create("mr", "", r, c(10, 20, 30, 40))
    msg <- 'options(gmulticontig.2d.mode = "full")'
    cg <- giterator.cartesian_grid(gintervals(c("chr1", "chr2"), 0, 100), c(-50, 50))

    # it returned no 2D intervals
    expect_error(gintervals.2d.band_intersect(.misha$ALLGENOME, c(-1e9, 1e9)), msg, fixed = TRUE)
    expect_error(gintervals.canonic(.misha$ALLGENOME), msg, fixed = TRUE)
    expect_error(gextract("mr", .misha$ALLGENOME, iterator = cg), msg, fixed = TRUE)
    expect_error(gtrack.create("cg_unset", "", "mr", iterator = cg), msg, fixed = TRUE)

    withr::with_options(list(gmulticontig.2d.mode = "full"), {
        # the 4 cis pairs within the band
        expect_equal(nrow(gintervals.2d.band_intersect(.misha$ALLGENOME, c(-1e9, 1e9))), 4)
        expect_equal(nrow(gintervals.canonic(.misha$ALLGENOME)), 16)
        # the cartesian grid of chr1 0-100 and chr2 0-100, with chr1-chr2 0-100 holding 20
        res <- gextract("mr", .misha$ALLGENOME, iterator = cg)
        expect_equal(res$mr[as.character(res$chrom1) == "chr1" & as.character(res$chrom2) == "chr2"], 20)
        gtrack.create("cg_full", "", "mr", iterator = cg)
        expect_equal(gextract("cg_full", gintervals.2d("chr1", 0, 100, "chr2", 0, 100))$cg_full, 20)
        # from the global environment, as a user calls it (the build is evaluated in misha's namespace)
        expect_equal(evalq(misha::gsummary("mr")[["Sum"]], globalenv()), 100)
    })
})

test_that("a rescan keeps the legacy flag that a lasting mark rewrote during its scan", {
    skip_on_os("windows")
    local_db_state()
    db <- review_db(withr::local_tempdir())
    legacy <- file.path(db, ".db.cache.dirty")
    .gdb.cache_mark_dirty(db)
    real_gcall <- misha:::.gcall
    other <- FALSE
    local_mocked_bindings(.gcall = function(...) {
        res <- real_gcall(...)
        # during the scan, another session marks the cache (a track it made is not in the scan)
        if (identical(..1, "gfind_tracks_n_intervals") && !other) {
            other <<- TRUE
            in_other_session(misha:::.gdb.cache_mark_dirty(db))
        }
        res
    })
    suppressMessages(gdb.reload())
    expect_true(other)
    expect_true(file.exists(legacy))
    expect_true(.gdb.cache_is_dirty(db))
})

test_that("a writer whose cache write fails leaves a lasting mark, the legacy flag included", {
    local_db_state()
    db <- review_db(withr::local_tempdir())
    local_mocked_bindings(.gdb.cache_write_file = function(...) stop("disk full"))
    expect_warning(gtrack.create_sparse("t10", "", gintervals("chr1", 0, 10), 1), "disk full")
    expect_true(file.exists(file.path(db, ".db.cache.dirty")))
    expect_true(dir.exists(file.path(db, "tracks", "t10.track")))
})
