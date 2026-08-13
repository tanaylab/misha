read_db_cache <- function(db) {
    p <- file.path(db, ".db.cache")
    if (!file.exists(p)) {
        return(NULL)
    }
    f <- file(p, "rb")
    on.exit(close(f))
    unserialize(f)
}

# Build two independent single-chromosome databases, each with one track.
make_two_dbs <- function(env = parent.frame()) {
    tmp <- tempfile("cacheprov_")
    dir.create(tmp)
    withr::defer(unlink(tmp, recursive = TRUE), envir = env)

    dbs <- lapply(c("A", "B"), function(tag) {
        db <- file.path(tmp, paste0("db", tag))
        fa <- file.path(tmp, paste0(tag, ".fasta"))
        cat(sprintf(">chr%s\nACGTACGTAC\n", tag), file = fa)
        suppressMessages(gdb.create(groot = db, fasta = fa, verbose = FALSE))
        suppressMessages(gsetroot(db))
        gtrack.create(paste0("track", tag), "t", "1",
            iterator = gintervals(paste0("chr", tag), 0, 10)
        )
        db
    })
    names(dbs) <- c("A", "B")
    dbs
}

test_that("track listing of one database is never persisted to another's .db.cache", {
    local_db_state()
    dbs <- make_two_dbs()

    expect_equal(read_db_cache(dbs$A)[[1]], "trackA")
    expect_equal(read_db_cache(dbs$B)[[1]], "trackB")

    # Reproduces the state left behind when gsetroot() is interrupted mid-reload,
    # or when a caller restores a memoized session snapshot (misha.ext::gset_genome
    # with force = FALSE): GROOT/GWD point at one database while GTRACKS still
    # holds the listing of another.
    suppressMessages(gsetroot(dbs$A))
    stale_tracks <- get("GTRACKS", envir = .misha)
    suppressMessages(gsetroot(dbs$B))
    assign("GTRACKS", stale_tracks, envir = .misha)
    assign("GTRACKS_SRC", dbs$A, envir = .misha)

    gtrack.create("newB", "t", "1", iterator = gintervals("chrB", 0, 10))

    # dbB's cache must not have been rewritten from dbA's listing
    expect_false("trackA" %in% read_db_cache(dbs$B)[[1]])
    expect_true(.gdb.cache_is_dirty(dbs$B))

    # ... and the next session rescans from disk and sees the truth
    suppressMessages(gsetroot(dbs$B))
    expect_setequal(gtrack.ls(), c("trackB", "newB"))
    expect_setequal(read_db_cache(dbs$B)[[1]], c("trackB", "newB"))
})

test_that("interrupting gsetroot() unloads the session instead of reporting success", {
    local_db_state()
    dbs <- make_two_dbs()

    suppressMessages(gsetroot(dbs$A))

    real_reload <- misha::gdb.reload
    assignInNamespace("gdb.reload", function(rescan = TRUE) {
        stop(structure(class = c("interrupt", "condition"), list(message = "", call = NULL)))
    }, ns = "misha")
    withr::defer(assignInNamespace("gdb.reload", real_reload, ns = "misha"))

    # Catch it here rather than with expect_condition(): testthat's runner
    # treats a real `interrupt` condition as Ctrl-C and aborts the file.
    caught <- FALSE
    tryCatch(suppressMessages(gsetroot(dbs$B)), interrupt = function(e) caught <<- TRUE)
    expect_true(caught)
    expect_null(get("GROOT", envir = .misha))

    assignInNamespace("gdb.reload", real_reload, ns = "misha")
    suppressMessages(gsetroot(dbs$B))
    expect_equal(gtrack.ls(), "trackB")
})

test_that("cache is still written incrementally in the normal single-database flow", {
    local_db_state()
    dbs <- make_two_dbs()

    suppressMessages(gsetroot(dbs$A))
    gtrack.create("extraA", "t", "1", iterator = gintervals("chrA", 0, 10))

    expect_setequal(read_db_cache(dbs$A)[[1]], c("trackA", "extraA"))
    expect_false(.gdb.cache_is_dirty(dbs$A))

    gintervals.save("ivA", gintervals("chrA", 0, 10))
    expect_setequal(read_db_cache(dbs$A)[[2]], "ivA")
    expect_false(.gdb.cache_is_dirty(dbs$A))
})
