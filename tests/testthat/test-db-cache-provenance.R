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

test_that("creating a track inserts its name into the track lists without sorting them again", {
    local_db_state()
    dbs <- make_two_dbs()
    suppressMessages(gsetroot(dbs$A))
    gdir.create("a1")
    for (nm in c("b1", "c_1", "a1.x", "zz", "B2")) {
        gtrack.create(nm, "t", "1", iterator = gintervals("chrA", 0, 10))
    }
    # a database can list ~10^5 tracks: sort() over them cost 0.36 s per list and track
    sorts_of_new <- 0L
    record <- function(x) if (is.character(x) && "m_new" %in% x) sorts_of_new <<- sorts_of_new + 1L
    suppressMessages(trace("sort.int", tracer = bquote(.(record)(x)), where = baseenv(), print = FALSE))
    withr::defer(suppressMessages(untrace("sort.int", where = baseenv())))
    gtrack.create("m_new", "t", "1", iterator = gintervals("chrA", 0, 10))
    suppressMessages(untrace("sort.int", where = baseenv()))

    expect_equal(sorts_of_new, 0L)
    want <- sort(c("trackA", "b1", "c_1", "a1.x", "zz", "B2", "m_new"))
    expect_identical(get("GTRACKS", envir = misha:::.misha), want)
    expect_identical(read_db_cache(dbs$A)[[1]], want)
})

test_that(".gdb.sorted_insert() is sort(unique(c(x, name))), or NULL", {
    set.seed(5)
    chars <- c(letters, LETTERS, 0:9, ".", "_")
    rand_name <- function() paste(sample(chars, sample(1:8, 1), replace = TRUE), collapse = "")
    for (i in 1:200) {
        x <- sort(unique(vapply(seq_len(sample(0:40, 1)), function(j) rand_name(), "")))
        name <- if (length(x) && i %% 4 == 0) sample(x, 1) else rand_name()
        got <- misha:::.gdb.sorted_insert(x, name)
        if (!is.null(got)) {
            expect_identical(got, sort(unique(c(x, name))))
        }
        expect_identical(misha:::.gdb.cache_list_add(x, name), misha:::.gdb.normalize_cache_list(c(x, name)))
    }
    expect_identical(misha:::.gdb.sorted_insert(character(0), "a"), "a")
    expect_identical(misha:::.gdb.sorted_insert(c("a", "c"), "b"), c("a", "b", "c"))
    expect_identical(misha:::.gdb.sorted_insert(c("a", "c"), "c"), c("a", "c"))
    # not sorted, NA: left to sort()
    expect_null(misha:::.gdb.sorted_insert(c("c", "a"), "b"))
    expect_null(misha:::.gdb.sorted_insert(c("a", NA), "b"))
    expect_identical(misha:::.gdb.cache_list_add(c("c", "a"), "b"), c("a", "b", "c"))
})
