# A track file that is not written completely must fail the call and leave no track. Writers
# used to close their files without looking at what fclose() returned, and the last buffer of a
# file is written there: a failure at that point (a full disk, a quota, an NFS write error) left
# a truncated track behind a call that reported success.
#
# RLIMIT_FSIZE makes a write fail on demand. Every track file here is smaller than the stdio
# buffer, so none of it reaches the file before fclose().

ensure_valid_groot()

with_file_size_limit <- function(bytes, code) {
    .Call("C_test_file_size_limit", as.numeric(bytes))
    on.exit(.Call("C_test_file_size_limit", -1))
    code
}

setup_write_db <- function(sizes) {
    create_test_db("write_db", chrom_sizes = data.frame(chrom = names(sizes), size = unname(sizes)))
    gdb.init("write_db")
}

# 400 bins of 1 bp: a 1604-byte dense file, over the 1000-byte limit; 100 bins: 404 bytes, under it
expect_write_failure <- function(sizes, create) {
    withr::with_tempdir({
        setup_write_db(sizes)
        expect_error(with_file_size_limit(1000, create()), "Failed to write .*/chr[0-9]+: ")
        expect_false(gtrack.exists("t"))
        create() # and the same call succeeds without the limit
        expect_true(gtrack.exists("t"))
    })
}

test_that("a failed write of the last chromosome's file fails gtrack.create", {
    local_db_state()
    withr::local_options(gmultitasking = FALSE)
    expect_write_failure(c(chr1 = 100, chr2 = 400), function() gtrack.create("t", "", "1", iterator = 1))
})

test_that("a failed write of a file closed for the next chromosome fails gtrack.create", {
    local_db_state()
    withr::local_options(gmultitasking = FALSE)
    expect_write_failure(c(chr1 = 400, chr2 = 100), function() gtrack.create("t", "", "1", iterator = 1))
})

test_that("a failed write fails gtrack.create_dense and gtrack.import", {
    local_db_state()
    expect_write_failure(c(chr1 = 400), function() {
        gtrack.create_dense("t", "", gintervals("chr1", 0, 400), 1, binsize = 1, defval = 0)
    })
    bg <- tempfile(fileext = ".bedGraph")
    writeLines("chr1\t0\t400\t1", bg)
    expect_write_failure(c(chr1 = 400), function() gtrack.import("t", "", bg, binsize = 1, defval = 0))
})

test_that("a failed write in a multitasking child fails the call", {
    local_db_state()
    withr::local_options(gmultitasking = TRUE, gmax.processes = 2, gmin.scope4process = 1)
    expect_write_failure(c(chr1 = 400, chr2 = 400), function() gtrack.create("t", "", "1", iterator = 1))
})

test_that("a call nested in a writer neither raises nor drops the writer's write error", {
    local_db_state()
    withr::local_options(gmultitasking = FALSE, gbuf.size = 10)
    withr::with_tempdir({
        setup_write_db(c(chr1 = 400, chr2 = 100))
        gtrack.create("src", "", "1", iterator = 1)
        # evaluated while gtrack.create holds its files; the inner call fails and is caught
        f <- function(x) {
            try(gtrack.info("no_such_track"), silent = TRUE)
            x
        }
        assign("f", f, envir = .GlobalEnv)
        withr::defer(rm("f", envir = .GlobalEnv))
        expect_error(with_file_size_limit(1000, gtrack.create("t", "", "f(src)", iterator = 1)), "Failed to write .*/chr1: ")
        expect_false(gtrack.exists("t"))
    })
})

test_that("gtrack.modify does not commit a chromosome whose staged file failed to write", {
    local_db_state()
    withr::local_options(gmultitasking = FALSE, gbuf.size = 10)
    withr::with_tempdir({
        setup_write_db(c(chr1 = 400, chr2 = 100))
        gtrack.create("t", "", "0", iterator = 1)
        # the cap comes after chr1 has been staged (copied in full), so only its new values fail
        calls <- 0
        g <- function(x) {
            calls <<- calls + 1
            if (calls == 5) .Call("C_test_file_size_limit", 1000)
            x + 1
        }
        assign("g", g, envir = .GlobalEnv)
        withr::defer(rm("g", envir = .GlobalEnv))
        withr::defer(.Call("C_test_file_size_limit", -1))
        expect_error(gtrack.modify("t", "g(t)", gintervals.all()), "Failed to write .*/chr1: ")
        .Call("C_test_file_size_limit", -1)
        expect_equal(unique(gextract("t", gintervals.all(), iterator = 1)$t), 0)
    })
})

test_that("a failed write fails a big intervals set and gintervals.save", {
    local_db_state()
    withr::local_options(gmultitasking = FALSE, gbig.intervals.size = 10)
    withr::with_tempdir({
        setup_write_db(c(chr1 = 240, chr2 = 240))
        starts <- seq(0, 238, by = 2)
        gtrack.create_sparse("sp", "", gintervals(rep(c("chr1", "chr2"), each = 120), starts, starts + 1), rep(1, 240))
        expect_error(with_file_size_limit(1000, gscreen("sp > 0", intervals.set.out = "out")), "Failed to write")
        expect_false(gintervals.exists("out"))

        expect_error(with_file_size_limit(1000, gintervals.save("s", gintervals("chr1", starts, starts + 1))), "Failed to write")
        expect_false(gintervals.exists("s"))
        gintervals.save("s", gintervals("chr1", starts, starts + 1))
        expect_equal(nrow(gintervals.load("s")), 120)

        # a set too large for the stdio buffer fails in serialize() rather than in close()
        many <- gintervals("chr1", 0:239, 1:240)
        expect_error(with_file_size_limit(1000, gintervals.save("s2", many)), "Failed to write")
        expect_false(file.exists(file.path(.misha$GWD, "s2.interv")))
    })
})

test_that("a failed write fails gtrack.array.import and leaves no track", {
    local_db_state()
    withr::with_tempdir({
        setup_write_db(c(chr1 = 1000))
        af <- file.path(getwd(), "arr.txt")
        cols <- sprintf("sample_column_%03d", 1:100) # .colnames is over 2000 bytes, the track data under it
        writeLines(c(
            paste(c("chrom", "start", "end", cols), collapse = "\t"),
            paste(c("chr1", 0, 10, rep(1, 100)), collapse = "\t")
        ), af)
        expect_error(with_file_size_limit(2000, gtrack.array.import("arr", "", af)), "Failed to write .*\\.colnames: ")
        expect_false(gtrack.exists("arr"))
        gtrack.array.import("arr", "", af)
        expect_equal(gtrack.array.get_colnames("arr"), cols)
    })
})

test_that("a failed gintervals.update drops the set it could not write", {
    local_db_state()
    withr::with_tempdir({
        setup_write_db(c(chr1 = 1000, chr2 = 1000))
        gintervals.save("sm", gintervals(c("chr1", "chr2"), 0, 10))
        expect_error(with_file_size_limit(500, gintervals.update("sm", gintervals("chr1", 0:299, 1:300), chrom = "chr1")), "Failed to write")
        expect_false(gintervals.exists("sm"))
    })
})

test_that("a track variable that fails to write midway is not left behind", {
    local_db_state()
    withr::with_tempdir({
        setup_write_db(c(chr1 = 100))
        gtrack.create("t", "", "1", iterator = 1)
        # 16 KB: larger than the stdio buffer, so serialize() itself fails, not close()
        expect_error(with_file_size_limit(1000, gtrack.var.set("t", "v", runif(2000))), "Failed to write .*/vars/v: ")
        expect_false("v" %in% gtrack.var.ls("t"))
        gtrack.var.set("t", "v", 1:3)
        expect_equal(gtrack.var.get("t", "v"), 1:3)
    })
})
