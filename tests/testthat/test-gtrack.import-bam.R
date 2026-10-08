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

    intervs <- gintervals("chr1", c(90, 100, 190, 200), c(100, 110, 200, 210))
    r <- gextract(tmptrack, intervs, iterator = 10, colnames = "v")
    # POS 100 / 200 are 1-based: the forward read covers [99, 109), the reverse
    # read (5' end at 208) covers [199, 209).
    expect_equal(as.numeric(r$v), c(0.1, 0.9, 0.1, 0.9), tolerance = 1e-6)
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

# Runs `cmd` in the background through the shell and returns a function that kills it.
# (Not system(intern = TRUE): a writer waiting on the pipe would hold R's output pipe open.)
start_background <- function(cmd) {
    pidfile <- tempfile()
    system(paste("(", cmd, ") > /dev/null 2>&1 & echo $! >", shQuote(pidfile)))
    pid <- as.integer(readLines(pidfile))
    unlink(pidfile)
    function() tools::pskill(pid, tools::SIGKILL)
}

# After `after` seconds, keeps opening the pipe for writing, so that an import stuck opening or
# reading it gets end-of-file and fails its expectations instead of hanging the test run.
# (No forked R child: misha reaps any child process, which confuses the parallel package.)
start_fifo_watchdog <- function(fifo, after = 30) {
    # stops by itself once the pipe is removed
    start_background(sprintf("sleep %d; while [ -p %s ]; do : > %s; sleep 0.2; done", after, shQuote(fifo), shQuote(fifo)))
}

# Starts `cat src > fifo` in the background; it blocks until a reader opens the pipe.
start_fifo_writer <- function(fifo, src) {
    start_background(paste("cat", shQuote(src), ">", shQuote(fifo)))
}

test_that("gtrack.import_mappedseq reads SAM from a named pipe", {
    skip_on_os("windows")
    skip_if(!nzchar(Sys.which("mkfifo")), "mkfifo not on PATH")

    sam <- tempfile(fileext = ".sam")
    writeLines(default_sam_text(), sam)
    sam_gz <- tempfile(fileext = ".sam.gz")
    con <- gzfile(sam_gz, "w")
    writeLines(default_sam_text(), con)
    close(con)

    for (gz in c(FALSE, TRUE)) {
        fifo <- tempfile()
        system2("mkfifo", fifo)
        stop_writer <- start_fifo_writer(fifo, if (gz) sam_gz else sam)
        stop_watchdog <- start_fifo_watchdog(fifo)
        withr::defer({
            stop_writer()
            stop_watchdog()
            unlink(fifo)
        })

        tmptrack <- random_track_name("test")
        withr::defer(gtrack.rm(tmptrack, force = TRUE))
        stats <- gtrack.import_mappedseq(tmptrack, "FIFO", fifo, cols.order = NULL)
        expect_equal(stats[[1]][["total.mapped"]], 2, info = paste("gzipped:", gz))
        expect_equal(stats[[1]][["total.unmapped"]], 1, info = paste("gzipped:", gz))
        expect_true(gtrack.exists(tmptrack))
    }
})

test_that("gtrack.import_mappedseq rejects BAM sent through a named pipe", {
    skip_on_os("windows")
    skip_if(!nzchar(Sys.which("mkfifo")), "mkfifo not on PATH")
    bam <- default_bam_path()
    fifo <- tempfile()
    system2("mkfifo", fifo)
    stop_writer <- start_fifo_writer(fifo, bam)
    stop_watchdog <- start_fifo_watchdog(fifo)
    withr::defer({
        stop_writer()
        stop_watchdog()
        unlink(fifo)
    })
    tmptrack <- random_track_name("test")
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    expect_error(gtrack.import_mappedseq(tmptrack, "FIFO", fifo, cols.order = NULL), "pipe carrying BAM")
})
