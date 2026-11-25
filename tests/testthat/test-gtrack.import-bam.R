# End-to-end tests of BAM input in gtrack.import_mappedseq. misha reads BAM itself; none of these
# needs samtools. See helper-bam.R for the fixtures.

local_bam_db()

# Imports file into a new track; returns its statistics, warnings and the bytes of its files
import_bam_track <- function(file, ...) {
    track <- random_track_name("test")
    gtrack.rm(track, force = TRUE)
    withr::defer(gtrack.rm(track, force = TRUE))
    warns <- character(0)
    stats <- withCallingHandlers(
        gtrack.import_mappedseq(track, "BAM test", file, ...),
        warning = function(w) {
            warns <<- c(warns, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    dir <- file.path(.misha$GROOT, "tracks", paste0(gsub(".", "/", track, fixed = TRUE), ".track"))
    files <- sort(setdiff(list.files(dir, all.files = TRUE, no.. = TRUE), ".attributes"))
    bytes <- lapply(stats::setNames(files, files), function(f) readBin(file.path(dir, f), "raw", file.size(file.path(dir, f))))
    list(stats = stats, warnings = warns, files = bytes)
}

bam_option_sets <- list(
    sparse = list(),
    sparse_all_dups = list(remove.dups = FALSE),
    dense = list(pileup = 100, binsize = 20),
    dense_all_dups = list(pileup = 100, binsize = 20, remove.dups = FALSE),
    min_mapq = list(min.mapq = 30),
    min_mapq_dense = list(min.mapq = 30, pileup = 200, binsize = 50),
    paired = list(paired = TRUE, binsize = 50),
    paired_all_dups = list(paired = TRUE, binsize = 50, remove.dups = FALSE),
    paired_min_mapq = list(paired = TRUE, binsize = 10, min.mapq = 30),
    paired_max_fraglen = list(paired = TRUE, binsize = 50, max.fraglen = 300),
    paired_no_max_fraglen = list(paired = TRUE, binsize = 50, max.fraglen = Inf)
)

test_that("a BAM gives the track and statistics of its samtools view text", {
    # split.bam holds single.bam's records in small BGZF blocks (records and header span blocks)
    fixtures <- c(default = "default", single = "single", paired = "paired", split = "single", edge = "edge", empty = "empty")
    for (bam in names(fixtures)) {
        sam <- bam_fixture(paste0(fixtures[[bam]], ".view.sam.gz"))
        for (opt in names(bam_option_sets)) {
            args <- bam_option_sets[[opt]]
            from_bam <- do.call(import_bam_track, c(list(bam_fixture(paste0(bam, ".bam"))), args))
            from_sam <- do.call(import_bam_track, c(list(sam, cols.order = NULL), args))
            expect_identical(from_bam, from_sam, info = paste(bam, opt))
        }
    }
})

test_that("BAM import statistics", {
    s <- import_bam_track(bam_fixture("single.bam"))$stats
    expect_equal(s[[1]], c(total = 2555, total.mapped = 1323, total.unmapped = 808, total.dups = 4, total.filtered = 424))
    expect_equal(s[[2]]$mapped, c(468, 408, 447))

    s <- import_bam_track(bam_fixture("paired.bam"), paired = TRUE, binsize = 50)$stats
    expect_equal(s[[1]], c(total = 1029, total.mapped = 466, total.unmapped = 255, total.dups = 1, total.filtered = 308))

    # a mapped record without a CIGAR or with an operation code SAM has no letter for, a record
    # placed without a position and a CIGAR spanning no reference base are unmapped
    s <- import_bam_track(bam_fixture("edge.bam"))$stats
    expect_equal(s[[1]], c(total = 13, total.mapped = 5, total.unmapped = 8, total.dups = 0, total.filtered = 0))

    res <- import_bam_track(bam_fixture("empty.bam"))
    expect_equal(res$stats[[1]][["total"]], 0)
    expect_match(res$warnings, "No read was imported")
})

test_that("the threads that inflate a BAM do not change the result", {
    withr::local_options(gmax.processes = 1)
    one <- import_bam_track(bam_fixture("split.bam"), pileup = 100, binsize = 20)
    withr::local_options(gmax.processes = 8)
    expect_identical(import_bam_track(bam_fixture("split.bam"), pileup = 100, binsize = 20), one)
})

test_that("gtrack.import_mappedseq imports BAM as sparse track", {
    stats <- import_bam_track(bam_fixture("default.bam"), cols.order = NULL, remove.dups = TRUE)$stats
    # 2 mapped records on chrom "chr1", 1 unmapped record (chrom "*" in the BAM).
    expect_equal(stats[[1]][["total.mapped"]], 2)
    expect_equal(stats[[1]][["total.unmapped"]], 1)
})

test_that("gtrack.import_mappedseq imports BAM as dense pileup", {
    tmptrack <- random_track_name("test")
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.rm(tmptrack, force = TRUE)

    stats <- gtrack.import_mappedseq(
        tmptrack, "BAM dense", bam_fixture("default.bam"),
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
    # Caller does NOT pass cols.order. The legacy default c(9, 11, 13, 14)
    # is the tab-delimited layout; BAM records are read as SAM (columns 10/3/4/2).
    # The wrapper must treat the default as "user didn't pass anything" and
    # switch to NULL; otherwise the import returns 0 mapped reads.
    stats <- import_bam_track(bam_fixture("default.bam"), remove.dups = TRUE)$stats
    expect_equal(stats[[1]][["total.mapped"]], 2)
})

test_that("gtrack.import_mappedseq errors on explicit cols.order with BAM", {
    tmptrack <- random_track_name("test")
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    expect_error(
        gtrack.import_mappedseq(
            tmptrack, "BAM explicit cols.order", bam_fixture("default.bam"),
            cols.order = c(9, 11, 13, 14), remove.dups = TRUE
        ),
        "BAM input forces SAM column layout"
    )
})

test_that("BAM import does not need samtools", {
    stats <- withr::with_envvar(
        c(PATH = "/nonexistent"),
        import_bam_track(bam_fixture("default.bam"))$stats
    )
    expect_equal(stats[[1]][["total.mapped"]], 2)
})

test_that("a truncated or corrupt BAM is an error, and leaves no track", {
    bam <- readBin(bam_fixture("single.bam"), "raw", file.size(bam_fixture("single.bam")))
    corrupt <- function(bytes) {
        f <- tempfile(fileext = ".bam")
        writeBin(bytes, f)
        f
    }
    # the second BGZF block: its size less 1 is in bytes 17-18 of each block header
    block_size <- function(at) as.integer(bam[at + 16]) + 256L * as.integer(bam[at + 17]) + 1L
    second <- 1L + block_size(1L)
    bad_crc <- bam
    crc_at <- second + block_size(second) - 8L
    bad_crc[crc_at] <- xor(bad_crc[crc_at], as.raw(0xff))

    cases <- list(
        list(corrupt(bam[seq_len(length(bam) %/% 2)]), "does not end with the BGZF end-of-file block"),
        list(corrupt(bam[seq_len(length(bam) - 28)]), "does not end with the BGZF end-of-file block"),
        list(corrupt(c(bam[seq_len(length(bam) %/% 2)], utils::tail(bam, 28))), "last BGZF block is incomplete"),
        list(corrupt(bad_crc), "CRC32 mismatch"),
        list(bam_fixture("bad-truncated-record.bam"), "last record is incomplete"),
        list(bam_fixture("bad-cigar-length.bam"), "CIGAR and query sequence lengths differ for bad"),
        list(bam_fixture("bad-refid.bam"), "refers to a reference that its header does not list"),
        list(bam_fixture("bad-n-ref.bam"), "truncated or corrupt: its header is incomplete")
    )
    for (case in cases) {
        track <- random_track_name("test")
        withr::defer(gtrack.rm(track, force = TRUE))
        expect_error(gtrack.import_mappedseq(track, "corrupt", case[[1]]), case[[2]])
        expect_false(gtrack.exists(track))
        # the error leaves the lifetime counters as it found them (C6 in src/rdbutils.h)
        expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    }
})

test_that("a BAM header with more references than the file holds is corrupt, not out of memory", {
    # bad-n-ref.bam says it lists 2^31 - 1 references. Reserving room for that many names asks for
    # some 80 GB, which fails as "Out of memory" on a machine with less, and succeeds through overcommit
    # on one with more: the address space is capped at 4 GB beyond what this process uses, so the
    # test sees the small machine's outcome on any machine.
    skip_if_not(file.exists("/proc/self/status"), "needs /proc to size the address space")
    status <- readLines("/proc/self/status")
    vm_kb <- as.numeric(sub("^VmSize:\\s*([0-9]+) kB$", "\\1", grep("^VmSize:", status, value = TRUE)))
    withr::defer(.Call("C_test_address_space_limit", -1))
    skip_if_not(.Call("C_test_address_space_limit", vm_kb * 1024 + 4 * 2^30), "cannot cap the address space")
    track <- random_track_name("test")
    withr::defer(gtrack.rm(track, force = TRUE))
    msg <- tryCatch(gtrack.import_mappedseq(track, "corrupt", bam_fixture("bad-n-ref.bam")), error = function(e) conditionMessage(e))
    .Call("C_test_address_space_limit", -1)
    expect_match(msg, "truncated or corrupt: its header is incomplete")
    expect_false(gtrack.exists(track))
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
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
        # local(): each iteration's deferred cleanup must see its own writer, watchdog, pipe and track
        local({
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
        })
    }
})

test_that("gtrack.import_mappedseq rejects BAM sent through a named pipe", {
    skip_on_os("windows")
    skip_if(!nzchar(Sys.which("mkfifo")), "mkfifo not on PATH")
    bam <- bam_fixture("default.bam")
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
