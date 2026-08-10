create_isolated_test_db()

# COVERAGE LIMIT, read before trusting this file: every test here asserts
# parallel == sequential. That is the right correctness property, but it is also
# satisfied when the multitasking path never runs - so if the gate silently stopped
# engaging (SEQ_MIN_INTERVALS4MT raised, the dynamic_cast regressing to NULL,
# prepare4multitasking_whole_intervals returning 0), every test below would still
# pass. There is currently no observable from R that says "this call forked", so
# that failure mode is unguarded here. It was checked by hand for this change
# (forced-parallel is ~1000x slower than sequential on warm data, and 41 live
# children were observed during a long extraction) but nothing re-checks it.

# gseq.extract decides whether to distribute by timing its first few reads, and the test DB is
# always in the page cache, so it would (correctly) stay sequential here. gseq.extract.probe.usec = 0
# means "always distribute" and is what lets these tests reach the multitasking path at all.
parallel_opts <- function(processes = 4) {
    list(gmultitasking = TRUE, gseq.extract.probe.usec = 0, gmax.processes = processes)
}

# Below 500 intervals gseq.extract never distributes, so every fixture here is comfortably above it.
rand_intervs <- function(n, size = 100, chroms = NULL) {
    all_chroms <- gintervals.all()
    if (!is.null(chroms)) {
        all_chroms <- all_chroms[as.character(all_chroms$chrom) %in% chroms, ]
    }
    i <- sample(nrow(all_chroms), n, replace = TRUE)
    starts <- floor(runif(n) * (all_chroms$end[i] - size))
    data.frame(
        chrom = as.character(all_chroms$chrom[i]),
        start = starts, end = starts + size,
        stringsAsFactors = FALSE
    )
}

extract_seq <- function(intervals) {
    withr::with_options(list(gmultitasking = FALSE), gseq.extract(intervals))
}

extract_par <- function(intervals, processes = 4) {
    withr::with_options(parallel_opts(processes), gseq.extract(intervals))
}

test_that("gseq.extract multitasking matches sequential across chromosomes", {
    set.seed(1)
    intervs <- rand_intervs(2000)
    expect_identical(extract_par(intervs), extract_seq(intervs))
})

test_that("gseq.extract multitasking splits within a single chromosome", {
    # The generic scope splitter only cuts on chromosome boundaries; gseq.extract needs whole
    # intervals, so it uses its own splitter, and a single-chromosome set is what exercises it.
    set.seed(2)
    chrom <- as.character(gintervals.all()$chrom[1])
    intervs <- rand_intervs(2000, chroms = chrom)
    expect_identical(extract_par(intervs), extract_seq(intervs))
})

test_that("gseq.extract multitasking preserves the original row order", {
    set.seed(3)
    intervs <- rand_intervs(1500)
    intervs <- intervs[sample(nrow(intervs)), ]
    rownames(intervs) <- NULL

    seqs <- extract_par(intervs)
    expect_identical(seqs, extract_seq(intervs))
    # independent check: a handful of rows extracted one at a time
    idx <- sample(nrow(intervs), 10)
    single <- vapply(idx, function(i) extract_seq(intervs[i, ]), character(1))
    expect_identical(single, seqs[idx])
})

test_that("gseq.extract multitasking reverse-complements per interval", {
    set.seed(4)
    intervs <- rand_intervs(1500)
    intervs$strand <- sample(c(-1L, 1L), nrow(intervs), replace = TRUE)
    expect_identical(extract_par(intervs), extract_seq(intervs))

    intervs$strand <- -1L
    expect_identical(extract_par(intervs), extract_seq(intervs))
})

test_that("gseq.extract multitasking handles duplicated and overlapping intervals", {
    set.seed(5)
    intervs <- rand_intervs(100)
    intervs <- intervs[rep(seq_len(nrow(intervs)), 15), ]
    rownames(intervs) <- NULL
    expect_identical(extract_par(intervs), extract_seq(intervs))
})

test_that("gseq.extract multitasking balances heavily skewed interval sizes", {
    set.seed(6)
    intervs <- rand_intervs(1500, size = 1)
    big <- rand_intervs(20, size = 100000)
    intervs <- rbind(intervs, big)
    expect_identical(extract_par(intervs), extract_seq(intervs))
})

test_that("gseq.extract result does not depend on the number of processes", {
    set.seed(7)
    intervs <- rand_intervs(2000)
    expected <- extract_seq(intervs)
    for (processes in c(1, 2, 8, 16)) {
        expect_identical(extract_par(intervs, processes), expected)
    }
})

test_that("gseq.extract falls back to sequential for big intervals sets", {
    set.seed(8)
    intervs <- rand_intervs(1500)
    withr::local_options(list(gbig.intervals.size = 100))
    set_name <- paste0("bigseq_", paste(sample(letters, 8, replace = TRUE), collapse = ""))
    gintervals.save(set_name, intervs)
    withr::defer(gintervals.rm(set_name, force = TRUE))

    expect_true(misha:::.gintervals.is_bigset(set_name))
    expect_identical(extract_par(set_name), extract_seq(set_name))
    expect_identical(extract_par(set_name), extract_seq(gintervals.load(set_name)))
})

test_that("gseq.extract timing probe does not disturb the sequential result", {
    # With default options on a cached DB the probe reads the first few intervals and the rest is
    # read sequentially; an off-by-one in that hand-off would corrupt or blank the leading rows.
    set.seed(10)
    for (n in c(499, 500, 501, 2000)) {
        intervs <- rand_intervs(n)
        expect_identical(
            withr::with_options(list(gmultitasking = TRUE), gseq.extract(intervs)),
            extract_seq(intervs)
        )
    }
})

test_that("gseq.extract multitasking still enforces gmax.data.size", {
    # The limit has to be above what the 16-interval timing probe reads (16 x 100 b) and below the
    # whole result, or the probe's own check throws and the multitasking path is never reached.
    set.seed(9)
    intervs <- rand_intervs(1500) # 150,000 b total
    limit <- 50000

    parallel_err <- tryCatch(
        withr::with_options(c(parallel_opts(), list(gmax.data.size = limit)), gseq.extract(intervs)),
        error = conditionMessage
    )
    sequential_err <- tryCatch(
        withr::with_options(list(gmultitasking = FALSE, gmax.data.size = limit), gseq.extract(intervs)),
        error = conditionMessage
    )

    # Both paths must refuse, and with the message that names gmax.data.size - the multitasking
    # path used to fail later and with the generic shared-memory wording instead.
    expect_match(parallel_err, "Result sequence size")
    expect_match(parallel_err, "gmax.data.size")
    expect_match(sequential_err, "Result sequence size")
})
