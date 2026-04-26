# Performance regression tests.
#
# These guard against the kind of subtle slowdowns that a C++ optimization
# audit accidentally introduced — they don't need to be precise benchmarks,
# only catch order-of-magnitude regressions in setup/per-track overhead.
#
# OPT-IN: skipped by default. They allocate tracks, do real I/O, and depend
# on wall-clock timing — running them inside a parallel `devtools::test()`
# under load would make timings noisy and slow the whole suite. To run them:
#
#     MISHA_PERF_TESTS=true R -e "devtools::test(filter='perf-regression')"
#
# or, from R:
#
#     Sys.setenv(MISHA_PERF_TESTS = "true")
#     devtools::test(filter = "perf-regression")

skip_unless_perf <- function() {
    if (!identical(tolower(Sys.getenv("MISHA_PERF_TESTS")), "true")) {
        skip("perf test (set MISHA_PERF_TESTS=true to run)")
    }
}

create_isolated_test_db()

test_that("gextract setup over many dense tracks stays fast (MAP_POPULATE regression)", {
    # The bug: a perf audit added MmapFile with MAP_POPULATE, which forced the
    # whole per-chrom track file into RAM at every mmap. With many tracks per
    # call, even a warm-cache extract on a tiny interval paid seconds of
    # page-table walking. v5.6.7 ran this in <1s; v5.6.11–v5.6.16 took 10–20s
    # with realistic motif track sizes (~38MB chr1 × 51 tracks).
    skip_unless_perf()

    n_tracks <- 30
    bin_size <- 50L
    chrom <- 1L
    chrom_len <- gintervals.all()[gintervals.all()$chrom == paste0("chr", chrom), "end"]
    if (length(chrom_len) == 0) skip("chr1 not found in chromkey")

    # Per-chrom file ~chrom_len/bin_size * 4 bytes ~= 16MB on chr1 — big enough
    # that MAP_POPULATE has real work to do, small enough that creating 30 of
    # them fits in the test budget.
    full_iv <- gintervals(chrom, 0, chrom_len)
    tracks <- vapply(seq_len(n_tracks), function(i) {
        nm <- random_track_name(prefix = "test")
        gtrack.create_dense(nm, "perf regression fixture", full_iv,
            values = runif(1), binsize = bin_size, defval = 0
        )
        nm
    }, character(1))
    on.exit(
        {
            for (t in tracks) try(gtrack.rm(t, force = TRUE), silent = TRUE)
        },
        add = TRUE
    )

    # Tiny interval — bug or no bug, the actual data we read is the same handful
    # of bins. The cost we're guarding against is per-track mmap setup.
    small <- gintervals(chrom, 1e6, 1e6 + 1000)

    invisible(gextract(tracks, small, iterator = bin_size)) # warm caches

    elapsed <- system.time(gextract(tracks, small, iterator = bin_size))[["elapsed"]]

    # Without the regression this is ~0.1s on a typical lab box; the regressed
    # code took ~5s on the same hardware (and ~14s with realistic motif sizes).
    # 2s catches a >5x slowdown without flaking on noisy hardware.
    expect_lt(elapsed, 2)
})
