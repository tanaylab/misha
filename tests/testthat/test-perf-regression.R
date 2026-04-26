# Performance regression tests.
#
# These guard against the kind of subtle slowdowns that a C++ optimization
# audit accidentally introduced — they don't need to be precise benchmarks,
# only catch order-of-magnitude regressions in setup/per-track overhead.
#
# Helper functions and the budget rule live in helper-perf.R. Each test
# embeds a single `baseline_ms` measured on the lab box. Budget for a test is
# `baseline_ms * 3`; the 3x factor keeps noisy NFS hardware from false-firing
# while still catching genuine >3x slowdowns.
#
# Pick `baseline_ms` as the FASTEST wall-clock you can reproduce for the
# operation across versions you care about. If a future fix makes an
# operation faster, lower the number (do not raise for convenience — the
# whole point is to keep the bar tight). If a planned change makes it
# legitimately slower, justify it before relaxing the baseline.
#
# OPT-IN: skipped by default. They allocate tracks, do real I/O, and depend
# on wall-clock timing — running them inside a parallel `devtools::test()`
# under load would make timings noisy and slow the whole suite.
#
#     MISHA_PERF_TESTS=true R -e "devtools::test(filter='perf-regression')"

create_isolated_test_db()

# Shared fixtures (created once, reused across tests via gtrack.* lookups).
.perf_setup <- function(n_dense_tracks = 30, dense_bin_size = 50L) {
    chrom_len <- gintervals.all()[gintervals.all()$chrom == "chr1", "end"]
    full_iv <- gintervals(1, 0, chrom_len)

    dense_tracks <- vapply(seq_len(n_dense_tracks), function(i) {
        nm <- random_track_name(prefix = "test")
        gtrack.create_dense(nm, "perf fixture", full_iv,
            values = runif(1), binsize = dense_bin_size, defval = 0
        )
        nm
    }, character(1))

    sparse_track <- random_track_name(prefix = "test")
    sparse_intervs <- gintervals(1, seq(0, 1e6, by = 100), seq(50, 1e6 + 50, by = 100))
    sparse_intervs <- sparse_intervs[sparse_intervs$end <= chrom_len, ]
    gtrack.create_sparse(sparse_track, "perf fixture", sparse_intervs,
        values = runif(nrow(sparse_intervs))
    )

    list(
        dense_tracks = dense_tracks,
        sparse_track = sparse_track,
        chrom1_len = chrom_len,
        full_chr1 = full_iv,
        small_chr1 = gintervals(1, 1e6, 1e6 + 1000)
    )
}

.perf_cleanup <- function(fix) {
    for (t in fix$dense_tracks) try(gtrack.rm(t, force = TRUE), silent = TRUE)
    try(gtrack.rm(fix$sparse_track, force = TRUE), silent = TRUE)
}

# --- Tests ---

test_that("gextract setup over many dense tracks stays fast (MAP_POPULATE regression)", {
    # Bug history: v5.6.11 added MmapFile with MAP_POPULATE which forced
    # eager paging-in of every per-chrom file at every mmap. v5.6.17 removed
    # MAP_POPULATE and switched the per-chrom-per-track validation loops to a
    # metadata-only path. This test would have caught both.
    skip_unless_perf()
    fix <- .perf_setup(n_dense_tracks = 30)
    on.exit(.perf_cleanup(fix), add = TRUE)

    measured <- time_op(gextract(fix$dense_tracks, fix$small_chr1, iterator = 50))
    expect_perf_baseline(measured, "many-dense setup", baseline_ms = 112)
})

test_that("gextract single dense track full chr1 stays fast (bin-scan path)", {
    # Guards the inner read loop in GenomeTrackFixedBin::read_interval — the
    # hot path that the perf audit's mmap zero-copy was supposed to speed up.
    # On chr1 (~250Mbp) at binsize=50 this scans ~5M bins.
    skip_unless_perf()
    fix <- .perf_setup(n_dense_tracks = 1)
    on.exit(.perf_cleanup(fix), add = TRUE)

    measured <- time_op(gextract(fix$dense_tracks[1], fix$full_chr1, iterator = 50))
    expect_perf_baseline(measured, "dense full-chr scan", baseline_ms = 2005)
})

test_that("gextract sparse track full chr1 stays fast", {
    # Guards GenomeTrackSparse::read_interval — uses a separate code path
    # (BufferedFile + load-into-memory) that didn't benefit from MmapFile but
    # could regress independently if the perf audit touched its hot loop.
    skip_unless_perf()
    fix <- .perf_setup(n_dense_tracks = 0)
    on.exit(.perf_cleanup(fix), add = TRUE)

    measured <- time_op(gextract(fix$sparse_track, fix$full_chr1))
    expect_perf_baseline(measured, "sparse full-chr", baseline_ms = 20)
})

test_that("gextract vtrack with avg + window stays fast", {
    # Realistic vtrack pattern: avg over a +/-N window. Touches the
    # iterator-modifier code path that the perf audit's hash-key changes
    # (BackendKey struct) sit on top of.
    skip_unless_perf()
    fix <- .perf_setup(n_dense_tracks = 1)
    on.exit(.perf_cleanup(fix), add = TRUE)
    on.exit(try(gvtrack.rm("v_avg"), silent = TRUE), add = TRUE)

    gvtrack.create("v_avg", fix$dense_tracks[1], func = "avg")
    gvtrack.iterator("v_avg", sshift = -250, eshift = 250)

    measured <- time_op(gextract("v_avg", fix$small_chr1, iterator = 50))
    expect_perf_baseline(measured, "vtrack avg+window", baseline_ms = 18)
})

test_that("gextract vtrack with LSE sliding window stays fast", {
    # LSE is the function Tamar uses for motif scoring. Worth its own test
    # because the LSE state machine in GenomeTrackFixedBin maintains a sliding
    # log-sum-exp across bins and is touched by several perf-audit changes.
    skip_unless_perf()
    fix <- .perf_setup(n_dense_tracks = 1)
    on.exit(.perf_cleanup(fix), add = TRUE)
    on.exit(try(gvtrack.rm("v_lse"), silent = TRUE), add = TRUE)

    gvtrack.create("v_lse", fix$dense_tracks[1], func = "lse")
    gvtrack.iterator("v_lse", sshift = -100, eshift = 100)

    measured <- time_op(gextract("v_lse", fix$small_chr1, iterator = 50))
    expect_perf_baseline(measured, "vtrack lse+window", baseline_ms = 18)
})

test_that("gscreen on many dense tracks stays fast", {
    # gscreen is the second hottest entry point after gextract for Tamar's
    # workloads. Same setup-cost characteristic as gextract — the validation
    # loops go through the same per-track path.
    skip_unless_perf()
    fix <- .perf_setup(n_dense_tracks = 30)
    on.exit(.perf_cleanup(fix), add = TRUE)

    expr <- paste(fix$dense_tracks, collapse = " > 0 & ")
    expr <- paste(expr, "> 0")

    measured <- time_op(gscreen(expr, intervals = fix$small_chr1, iterator = 50))
    expect_perf_baseline(measured, "gscreen many-dense", baseline_ms = 109)
})

test_that("gsummary on dense full chr1 stays fast", {
    # Reduction path — single-function fast-path territory in
    # GenomeTrackFixedBin. The audit added a "mode 3" single-function fast
    # path; this test would catch a regression there.
    skip_unless_perf()
    fix <- .perf_setup(n_dense_tracks = 1)
    on.exit(.perf_cleanup(fix), add = TRUE)

    measured <- time_op(gsummary(fix$dense_tracks[1], intervals = fix$full_chr1))
    expect_perf_baseline(measured, "gsummary dense full-chr", baseline_ms = 268)
})

test_that("gquantiles on dense full chr1 stays fast", {
    # Streaming-percentile path (StreamPercentiler) — the audit templated its
    # comparator for inlining. Different code path than the reducer/agg paths.
    skip_unless_perf()
    fix <- .perf_setup(n_dense_tracks = 1)
    on.exit(.perf_cleanup(fix), add = TRUE)

    measured <- time_op(gquantiles(fix$dense_tracks[1], percentiles = c(0.1, 0.5, 0.9), intervals = fix$full_chr1))
    expect_perf_baseline(measured, "gquantiles dense full-chr", baseline_ms = 429)
})

test_that("2D gextract on rect track full chr1xchr1 stays fast", {
    # 2D path is independent of the 1D hot paths above. The perf audit
    # touched StatQuadTree (uint8_t instead of bool) — that would show up here.
    skip_unless_perf()
    rect_tracks <- gtrack.ls("test\\.rects$")
    if (length(rect_tracks) == 0) skip("no rects track in test_db")

    iv2d <- gintervals.2d("chr1", 0, 1e7, "chr1", 0, 1e7)

    measured <- time_op(gextract(rect_tracks[1], intervals = iv2d, iterator = c(1e5, 1e5)))
    expect_perf_baseline(measured, "2D rect extract", baseline_ms = 21)
})
