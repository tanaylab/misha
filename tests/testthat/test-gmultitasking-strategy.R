# Correctness + heuristic tests for getOption("gmultitasking.strategy").
# Track-parallel must produce IDENTICAL output to tile-parallel for any
# query it accepts; the heuristic must keep small queries on tiles so we
# don't regress them.

create_isolated_test_db()

with_strategy <- function(strategy, expr) {
    old <- options(gmultitasking.strategy = strategy)
    on.exit(options(old), add = TRUE)
    eval.parent(substitute(expr))
}

# Build a multi-track fixture once for all correctness tests.
.mt_fixture <- function() {
    chrom_len <- gintervals.all()[gintervals.all()$chrom == "chr1", "end"]
    full_iv <- gintervals(1, 0, chrom_len)
    tracks <- vapply(seq_len(6L), function(i) {
        nm <- random_track_name(prefix = "test")
        if (i %% 2L == 0L) {
            gtrack.create_dense(nm, "strategy fixture", full_iv,
                values = runif(1), binsize = 50L, defval = 0
            )
        } else {
            ints <- gintervals(
                1,
                seq(0, 5e5, by = 200),
                seq(100, 5e5 + 100, by = 200)
            )
            ints <- ints[ints$end <= chrom_len, ]
            gtrack.create_sparse(nm, "strategy fixture", ints,
                values = runif(nrow(ints))
            )
        }
        nm
    }, character(1))
    list(
        tracks = tracks,
        scope_small = gintervals(1, 1e6, 1e6 + 5000L),
        scope_medium = gintervals(1, 1e6, 1e6 + 1e5)
    )
}

test_that("track-parallel returns identical data.frame to tile-parallel", {
    fix <- .mt_fixture()
    on.exit(for (t in fix$tracks) try(gtrack.rm(t, force = TRUE), silent = TRUE),
        add = TRUE
    )

    res_tiles <- with_strategy(
        "tiles",
        gextract(fix$tracks, intervals = fix$scope_medium, iterator = 50)
    )
    res_tracks <- with_strategy(
        "tracks",
        gextract(fix$tracks, intervals = fix$scope_medium, iterator = 50)
    )

    expect_equal(nrow(res_tracks), nrow(res_tiles))
    expect_equal(names(res_tracks), names(res_tiles))
    expect_equal(
        res_tracks[, c("chrom", "start", "end")],
        res_tiles[, c("chrom", "start", "end")]
    )
    for (col in fix$tracks) {
        expect_equal(res_tracks[[col]], res_tiles[[col]],
            tolerance = 1e-12, label = col
        )
    }
})

test_that("track-parallel respects custom colnames", {
    fix <- .mt_fixture()
    on.exit(for (t in fix$tracks) try(gtrack.rm(t, force = TRUE), silent = TRUE),
        add = TRUE
    )
    custom <- paste0("col_", seq_along(fix$tracks))

    res_tracks <- with_strategy(
        "tracks",
        gextract(fix$tracks,
            intervals = fix$scope_small, iterator = 50,
            colnames = custom
        )
    )
    expect_true(all(custom %in% names(res_tracks)))
    # Column order: interval cols, then values, then intervalID.
    n <- ncol(res_tracks)
    expect_equal(names(res_tracks)[seq.int(4, 4 + length(custom) - 1)], custom)
})

test_that("auto strategy picks 'tiles' for tiny queries", {
    skip_if_not(.ggetOption("gmultitasking"))
    iv <- gintervals(1, 0, 1000L)
    expect_equal(.gmultitasking_strategy(c("a", "b"), iv, iterator = iv), "tiles")
})

test_that("auto strategy picks 'tracks' for many-track interval-based queries", {
    skip_if_not(.ggetOption("gmultitasking"))
    iv <- data.frame(chrom = "chr1", start = 0:9999, end = 1:10000L)
    # 200 tracks × 10k intervals — interval-based iterator → tracks.
    expect_equal(
        .gmultitasking_strategy(rep("trk", 200), iv, iterator = iv),
        "tracks"
    )
})

test_that("auto strategy stays on 'tiles' for streaming (numeric) iterators", {
    # Critical: streaming iterators are catastrophic for track-parallel
    # (measured 6-17x slower because per-worker scan is unsplit). The
    # heuristic must NEVER pick 'tracks' for them.
    iv <- data.frame(chrom = "chr1", start = 0:9999, end = 1:10000L)
    expect_equal(
        .gmultitasking_strategy(rep("trk", 200), iv, iterator = 50),
        "tiles"
    )
    expect_equal(
        .gmultitasking_strategy(rep("trk", 200), iv, iterator = NULL),
        "tiles"
    )
    # Numeric vector (2D rect iterator) must also stay on tiles.
    expect_equal(
        .gmultitasking_strategy(rep("trk", 200), iv, iterator = c(50, 50)),
        "tiles"
    )
})

test_that("auto strategy picks 'tracks' for saved-intervals-set iterator", {
    # A character iterator that names a saved intervals set is interval-based
    # (same row geometry as a data.frame iterator); track-parallel benefits
    # apply.
    iv <- data.frame(chrom = "chr1", start = 0:9999, end = 1:10000L)
    set_name <- random_track_name(prefix = "test")
    gintervals.save(set_name, iv)
    on.exit(try(gintervals.rm(set_name, force = TRUE), silent = TRUE), add = TRUE)
    expect_equal(
        .gmultitasking_strategy(rep("trk", 200), iv, iterator = set_name),
        "tracks"
    )
})

test_that("auto strategy stays on 'tiles' when character iterator is not a saved set", {
    # Unknown character iterator (could be a track name → dense streaming).
    # Without confirming it's an intervals set, stay conservative.
    iv <- data.frame(chrom = "chr1", start = 0:9999, end = 1:10000L)
    expect_equal(
        .gmultitasking_strategy(rep("trk", 200), iv,
            iterator = "definitely_not_a_set_xyz"
        ),
        "tiles"
    )
})

test_that("auto strategy stays on 'tiles' for fewer than 8 tracks", {
    # Below 8 tracks the parallelism is too capped and warm-cache penalty
    # outweighs the cold-cache benefit (matrix bench: borderline at n=5).
    iv <- data.frame(chrom = "chr1", start = 0:9999, end = 1:10000L)
    for (n in c(2, 5, 7)) {
        expect_equal(
            .gmultitasking_strategy(rep("trk", n), iv, iterator = iv),
            "tiles"
        )
    }
    # 8 tracks crosses the threshold (large enough query).
    expect_equal(
        .gmultitasking_strategy(rep("trk", 8), iv, iterator = iv),
        "tracks"
    )
})

test_that("auto strategy stays on 'tiles' when output is a file", {
    skip_if_not(.ggetOption("gmultitasking"))
    iv <- data.frame(chrom = "chr1", start = 0:9999, end = 1:10000L)
    expect_equal(
        .gmultitasking_strategy(rep("trk", 200), iv,
            iterator = iv,
            file = "/tmp/x.tab"
        ),
        "tiles"
    )
    expect_equal(
        .gmultitasking_strategy(rep("trk", 200), iv,
            iterator = iv,
            intervals.set.out = "set"
        ),
        "tiles"
    )
})

test_that("auto strategy stays on 'tiles' for single-track queries", {
    iv <- data.frame(chrom = "chr1", start = 0:9999, end = 1:10000L)
    expect_equal(.gmultitasking_strategy("just_one", iv, iterator = iv), "tiles")
})

test_that("auto strategy stays on 'tiles' for 2D band queries", {
    iv <- data.frame(chrom = "chr1", start = 0:9999, end = 1:10000L)
    expect_equal(
        .gmultitasking_strategy(rep("trk", 200), iv,
            iterator = iv,
            band = c(0, 1e6)
        ),
        "tiles"
    )
})

test_that("explicit strategy override wins over heuristic", {
    iv <- gintervals(1, 0, 1000L)
    with_strategy(
        "tracks",
        expect_equal(
            .gmultitasking_strategy(c("a", "b"), iv, iterator = iv),
            "tracks"
        )
    )
    with_strategy(
        "tiles",
        expect_equal(
            .gmultitasking_strategy(c("a", "b"), iv, iterator = iv),
            "tiles"
        )
    )
})

test_that("unknown strategy falls back to auto with warning", {
    iv <- gintervals(1, 0, 1000L)
    with_strategy("nonsense", {
        expect_warning(
            s <- .gmultitasking_strategy(c("a", "b"), iv,
                iterator = iv
            ),
            "Unknown gmultitasking.strategy"
        )
        expect_equal(s, "tiles") # auto picks tiles for 2-track tiny query
    })
})

test_that("strategy='tracks' and 'tiles' return identical row counts for tiny iterators", {
    # Sanity check at the small end: a few-row scope should produce the same
    # output structure regardless of strategy.
    fix <- .mt_fixture()
    on.exit(for (t in fix$tracks) try(gtrack.rm(t, force = TRUE), silent = TRUE),
        add = TRUE
    )
    tiny <- gintervals(1, 1e6, 1e6 + 250L)
    res_tiles <- with_strategy(
        "tiles",
        gextract(fix$tracks, intervals = tiny, iterator = 50)
    )
    res_tracks <- with_strategy(
        "tracks",
        gextract(fix$tracks, intervals = tiny, iterator = 50)
    )
    expect_equal(nrow(res_tracks), nrow(res_tiles))
    expect_equal(names(res_tracks), names(res_tiles))
})
