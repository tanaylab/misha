# Regression tests for the 2026-06-28 full-audit fixes. See
# dev/notes/2026-06-28_full-audit.md. Uses the hg-scale isolated test DB
# (real sequence symlinked) so large coordinates and PWM scoring are available.

create_isolated_test_db()

test_that("H2: gintervals.quantiles bigset fills un-scanned scope chromosomes with NaN (not stale values)", {
    # The buggy path is the NON-multitask streaming-bigset save: it sized the
    # finish-loop buffer `size` instead of `size*num_percentiles`, so a skipped
    # scope chromosome's read ran out of bounds (verified via instrumentation:
    # idx 10..29 into a size-10 vector). Force that path: gmultitasking=FALSE and
    # a low big-intervals threshold so the result streams to a bigset. chr3 has
    # scope intervals but no iterator coverage -> it is the skipped chromosome.
    withr::local_options(list(gmultitasking = FALSE, gbig.intervals.size = 10))
    temp_track_name <- paste0("test.tmpq_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    # chr2: many wide scanned intervals (finite quantiles, fills the buffer);
    # chr3: small skipped intervals (read into the buffer's tail).
    s2 <- gintervals(2, seq(0, by = 1000, length.out = 100), seq(1000, by = 1000, length.out = 100))
    s3 <- gintervals(3, seq(0, by = 20, length.out = 10), seq(10, by = 20, length.out = 10))
    scope <- rbind(s2, s3)
    iter <- s2 # iterator covers chr2 only -> chr3 is the un-scanned scope chrom
    pct <- c(0.2, 0.5, 0.9)

    gintervals.quantiles("test.fixedbin", pct,
        intervals = scope, iterator = iter, intervals.set.out = temp_track_name
    )
    in_memory <- gintervals.quantiles("test.fixedbin", pct, intervals = scope, iterator = iter)
    from_disk <- gintervals.load(temp_track_name)

    expect_equal(
        from_disk[order(from_disk$chrom, from_disk$start), ],
        in_memory[order(in_memory$chrom, in_memory$start), ],
        ignore_attr = TRUE
    )
    # chr3 (skipped) quantiles must all be NaN, never stale finite values.
    chr3 <- from_disk[from_disk$chrom == "chr3", 4:6]
    expect_true(all(is.na(unlist(chr3))))
})

test_that("H6: bedGraph export writes large coordinates as plain integers", {
    # gextract returns coords as doubles; paste() rendered round large values
    # (1e8, 2e8) in scientific notation, which bedGraphToBigWig mis-parses.
    track <- paste0("test.tmpexport_", sample(1:1e9, 1))
    gtrack.rm(track, force = TRUE)
    withr::defer(gtrack.rm(track, force = TRUE))

    ivs <- gintervals(c(1, 1, 1), c(1e8, 15e7, 2e8), c(1e8 + 50, 15e7 + 50, 2e8 + 50))
    gtrack.create_sparse(track, "audit H6", ivs, c(1, 2, 3))

    outfile <- tempfile(fileext = ".bedgraph")
    withr::defer(unlink(outfile))
    gtrack.export_bedgraph(track, outfile)

    cols <- strsplit(readLines(outfile)[-1], "\t")
    starts <- vapply(cols, `[`, character(1), 2)
    ends <- vapply(cols, `[`, character(1), 3)
    expect_false(any(grepl("e", c(starts, ends), ignore.case = TRUE)))
    expect_true(all(c("100000000", "200000000") %in% starts))
})

test_that("H7: track-parallel gextract returns value columns in requested track order", {
    # Round-robin split across workers then merge-in-worker-order scrambled the
    # positional column order (names stayed correct). Force the track-parallel
    # path with more tracks than workers so the scramble would occur.
    withr::local_options(list(
        gmultitasking = TRUE,
        gmultitasking.strategy = "tracks",
        gmax.processes = 3
    ))
    exprs <- c(
        "test.fixedbin", "test.fixedbin+10", "test.fixedbin+20",
        "test.fixedbin+30", "test.fixedbin+40", "test.fixedbin+50"
    )
    res <- gextract(exprs, gintervals(1, 0, 50000))
    expect_false(is.null(res))

    val_names <- setdiff(names(res), c("chrom", "start", "end", "intervalID"))
    # Positional order must match requested order (the bug scrambled this).
    expect_identical(val_names, exprs)
    # And names stay attached to the right data: col k == base + offset.
    base <- res[["test.fixedbin"]]
    for (off in c(10, 20, 30, 40, 50)) {
        expect_equal(res[[paste0("test.fixedbin+", off)]], base + off)
    }
})

test_that("H8: updating the cache from a database subdirectory does not write a truncated root cache", {
    # In a subdirectory, gdb.reload scans only the subtree, so GTRACKS is a
    # partial, subdir-relative list. The unfixed .gdb.cache_update_lists wrote
    # that partial list to the ROOT .db.cache and cleared the dirty flag, so the
    # next session loaded a truncated track listing. The fix marks the cache
    # dirty instead (forcing a rescan).
    iv <- gintervals(1, 0, 1000)
    gdir.create("auditsub", showWarnings = FALSE)
    gtrack.create_sparse("auditsub.tA", "audit", iv, 1)
    gdb.reload(rescan = TRUE)
    withr::defer(suppressMessages(gdb.init(.misha$GROOT))) # restore GWD to root
    full_n <- length(gtrack.ls())
    expect_gt(full_n, 1)

    gdir.cd("auditsub")
    expect_lt(length(get("GTRACKS", envir = .misha)), full_n) # GTRACKS now partial
    misha:::.gdb.cache_update_lists(get("GROOT", envir = .misha))
    expect_true(misha:::.gdb.cache_is_dirty(get("GROOT", envir = .misha)))
})

test_that("H4: pwm.edit_distance vtrack with a >64bp motif does not read out of bounds", {
    # The subs-only solvers indexed fixed [64] tables over the full motif length;
    # a motif longer than 64bp read past them. Just exercise the path on real
    # sequence and require a finite/NaN result with no crash.
    remove_all_vtracks()
    set.seed(1)
    L <- 70
    pssm <- matrix(runif(L * 4, 0.1, 0.9), nrow = L, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("edist_long", NULL, "pwm.edit_distance",
        pssm = pssm, score.thresh = -50, score.min = -120, score.max = -40,
        max_edits = 3L, max_indels = 0L, bidirect = TRUE
    )
    interv <- gintervals(1, 100000, 100200)
    val <- gextract("edist_long", interv, iterator = interv)
    expect_true(is.na(val$edist_long) || is.finite(val$edist_long))
})
