# Second adversarial pass for the 2026-06-28 audit fixes: targets fixes the first
# pass did not cover (L4, H1, M3, M8, M11, L2) and probes their boundaries.

create_isolated_test_db()

test_that("ADV L4: a 2D bigset written under the indexed format round-trips", {
    # Exercises gbigintervs_2d_indexed_finalize, the streaming writer fixed in L4
    # (the convert path the existing suite covers is a different code path). With the
    # broken 36-byte/24-byte layout the load fails ("Failed to read entry 1").
    withr::local_options(list(gmulticontig.indexed_format = TRUE, gbig.intervals.size = 10))
    n <- 40
    ivs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = n), starts1 = rep(seq(0, by = 1000, length.out = n), 2),
        ends1 = rep(seq(500, by = 1000, length.out = n), 2),
        chroms2 = rep(c(1, 2), each = n), starts2 = rep(seq(0, by = 1000, length.out = n), 2),
        ends2 = rep(seq(500, by = 1000, length.out = n), 2)
    )
    nm <- paste0("test.idx2d_", sample(1:1e9, 1))
    gintervals.rm(nm, force = TRUE)
    withr::defer(gintervals.rm(nm, force = TRUE))
    gintervals.save(intervals = ivs, intervals.set.out = nm)
    back <- gintervals.load(nm)
    ord <- function(d) d[order(d$chrom1, d$start1, d$chrom2, d$start2), ]
    expect_equal(nrow(back), nrow(ivs))
    expect_equal(ord(back), ord(ivs), ignore_attr = TRUE)
})

test_that("ADV H1: track-expression errors with '%' propagate verbatim and never crash", {
    msgs <- c("plain error", "value is 50%", "fmt %d %s %n marker", "literal %% sign", "100% done at /a%2Fb")
    for (m in msgs) {
        e <- tryCatch(
            gextract(sprintf('stop("%s")', m), gintervals(1, 0, 1000), iterator = 1000),
            error = function(e) conditionMessage(e)
        )
        expect_true(grepl(m, e, fixed = TRUE), info = m)
    }
})

test_that("ADV M11: a failed gsetroot leaves the current session fully usable", {
    cur <- get("GROOT", envir = .misha)
    n_tracks <- length(gtrack.ls())
    bad <- tempfile(pattern = "baddb_")
    dir.create(file.path(bad, "tracks"), recursive = TRUE)
    dir.create(file.path(bad, "seq"))
    withr::defer(unlink(bad, recursive = TRUE))

    expect_error(gsetroot(bad), "chrom_sizes")
    # session intact: same root, genome loaded, tracks listable and extractable
    expect_identical(get("GROOT", envir = .misha), cur)
    expect_false(is.null(get("ALLGENOME", envir = .misha)))
    expect_equal(length(gtrack.ls()), n_tracks)
    expect_true(is.data.frame(gextract("test.fixedbin", gintervals(1, 0, 1000))))
})

test_that("ADV M8: gtrack.var.ls usage - rejects a missing track, accepts name and expression", {
    gtrack.create_sparse("m8btrk", "t", gintervals(1, 0, 10000), 1)
    withr::defer(gtrack.rm("m8btrk", force = TRUE))
    gtrack.var.set("m8btrk", "vv", 7)
    expect_error(gtrack.var.ls(), "Usage")
    expect_true("vv" %in% gtrack.var.ls("m8btrk"))
    trks <- c("m8btrk")
    expect_true("vv" %in% gtrack.var.ls(trks[1]))
})

test_that("ADV M3: gcompute_strands_autocorr handles reverse reads at a contig end without crashing", {
    chromsize <- gintervals.all()$end[gintervals.all()$chrom == "chr1"]
    mk_reads <- function(path, df) {
        lines <- vapply(seq_len(nrow(df)), function(i) {
            cols <- as.character(rep("x", 14))
            cols[9] <- df$seq[i]
            cols[11] <- df$chrom[i]
            cols[13] <- as.character(df$coord[i])
            cols[14] <- df$strand[i]
            paste(cols, collapse = "\t")
        }, character(1))
        writeLines(lines, path)
    }
    rf <- tempfile()
    withr::defer(unlink(rf))
    s50 <- paste(rep("A", 50), collapse = "")
    reads <- data.frame(
        chrom = "chr1",
        coord = c(1000, 2000, 3000, chromsize - 10), # last: reverse read past the end
        strand = c("+", "-", "+", "-"),
        seq = s50, stringsAsFactors = FALSE
    )
    mk_reads(rf, reads)
    # must not crash / corrupt (the M3 clamp); returns without error
    expect_error(
        gcompute_strands_autocorr(rf, "chr1", 1000, maxread = 100, max.coord = 3e8),
        NA
    )
})

test_that("ADV L2: ggenome.implant succeeds in-bounds and errors cleanly out-of-bounds", {
    ref <- tempfile(fileext = ".fa")
    out <- tempfile(fileext = ".fa")
    withr::defer(unlink(c(ref, out, paste0(out, ".fai"))))
    cat(">chrA\nACGTACGTAC\n", file = ref) # 10 bp

    # success path
    res <- ggenome.implant(data.frame(chrom = "chrA", start = 2, end = 6), "NNNN", out,
        genome_fasta = ref, create_trackdb = FALSE, overwrite = TRUE
    )
    expect_equal(res, out)

    # out-of-bounds interval -> clean error (the verror path that now runs destructors)
    expect_error(
        ggenome.implant(data.frame(chrom = "chrA", start = 100, end = 104), "NNNN", tempfile(fileext = ".fa"),
            genome_fasta = ref, create_trackdb = FALSE, overwrite = TRUE
        )
    )
})
