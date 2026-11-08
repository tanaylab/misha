restore_groot_on_exit()

# vtrack func = "weighted.sum" on 1D sources: the sum of value x the number of bp of each bin or interval that
# falls inside the vtrack window (after sshift / eshift, clipped to the chromosome). "sum" counts every bin the
# window touches in full.
#
# Every value below is a multiple of 1/1024 under 4, so it is stored exactly as a float, and every sum of products
# with bp counts is exact in double: the references are compared without tolerance.

.wsum_chrom_sizes <- c(chr1 = 20013, chr2 = 9991)

.wsum_random_seq <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")

# Random non-overlapping intervals with gaps, some values NaN
.wsum_random_intervals <- function(chrom, len, nan_frac = 0.1) {
    starts <- c()
    ends <- c()
    pos <- sample(0:40, 1)
    while (TRUE) {
        size <- sample(1:120, 1)
        if (pos + size > len) {
            break
        }
        starts <- c(starts, pos)
        ends <- c(ends, pos + size)
        pos <- pos + size + sample(0:80, 1)
    }
    vals <- sample(1:4000, length(starts), replace = TRUE) / 1024
    vals[runif(length(vals)) < nan_frac] <- NaN
    data.frame(chrom = chrom, start = starts, end = ends, value = vals)
}

# A small database (per-chromosome or indexed) with dense tracks at several bin sizes, a constant track, a sparse
# track and an array track. The chromosome sizes do not divide by the bin sizes, so the last bin is short.
.wsum_make_db <- function(indexed) {
    set.seed(17)
    root <- tempfile(pattern = "misha_wsum_")
    dir.create(root)
    fastas <- file.path(root, paste0(names(.wsum_chrom_sizes), ".fa"))
    for (i in seq_along(fastas)) {
        cat(">", names(.wsum_chrom_sizes)[i], "\n", .wsum_random_seq(.wsum_chrom_sizes[[i]]), "\n",
            sep = "", file = fastas[i]
        )
    }
    db <- file.path(root, "db")
    withr::with_options(list(gmulticontig.indexed_format = indexed), {
        suppressMessages(gdb.create(groot = db, fasta = fastas, verbose = FALSE))
    })
    suppressMessages(gsetroot(db))

    for (binsize in c(1, 20, 32, 50)) {
        bins <- do.call(rbind, lapply(names(.wsum_chrom_sizes), function(chrom) {
            len <- .wsum_chrom_sizes[[chrom]]
            starts <- seq(0, len - 1, by = binsize)
            data.frame(chrom = chrom, start = starts, end = pmin(starts + binsize, len))
        }))
        vals <- sample(1:4000, nrow(bins), replace = TRUE) / 1024
        # runs of bins without a value (defval NaN): a few long ones and scattered single bins
        keep <- runif(nrow(bins)) > 0.05
        for (run in sample(seq_len(nrow(bins) - 40), 4)) {
            keep[run:(run + 39)] <- FALSE
        }
        gtrack.create_dense(paste0("d", binsize), "", bins[keep, ], vals[keep], binsize, defval = NaN)
    }

    const_bins <- gintervals(names(.wsum_chrom_sizes), 0, .wsum_chrom_sizes)
    gtrack.create_dense("const20", "", const_bins, rep(1, nrow(const_bins)), 20, defval = NaN)

    sp <- do.call(rbind, lapply(names(.wsum_chrom_sizes), function(chrom) .wsum_random_intervals(chrom, .wsum_chrom_sizes[[chrom]])))
    gtrack.create_sparse("sp", "", sp[, c("chrom", "start", "end")], sp$value)

    arr <- .wsum_random_intervals("chr1", .wsum_chrom_sizes[["chr1"]], nan_frac = 0)
    arr <- rbind(arr, .wsum_random_intervals("chr2", .wsum_chrom_sizes[["chr2"]], nan_frac = 0))
    arr$col1 <- arr$value
    arr$col2 <- sample(1:4000, nrow(arr), replace = TRUE) / 1024
    arr$col2[runif(nrow(arr)) < 0.3] <- NaN
    arr$col3 <- sample(1:4000, nrow(arr), replace = TRUE) / 1024
    arr_file <- file.path(root, "arr.tsv")
    write.table(arr[, c("chrom", "start", "end", "col1", "col2", "col3")], arr_file,
        sep = "\t", quote = FALSE, row.names = FALSE, na = ""
    )
    suppressMessages(gtrack.array.import("arr", "", arr_file))

    if (indexed) {
        for (track in c("d1", "d20", "d32", "d50", "const20", "sp", "arr")) {
            suppressMessages(gtrack.convert_to_indexed(track))
        }
        stopifnot(file.exists(file.path(db, "tracks", "d20.track", "track.idx")))
    }

    list(db = db, sp = sp, arr = arr)
}

# Per-bp values of a chromosome (NA where no value), from intervals and values
.wsum_bp_values <- function(len, starts, ends, vals) {
    bp <- rep(NA_real_, len)
    for (i in seq_along(starts)) {
        if (!is.na(vals[i])) {
            bp[(starts[i] + 1):ends[i]] <- vals[i]
        }
    }
    bp
}

# The raw bins of a dense track (or the items of a sparse track) as per-bp values, one vector per chromosome
.wsum_bp_from_track <- function(track) {
    raw <- gextract(track, gintervals.all(), iterator = track, colnames = "v")
    raw$chrom <- as.character(raw$chrom)
    lapply(setNames(names(.wsum_chrom_sizes), names(.wsum_chrom_sizes)), function(chrom) {
        x <- raw[raw$chrom == chrom, ]
        .wsum_bp_values(.wsum_chrom_sizes[[chrom]], x$start, x$end, x$v)
    })
}

# Reference weighted sum of the windows [start + sshift, end + eshift) clipped to the chromosome, optionally with
# a mask (data frame of chrom/start/end) whose bp do not count
.wsum_reference <- function(bp, iv, sshift = 0, eshift = 0, mask = NULL) {
    chroms <- as.character(iv$chrom)
    if (!is.null(mask)) {
        for (i in seq_len(nrow(mask))) {
            chrom <- as.character(mask$chrom[i])
            bp[[chrom]][(mask$start[i] + 1):mask$end[i]] <- NA
        }
    }
    vapply(seq_len(nrow(iv)), function(i) {
        chrom <- chroms[i]
        s <- max(0, iv$start[i] + sshift)
        e <- min(.wsum_chrom_sizes[[chrom]], iv$end[i] + eshift)
        if (s >= e) {
            return(NaN)
        }
        x <- bp[[chrom]][(s + 1):e]
        if (all(is.na(x))) NaN else sum(x, na.rm = TRUE)
    }, numeric(1))
}

.wsum_extract <- function(vtrack, iterator) {
    gextract(vtrack, gintervals.all(), iterator = iterator, colnames = "v")
}

.wsum_dbs <- list()
.wsum_db <- function(indexed = FALSE) {
    key <- if (indexed) "indexed" else "per_chrom"
    if (is.null(.wsum_dbs[[key]])) {
        .wsum_dbs[[key]] <<- .wsum_make_db(indexed)
    }
    suppressMessages(gsetroot(.wsum_dbs[[key]]$db))
    gvtrack.clear()
    .wsum_dbs[[key]]
}
withr::defer(for (info in .wsum_dbs) unlink(dirname(info$db), recursive = TRUE), teardown_env())

test_that("weighted.sum on a constant track is the window length; sum counts the touched bins", {
    .wsum_db()
    gvtrack.create("ws", "const20", func = "weighted.sum", sshift = -140, eshift = 140)
    gvtrack.create("s", "const20", func = "sum", sshift = -140, eshift = 140)
    res <- gextract(c("ws", "s"), gintervals.all(), iterator = 32)
    len <- unname(.wsum_chrom_sizes[as.character(res$chrom)])
    expected <- pmin(len, res$end + 140) - pmax(0, res$start - 140)
    expect_equal(res$ws, expected)
    inside <- res$start >= 140 & res$end + 140 <= len
    expect_true(all(res$ws[inside] == 312))
    # the 312 bp window touches 16 or 17 of the 20 bp bins, by its phase
    expect_setequal(unique(res$s[inside]), c(16, 17))
    expect_equal(res$ws[1], 32 + 140)
})

test_that("weighted.sum of dense tracks matches the raw bins at bin sizes 1, 20, 32 and 50", {
    for (indexed in c(FALSE, TRUE)) {
        .wsum_db(indexed)
        for (binsize in c(1, 20, 32, 50)) {
            track <- paste0("d", binsize)
            bp <- .wsum_bp_from_track(track)
            for (shifts in list(c(-140, 140), c(-7, 13), c(0, 0), c(25, 61))) {
                gvtrack.create("ws", track, func = "weighted.sum", sshift = shifts[1], eshift = shifts[2])
                for (iterator in c(32, 25)) {
                    res <- .wsum_extract("ws", iterator)
                    expect_identical(res$v, .wsum_reference(bp, res, shifts[1], shifts[2]),
                        info = sprintf("indexed %s, bin size %d, shifts %d %d, iterator %d", indexed, binsize, shifts[1], shifts[2], iterator)
                    )
                }
            }
            # windows inside a run of bins without a value are NaN
            res <- .wsum_extract("ws", 32)
            expect_true(any(is.na(res$v)))
        }
    }
})

test_that("weighted.sum over intervals that are not bin-aligned and of varying width matches the raw bins", {
    .wsum_db()
    bp <- .wsum_bp_from_track("d20")
    iv <- gintervals(c("chr1", "chr1", "chr1", "chr2", "chr2"), c(0, 3, 1017, 9000, 9970), c(1, 1500, 1019, 9991, 9991))
    gvtrack.create("ws", "d20", func = "weighted.sum")
    res <- gextract("ws", iv, iterator = iv, colnames = "v")
    expect_identical(res$v, .wsum_reference(bp, res))
})

test_that("weighted.sum over a bin-aligned window is sum x bin size", {
    .wsum_db()
    for (binsize in c(20, 50)) {
        track <- paste0("d", binsize)
        gvtrack.create("ws", track, func = "weighted.sum", sshift = -2 * binsize, eshift = 3 * binsize)
        gvtrack.create("s", track, func = "sum", sshift = -2 * binsize, eshift = 3 * binsize)
        res <- gextract(c("ws", "s"), gintervals.all(), iterator = 2 * binsize)
        # windows that end in the short last bin of a chromosome are not bin-aligned
        len <- .wsum_chrom_sizes[as.character(res$chrom)]
        aligned <- res$end + 3 * binsize <= floor(len / binsize) * binsize
        expect_true(sum(aligned & !is.na(res$s)) > 100)
        expect_identical(res$ws[aligned], res$s[aligned] * binsize)
    }
})

test_that("weighted.sum of a sparse track and of a value-based vtrack matches the intervals", {
    info <- .wsum_db()
    sp <- info$sp
    bp <- lapply(setNames(names(.wsum_chrom_sizes), names(.wsum_chrom_sizes)), function(chrom) {
        x <- sp[sp$chrom == chrom, ]
        .wsum_bp_values(.wsum_chrom_sizes[[chrom]], x$start, x$end, x$value)
    })
    for (shifts in list(c(-140, 140), c(0, 0), c(-3, 50))) {
        gvtrack.create("ws_sparse", "sp", func = "weighted.sum", sshift = shifts[1], eshift = shifts[2])
        gvtrack.create("ws_values", sp, func = "weighted.sum", sshift = shifts[1], eshift = shifts[2])
        for (iterator in c(32, 7)) {
            res <- gextract(c("ws_sparse", "ws_values"), gintervals.all(), iterator = iterator)
            ref <- .wsum_reference(bp, res, shifts[1], shifts[2])
            expect_identical(res$ws_sparse, ref)
            expect_identical(res$ws_values, ref)
        }
    }
})

test_that("weighted.sum of an array track weighs the slice value of each interval", {
    info <- .wsum_db()
    arr <- info$arr
    item_bp <- function(vals) {
        lapply(setNames(names(.wsum_chrom_sizes), names(.wsum_chrom_sizes)), function(chrom) {
            sel <- arr$chrom == chrom
            .wsum_bp_values(.wsum_chrom_sizes[[chrom]], arr$start[sel], arr$end[sel], vals[sel])
        })
    }
    cols <- as.matrix(arr[, c("col1", "col2", "col3")])

    gvtrack.create("ws", "arr", func = "weighted.sum", sshift = -140, eshift = 140)
    res <- .wsum_extract("ws", 32)
    expect_equal(res$v, .wsum_reference(item_bp(rowMeans(cols, na.rm = TRUE)), res, -140, 140))

    gvtrack.array.slice("ws", c("col1", "col3"), "sum")
    res <- .wsum_extract("ws", 32)
    expect_identical(res$v, .wsum_reference(item_bp(arr$col1 + arr$col3), res, -140, 140))
})

test_that("with a filter weighted.sum counts the unmasked bp of each bin or interval once", {
    info <- .wsum_db()
    mask <- gintervals(c("chr1", "chr1", "chr1", "chr2"), c(105, 333, 5000, 41), c(117, 600, 5013, 9000))
    sp <- info$sp
    sp_bp <- lapply(setNames(names(.wsum_chrom_sizes), names(.wsum_chrom_sizes)), function(chrom) {
        x <- sp[sp$chrom == chrom, ]
        .wsum_bp_values(.wsum_chrom_sizes[[chrom]], x$start, x$end, x$value)
    })
    for (track in c("d20", "d50", "sp")) {
        bp <- if (track == "sp") sp_bp else .wsum_bp_from_track(track)
        gvtrack.create("ws", track, func = "weighted.sum", sshift = -140, eshift = 140, filter = mask)
        res <- .wsum_extract("ws", 32)
        expect_identical(res$v, .wsum_reference(bp, res, -140, 140, mask), info = track)
    }
    gvtrack.create("ws", sp, func = "weighted.sum", sshift = -140, eshift = 140, filter = mask)
    res <- .wsum_extract("ws", 32)
    expect_identical(res$v, .wsum_reference(sp_bp, res, -140, 140, mask))
})

test_that("weighted.sum is the same with multitasking, on an indexed database and with a shared backend", {
    results <- list()
    for (indexed in c(FALSE, TRUE)) {
        for (shared in c(FALSE, TRUE)) {
            for (multitask in c(FALSE, TRUE)) {
                .wsum_db(indexed)
                withr::local_options(list(
                    gshared.track.backend = shared, gmultitasking = multitask,
                    gmax.processes = 4L, gmax.processes2core = 4L, gmin.scope4process = 1L
                ))
                # the same track under two shifts and next to other functions, so the track object is shared
                gvtrack.create("ws1", "d20", func = "weighted.sum", sshift = -140, eshift = 140)
                gvtrack.create("ws2", "d20", func = "weighted.sum", sshift = -60, eshift = 10)
                gvtrack.create("s1", "d20", func = "sum", sshift = -140, eshift = 140)
                gvtrack.create("a2", "d20", func = "avg", sshift = -60, eshift = 10)
                gvtrack.create("wsp", "sp", func = "weighted.sum", sshift = -140, eshift = 140)
                res <- gextract(c("ws1", "ws2", "s1", "a2", "wsp"), gintervals.all(), iterator = 32)
                res <- res[order(res$chrom, res$start), ]
                rownames(res) <- NULL
                results[[sprintf("%s/%s/%s", indexed, shared, multitask)]] <- res[, c("chrom", "start", "end", "ws1", "ws2", "s1", "a2", "wsp")]
            }
        }
    }
    for (key in names(results)[-1]) {
        expect_identical(results[[key]], results[[1]], info = key)
    }
    bp <- {
        .wsum_db()
        .wsum_bp_from_track("d20")
    }
    expect_identical(results[[1]]$ws1, .wsum_reference(bp, results[[1]], -140, 140))
    expect_identical(results[[1]]$ws2, .wsum_reference(bp, results[[1]], -60, 10))
})

test_that("weighted.sum does not change sum and accepts no parameters", {
    .wsum_db()
    gvtrack.create("s", "d20", func = "sum", sshift = -140, eshift = 140)
    alone <- .wsum_extract("s", 32)
    gvtrack.create("ws", "d20", func = "weighted.sum", sshift = -140, eshift = 140)
    both <- gextract(c("s", "ws"), gintervals.all(), iterator = 32)
    expect_identical(both$s, alone$v)
    expect_error(gvtrack.create("bad", "d20", func = "weighted.sum", params = 1))
})
