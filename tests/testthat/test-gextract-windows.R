# gseq.onehot and gextract.windows, on the example database (it ships with the package: its
# sequence has lowercase runs and N runs, and pymisha ships the same one, see
# tests/test_data/windows in pymisha).

# Re-roots into a database of its own; hand this process's overlay back to the next
# file in the parallel worker (helper-test_db.R).
restore_groot_on_exit()

gdb.init_examples(withr::local_tempdir(.local_envir = teardown_env()))

# the one-hot encoding of strings, in R: x[n, c, p] = 1 if toupper(base p of string n) is
# channel c
onehot_ref <- function(seqs) {
    n <- length(seqs)
    w <- nchar(seqs[1])
    res <- array(0L, dim = c(n, 4, w))
    for (i in seq_len(n)) {
        ch <- match(strsplit(toupper(seqs[i]), "")[[1]], c("A", "C", "G", "T"))
        for (p in which(!is.na(ch))) {
            res[i, ch[p], p] <- 1L
        }
    }
    res
}

revcomp_str <- function(x) {
    vapply(x, function(s) {
        paste(rev(strsplit(chartr("ACGTacgt", "TGCAtgca", s), "")[[1]]), collapse = "")
    }, character(1), USE.NAMES = FALSE)
}

windows_iv <- function(num = 60, width = 300, seed = 3) {
    set.seed(seed)
    sizes <- gintervals.all()
    chroms <- sample(as.character(sizes$chrom), num, replace = TRUE)
    lens <- sizes$end[match(chroms, sizes$chrom)]
    starts <- floor(runif(num) * (lens - width))
    data.frame(chrom = chroms, start = starts, end = starts + width, stringsAsFactors = FALSE)
}

test_that("gseq.onehot of intervals equals the encoding of gseq.extract, lowercase and N included", {
    iv <- windows_iv()
    seqs <- gseq.extract(iv)
    # the windows do hold lowercase bases and N
    expect_true(any(grepl("[acgt]", seqs)))
    expect_true(any(grepl("N", seqs)))

    x <- gseq.onehot(iv)
    expect_identical(typeof(x), "raw")
    expect_identical(dim(x), c(nrow(iv), 4L, 300L))
    expect_identical(dimnames(x)[[2]], c("A", "C", "G", "T"))
    ref <- onehot_ref(seqs)
    expect_identical(array(as.integer(x), dim(x)), ref)
    # N is all zero, every other position has exactly one 1
    sums <- apply(array(as.integer(x), dim(x)), c(1, 3), sum)
    expect_identical(sums == 0, matrix(strsplit(paste(seqs, collapse = ""), "")[[1]] == "N", nrow = nrow(iv), byrow = TRUE))

    expect_identical(gseq.onehot(iv, dtype = "integer"), {
        r <- ref
        dimnames(r) <- list(NULL, c("A", "C", "G", "T"), NULL)
        r
    })
    xd <- gseq.onehot(iv, dtype = "double")
    expect_identical(typeof(xd), "double")
    expect_equal(array(as.numeric(xd), dim(xd)), array(as.numeric(ref), dim(ref)))
})

test_that("gseq.onehot reverse complements strand -1 rows", {
    iv <- windows_iv(40, 101)
    iv$strand <- rep(c(1, -1), 20)
    seqs <- gseq.extract(iv[, 1:3])
    x <- gseq.onehot(iv)
    seqs[iv$strand == -1] <- revcomp_str(seqs[iv$strand == -1])
    expect_identical(array(as.integer(x), dim(x)), onehot_ref(seqs))
    # = gseq.extract, which reverse complements as well
    expect_identical(gseq.onehot(gseq.extract(iv)), x)
    # character strand
    iv$strand <- ifelse(iv$strand == 1, "+", "-")
    expect_identical(gseq.onehot(iv), x)
})

test_that("gseq.onehot of strings: lowercase, N, IUPAC codes and gaps", {
    seqs <- c("ACGTNRYacgtn-.", "tttTTTgggGGa*a")
    x <- gseq.onehot(seqs, dtype = "integer")
    expect_identical(unname(x), onehot_ref(seqs))
    expect_error(gseq.onehot(c("ACGT", "ACG")), "same length")
    expect_error(gseq.onehot(NA_character_), "NA")
})

test_that("gseq.onehot pad = TRUE gives zeros outside the chromosome, as gseq.extract(pad = TRUE) gives N", {
    len1 <- gintervals.all()$end[1]
    iv <- data.frame(
        chrom = "chr1", start = c(-50, len1 - 70, -10, len1 + 5),
        end = c(50, len1 + 30, 90, len1 + 105), strand = c(1, -1, -1, 1)
    )
    expect_error(gseq.onehot(iv), "runs past")
    x <- gseq.onehot(iv, pad = TRUE)
    expect_identical(array(as.integer(x), dim(x)), onehot_ref(gseq.extract(iv, pad = TRUE)))
})

test_that("gseq.onehot checks its input", {
    expect_error(gseq.onehot(gintervals(1, c(0, 100), c(10, 120))), "same positive width")
    expect_error(gseq.onehot(gintervals(1, 0, 10), dtype = "float"))
    expect_error(gseq.onehot(data.frame(chrom = "chrZ", start = 0, end = 10)))
})

test_that("gextract.windows centres the windows and matches gseq.onehot and gextract.matrix", {
    iv <- windows_iv(50, 1, seed = 11)
    iv$end <- iv$start + sample(1:400, nrow(iv), replace = TRUE)
    iv$name <- paste0("r", seq_len(nrow(iv)))
    w <- gextract.windows(iv, tracks = c("dense_track", sp = "sparse_track"), width = 301, signal_width = 200, binsize = 20)

    kept <- w$intervals$index
    expect_true(length(kept) > 30)
    mids <- iv$start[kept] + (iv$end[kept] - iv$start[kept]) %/% 2
    expect_identical(w$intervals$start, mids - 150)
    expect_identical(w$intervals$end, mids + 151)
    expect_identical(w$intervals$signal_start, mids - 100)
    expect_identical(w$intervals$signal_end, mids + 100)
    expect_identical(w$intervals$name, iv$name[kept])
    expect_identical(sort(c(kept, w$dropped$index)), seq_len(nrow(iv)))

    expect_identical(w$x, gseq.onehot(w$intervals[, c("chrom", "start", "end")]))
    m <- gextract.matrix("dense_track", "sparse_track",
        intervals = data.frame(chrom = w$intervals$chrom, start = w$intervals$signal_start, end = w$intervals$signal_end),
        iterator = 20
    )
    expect_identical(dim(w$y), c(length(kept), 2L, 10L))
    expect_identical(dimnames(w$y)[[2]], c("dense_track", "sp"))
    expect_identical(unname(w$y[, 1, ]), unname(m[, , 1]))
    expect_identical(unname(w$y[, 2, ]), unname(m[, , 2]))

    # sequence only
    w2 <- gextract.windows(iv, width = 301)
    expect_null(w2$y)
    expect_identical(w2$x, gseq.onehot(w2$intervals[, 1:3]))
    expect_null(gextract.windows(iv, tracks = "dense_track", width = 301, sequence = FALSE)$x)
})

test_that("gextract.windows drops windows off the chromosome, and jitter widens both windows", {
    len1 <- gintervals.all()$end[1]
    iv <- data.frame(chrom = "chr1", start = c(100, 1000, len1 - 1000, len1 - 140, 5000), end = c(101, 1001, len1 - 999, len1 - 139, 5001))
    w <- gextract.windows(iv, tracks = "dense_track", width = 200, signal_width = 100, jitter = 50)
    expect_identical(w$intervals$index, c(2L, 3L, 5L))
    expect_identical(w$dropped, data.frame(index = c(1L, 4L), reason = "off_chrom", stringsAsFactors = FALSE))
    expect_identical(dim(w$x), c(3L, 4L, 300L))
    expect_identical(dim(w$y), c(3L, 1L, 200L))
    expect_identical(w$intervals$start, iv$start[c(2, 3, 5)] - 150)
    expect_identical(w$intervals$signal_start, iv$start[c(2, 3, 5)] - 100)

    # a signal window wider than the sequence window decides by itself
    w <- gextract.windows(iv, tracks = "dense_track", width = 10, signal_width = 300)
    expect_identical(w$dropped$index, c(1L, 4L))
    expect_identical(gextract.windows(iv, width = 10)$dropped$index, integer(0))
})

test_that("gextract.windows reverse complements strand -1 rows and reverses their signals", {
    iv <- windows_iv(30, 1, seed = 5)
    iv$strand <- rep(c(-1, 1, 0), 10)
    w <- gextract.windows(iv, tracks = "dense_track", width = 100, binsize = 10)
    wf <- gextract.windows(iv, tracks = "dense_track", width = 100, binsize = 10, use_strand = FALSE)
    neg <- w$intervals$strand == -1
    expect_true(any(neg))
    expect_identical(w$x[!neg, , , drop = FALSE], wf$x[!neg, , , drop = FALSE])
    expect_identical(unname(w$x[neg, , , drop = FALSE]), unname(wf$x[neg, 4:1, 100:1, drop = FALSE]))
    expect_identical(w$y[neg, , , drop = FALSE], wf$y[neg, , 10:1, drop = FALSE])
    expect_identical(w$y[!neg, , , drop = FALSE], wf$y[!neg, , , drop = FALSE])
    seqs <- gseq.extract(w$intervals[, c("chrom", "start", "end", "strand")])
    expect_identical(array(as.integer(w$x), dim(w$x)), onehot_ref(seqs))
})

test_that("gextract.windows drops windows that overlap the exclusion intervals", {
    iv <- data.frame(chrom = "chr1", start = c(1000, 2000, 3000, 4000), end = c(1001, 2001, 3001, 4001))
    # windows [950, 1050), [1950, 2050), [2950, 3050), [3950, 4050)
    excl <- data.frame(
        chrom = c("chr1", "chr1", "chr1", "chrNotThere", "chr1"),
        start = c(1050, 2049, 2900, 3960, 3990),
        end = c(1100, 2060, 2950, 3970, 3990)
    )
    w <- gextract.windows(iv, width = 100, exclude = excl)
    # [1050, 1100) touches but does not overlap; [2049, 2060) overlaps by one base; [2900, 2950)
    # ends where the window starts; the empty [3990, 3990) and the unknown chromosome count for
    # nothing
    expect_identical(w$dropped, data.frame(index = 2L, reason = "excluded", stringsAsFactors = FALSE))
    # a long exclusion interval that starts well before the window
    excl2 <- data.frame(chrom = "chr1", start = c(10, 500), end = c(3500, 600))
    expect_identical(gextract.windows(iv, width = 100, exclude = excl2)$dropped$index, 1:3)
    # off_chrom wins over excluded
    w <- gextract.windows(rbind(data.frame(chrom = "chr1", start = 10, end = 11), iv), width = 100, exclude = data.frame(chrom = "chr1", start = 0, end = 100))
    expect_identical(w$dropped$reason, "off_chrom")
})

test_that("gextract.windows count filters", {
    iv <- windows_iv(40, 1, seed = 8)
    w <- gextract.windows(iv, tracks = c(d = "dense_track", s = "sparse_track"), width = 200, binsize = 10)
    counts <- rowSums(w$y[, "s", ], na.rm = TRUE)
    lo <- quantile(counts, 0.25)
    hi <- quantile(counts, 0.75)
    wf <- gextract.windows(iv, tracks = c(d = "dense_track", s = "sparse_track"), width = 200, binsize = 10, min_counts = lo, max_counts = hi, count_track = "s")
    keep <- counts >= lo & counts <= hi
    expect_identical(wf$intervals$index, w$intervals$index[keep])
    expect_identical(wf$x, w$x[keep, , , drop = FALSE])
    expect_identical(wf$y, w$y[keep, , , drop = FALSE])
    dropped <- wf$dropped[wf$dropped$reason != "off_chrom", ]
    expect_identical(dropped$index, w$intervals$index[!keep])
    expect_identical(dropped$reason, ifelse(counts[!keep] < lo, "min_counts", "max_counts"))
    # by index
    expect_identical(gextract.windows(iv, tracks = c(d = "dense_track", s = "sparse_track"), width = 200, binsize = 10, min_counts = lo, count_track = 2)$intervals$index, w$intervals$index[counts >= lo])
    expect_error(gextract.windows(iv, width = 200, min_counts = 1), "tracks")
    expect_error(gextract.windows(iv, tracks = "dense_track", width = 200, min_counts = 1, count_track = 2), "count_track")
})

test_that("gextract.windows labels the sets of a list and keeps their order", {
    pos <- data.frame(chrom = "chr1", start = c(1000, 2000), end = c(1001, 2001), score = c(1, 2))
    neg <- data.frame(chrom = "chr2", start = c(5000, 10, 7000), end = c(5001, 11, 7001), score = c(0, 0, 0))
    w <- gextract.windows(list(pos = pos, neg = neg), width = 100)
    expect_identical(w$intervals$set, c("pos", "pos", "neg", "neg"))
    expect_identical(w$intervals$index, c(1L, 2L, 3L, 5L))
    expect_identical(w$intervals$score, c(1, 2, 0, 0))
    expect_identical(w$dropped$index, 4L)
    expect_identical(gextract.windows(list(pos, neg), width = 100)$intervals$set, c(1L, 1L, 2L, 2L))
})

test_that("gextract.windows checks its arguments and handles no windows", {
    iv <- gintervals(1, 1000, 1001)
    expect_error(gextract.windows(iv), "Usage")
    expect_error(gextract.windows(iv, width = 0), "width")
    expect_error(gextract.windows(iv, tracks = "dense_track", width = 100, binsize = 30), "binsize")
    w <- gextract.windows(gintervals(1, 10, 11), tracks = "dense_track", width = 100)
    expect_identical(dim(w$x), c(0L, 4L, 100L))
    expect_identical(dim(w$y), c(0L, 1L, 100L))
    expect_identical(nrow(w$intervals), 0L)
})

test_that("gextract.windows works with virtual tracks and track expressions", {
    gvtrack.create("dsum", "dense_track", func = "sum")
    withr::defer(gvtrack.rm("dsum"))
    iv <- windows_iv(20, 1, seed = 2)
    w <- gextract.windows(iv, tracks = c("dsum", "dense_track * 2"), width = 100, binsize = 25)
    m <- gextract.matrix("dsum", "dense_track * 2", intervals = w$intervals[, 1:3], iterator = 25)
    expect_identical(unname(w$y[, 1, ]), unname(m[, , 1]))
    expect_identical(unname(w$y[, 2, ]), unname(m[, , 2]))
})

test_that("gseq.onehot leaves the lifetime counters balanced, on errors too", {
    x <- gseq.onehot(gintervals(1, 0, 10))
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    expect_error(gseq.onehot(data.frame(chrom = "chr1", start = -5, end = 5)), "runs past")
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    expect_error(gseq.onehot(c("ACGT", "AC")), "same length")
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    expect_error(.gcall("C_gseq_onehot", "chrNope", 0, 10, FALSE, FALSE, 0L, .misha_env()), "does not exist")
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
})
