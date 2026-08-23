# Regression tests for the uninitialised struct padding that the 2D quad-tree
# serializer used to write into track files.
#
# Rectangle_val<float> is 40 bytes long but holds only 36 bytes of members, and
# Point_val<float> is 24 bytes long for 20 bytes of members. StatQuadTreeCached
# writes those records to the track file verbatim, so the tail padding used to
# carry whatever happened to be in the writing process's memory. Consequences:
# a 2D track was not byte-reproducible (the same input produced files differing
# by 4 bytes per record), heap bytes leaked into files that get copied between
# users, and valgrind reported "write(buf) points to uninitialised byte(s)".
#
# What actually guards the fix: writing the same track twice *inside one R
# session* is NOT sufficient - both calls leaked the same stale bytes from the
# same stack slot, so the two files matched even on the broken build. The
# load-bearing assertion is therefore that the padding bytes are zero; on the
# broken build they are not (e.g. 2f 05 00 00). The two-writes comparison is
# kept as a cheap guard against nondeterminism from other sources.
#
# The record sizes below double as a pin on the on-disk layout: the fix must not
# change it, since existing tracks have to stay readable and new tracks have to
# stay readable by older misha.

# Layout of a 2D quad-tree track file that holds a single leaf (see
# StatQuadTreeCached::serialize / serialize_subtree).
QT_HEADER <- 4 + 8 + 8 # format signature, number of objects, root chunk position
QT_CHUNK_HEADER <- 8 + 8 # chunk size, offset of the chunk's top node
QT_LEAF <- 80 # is_leaf + padding + stat + arena + num_objs + padding
QT_RECT_OBJ <- 48 # id + x1,y1,x2,y2 + float value + 4 padding bytes
QT_POINT_OBJ <- 32 # id + x,y + float value + 4 padding bytes
QT_OBJS_START <- QT_HEADER + QT_CHUNK_HEADER + QT_LEAF
QT_PAD <- 4 # tail padding of every object record

# `off` is a 0-based offset into the raw vector `r`.
qt_u32 <- function(r, off) {
    v <- readBin(r[off + 1:4], "integer", 1, size = 4, endian = "little")
    if (v < 0) v + 2^32 else as.numeric(v)
}

qt_i64 <- function(r, off) {
    lo <- qt_u32(r, off)
    hi <- readBin(r[off + 5:8], "integer", 1, size = 4, endian = "little")
    hi * 2^32 + lo
}

qt_f32 <- function(r, off) readBin(r[off + 1:4], "numeric", 1, size = 4, endian = "little")

qt_pad <- function(r, off) as.integer(r[off + 1:QT_PAD])

read_track_pair_file <- function(track, pair = "chr1-chr2") {
    f <- file.path(.misha$GROOT, "tracks", paste0(gsub("\\.", "/", track), ".track"), pair)
    expect_true(file.exists(f))
    readBin(f, "raw", file.size(f))
}

pad_test_db <- function() {
    # Force the per-chromosome layout: an indexed database moves the same bytes
    # into a single .dat, which would only make them harder to address here.
    withr::local_options(list(gmulticontig.indexed_format = FALSE), .local_envir = parent.frame())
    seq <- paste(rep("A", 2000), collapse = "")
    setup_db(list(sprintf(">chr1\n%s\n", seq), sprintf(">chr2\n%s\n", seq)), return_db = TRUE)
}

test_that("rectangle records of a 2D track carry no uninitialised padding", {
    local_db_state()
    pad_test_db()

    x1 <- c(100, 300, 500, 700, 900)
    y1 <- c(150, 350, 550, 750, 950)
    vals <- c(1, 2, 3, 4, 5)
    intervs <- gintervals.2d("chr1", x1, x1 + 100, "chr2", y1, y1 + 100)
    gtrack.2d.create("pad_rects", "padding regression", intervs, vals)

    r <- read_track_pair_file("pad_rects")

    # A handful of rectangles fit in the root leaf, so the file layout is fully
    # determined. If this fails the on-disk format has changed.
    expect_equal(length(r), QT_OBJS_START + length(x1) * QT_RECT_OBJ)
    expect_equal(qt_i64(r, 4), as.numeric(length(x1))) # number of objects

    # The leaf header's own padding (after the 1-byte is_leaf, and after the
    # 4-byte object count) must be zero as well.
    expect_equal(as.integer(r[38:44]), rep(0L, 7))
    expect_equal(as.integer(r[113:116]), rep(0L, 4))

    # Decode the records so that the padding assertion below is anchored to real
    # record data rather than to arbitrary zeroed bytes.
    recs <- do.call(rbind, lapply(seq_along(x1), function(i) {
        off <- QT_OBJS_START + (i - 1) * QT_RECT_OBJ
        data.frame(
            x1 = qt_i64(r, off + 8), y1 = qt_i64(r, off + 16),
            x2 = qt_i64(r, off + 24), y2 = qt_i64(r, off + 32),
            v = qt_f32(r, off + 40)
        )
    }))
    recs <- recs[order(recs$x1), ]
    expect_equal(recs$x1, x1)
    expect_equal(recs$x2, x1 + 100)
    expect_equal(recs$y1, y1)
    expect_equal(recs$y2, y1 + 100)
    expect_equal(recs$v, vals)

    # The four bytes that follow the float value are structure padding. They are
    # written to the file, so they must be defined.
    for (i in seq_along(x1)) {
        off <- QT_OBJS_START + (i - 1) * QT_RECT_OBJ
        expect_equal(qt_pad(r, off + 44), rep(0L, QT_PAD),
            info = sprintf("padding of rectangle record %d", i)
        )
    }
})

test_that("point records of a 2D track carry no uninitialised padding", {
    local_db_state()
    pad_test_db()

    x <- c(100, 300, 500)
    y <- c(150, 350, 550)
    vals <- c(1, 2, 3)
    # All-1bp intervals make gtrack.2d.import produce a POINTS track, which
    # serializes Point_val<float> instead of Rectangle_val<float>.
    src <- tempfile(fileext = ".txt")
    withr::defer(unlink(src))
    utils::write.table(
        data.frame(
            chrom1 = "chr1", start1 = x, end1 = x + 1,
            chrom2 = "chr2", start2 = y, end2 = y + 1, value = vals
        ),
        src,
        sep = "\t", row.names = FALSE, quote = FALSE
    )
    gtrack.2d.import("pad_points", "padding regression", src)
    expect_equal(gtrack.info("pad_points")$type, "points")

    r <- read_track_pair_file("pad_points")
    expect_equal(length(r), QT_OBJS_START + length(x) * QT_POINT_OBJ)
    expect_equal(qt_i64(r, 4), as.numeric(length(x)))

    recs <- do.call(rbind, lapply(seq_along(x), function(i) {
        off <- QT_OBJS_START + (i - 1) * QT_POINT_OBJ
        data.frame(x = qt_i64(r, off + 8), y = qt_i64(r, off + 16), v = qt_f32(r, off + 24))
    }))
    recs <- recs[order(recs$x), ]
    expect_equal(recs$x, x)
    expect_equal(recs$y, y)
    expect_equal(recs$v, vals)

    for (i in seq_along(x)) {
        off <- QT_OBJS_START + (i - 1) * QT_POINT_OBJ
        expect_equal(qt_pad(r, off + 28), rep(0L, QT_PAD),
            info = sprintf("padding of point record %d", i)
        )
    }
})

test_that("writing the same 2D track twice produces byte-identical files", {
    local_db_state()
    pad_test_db()

    set.seed(17)
    x1 <- (seq_len(500) - 1) * 4
    y1 <- sample(seq(0, 1900, by = 4), 500, replace = TRUE)
    intervs <- gintervals.2d("chr1", x1, x1 + 2, "chr2", y1, y1 + 2)
    vals <- seq_len(500) / 3

    gtrack.2d.create("pad_twice1", "padding regression", intervs, vals)
    gtrack.2d.create("pad_twice2", "padding regression", intervs, vals)

    expect_identical(
        read_track_pair_file("pad_twice1"),
        read_track_pair_file("pad_twice2")
    )

    # ... and the tree that spans several leaves must be free of padding leftovers
    # too: every byte of the file is either a decoded field or a defined pad byte,
    # so a whole-file scan for the tell-tale stale values is not possible. Instead
    # check the values survive the round trip unchanged.
    extracted <- gextract("pad_twice1", gintervals.2d.all(), colnames = "v")
    extracted <- extracted[order(extracted$start1), ]
    expect_equal(nrow(extracted), length(x1))
    expect_equal(extracted$v, vals, tolerance = 1e-6)
})
