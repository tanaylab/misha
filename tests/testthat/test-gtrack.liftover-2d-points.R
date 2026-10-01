# gtrack.liftover of 2D POINTS tracks (Hi-C contacts): the target stays a POINTS track,
# each point is lifted by both ends as gintervals.liftover lifts a 1bp interval, a point
# with an unmapped end is dropped, and points landing on the same target point are
# merged by multi_target_agg.

# Creates a DB with the given chromosomes and opens it. gdb.create builds an indexed DB only
# from a single multi-FASTA file; from one FASTA per chromosome it builds a per-chromosome DB
# whatever 'format' says.
mk_points_db <- function(chroms, size, format = "per-chromosome") {
    d <- tempfile("lift2dpts_")
    dir.create(d)
    seqs <- sprintf(">%s\n%s\n", chroms, paste(rep("A", size), collapse = ""))
    if (format == "indexed") {
        fas <- file.path(d, "genome.fasta")
        cat(seqs, sep = "", file = fas)
    } else {
        fas <- file.path(d, paste0(chroms, ".fasta")) # chrom name derives from the file name
        for (i in seq_along(chroms)) cat(seqs[i], file = fas[i])
    }
    db <- tempfile("lift2dpts_db_")
    suppressMessages(gdb.create(groot = db, fasta = fas, format = format))
    withr::defer(
        {
            unlink(db, recursive = TRUE)
            unlink(d, recursive = TRUE)
        },
        envir = testthat::teardown_env()
    )
    gdb.init(db)
    expect_equal(.gdb.is_indexed(), format == "indexed")
    db
}

# Imports contacts (1bp ends) as a POINTS track in the current DB; cis pairs are mirrored
# and trans pairs stored in both orientations, as for any Hi-C track.
import_points <- function(track, contacts) {
    f <- tempfile(fileext = ".tsv")
    withr::defer(unlink(f))
    write.table(data.frame(
        chrom1 = contacts$chrom1, start1 = contacts$start1, end1 = contacts$start1 + 1,
        chrom2 = contacts$chrom2, start2 = contacts$start2, end2 = contacts$start2 + 1, value = contacts$v
    ), f, sep = "\t", quote = FALSE, row.names = FALSE)
    suppressMessages(gtrack.2d.import_contacts(track, "x", f, fends = NULL))
}

extract_points <- function(track) {
    x <- gextract(track, gintervals.2d.all(), colnames = "v")
    if (is.null(x)) {
        return(data.frame(chrom1 = character(0), start1 = numeric(0), chrom2 = character(0), start2 = numeric(0), v = numeric(0)))
    }
    expect_true(all(x$end1 - x$start1 == 1 & x$end2 - x$start2 == 1))
    x <- data.frame(chrom1 = as.character(x$chrom1), start1 = x$start1, chrom2 = as.character(x$chrom2), start2 = x$start2, v = x$v)
    x <- x[do.call(order, x), , drop = FALSE]
    rownames(x) <- NULL
    x
}

# every point has its mirror image (y, x) with the same value
expect_mirrored <- function(pts) {
    mirror <- data.frame(chrom1 = pts$chrom2, start1 = pts$start2, chrom2 = pts$chrom1, start2 = pts$start1, v = pts$v)
    mirror <- mirror[do.call(order, mirror), , drop = FALSE]
    rownames(mirror) <- NULL
    expect_equal(mirror, pts)
}

# source chrS1 (1000bp): [0,300) -> chrT1 [100,400) +, [300,400) unmapped,
# [400,700) -> chrT1 [1500,1800) -, [700,1000) -> chrT2 [0,300) +;
# chrS2 (1000bp): [0,1000) -> chrT2 [1000,2000) +
write_points_chain <- function() {
    chain <- new_chain_file()
    # the "-" block is given on the target's reverse strand: [1500,1800) is [200,500) there
    write_chain_entry(chain, "chrS1", 1000, "+", 0, 300, "chrT1", 2000, "+", 100, 400, 1)
    write_chain_entry(chain, "chrS1", 1000, "+", 400, 700, "chrT1", 2000, "-", 200, 500, 2)
    write_chain_entry(chain, "chrS1", 1000, "+", 700, 1000, "chrT2", 2000, "+", 0, 300, 3)
    write_chain_entry(chain, "chrS2", 1000, "+", 0, 1000, "chrT2", 2000, "+", 1000, 2000, 4)
    chain
}

test_that("gtrack.liftover of a POINTS track gives a POINTS track lifted end by end", {
    local_db_state()

    src_db <- mk_points_db(c("chrS1", "chrS2"), 1000)
    set.seed(17)
    n <- 300
    contacts <- data.frame(
        chrom1 = sample(c("chrS1", "chrS2"), n, replace = TRUE, prob = c(3, 1)),
        start1 = sample(0:999, n, replace = TRUE),
        chrom2 = sample(c("chrS1", "chrS2"), n, replace = TRUE, prob = c(3, 1)),
        start2 = sample(0:999, n, replace = TRUE),
        v = sample(1:9, n, replace = TRUE)
    )
    import_points("src", contacts)
    expect_equal(gtrack.info("src")$type, "points")
    src <- extract_points("src")
    src_dir <- file.path(src_db, "tracks", "src.track")

    # an indexed target: gtrack.liftover converts the lifted track to indexed at the end
    tgt_db <- mk_points_db(c("chrT1", "chrT2"), 2000, format = "indexed")
    chain <- write_points_chain()
    gtrack.liftover("lifted", "x", src_dir, chain, multi_target_agg = "sum")
    expect_equal(gtrack.info("lifted")$type, "points")
    expect_true(file.exists(file.path(tgt_db, "tracks", "lifted.track", "track.idx")))
    got <- extract_points("lifted")

    # expected: both ends lifted as 1bp intervals; the chain is 1:1, so nothing collides
    lift_pos <- function(chrom, start) {
        l <- gintervals.liftover(data.frame(chrom = chrom, start = start, end = start + 1), chain)
        out <- data.frame(chrom = rep(NA_character_, length(chrom)), start = NA_real_)
        out$chrom[l$intervalID] <- as.character(l$chrom)
        out$start[l$intervalID] <- l$start
        out
    }
    e1 <- lift_pos(src$chrom1, src$start1)
    e2 <- lift_pos(src$chrom2, src$start2)
    exp <- data.frame(chrom1 = e1$chrom, start1 = e1$start, chrom2 = e2$chrom, start2 = e2$start, v = src$v)
    exp <- exp[!is.na(exp$chrom1) & !is.na(exp$chrom2), ]
    exp <- exp[do.call(order, exp), ]
    rownames(exp) <- NULL
    expect_equal(got, exp)

    # the data exercises every case: an unmapped end, the inverted block, cis -> trans
    in_gap <- function(chrom, start) chrom == "chrS1" & start >= 300 & start < 400
    expect_true(any(in_gap(src$chrom1, src$start1) | in_gap(src$chrom2, src$start2)))
    expect_equal(nrow(got), sum(!(in_gap(src$chrom1, src$start1) | in_gap(src$chrom2, src$start2))))
    expect_true(any(got$chrom1 == "chrT1" & got$start1 >= 1500))
    expect_true(any(got$chrom1 == "chrT1" & got$chrom2 == "chrT2" & got$start2 < 300))

    # mirrored cis source -> mirrored target, including the cis pairs that became trans
    expect_mirrored(got)
})

test_that("gtrack.liftover of POINTS: explicit points on each kind of block", {
    local_db_state()

    src_db <- mk_points_db(c("chrS1", "chrS2"), 1000)
    import_points("src", data.frame(
        chrom1 = c("chrS1", "chrS1", "chrS1", "chrS1", "chrS1", "chrS1"),
        start1 = c(10, 450, 50, 20, 350, 55),
        chrom2 = c("chrS1", "chrS1", "chrS1", "chrS2", "chrS1", "chrS1"),
        start2 = c(20, 460, 750, 5, 30, 55),
        v = c(1, 2, 3, 4, 5, 6)
    ))
    src_dir <- file.path(src_db, "tracks", "src.track")

    tgt_db <- mk_points_db(c("chrT1", "chrT2"), 2000)
    gtrack.liftover("lifted", "x", src_dir, write_points_chain())

    # (10, 20) -> (110, 120) and its mirror; (450, 460) is in the inverted block,
    # x -> 1500 + (699 - x); (50, 750) has 750 in the block to chrT2, so the cis point
    # becomes the trans pair chrT1-chrT2 (150, 50) plus chrT2-chrT1 (50, 150); the trans
    # (chrS1 20, chrS2 5) gives both orientations; (350, 30) has an end in the unmapped
    # [300, 400) and is dropped; the diagonal (55, 55) stays a single point.
    exp <- data.frame(
        chrom1 = c("chrT1", "chrT1", "chrT1", "chrT1", "chrT1", "chrT1", "chrT1", "chrT2", "chrT2"),
        start1 = c(110, 120, 120, 150, 155, 1739, 1749, 50, 1005),
        chrom2 = c("chrT1", "chrT1", "chrT2", "chrT2", "chrT1", "chrT1", "chrT1", "chrT1", "chrT1"),
        start2 = c(120, 110, 1005, 50, 155, 1749, 1739, 150, 120),
        v = c(1, 1, 4, 3, 6, 2, 2, 3, 4)
    )
    expect_equal(extract_points("lifted"), exp)
})

test_that("gtrack.liftover of POINTS merges points landing on one target point by multi_target_agg", {
    local_db_state()

    src_db <- mk_points_db("chrS1", 1000)
    import_points("src", data.frame(
        chrom1 = "chrS1", start1 = c(10, 510, 30, 45, 545), chrom2 = "chrS1", start2 = c(20, 520, 40, 45, 545),
        v = c(2, 5, 7, 1, 4)
    ))
    src_dir <- file.path(src_db, "tracks", "src.track")

    tgt_db <- mk_points_db("chrT1", 1000)
    # two source blocks onto the same target block (kept by tgt_overlap_policy = "keep"):
    # (10, 20) and (510, 520) both land on (10, 20); diagonal (45, 45) and (545, 545) on (45, 45)
    chain <- new_chain_file()
    write_chain_entry(chain, "chrS1", 1000, "+", 0, 100, "chrT1", 1000, "+", 0, 100, 1)
    write_chain_entry(chain, "chrS1", 1000, "+", 500, 600, "chrT1", 1000, "+", 0, 100, 2)

    lift <- function(agg, ...) {
        if (gtrack.exists("lifted")) gtrack.rm("lifted", force = TRUE)
        gtrack.liftover("lifted", "x", src_dir, chain, tgt_overlap_policy = "keep", multi_target_agg = agg, ...)
        expect_equal(gtrack.info("lifted")$type, "points")
        extract_points("lifted")
    }
    value_at <- function(got, x, y) got$v[got$start1 == x & got$start2 == y]

    expected <- list(sum = c(7, 7, 5), mean = c(3.5, 7, 2.5), count = c(2, 1, 2), max = c(5, 7, 4), min = c(2, 7, 1))
    for (agg in names(expected)) {
        got <- lift(agg)
        expect_mirrored(got)
        expect_equal(nrow(got), 5, info = agg) # (10,20), (20,10), (30,40), (40,30), (45,45)
        expect_equal(c(value_at(got, 10, 20), value_at(got, 30, 40), value_at(got, 45, 45)), expected[[agg]], info = agg)
    }
    # first and last pick different contributors
    first <- value_at(lift("first"), 10, 20)
    last <- value_at(lift("last"), 10, 20)
    expect_setequal(c(first, last), c(2, 5))
    got <- lift("nth", params = 2)
    expect_equal(value_at(got, 10, 20), last)
    expect_true(is.na(value_at(got, 30, 40)))
    got <- lift("mean", min_n = 2)
    expect_equal(value_at(got, 10, 20), 3.5)
    expect_true(is.na(value_at(got, 30, 40)))
})

test_that("gtrack.liftover of POINTS turns an infinite value into NaN, as for 1D and rects", {
    local_db_state()

    src_db <- mk_points_db("chrS1", 1000)
    import_points("src", data.frame(
        chrom1 = "chrS1", start1 = c(10, 30, 510), chrom2 = "chrS1", start2 = c(20, 40, 520),
        v = c(Inf, -Inf, 3)
    ))
    expect_equal(extract_points("src")$v, c(Inf, Inf, -Inf, -Inf, 3, 3))
    src_dir <- file.path(src_db, "tracks", "src.track")

    tgt_db <- mk_points_db("chrT1", 1000)
    # (10, 20) and (510, 520) land on (10, 20)
    chain <- new_chain_file()
    write_chain_entry(chain, "chrS1", 1000, "+", 0, 100, "chrT1", 1000, "+", 0, 100, 1)
    write_chain_entry(chain, "chrS1", 1000, "+", 500, 600, "chrT1", 1000, "+", 0, 100, 2)
    gtrack.liftover("lifted", "x", src_dir, chain, tgt_overlap_policy = "keep", multi_target_agg = "sum")
    got <- extract_points("lifted")
    expect_equal(got$start1, c(10, 20, 30, 40))
    # the NaN from Inf is dropped by na.rm; -Inf alone gives NaN
    expect_equal(got$v, c(3, 3, NaN, NaN))
})

test_that("gtrack.liftover of POINTS gives the same track when the pair is split into subtrees", {
    local_db_state()

    src_db <- mk_points_db(c("chrS1", "chrS2"), 1000)
    set.seed(3)
    n <- 500
    import_points("src", data.frame(
        chrom1 = sample(c("chrS1", "chrS2"), n, replace = TRUE), start1 = sample(0:999, n, replace = TRUE),
        chrom2 = sample(c("chrS1", "chrS2"), n, replace = TRUE), start2 = sample(0:999, n, replace = TRUE),
        v = runif(n)
    ))
    src_dir <- file.path(src_db, "tracks", "src.track")

    tgt_db <- mk_points_db(c("chrT1", "chrT2"), 2000)
    chain <- write_points_chain()
    gtrack.liftover("whole", "x", src_dir, chain)
    withr::with_options(list(gmax.data.size = 10), gtrack.liftover("split", "x", src_dir, chain))
    expect_equal(gtrack.info("split")$type, "points")
    whole <- extract_points("whole")
    expect_gt(nrow(whole), 400)
    expect_equal(extract_points("split"), whole)
    # no buffer or subtree files left behind
    expect_setequal(list.files(file.path(tgt_db, "tracks", "split.track")), list.files(file.path(tgt_db, "tracks", "whole.track")))
})

test_that("gtrack.liftover of an indexed 2D track reads each chromosome pair by the source genome's ids", {
    local_db_state()

    # the same track as per-pair files and in an indexed DB, where the import converts it to
    # indexed; chrS1 has id 0 and chrS2 id 1 in both genomes
    contacts <- data.frame(chrom1 = c("chrS1", "chrS2"), start1 = c(10, 40), chrom2 = c("chrS1", "chrS2"), start2 = c(100, 300), v = c(1, 4))
    src_db <- mk_points_db(c("chrS1", "chrS2"), 1000)
    import_points("src", contacts)
    expect_false(file.exists(file.path(src_db, "tracks", "src.track", "track.idx")))
    src_idx_db <- mk_points_db(c("chrS1", "chrS2"), 1000, format = "indexed")
    import_points("src_idx", contacts)
    expect_true(file.exists(file.path(src_idx_db, "tracks", "src_idx.track", "track.idx")))

    tgt_db <- mk_points_db("chrT2", 2000)
    # the chain covers chrS2 only, so its chrom ids differ from the source genome's
    chain <- new_chain_file()
    write_chain_entry(chain, "chrS2", 1000, "+", 0, 1000, "chrT2", 2000, "+", 0, 1000, 1)
    gtrack.liftover("from_files", "x", file.path(src_db, "tracks", "src.track"), chain)
    gtrack.liftover("from_idx", "x", file.path(src_idx_db, "tracks", "src_idx.track"), chain)
    exp <- extract_points("from_files")
    expect_equal(exp$start1, c(40, 300))
    expect_equal(extract_points("from_idx"), exp)
})
