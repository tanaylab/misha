# On-disk safety of the three writers that could leave a corrupt or orphaned
# artefact behind (2026-08-23 lifetime-contract audit, section C):
#
#  1. ~GenomeTrackArrays() wrote the array-track footer only from the
#     destructor, and threw on a short write - std::terminate, so R never
#     regained control and the staging directory was orphaned.
#  2. ~StatQuadTreeCachedSerializer() had the same shape, and end() ignored
#     every write and seek return, so a full disk stamped a wrong 2D chunk
#     header over a truncated tree and the import reported success.
#  3. gtrack.modify overwrote the live track in place, so an interrupt left it
#     durably half-old/half-new under its real name.
#
# The out-of-disk halves cannot be provoked from the suite (they need a write
# to fail); they were verified by hand against a simulated full disk. What is
# pinned here is everything that can be: the staging contract, the sweep that
# removes what a killed process leaves behind, and the round trips through the
# two footer writers that are now called explicitly.

# md5 of a track's data files (everything but the dotfiles misha keeps
# attributes in, which a successful modify is supposed to update).
track_data_digest <- function(track) {
    dir <- file.path(.misha$GROOT, "tracks", paste0(gsub("\\.", "/", track), ".track"))
    files <- sort(list.files(dir, full.names = TRUE))
    vapply(files, function(f) tools::md5sum(f)[[1]], character(1))
}

staging_dirs <- function(track) {
    dir <- file.path(.misha$GROOT, "tracks", paste0(gsub("\\.", "/", track), ".track"))
    list.files(dirname(dir), all.files = TRUE, pattern = "\\.tmp\\.")
}

# A small genome with a 1bp-bin dense track. 4000 bins is four passes of the
# scanner's 1000-interval evaluation buffer, so an expression that fails on the
# third pass has already produced two buffers of new values.
ondisk_test_track <- function(env = parent.frame()) {
    withr::local_options(list(gmulticontig.indexed_format = FALSE), .local_envir = env)
    seq <- paste(rep("A", 4000), collapse = "")
    setup_db(list(sprintf(">chr1\n%s\n", seq)), return_db = TRUE)
    gtrack.create_sparse(
        "src", "source",
        gintervals("chr1", seq(0, 3900, 100), seq(100, 4000, 100)),
        seq_len(40) / 10
    )
    gtrack.create("dense1", "1bp dense", "src", iterator = 1)
}

test_that("gtrack.modify writes the new values and leaves no staging directory", {
    local_db_state()
    ondisk_test_track()

    before <- gextract("dense1", gintervals.all(), colnames = "v")
    before <- before[order(before$start), ]

    gtrack.modify("dense1", "dense1 * 2", gintervals("chr1", 1000, 2000))

    after <- gextract("dense1", gintervals.all(), colnames = "v")
    after <- after[order(after$start), ]
    touched <- after$start >= 1000 & after$start < 2000

    expect_equal(sum(touched), 1000)
    expect_identical(after$v[!touched], before$v[!touched])
    expect_equal(after$v[touched], before$v[touched] * 2)
    expect_length(staging_dirs("dense1"), 0)
})

test_that("gtrack.modify leaves the track byte-identical when the expression fails mid-scan", {
    local_db_state()
    ondisk_test_track()

    before <- track_data_digest("dense1")

    # The scanner evaluates the track expression one buffer of iterator
    # intervals at a time, so failing on the third call means two buffers of
    # replacement values have already reached the writer. Written in place -
    # as gtrack.modify used to be - that is a durably half-modified track.
    state <- new.env(parent = emptyenv())
    state$n <- 0
    assign(".misha_ondisk_boom", function(x) {
        state$n <- state$n + 1
        if (state$n >= 3) stop("boom")
        x * 2
    }, envir = globalenv())
    withr::defer(rm(".misha_ondisk_boom", envir = globalenv()))

    expect_error(
        gtrack.modify("dense1", ".misha_ondisk_boom(dense1)", gintervals.all()),
        "boom"
    )
    expect_gte(state$n, 3)
    expect_identical(track_data_digest("dense1"), before)
    expect_length(staging_dirs("dense1"), 0)
})

test_that("an array track whose data ends on the last chromosome round-trips", {
    # The footer of the last chromosome used to be written only by
    # ~GenomeTrackArrays(); it is now written by an explicit finish_writing().
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    seq <- paste(rep("A", 2000), collapse = "")
    setup_db(list(sprintf(">chr1\n%s\n", seq), sprintf(">chr2\n%s\n", seq)), return_db = TRUE)

    src <- withr::local_tempfile()
    n <- 100
    d <- data.frame(
        chrom = "chr2",
        start = as.integer(seq(0L, by = 20L, length.out = n)),
        end = as.integer(seq(0L, by = 20L, length.out = n)) + 20L,
        col1 = seq_len(n) / 10,
        col2 = seq_len(n) / 5
    )
    write.table(d, src, sep = "\t", quote = FALSE, row.names = FALSE)

    gtrack.array.import("arr", "last-chrom array", src)
    withr::defer(try(gtrack.rm("arr", TRUE), silent = TRUE))

    r <- gtrack.array.extract("arr", NULL, gintervals.all())
    r <- r[order(r$chrom, r$start), ]
    expect_equal(nrow(r), n)
    # track values are float32; compare at that precision
    expect_equal(r$col1, d$col1, tolerance = 1e-6)
    expect_equal(r$col2, d$col2, tolerance = 1e-6)
})

test_that("a multi-subtree 2D track round-trips", {
    # Exercises StatQuadTreeCachedSerializer::end()'s top-tree header, whose
    # writes and seeks are now checked.
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    seq <- paste(rep("A", 20000), collapse = "")
    setup_db(list(sprintf(">chr1\n%s\n", seq)), return_db = TRUE)

    k <- 20
    g <- expand.grid(i = 0:(k - 1), j = 0:(k - 1))
    src <- withr::local_tempfile()
    d <- data.frame(
        chrom1 = "chr1", start1 = as.integer(g$i * 1000L), end1 = as.integer(g$i * 1000L + 500L),
        chrom2 = "chr1", start2 = as.integer(g$j * 1000L), end2 = as.integer(g$j * 1000L + 500L),
        value = seq_len(nrow(g)) / 10
    )
    write.table(d, src, sep = "\t", quote = FALSE, row.names = FALSE)

    # force more than one subtree, i.e. the serializer's chunked path
    withr::local_options(list(gmax.data.size = 20))
    gtrack.2d.import("t2d", "multi subtree", src)
    withr::defer(try(gtrack.rm("t2d", TRUE), silent = TRUE))

    withr::local_options(list(gmax.data.size = 1e6))
    r <- gextract("t2d", gintervals.2d("chr1", 0, -1, "chr1", 0, -1))
    expect_equal(nrow(r), nrow(d))
    expect_equal(sum(r$t2d), sum(d$value), tolerance = 1e-6)
})
