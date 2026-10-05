# A per-pair 2D track whose pair files are named by an alias of the chromosomes ("1-2" or
# "chr1-2" for chr1 and chr2) is read as if its files had the canonical names, as the 1D
# per-chromosome files of a track are (GenomeTrack::find_existing_1d_filename). Older pymisha
# named pair files by chrom_sizes.txt names in per-chromosome databases.

restore_groot_on_exit()

test_that("a 2D track whose pair files are named by chromosome aliases reads as with canonical names", {
    local_db_state()
    td <- tempfile("pair_names_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))

    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    gsetroot(db)

    objs <- data.frame(
        chrom1 = c("chr1", "chr2", "chr1", "chrX"), start1 = c(10, 20, 50, 5), end1 = c(11, 21, 51, 6),
        chrom2 = c("chr1", "chr2", "chr2", "chr1_KI270706v1_random"), start2 = c(30, 40, 60, 7), end2 = c(31, 41, 61, 8)
    )
    rects <- objs
    rects$end1 <- rects$end1 + 4
    rects$end2 <- rects$end2 + 4
    gtrack.2d.create("rects", "x", rects, c(1, 2, 3, 4))
    contacts <- file.path(td, "contacts.tsv")
    write.table(cbind(objs, value = c(1, 2, 3, 4)), contacts, sep = "\t", quote = FALSE, row.names = FALSE)
    suppressMessages(gtrack.2d.import_contacts("pts", "x", contacts, fends = NULL))
    expect_equal(gtrack.info("rects")$type, "rectangles")
    expect_equal(gtrack.info("pts")$type, "points")

    scope <- gintervals.2d.all()
    read_all <- function(track) {
        list(
            # the track's own objects
            objs = gextract(track, scope),
            # evaluated over another iterator
            iter = gextract(track, scope, iterator = scope),
            # loaded as an intervals set, per chromosome pair
            load = gintervals.load(track, chrom1 = "chr1", chrom2 = "chr2"),
            # used as an intervals set
            as_intervals = gextract(track, track)
        )
    }

    # "chr1-chr2" becomes "1-2", "1-chr2" or "chr1-2" in turn
    alias_name <- function(pair_file, i) {
        chroms <- strsplit(pair_file, "-", fixed = TRUE)[[1]]
        strip <- list(c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, TRUE))[[(i - 1) %% 3 + 1]]
        chroms[strip] <- sub("^chr", "", chroms[strip])
        paste(chroms, collapse = "-")
    }
    for (track in c("rects", "pts")) {
        expected <- suppressMessages(read_all(track))
        expect_gte(nrow(expected$objs), 4)

        track_dir <- file.path(db, "tracks", paste0(track, ".track"))
        pair_files <- grep("^chr[^-]*-chr", list.files(track_dir), value = TRUE)
        expect_gte(length(pair_files), 4)
        aliases <- mapply(alias_name, pair_files, seq_along(pair_files))
        expect_true(all(file.rename(file.path(track_dir, pair_files), file.path(track_dir, aliases))))
        # stats cached by gintervals.load() from the canonical files
        unlink(file.path(track_dir, ".meta"))

        expect_equal(suppressMessages(read_all(track)), expected, info = track)
    }
})
