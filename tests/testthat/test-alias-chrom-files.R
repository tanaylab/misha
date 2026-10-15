# Per-chromosome files of tracks and big intervals sets that are named by an alias of their
# chromosomes ("1" for chr1, "MT" for chrM; "1-2", "1-chr2" or "chr1-2" for chr1 and chr2), as earlier
# pymisha named them in databases whose chrom_sizes.txt names lack "chr". One listing of the
# directory (GenomeTrack::ChromFiles) maps them to their chromosomes for every reader and converter,
# so each path below reads the same as with the canonical names.

restore_groot_on_exit()

# A per-chromosome database as create_db_with_unsorted_chrom_sizes(), with a chrM, which has an alias
# other than its name without "chr" (MT)
create_alias_test_db <- function(path) {
    chroms <- c("2", "10", "1", "X", "M")
    sizes <- c(2000, 1500, 1000, 1200, 500)
    dir.create(file.path(path, "tracks"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(path, "seq"), showWarnings = FALSE)
    for (i in seq_along(chroms)) {
        writeBin(charToRaw(strrep("A", sizes[i])), file.path(path, "seq", paste0("chr", chroms[i], ".seq")))
    }
    writeLines(paste(chroms, sizes, sep = "\t"), file.path(path, "chrom_sizes.txt"))
    invisible(normalizePath(path))
}

# "chr1" becomes "1" and "chrM" "MT"; "chr1-chr2" becomes "1-2", "1-chr2" or "chr1-2" in turn
alias_names <- function(files) {
    vapply(seq_along(files), function(i) {
        chroms <- strsplit(files[i], "-", fixed = TRUE)[[1]]
        strip <- if (length(chroms) == 1) TRUE else list(c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, TRUE))[[(i - 1) %% 3 + 1]]
        chroms[strip] <- ifelse(chroms[strip] == "chrM", "MT", sub("^chr", "", chroms[strip]))
        paste(chroms, collapse = "-")
    }, "")
}

alias_test_objects <- c(
    "sp.track", "dn.track", "arr.track", "rects.track", "pts.track", "bigset.interv", "bigset2d.interv"
)

# In the database at db: a sparse, a dense and an array track, a rectangles and a points track, and a
# 1D and a 2D per-chromosome big intervals set, all with files on chrM
create_alias_test_objects <- function(db) {
    gsetroot(db)
    intervs <- gintervals(c("chr1", "chr1", "chr2", "chrX", "chrM"), c(0, 100, 0, 10, 5), c(50, 150, 100, 20, 15))
    intervs2d <- gintervals.2d(
        c("chr1", "chr2", "chr1", "chrM"), c(0, 0, 10, 0), c(50, 100, 20, 30),
        c("chr1", "chr2", "chrX", "chr1"), c(0, 0, 10, 40), c(50, 100, 20, 60)
    )
    gtrack.create_sparse("sp", "x", intervs, seq_len(nrow(intervs)))
    gtrack.create("dn", "x", "sp", iterator = 10)
    arr_file <- tempfile("arr_", fileext = ".tsv")
    write.table(cbind(intervs, a = 1:5, b = 6:10), arr_file, sep = "\t", quote = FALSE, row.names = FALSE)
    suppressMessages(gtrack.array.import("arr", "x", arr_file))
    gtrack.2d.create("rects", "x", intervs2d, seq_len(nrow(intervs2d)))
    contacts <- tempfile("contacts_", fileext = ".tsv")
    pts <- intervs2d
    pts$end1 <- pts$start1 + 1
    pts$end2 <- pts$start2 + 1
    write.table(cbind(pts, value = seq_len(nrow(pts))), contacts, sep = "\t", quote = FALSE, row.names = FALSE)
    suppressMessages(gtrack.2d.import_contacts("pts", "x", contacts, fends = NULL))
    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("bigset", intervs)
        gintervals.save("bigset2d", intervs2d)
    })
    unlink(c(arr_file, contacts))
}

# Renames the per-chromosome files of the objects to aliases and removes the tracks' .meta (the stats
# of a track used as an intervals set, which gintervals.load() saved from the canonical files)
rename_alias_test_objects <- function(db) {
    for (obj in alias_test_objects) {
        dir <- file.path(db, "tracks", obj)
        files <- grep("^chr", list.files(dir), value = TRUE)
        expect_gte(length(files), 3)
        expect_true(all(file.rename(file.path(dir, files), file.path(dir, alias_names(files)))))
        if (endsWith(obj, ".track")) {
            unlink(file.path(dir, ".meta"))
        }
    }
}

# Every path that reads the objects
read_alias_test_objects <- function() {
    scope <- gintervals.all()
    scope2d <- gintervals.2d.all()
    gvtrack.create("v_dn", "dn", "avg")
    gvtrack.create("v_sp1", "sp", "sum")
    gvtrack.iterator("v_sp1", dim = 1)
    withr::defer({
        gvtrack.rm("v_dn")
        gvtrack.rm("v_sp1")
    })
    suppressMessages(list(
        # 1D tracks: by their own iterator and by others, in expressions and virtual tracks
        sp = gextract("sp", scope),
        sp_bins = gextract("sp", scope, iterator = 50),
        dn = gextract("dn", scope),
        dn_vtrack = gextract("v_dn * 2", scope, iterator = 100),
        arr = gextract("arr", scope),
        arr_cols = gtrack.array.extract("arr", "b", scope),
        screen = gscreen("sp > 1", scope),
        projected = gextract("v_sp1", scope2d, iterator = "rects"),
        info = lapply(c("sp", "dn", "arr", "rects", "pts"), function(track) gtrack.info(track)[c("type", "dimensions")]),
        # a 1D track used as an intervals set
        sp_load = gintervals.load("sp"),
        sp_as_scope = gextract("dn", "sp"),
        # 2D tracks
        rects = gextract("rects", scope2d),
        rects_bins = gextract("rects", scope2d, iterator = scope2d),
        pts = gextract("pts", scope2d),
        pts_bins = gextract("pts", scope2d, iterator = scope2d),
        # 2D tracks used as intervals sets
        rects_load = gintervals.load("rects"),
        pts_load = gintervals.load("pts", chrom1 = "chrM", chrom2 = "chr1"),
        rects_as_scope = gextract("rects", "rects"),
        pts_as_scope = gextract("pts", "pts"),
        # big intervals sets
        bigset = gintervals.load("bigset"),
        bigset_chrom = gintervals.load("bigset", chrom = "chrM"),
        bigset_as_scope = gextract("sp", "bigset"),
        bigset_as_iterator = gextract("sp", scope, iterator = "bigset"),
        bigset2d = gintervals.load("bigset2d"),
        bigset2d_pair = gintervals.load("bigset2d", chrom1 = "chrM", chrom2 = "chr1"),
        bigset2d_as_scope = gextract("rects", "bigset2d"),
        bigset2d_as_iterator = gextract("rects", scope2d, iterator = "bigset2d")
    ))
}

test_that("tracks and intervals sets whose files are named by chromosome aliases read as with canonical names", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_read_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    create_alias_test_objects(db)

    expected <- read_alias_test_objects()
    expect_equal(nrow(expected$sp), 5)
    expect_equal(nrow(expected$rects), 4)
    expect_equal(nrow(expected$bigset2d), 4)
    expect_equal(nrow(expected$bigset_chrom), 1)

    rename_alias_test_objects(db)
    expect_true("MT" %in% list.files(file.path(db, "tracks", "sp.track")))
    expect_equal(read_alias_test_objects(), expected)
})

test_that("tracks and intervals sets whose files are named by chromosome aliases convert to the indexed format", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_convert_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    create_alias_test_objects(db)
    expected <- read_alias_test_objects()
    rename_alias_test_objects(db)

    suppressMessages({
        for (track in c("sp", "dn", "arr", "rects", "pts")) {
            gtrack.convert_to_indexed(track)
        }
        gintervals.convert_to_indexed("bigset", remove.old = TRUE)
        gintervals.2d.convert_to_indexed("bigset2d", remove.old = TRUE)
    })
    for (obj in alias_test_objects) {
        files <- list.files(file.path(db, "tracks", obj))
        expect_true(any(c("track.idx", "intervals.idx", "intervals2d.idx") %in% files), info = obj)
        # no per-chromosome file is left
        expect_equal(files[!grepl("^(track|intervals|intervals2d)\\.(idx|dat)$", files)], character(0), info = obj)
    }
    expect_equal(read_alias_test_objects(), expected)
})

test_that("tracks and intervals sets whose files are named by chromosome aliases are copied, modified and updated", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_write_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    dest_perchrom <- create_alias_test_db(file.path(td, "dest_perchrom"))
    dest_indexed <- normalizePath(create_test_db(file.path(td, "dest_indexed"), chrom_sizes = data.frame(
        chrom = c("chr2", "chr10", "chr1", "chrX", "chrM"),
        size = c(2000, 1500, 1000, 1200, 500)
    )))
    gdb.init(dest_indexed)
    gdb.convert_to_indexed(force = TRUE, verbose = FALSE)

    create_alias_test_objects(db)
    expected <- read_alias_test_objects()
    expected_arr <- gtrack.array.extract("arr", NULL, gintervals.all())
    rename_alias_test_objects(db)

    # tracks written from alias-named ones have files with canonical names
    gtrack.create("sp_copy", "x", "sp", iterator = "sp")
    suppressMessages(gtrack.array.import("arr_copy", "x", "arr"))
    expect_equal(gextract("sp_copy", gintervals.all())$sp_copy, expected$sp$sp)
    expect_equal(gtrack.array.extract("arr_copy", NULL, gintervals.all()), expected_arr)
    for (track in c("sp_copy", "arr_copy")) {
        files <- list.files(file.path(db, "tracks", paste0(track, ".track")))
        expect_true(all(c("chr1", "chrM") %in% files))
        expect_true(all(files %in% c(as.character(gintervals.all()$chrom), "vars")))
    }

    # gtrack.modify() changes the alias-named file, as earlier pymisha, which reads it too, does
    gtrack.modify("dn", "dn * 2", gintervals.all())
    expect_equal(gextract("dn", gintervals.all())$dn, expected$dn$dn * 2)
    dn_files <- list.files(file.path(db, "tracks", "dn.track"))
    expect_true("MT" %in% dn_files)
    expect_false("chrM" %in% dn_files)

    # gtrack.copy() into a per-chromosome and into an indexed database (which orders the chromosomes
    # differently)
    by_chrom <- function(res) {
        res <- res[order(as.character(res$chrom), res$start), ]
        paste(res$chrom, res$start, res$end, res$sp)
    }
    for (dest in c(dest_perchrom, dest_indexed)) {
        gsetroot(db)
        gtrack.copy("sp", "sp", db = dest)
        gsetroot(dest)
        expect_equal(by_chrom(gextract("sp", gintervals.all())), by_chrom(expected$sp), info = dest)
    }
    gsetroot(db)

    # gintervals.update() of a chromosome leaves no file named by its alias behind
    withr::with_options(list(gmax.data.size = 1), {
        gintervals.update("bigset", NULL, chrom = "chrM")
        gintervals.update("bigset", gintervals("chr1", 300, 400), chrom = "chr1")
    })
    files <- list.files(file.path(db, "tracks", "bigset.interv"))
    expect_false(any(c("MT", "1") %in% files))
    expect_true("chr1" %in% files)
    res <- gintervals.load("bigset")
    expect_equal(paste(res$chrom, res$start, res$end), c("chr1 300 400", "chr2 0 100", "chrX 10 20"))
})

test_that("a canonically named file wins over a file named by an alias, also when converted", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_canonical_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    gtrack.create_sparse("sp", "x", gintervals("chr1", 0, 100), 1)
    gtrack.create_sparse("other", "x", gintervals("chr1", 200, 300), 2)
    rects <- gintervals.2d("chr1", 10, 15, "chr2", 30, 35)
    gtrack.2d.create("rects", "x", rects, 1)
    gtrack.2d.create("other2d", "x", gintervals.2d("chr1", 100, 105, "chr2", 300, 305), 2)
    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("set2d", gintervals.2d(c("chr1", "chr2"), 10, 15, "chr2", 30, 35))
        gintervals.save("other2dset", gintervals.2d(c("chr1", "chr2"), 100, 105, "chr2", 300, 305))
    })
    expected <- list(gextract("sp", gintervals.all()), gextract("rects", gintervals.2d.all()), gintervals.load("set2d"))

    # a file named by an alias next to the canonical one
    tracks <- file.path(db, "tracks")
    expect_true(file.copy(file.path(tracks, "other.track", "chr1"), file.path(tracks, "sp.track", "1")))
    expect_true(file.copy(file.path(tracks, "other2d.track", "chr1-chr2"), file.path(tracks, "rects.track", "1-2")))
    expect_true(file.copy(file.path(tracks, "other2dset.interv", "chr1-chr2"), file.path(tracks, "set2d.interv", "1-chr2")))
    unlink(file.path(tracks, "rects.track", ".meta"))
    read <- function() list(gextract("sp", gintervals.all()), gextract("rects", gintervals.2d.all()), gintervals.load("set2d"))
    expect_equal(read(), expected)

    # converted, the canonical file is packed and remove.old removes both
    suppressMessages({
        gtrack.convert_to_indexed("sp")
        gtrack.2d.convert_to_indexed("rects", remove.old = TRUE)
        gintervals.2d.convert_to_indexed("set2d", remove.old = TRUE)
    })
    expect_false(any(c("chr1", "1") %in% list.files(file.path(tracks, "sp.track"))))
    expect_false(any(c("chr1-chr2", "1-2") %in% list.files(file.path(tracks, "rects.track"))))
    expect_false(any(c("chr1-chr2", "1-chr2") %in% list.files(file.path(tracks, "set2d.interv"))))
    expect_equal(read(), expected)
})

test_that("without a canonical name, the \"chr\"-toggled one wins, and two other alias-named files are an error", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_twice_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    gtrack.create_sparse("sp", "x", gintervals(c("chr1", "chrM"), 0, 100), c(1, 2))
    gtrack.create_sparse("other", "x", gintervals("chrM", 200, 300), 3)
    gtrack.2d.create("rects", "x", gintervals.2d("chr1", 10, 15, "chr2", 30, 35), 1)
    gtrack.2d.create("other2d", "x", gintervals.2d("chr1", 100, 105, "chr2", 300, 305), 2)
    withr::with_options(list(gmax.data.size = 1), gintervals.save("set", gintervals(c("chr1", "chrM"), 0, 100)))
    expected <- gextract("sp", gintervals.all())
    expected2d <- gextract("rects", gintervals.2d.all())

    # chrM's chrom_sizes.txt name is "M": "M" wins over the other aliases ("MT", "chrMT"), as earlier
    # pymisha writes it and master looked it up
    tracks <- file.path(db, "tracks")
    sp_dir <- file.path(tracks, "sp.track")
    expect_true(file.rename(file.path(sp_dir, "chrM"), file.path(sp_dir, "M")))
    expect_true(file.copy(file.path(tracks, "other.track", "chrM"), file.path(sp_dir, "MT")))
    expect_equal(gextract("sp", gintervals.all()), expected)
    # without it, two other aliases
    expect_true(file.rename(file.path(sp_dir, "M"), file.path(sp_dir, "chrMT")))
    expect_error(gextract("sp", gintervals.all()), "two files of chromosome chrM, MT and chrMT, and none named chrM or M")
    expect_error(gextract("other", "sp"), "two files of chromosome chrM")
    expect_error(suppressMessages(gtrack.convert_to_indexed("sp")), "two files of chromosome chrM")
    # two names of one file are not two files
    unlink(file.path(sp_dir, "MT"))
    expect_true(file.symlink(file.path(sp_dir, "chrMT"), file.path(sp_dir, "MT")))
    expect_equal(gextract("sp", gintervals.all()), expected)

    rects_dir <- file.path(tracks, "rects.track")
    expect_true(file.rename(file.path(rects_dir, "chr1-chr2"), file.path(rects_dir, "1-2")))
    expect_true(file.copy(file.path(tracks, "other2d.track", "chr1-chr2"), file.path(rects_dir, "chr1-2")))
    unlink(file.path(rects_dir, ".meta"))
    expect_equal(gextract("rects", gintervals.2d.all()), expected2d)
    expect_true(file.rename(file.path(rects_dir, "1-2"), file.path(rects_dir, "1-chr2")))
    expect_error(gextract("rects", gintervals.2d.all()), "two files of chromosome pair \\(chr1, chr2\\), 1-chr2 and chr1-2")
    expect_error(gintervals.load("rects", chrom1 = "chr1", chrom2 = "chr2"), "two files of chromosome pair")
    expect_error(suppressMessages(gtrack.2d.convert_to_indexed("rects")), "two files of chromosome pair")

    set_dir <- file.path(tracks, "set.interv")
    expect_true(file.rename(file.path(set_dir, "chrM"), file.path(set_dir, "MT")))
    expect_true(file.copy(file.path(set_dir, "chr1"), file.path(set_dir, "chrMT")))
    expect_error(gintervals.load("set"), "two files of chromosome chrM")
    expect_error(suppressMessages(gintervals.convert_to_indexed("set")), "two files of chromosome chrM")
})

test_that("a conversion stops at a file whose name is not a chromosome's", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_unparsed_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    gtrack.create_sparse("sp", "x", gintervals("chr1", 0, 100), 1)
    gtrack.2d.create("rects", "x", gintervals.2d("chr1", 10, 15, "chr2", 30, 35), 1)
    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("set", gintervals(c("chr1", "chr2"), 0, 100))
        gintervals.save("set2d", gintervals.2d(c("chr1", "chr2"), 0, 100, "chr2", 0, 100))
    })
    tracks <- file.path(db, "tracks")
    for (file in c("sp.track/chr1.bak", "rects.track/chr1-chr7", "set.interv/notes", "set2d.interv/notes")) {
        writeLines("notes", file.path(tracks, file))
    }

    expect_error(suppressMessages(gtrack.convert_to_indexed("sp")), "file chr1.bak is not named by a chromosome")
    expect_error(suppressMessages(gtrack.2d.convert_to_indexed("rects")), "file chr1-chr7 is not named by a chromosome pair")
    expect_error(suppressMessages(gintervals.convert_to_indexed("set")), "file notes is not named by a chromosome")
    expect_error(suppressMessages(gintervals.2d.convert_to_indexed("set2d")), "file notes is not named by a chromosome pair")
    for (obj in c("sp.track", "rects.track", "set.interv", "set2d.interv")) {
        expect_false(any(grepl("\\.(idx|dat)", list.files(file.path(tracks, obj)))), info = obj)
    }
})

test_that("dangling symlinks are ignored, and a symlink counts only when named by a chromosome", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_links_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    dest <- create_alias_test_db(file.path(td, "dest"))
    gsetroot(db)

    gtrack.create_sparse("sp", "x", gintervals(c("chr1", "chr2"), 0, 100), c(1, 2))
    gtrack.2d.create("rects", "x", gintervals.2d(c("chr1", "chr2"), 10, 15, "chr2", 30, 35), c(1, 2))
    withr::with_options(list(gmax.data.size = 1), gintervals.save("set2d", gintervals.2d(c("chr1", "chr2"), 0, 100, "chr2", 0, 100)))
    expected <- list(gextract("sp", gintervals.all()), gextract("rects", gintervals.2d.all()), gintervals.load("set2d"))
    read <- function() list(gextract("sp", gintervals.all()), gextract("rects", gintervals.2d.all()), gintervals.load("set2d"))

    # each canonical file moves away and is read through a symlink named by an alias; a dangling
    # symlink keeps its name, and a link to a notes file is not named by a chromosome
    writeLines("Notes on this track", file.path(td, "notes.txt"))
    tracks <- file.path(db, "tracks")
    for (file in c("sp.track/chr1", "rects.track/chr1-chr2", "set2d.interv/chr1-chr2")) {
        real <- file.path(td, gsub("/", "_", file))
        expect_true(file.rename(file.path(tracks, file), real))
        expect_true(file.symlink(real, file.path(tracks, sub("chr1", "1", file))))
        expect_true(file.symlink(file.path(td, "missing"), file.path(tracks, file)))
        expect_true(file.symlink(file.path(td, "notes.txt"), file.path(tracks, dirname(file), "README")))
    }
    unlink(file.path(tracks, "rects.track", ".meta"))
    gdb.reload()
    expect_equal(gtrack.info("sp")$type, "sparse")
    expect_equal(gtrack.info("rects")$type, "rectangles")
    expect_equal(read(), expected)

    unlink(file.path(tracks, c("sp.track", "rects.track", "set2d.interv"), "README"))
    # a copy stops at a dangling symlink (its target may come back), as on master
    expect_error(suppressWarnings(gtrack.copy("rects", "rects", db = dest)), "Failed to copy")
    unlink(file.path(tracks, c("sp.track/chr1", "rects.track/chr1-chr2", "set2d.interv/chr1-chr2")))
    gtrack.copy("sp", "sp", db = dest)
    gtrack.copy("rects", "rects", db = dest)
    suppressMessages({
        gtrack.convert_to_indexed("sp")
        gtrack.2d.convert_to_indexed("rects", remove.old = TRUE)
        gintervals.2d.convert_to_indexed("set2d", remove.old = TRUE)
    })
    expect_equal(read(), expected)
    gsetroot(dest)
    expect_equal(gextract("sp", gintervals.all()), expected[[1]])
    expect_equal(gextract("rects", gintervals.2d.all()), expected[[2]])
})

test_that("a conversion stops at a file it cannot read, and removes nothing", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_unreadable_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    gtrack.create_sparse("sp", "x", gintervals(c("chr1", "chr2"), 0, 100), c(1, 2))
    gtrack.2d.create("rects", "x", gintervals.2d(c("chr1", "chr2"), 10, 15, "chr2", 30, 35), c(1, 2))
    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("set", gintervals(c("chr1", "chr2"), 0, 100))
        gintervals.save("set2d", gintervals.2d(c("chr1", "chr2"), 0, 100, "chr2", 0, 100))
    })
    tracks <- file.path(db, "tracks")
    files <- file.path(tracks, c("sp.track/chr2", "rects.track/chr2-chr2", "set.interv/chr2", "set2d.interv/chr2-chr2"))
    Sys.chmod(files, "000")
    skip_if(any(file.access(files, 4) == 0), "the files stay readable to this user")

    expect_error(suppressMessages(gtrack.convert_to_indexed("sp")), "Failed to open")
    expect_error(suppressMessages(gtrack.2d.convert_to_indexed("rects", remove.old = TRUE)), "Failed to open")
    expect_error(suppressMessages(gintervals.convert_to_indexed("set", remove.old = TRUE)), "Failed to open")
    expect_error(suppressMessages(gintervals.2d.convert_to_indexed("set2d", remove.old = TRUE)), "Failed to open")
    expect_true(all(file.exists(files)))
    for (obj in c("sp.track", "rects.track", "set.interv", "set2d.interv")) {
        expect_false(any(grepl("\\.(idx|dat)", list.files(file.path(tracks, obj)))), info = obj)
    }
})

test_that("a forced conversion of an indexed 2D track keeps its index when the per-pair files fail the checks", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_force2d_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    gtrack.2d.create("rects", "x", gintervals.2d(c("chr1", "chr2"), 10, 15, c("chr1", "chrX"), 30, 35), c(1, 2))
    suppressMessages(gtrack.2d.convert_to_indexed("rects", remove.old = TRUE))
    expected <- gextract("rects", gintervals.2d.all())
    dir <- file.path(db, "tracks", "rects.track")
    writeLines("notes", file.path(dir, "README"))

    expect_error(
        suppressMessages(gtrack.2d.convert_to_indexed("rects", remove.old = TRUE, force = TRUE)),
        "file README is not named by a chromosome pair"
    )
    expect_true(all(file.exists(file.path(dir, c("track.idx", "track.dat")))))
    expect_equal(gextract("rects", gintervals.2d.all()), expected)

    # a pair file it cannot open stops it later, after the new files are begun: the old ones come back too
    unlink(file.path(dir, "README"))
    writeLines("x", file.path(dir, "chr1-chr2"))
    Sys.chmod(file.path(dir, "chr1-chr2"), "000")
    if (file.access(file.path(dir, "chr1-chr2"), 4) != 0) {
        expect_error(suppressMessages(gtrack.2d.convert_to_indexed("rects", remove.old = TRUE, force = TRUE)), "Failed to open")
        expect_true(all(file.exists(file.path(dir, c("track.idx", "track.dat")))))
        expect_equal(gextract("rects", gintervals.2d.all()), expected)
    }

    # once it succeeds nothing is left aside
    unlink(file.path(dir, "chr1-chr2"))
    suppressMessages(gtrack.2d.convert_to_indexed("rects", remove.old = TRUE, force = TRUE))
    expect_equal(gextract("rects", gintervals.2d.all()), expected)
    expect_equal(sort(list.files(dir, all.files = TRUE, no.. = TRUE)), c(".attributes", "track.dat", "track.idx"))
})

test_that("files named by the \"chr\" toggle convert in a database of over 1000 contigs without \"chr\"", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_toggle_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    # a per-chromosome database names its files chr<contig>; the indexed one names its chromosomes as
    # chrom_sizes.txt does, and with over 1000 such contigs its aliases leave out the toggle (.compute_chrom_aliases)
    chroms <- paste0("s", 1:1001)
    pc <- file.path(td, "pc")
    dir.create(file.path(pc, "tracks"), recursive = TRUE)
    dir.create(file.path(pc, "seq"))
    for (chrom in chroms) writeBin(charToRaw(strrep("A", 100)), file.path(pc, "seq", paste0("chr", chrom, ".seq")))
    writeLines(paste(chroms, 100, sep = "\t"), file.path(pc, "chrom_sizes.txt"))
    fasta <- file.path(td, "genome.fa")
    writeLines(rbind(paste0(">", chroms), strrep("A", 100)), fasta)
    ix <- file.path(td, "ix")
    suppressMessages(gdb.create(ix, fasta, format = "indexed"))

    gsetroot(pc)
    intervs <- gintervals(c("chrs1", "chrs2", "chrs500"), c(0, 0, 10), c(50, 60, 20))
    gtrack.create_sparse("sp", "x", intervs, 1:3)
    withr::with_options(list(gmax.data.size = 1), gintervals.save("set", data.frame(intervs, v = 1:3)))
    expect_true(file.exists(file.path(pc, "tracks", "sp.track", "chrs500")))
    expect_true(all(file.copy(file.path(pc, "tracks", c("sp.track", "set.interv")), file.path(ix, "tracks"), recursive = TRUE)))

    gsetroot(ix)
    intervs <- gintervals(c("s1", "s2", "s500"), c(0, 0, 10), c(50, 60, 20))
    suppressMessages({
        gtrack.convert_to_indexed("sp")
        gintervals.convert_to_indexed("set", remove.old = TRUE)
    })
    expect_equal(gextract("sp", intervs, colnames = "v")$v, 1:3)
    expect_equal(gintervals.load("set")$v, 1:3)
})

test_that("loading one chromosome (pair) of a big set lists its directory only when the canonical file is missing", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_load_one_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("set", gintervals(c("chr1", "chr2"), 0, 100))
        gintervals.save("set2d", gintervals.2d(c("chr1", "chr2"), 0, 100, "chr2", 0, 100))
    })
    listed <- 0
    gchrom_files <- .gchrom_files
    local_mocked_bindings(.gchrom_files = function(...) {
        listed <<- listed + 1
        gchrom_files(...)
    })
    expect_equal(nrow(gintervals.load("set", chrom = "chr2")), 1)
    expect_equal(nrow(gintervals.load("set2d", chrom1 = "chr1", chrom2 = "chr2")), 1)
    expect_equal(listed, 0)

    tracks <- file.path(db, "tracks")
    file.rename(file.path(tracks, "set.interv", "chr2"), file.path(tracks, "set.interv", "2"))
    file.rename(file.path(tracks, "set2d.interv", "chr1-chr2"), file.path(tracks, "set2d.interv", "1-2"))
    expect_equal(nrow(gintervals.load("set", chrom = "chr2")), 1)
    expect_equal(nrow(gintervals.load("set2d", chrom1 = "chr1", chrom2 = "chr2")), 1)
    expect_equal(listed, 2)
})

test_that("next to an intervals set's index a file named by an alias is ignored, as older than the index", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_stale_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("old", data.frame(gintervals(c("chr1", "chr2"), 0, 10), v = 1))
        gintervals.save("s", data.frame(gintervals(c("chr1", "chr2"), c(500, 0), c(600, 10)), v = c(99, 2)))
        gintervals.save("old2d", data.frame(gintervals.2d(c("chr1", "chr2"), 0, 10, "chr2", 0, 10), v = 1))
        gintervals.save("s2d", data.frame(gintervals.2d(c("chr1", "chr2"), c(500, 0), c(600, 10), "chr2", c(500, 0), c(600, 10)), v = c(99, 2)))
    })
    suppressMessages({
        gintervals.convert_to_indexed("s", remove.old = TRUE)
        gintervals.2d.convert_to_indexed("s2d", remove.old = TRUE)
    })
    expected <- list(gintervals.load("s"), gintervals.load("s2d"))
    # as an earlier conversion left them: the canonical file packed, an older one named by an alias kept
    tracks <- file.path(db, "tracks")
    expect_true(file.copy(file.path(tracks, "old.interv", "chr1"), file.path(tracks, "s.interv", "1")))
    expect_true(file.copy(file.path(tracks, "old2d.interv", "chr1-chr2"), file.path(tracks, "s2d.interv", "1-2")))

    expect_equal(list(gintervals.load("s"), gintervals.load("s2d")), expected)
    expect_equal(gintervals.load("s", chrom = "chr1")$v, 99)
    expect_equal(gintervals.load("s2d", chrom1 = "chr1", chrom2 = "chr2")$v, 99)
    gtrack.create_sparse("sp", "x", gintervals(c("chr1", "chr2"), 0, 1000), 1:2)
    expect_equal(gextract("sp", intervals = "s", iterator = "s")$start, c(500, 0))
})

test_that("a 2D copy that is converted stops before the destination at a file not named by a pair", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_copy2d_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    p <- create_alias_test_db(file.path(td, "p"))
    fasta <- file.path(td, "genome.fa")
    writeLines(c(">chr1", strrep("A", 1000), ">chr2", strrep("A", 2000)), fasta)
    w <- file.path(td, "w")
    suppressMessages(gdb.create(w, fasta, format = "indexed"))

    gsetroot(p)
    intervs <- gintervals.2d(c("chr1", "chr2"), 10, 15, "chr2", 30, 35)
    gtrack.2d.create("src", "x", intervs, c(1, 2))
    gtrack.2d.create("dst", "x", intervs, c(7, 8))
    expect_true(all(file.copy(file.path(p, "tracks", c("src.track", "dst.track")), file.path(w, "tracks"), recursive = TRUE)))
    writeLines("notes", file.path(w, "tracks", "src.track", "README"))
    gsetroot(w)
    gdb.reload()
    suppressMessages(gtrack.2d.convert_to_indexed("dst", remove.old = TRUE))
    expected <- gextract("dst", gintervals.2d.all())

    expect_error(gtrack.copy("src", "dst", overwrite = TRUE), "file README is not named by a chromosome pair")
    expect_equal(gextract("dst", gintervals.2d.all()), expected)
    # the copy makes a file of a symlink, which the conversion stops at when it is not named by a pair
    src_dir <- file.path(w, "tracks", "src.track")
    unlink(file.path(src_dir, "README"))
    expect_true(file.symlink(file.path(src_dir, "chr1-chr2"), file.path(src_dir, "chr1")))
    expect_error(gtrack.copy("src", "dst", overwrite = TRUE), "file chr1 is not named by a chromosome pair")
    expect_equal(gextract("dst", gintervals.2d.all()), expected)
    # nor at a file it cannot read, hidden ones too
    unlink(file.path(src_dir, "chr1"))
    Sys.chmod(file.path(src_dir, ".attributes"), "000")
    if (file.access(file.path(src_dir, ".attributes"), 4) != 0) {
        expect_error(gtrack.copy("src", "dst", overwrite = TRUE), "file .attributes cannot be read")
        expect_equal(gextract("dst", gintervals.2d.all()), expected)
    }
    Sys.chmod(file.path(src_dir, ".attributes"), "644")
    # a file the conversion skips does not stop it
    writeLines("x", file.path(src_dir, "track.dat.tmp"))
    gtrack.copy("src", "dst", overwrite = TRUE)
    expect_equal(gextract("dst", gintervals.2d.all())$dst, c(1, 2))
})

test_that("a forced re-conversion of an indexed set packs its files named by aliases, as before", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_force_set_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("s", data.frame(gintervals(c("chr1", "chr2", "chrX"), 0, c(10, 20, 30)), v = 1:3))
        gintervals.save("s2d", data.frame(gintervals.2d(c("chr1", "chr2", "chrX"), 0, 10, "chr2", 0, 10), v = 1:3))
    })
    # named by aliases, as earlier pymisha wrote them, and converted keeping the files
    tracks <- file.path(db, "tracks")
    for (set in c("s.interv", "s2d.interv")) {
        files <- list.files(file.path(tracks, set))
        expect_true(all(file.rename(file.path(tracks, set, files), file.path(tracks, set, alias_names(files)))))
    }
    suppressMessages({
        gintervals.convert_to_indexed("s")
        gintervals.2d.convert_to_indexed("s2d")
    })
    withr::with_options(list(gmax.data.size = 1), {
        gintervals.update("s", data.frame(gintervals("chr1", 500, 600), v = 99L), chrom = "chr1")
        gintervals.update("s2d", data.frame(gintervals.2d("chr1", 500, 600, "chr2", 500, 600), v = 99L), chrom1 = "chr1", chrom2 = "chr2")
        # a deleted chromosome stays deleted
        gintervals.update("s", NULL, chrom = "chrX")
        gintervals.update("s2d", NULL, chrom1 = "chrX", chrom2 = "chr2")
    })
    expected <- list(gintervals.load("s"), gintervals.load("s2d"))
    expect_equal(expected[[1]]$v, c(99, 2))
    expect_equal(expected[[2]]$v, c(99, 2))

    repack <- function() {
        suppressMessages({
            gintervals.convert_to_indexed("s", remove.old = TRUE, force = TRUE)
            gintervals.2d.convert_to_indexed("s2d", remove.old = TRUE, force = TRUE)
        })
        list(gintervals.load("s"), gintervals.load("s2d"))
    }
    expect_equal(repack(), expected)
    expect_equal(list.files(file.path(tracks, "s.interv")), c("intervals.dat", "intervals.idx"))
    expect_equal(list.files(file.path(tracks, "s2d.interv")), c("intervals2d.dat", "intervals2d.idx"))
})

test_that("a symlink named only by the other dimension's names is ignored", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_link_dim_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    gtrack.create_sparse("sp", "x", gintervals(c("chr1", "chr2"), 0, 100), c(1, 2))
    gtrack.2d.create("rects", "x", gintervals.2d(c("chr1", "chr2"), 10, 15, "chr2", 30, 35), c(1, 2))
    expected <- list(gextract("sp", gintervals.all()), gextract("rects", gintervals.2d.all()))
    tracks <- file.path(db, "tracks")
    expect_true(file.symlink(file.path(tracks, "rects.track", "chr1-chr2"), file.path(tracks, "sp.track", "chr1-chr2")))
    expect_true(file.symlink(file.path(tracks, "sp.track", "chr1"), file.path(tracks, "rects.track", "chr1")))
    unlink(file.path(tracks, "rects.track", ".meta"))
    gdb.reload()

    expect_equal(gtrack.info("sp")$type, "sparse")
    expect_equal(gtrack.info("rects")$type, "rectangles")
    expect_equal(list(gextract("sp", gintervals.all()), gextract("rects", gintervals.2d.all())), expected)
    suppressMessages({
        gtrack.convert_to_indexed("sp")
        gtrack.2d.convert_to_indexed("rects", remove.old = TRUE)
    })
    expect_equal(list(gextract("sp", gintervals.all()), gextract("rects", gintervals.2d.all())), expected)
    # the conversions' remove.old keeps them (Sys.readlink() gives NA when there is no such file)
    links <- Sys.readlink(file.path(tracks, c("sp.track/chr1-chr2", "rects.track/chr1")))
    expect_true(all(!is.na(links) & nzchar(links)))
})

test_that("a file named by the destination's mitochondrial name is copied to it across genomes", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_copy_mito_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    mkdb <- function(path, chroms) {
        dir.create(file.path(path, "tracks"), recursive = TRUE)
        dir.create(file.path(path, "seq"))
        for (chrom in chroms) writeBin(charToRaw(strrep("A", 1000)), file.path(path, "seq", paste0("chr", chrom, ".seq")))
        writeLines(paste(chroms, 1000, sep = "\t"), file.path(path, "chrom_sizes.txt"))
        normalizePath(path)
    }
    src <- mkdb(file.path(td, "src"), c("1", "MT"))
    dest <- mkdb(file.path(td, "dest"), c("1", "M"))
    gsetroot(src)
    gtrack.create_sparse("sp", "x", gintervals(c("chr1", "chrMT"), 0, 100), c(1, 2))
    # chrMT's file named "M", an alias of chrMT here and the name of chrM in dest
    expect_true(file.rename(file.path(src, "tracks", "sp.track", "chrMT"), file.path(src, "tracks", "sp.track", "M")))
    gdb.reload()
    expect_equal(gextract("sp", gintervals("chrMT", 0, 1000))$sp, 2)

    gtrack.copy("sp", "sp", db = dest)
    gsetroot(dest)
    expect_equal(gextract("sp", gintervals(c("chr1", "chrM"), 0, 1000))$sp, c(1, 2))
})

test_that("a copy of a file with two names (a symlink, a hard link), neither canonical, reads", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_two_names_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    dest <- create_alias_test_db(file.path(td, "dest"))
    gsetroot(db)

    gtrack.create_sparse("sp", "x", gintervals(c("chr1", "chrM"), 0, 100), c(1, 2))
    dir <- file.path(db, "tracks", "sp.track")
    expect_true(file.rename(file.path(dir, "chrM"), file.path(dir, "MT")))
    expect_true(file.symlink("MT", file.path(dir, "chrMT")))
    gdb.reload()
    expected <- gextract("sp", gintervals.all())
    expect_equal(expected$sp, c(1, 2))

    # each copy makes a file of each name: the same bytes are not two files of chrM
    gtrack.copy("sp", "sp2")
    expect_equal(gextract("sp2", gintervals.all())$sp2, c(1, 2))
    gtrack.copy("sp", "sp", db = dest)
    gsetroot(dest)
    expect_equal(gextract("sp", gintervals.all())$sp, c(1, 2))
    # different bytes are
    gsetroot(db)
    dir2 <- file.path(db, "tracks", "sp2.track")
    expect_true(file.copy(file.path(dir2, "chr1"), file.path(dir2, "chrMT"), overwrite = TRUE))
    gdb.reload()
    expect_error(gextract("sp2", gintervals.all()), "two files of chromosome chrM")
})

test_that("gtrack.modify of a chromosome with two names, neither canonical, leaves both naming the new data", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_modify_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    two_names <- list(
        symlink = function(dir) file.symlink("MT", file.path(dir, "M")),
        hardlink = function(dir) file.link(file.path(dir, "MT"), file.path(dir, "M")),
        copy = function(dir) file.copy(file.path(dir, "MT"), file.path(dir, "M"))
    )
    for (how in names(two_names)) {
        gtrack.create_dense("d", "x", gintervals(c("chr1", "chrM"), 0, 30), c(1, 1), 10, 0)
        dir <- file.path(db, "tracks", "d.track")
        expect_true(file.rename(file.path(dir, "chrM"), file.path(dir, "MT")), info = how)
        expect_true(two_names[[how]](dir), info = how)
        gdb.reload()
        gtrack.modify("d", "d * 7", gintervals("chrM", 0, 30))
        expect_equal(gextract("d", gintervals("chrM", 0, 30))$d, c(7, 7, 7), info = how)
        expect_equal(gextract("d", gintervals("chr1", 0, 30))$d, c(1, 1, 1), info = how)
        # both name the modified file, so what links to them reads it too; no canonical name is added
        sums <- tools::md5sum(file.path(dir, c("M", "MT")))
        expect_true(all(sums == sums[1]), info = how)
        expect_false("chrM" %in% list.files(dir), info = how)
        gtrack.rm("d", force = TRUE)
    }
})

test_that("a copy across genomes takes the name of a chromosome's file the destination has", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_copy_dups_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    mkdb <- function(path, chroms) {
        dir.create(file.path(path, "tracks"), recursive = TRUE)
        dir.create(file.path(path, "seq"))
        for (chrom in chroms) writeBin(charToRaw(strrep("A", 1000)), file.path(path, "seq", paste0("chr", chrom, ".seq")))
        writeLines(paste(chroms, 1000, sep = "\t"), file.path(path, "chrom_sizes.txt"))
        normalizePath(path)
    }
    src <- mkdb(file.path(td, "src"), c("1", "M"))
    dest <- mkdb(file.path(td, "dest"), c("1", "MT"))
    gsetroot(src)
    gtrack.create_sparse("sp", "x", gintervals(c("chr1", "chrM"), 0, 100), c(1, 2))
    # chrM in "M" and the same bytes in "MT"; dest knows chrMT
    dir <- file.path(src, "tracks", "sp.track")
    expect_true(file.rename(file.path(dir, "chrM"), file.path(dir, "M")))
    expect_true(file.copy(file.path(dir, "M"), file.path(dir, "MT")))
    gdb.reload()

    suppressWarnings(gtrack.copy("sp", "sp", db = dest))
    gsetroot(dest)
    expect_equal(gextract("sp", gintervals(c("chr1", "chrMT"), 0, 1000))$sp, c(1, 2))
})

test_that("a symlink to a file that gtrack.modify replaced reads the new data", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_modify_link_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    gtrack.create_dense("d", "x", gintervals(c("chr1", "chr2"), 0, 30), c(1, 2), 10, 0)
    tracks <- file.path(db, "tracks")
    expect_true(file.rename(file.path(tracks, "d.track", "chr1"), file.path(tracks, "d.track", "1")))
    # another track's chr1 is d's file
    dir.create(file.path(tracks, "e.track"))
    files <- setdiff(list.files(file.path(tracks, "d.track"), all.files = TRUE, no.. = TRUE), "1")
    expect_true(all(file.copy(file.path(tracks, "d.track", files), file.path(tracks, "e.track"))))
    expect_true(file.symlink(file.path(tracks, "d.track", "1"), file.path(tracks, "e.track", "chr1")))
    gdb.reload()
    expect_equal(gextract("e", gintervals("chr1", 0, 30))$e, c(1, 1, 1))

    gtrack.modify("d", "d * 7", gintervals("chr1", 0, 30))
    expect_equal(gextract("d", gintervals("chr1", 0, 30))$d, c(7, 7, 7))
    expect_equal(gextract("e", gintervals("chr1", 0, 30))$e, c(7, 7, 7))
})

test_that("an indexed big set does not return a chromosome deleted with gintervals.update(NULL)", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_deleted_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    withr::with_options(list(gmax.data.size = 1, gbig.intervals.size = 1), {
        gintervals.save("s", data.frame(gintervals(c("chr1", "chr2", "chrX"), 0, 10), v = 1:3))
        gintervals.save("s2d", data.frame(gintervals.2d(c("chr1", "chr2", "chrX"), 0, 10, "chr2", 0, 10), v = 1:3))
        suppressMessages({
            gintervals.convert_to_indexed("s", remove.old = TRUE)
            gintervals.2d.convert_to_indexed("s2d", remove.old = TRUE)
        })
        gintervals.update("s", NULL, chrom = "chr1")
        gintervals.update("s2d", NULL, chrom1 = "chr1", chrom2 = "chr2")
    })
    expect_equal(gintervals.load("s")$v, 2:3)
    expect_equal(gintervals.load("s2d")$v, 2:3)
})

test_that("gintervals.update over a canonically named file does not list the set's directory", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_update_list_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    listed <- 0
    gchrom_files <- .gchrom_files
    local_mocked_bindings(.gchrom_files = function(...) {
        listed <<- listed + 1
        gchrom_files(...)
    })
    withr::with_options(list(gmax.data.size = 1, gbig.intervals.size = 1), {
        gintervals.save("s", data.frame(gintervals(c("chr1", "chr2", "chrX"), 0, 10), v = 1:3))
        gintervals.update("s", data.frame(gintervals("chr1", 0, 20), v = 5L), chrom = "chr1")
        expect_equal(listed, 0)
        # a deletion lists it: a file named by an alias would come back
        writeLines("x", file.path(db, "tracks", "s.interv", "2"))
        gintervals.update("s", NULL, chrom = "chr2")
        expect_equal(listed, 1)
    })
    expect_false(file.exists(file.path(db, "tracks", "s.interv", "2")))
    expect_equal(gintervals.load("s")$v, c(5, 3))
})

test_that("gtrack.modify keeps a file named by an alias where earlier pymisha also writes it", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_modify_pymisha_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    gtrack.create_dense("d", "x", gintervals(c("chr1", "chr2"), 0, 30), c(1, 2), 10, 0)
    dir <- file.path(db, "tracks", "d.track")
    files <- list.files(dir)
    expect_true(all(file.rename(file.path(dir, files), file.path(dir, alias_names(files)))))
    gdb.reload()
    gtrack.modify("d", "d + 100", gintervals("chr1", 0, 30))
    expect_equal(sort(list.files(dir)), sort(alias_names(files)))

    # earlier pymisha writes by renaming a new file over "1": misha reads what it wrote
    saved <- file.path(td, "chr1_101")
    expect_true(file.copy(file.path(dir, "1"), saved))
    gtrack.modify("d", "d + 1000", gintervals("chr1", 0, 30))
    expect_equal(gextract("d", gintervals("chr1", 0, 30))$d, c(1101, 1101, 1101))
    expect_true(file.rename(saved, file.path(dir, "1")))
    gdb.reload()
    expect_equal(gextract("d", gintervals("chr1", 0, 30))$d, c(101, 101, 101))
})

test_that("after gtrack.modify, earlier pymisha's write of the \"chr\"-toggled name is read", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_modify_toggled_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    second_name <- list(
        symlink = function(dir) file.symlink("M", file.path(dir, "MT")),
        hardlink = function(dir) file.link(file.path(dir, "M"), file.path(dir, "MT")),
        copy = function(dir) file.copy(file.path(dir, "M"), file.path(dir, "MT"))
    )
    for (how in names(second_name)) {
        gtrack.create_dense("d", "x", gintervals(c("chr1", "chrM"), 0, 30), c(1, 5), 10, 0)
        dir <- file.path(db, "tracks", "d.track")
        expect_true(file.rename(file.path(dir, "chrM"), file.path(dir, "M")), info = how)
        expect_true(second_name[[how]](dir), info = how)
        gdb.reload()
        saved <- file.path(td, "M_5")
        expect_true(file.copy(file.path(dir, "M"), saved, overwrite = TRUE), info = how)
        gtrack.modify("d", "d + 100", gintervals("chrM", 0, 30))
        expect_equal(gextract("d", gintervals("chrM", 0, 30))$d, c(105, 105, 105), info = how)
        # earlier pymisha writes chrM's data by renaming a new file over "M", its name for it
        expect_true(file.rename(saved, file.path(dir, "M")), info = how)
        gdb.reload()
        expect_equal(gextract("d", gintervals("chrM", 0, 30))$d, c(5, 5, 5), info = how)
        expect_equal(gextract("d", gintervals("chr1", 0, 30))$d, c(1, 1, 1), info = how)
        gtrack.rm("d", force = TRUE)
    }
})

test_that("a file read through a symlink that cannot be reached is an error, not a missing chromosome", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_unreachable_")
    dir.create(td)
    withr::defer({
        Sys.chmod(file.path(td, "closed"), "755")
        unlink(td, recursive = TRUE)
    })
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    gtrack.create_sparse("sp", "x", gintervals(c("chr1", "chr2"), 0, 100), c(1, 2))
    gtrack.2d.create("rects", "x", gintervals.2d(c("chr1", "chr2"), 10, 15, "chr2", 30, 35), c(1, 2))
    expected <- gextract("sp", gintervals.all())
    tracks <- file.path(db, "tracks")
    closed <- file.path(td, "closed")
    dir.create(closed)
    for (file in c("sp.track/chr2", "rects.track/chr1-chr2")) {
        target <- file.path(closed, basename(file))
        expect_true(file.rename(file.path(tracks, file), target))
        expect_true(file.symlink(target, file.path(tracks, file)))
    }
    unlink(file.path(tracks, "rects.track", ".meta"))
    Sys.chmod(closed, "000")
    skip_if(file.access(file.path(closed, "chr2"), 0) == 0, "the directory stays searchable to this user")
    gdb.reload()

    expect_error(gextract("sp", gintervals.all()), "Permission denied")
    expect_error(gextract("rects", gintervals.2d.all()), "Permission denied")
    expect_error(gintervals.load("sp"), "Permission denied")
    expect_error(gtrack.copy("rects", "rects2"))
    expect_false(gtrack.exists("rects2"))
    # a canonically named file wins over it, as master read
    Sys.chmod(closed, "755")
    unlink(file.path(tracks, "sp.track", "chr2"))
    expect_true(file.copy(file.path(closed, "chr2"), file.path(tracks, "sp.track", "chr2")))
    expect_true(file.symlink(file.path(closed, "chr2"), file.path(tracks, "sp.track", "2")))
    Sys.chmod(closed, "000")
    gdb.reload()
    expect_equal(gextract("sp", gintervals.all()), expected)
    expect_equal(gtrack.info("sp")$type, "sparse")

    # a copy stops at an entry it cannot follow (vars), where it skips a dangling symlink
    Sys.chmod(closed, "755")
    gtrack.var.set("sp", "v", 1)
    expect_true(file.rename(file.path(tracks, "sp.track", "vars"), file.path(closed, "vars")))
    expect_true(file.symlink(file.path(closed, "vars"), file.path(tracks, "sp.track", "vars")))
    Sys.chmod(closed, "000")
    expect_error(suppressWarnings(gtrack.copy("sp", "sp_copy")))
    Sys.chmod(closed, "755")
})

test_that("gintervals.update removes a chromosome's alias-named files only once its write succeeded", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_update_fail_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_alias_test_db(file.path(td, "db"))
    gsetroot(db)

    withr::with_options(list(gmax.data.size = 1, gbig.intervals.size = 1), {
        gintervals.save("s", data.frame(gintervals(c("chr1", "chr2", "chrX"), 0, 10), v = 1:3))
        dir <- file.path(db, "tracks", "s.interv")
        expect_true(file.rename(file.path(dir, "chr1"), file.path(dir, "1")))
        # the write of chr1 fails (a directory in its place)
        dir.create(file.path(dir, "chr1"))
        expect_error(suppressWarnings(gintervals.update("s", data.frame(gintervals("chr1", 0, 20), v = 5L), chrom = "chr1")))
        expect_true(file.exists(file.path(dir, "1")))
        unlink(file.path(dir, "chr1"), recursive = TRUE)
    })
    expect_equal(gintervals.load("s")$v, 1:3)
})

test_that("an overwriting copy of a computed track into a database without its chromosomes keeps the destination", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("alias_copy_computed_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    src <- create_alias_test_db(file.path(td, "src"))
    dest <- file.path(td, "dest")
    dir.create(file.path(dest, "tracks"), recursive = TRUE)
    dir.create(file.path(dest, "seq"))
    for (chrom in c("A", "B")) writeBin(charToRaw(strrep("A", 1000)), file.path(dest, "seq", paste0("chr", chrom, ".seq")))
    writeLines(paste(c("A", "B"), 1000, sep = "\t"), file.path(dest, "chrom_sizes.txt"))
    gsetroot(dest)
    gtrack.create_sparse("victim", "x", gintervals(c("chrA", "chrB"), 0, 10), c(7, 8))

    gsetroot(src)
    .gtrack.create_test_computer2d("cmp", 0, 3, 10)
    expect_error(suppressWarnings(gtrack.copy("cmp", "victim", db = dest, overwrite = TRUE)), "no chromosomes from source database")
    gsetroot(dest)
    expect_equal(gextract("victim", gintervals.all())$victim, c(7, 8))
})
