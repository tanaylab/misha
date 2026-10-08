# A big intervals set whose per-chromosome or per-pair files are named by an alias of the
# chromosomes ("1" for chr1, "1-chr2" for chr1 and chr2) is read as if its files had the
# canonical names, as a track's files are.

restore_groot_on_exit()

test_that("a big intervals set whose files are named by chromosome aliases reads as with canonical names", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("bigset_alias_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))

    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    gsetroot(db)

    intervs <- gintervals(
        c("chr1", "chr1", "chr2", "chrX", "chr1_KI270706v1_random"),
        c(0, 100, 0, 10, 5), c(50, 150, 100, 20, 15)
    )
    intervs2d <- gintervals.2d(
        c("chr1", "chr2", "chr1"), c(0, 0, 10), c(50, 100, 20),
        c("chr1", "chr2", "chrX"), c(0, 0, 10), c(50, 100, 20)
    )
    # per-chromosome big sets
    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("bigset", intervs)
        gintervals.save("bigset2d", intervs2d)
    })
    expect_true(all(c("chr1", "chr2") %in% list.files(file.path(db, "tracks", "bigset.interv"))))
    expect_true("chr1-chrX" %in% list.files(file.path(db, "tracks", "bigset2d.interv")))
    gtrack.create_sparse("sp", "x", intervs, seq_len(nrow(intervs)))
    gtrack.2d.create("rects", "x", intervs2d, seq_len(nrow(intervs2d)))

    read_all <- function() {
        list(
            load = gintervals.load("bigset"),
            load_chrom = gintervals.load("bigset", chrom = "chr1"),
            as_scope = gextract("sp", "bigset"),
            load2d = gintervals.load("bigset2d"),
            load_pair = gintervals.load("bigset2d", chrom1 = "chr1", chrom2 = "chrX"),
            as_scope2d = gextract("rects", "bigset2d")
        )
    }
    expected <- read_all()
    expect_equal(nrow(expected$load), 5)
    expect_equal(nrow(expected$load2d), 3)

    # "chr1" becomes "1"; "chr1-chr2" becomes "1-2", "1-chr2" or "chr1-2" in turn
    alias_name <- function(file, i) {
        chroms <- strsplit(file, "-", fixed = TRUE)[[1]]
        strip <- if (length(chroms) == 1) TRUE else list(c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, TRUE))[[(i - 1) %% 3 + 1]]
        chroms[strip] <- sub("^chr", "", chroms[strip])
        paste(chroms, collapse = "-")
    }
    for (set in c("bigset", "bigset2d")) {
        set_dir <- file.path(db, "tracks", paste0(set, ".interv"))
        files <- grep("^chr", list.files(set_dir), value = TRUE)
        expect_true(all(file.rename(file.path(set_dir, files), file.path(set_dir, mapply(alias_name, files, seq_along(files))))))
    }

    expect_equal(read_all(), expected)
})

test_that("gintervals.2d.convert_to_indexed packs the file the readers use when aliases name a pair twice", {
    local_db_state()
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    td <- tempfile("bigset_dup_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    gsetroot(db)

    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("dupset", gintervals.2d(c("chr1", "chr2"), c(0, 0), c(50, 100), c("chr2", "chr2"), c(0, 0), c(50, 100)))
        gintervals.save("otherset", gintervals.2d(c("chr1", "chr2"), c(200, 200), c(300, 300), c("chr2", "chr2"), c(200, 200), c(300, 300)))
    })
    # an alias-named file of a pair that also has its canonical file, which the readers use
    set_dir <- file.path(db, "tracks", "dupset.interv")
    expect_true(file.copy(file.path(db, "tracks", "otherset.interv", "chr1-chr2"), file.path(set_dir, "1-2")))
    expected <- gintervals.load("dupset")
    expect_equal(nrow(expected), 2)

    suppressMessages(gintervals.2d.convert_to_indexed("dupset", remove.old = TRUE))
    expect_true(file.exists(file.path(set_dir, "intervals2d.idx")))
    # remove.old removes every file of a packed pair
    expect_false(any(c("chr1-chr2", "1-2") %in% list.files(set_dir)))
    expect_equal(gintervals.load("dupset"), expected)
})
