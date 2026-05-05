# Regression test for gintervals.save / gintervals.load on a bigset whose
# input has chrom as character (e.g., a tibble from dplyr). Previously this
# failed with "invalid columns definition" because the on-disk per-chrom
# files store $chrom as factor while the persisted .meta zeroline kept the
# original character class.

test_that("bigset save then load works when input has character chrom", {
    skip_if_not_installed("withr")
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    old_gwd <- get("GWD", envir = misha:::.misha)
    old_root <- dirname(old_gwd)
    tmp_root <- withr::local_tempdir()

    utils::write.table(
        data.frame(chrom = c("chr1", "chr2", "chr3"), size = c(1e6, 1e6, 1e6)),
        file.path(tmp_root, "chrom_sizes.txt"),
        quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
    )
    dir.create(file.path(tmp_root, "tracks", "intervs", "global"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(tmp_root, "seq"), recursive = TRUE, showWarnings = FALSE)
    gdb.init(tmp_root)

    withr::defer(
        {
            if (dir.exists(old_root)) {
                try(gsetroot(old_root), silent = TRUE)
                try(gdb.reload(), silent = TRUE)
            }
        },
        envir = parent.frame()
    )
    gsetroot(tmp_root)
    gdb.reload()

    withr::defer(gintervals.rm("intervs.global.rmsk_test", force = TRUE), envir = parent.frame())

    # Plain data.frame with character chrom and extra metadata columns
    intervs <- data.frame(
        chrom = c("chr1", "chr1", "chr2", "chr2", "chr3"),
        start = c(0L, 100L, 0L, 200L, 0L),
        end = c(50L, 150L, 100L, 300L, 50L),
        strand = c(1L, -1L, 1L, 1L, -1L),
        name = c("a", "b", "c", "d", "e"),
        stringsAsFactors = FALSE
    )
    expect_identical(class(intervs$chrom), "character")

    # Force the bigset (per-chromosome) save path
    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("intervs.global.rmsk_test", intervs)
    })

    # Loading must not error out
    res <- gintervals.load("intervs.global.rmsk_test")
    expect_equal(nrow(res), nrow(intervs))
    expect_setequal(as.character(res$chrom), as.character(intervs$chrom))
    # On-disk format normalizes chrom to factor; load returns it as factor
    expect_s3_class(res$chrom, "factor")
    expect_true(all(c("strand", "name") %in% colnames(res)))
})

test_that("bigset save then load works for tibble input with character chrom", {
    skip_if_not_installed("withr")
    skip_if_not_installed("tibble")
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    old_gwd <- get("GWD", envir = misha:::.misha)
    old_root <- dirname(old_gwd)
    tmp_root <- withr::local_tempdir()

    utils::write.table(
        data.frame(chrom = c("chr1", "chr2", "chr3"), size = c(1e6, 1e6, 1e6)),
        file.path(tmp_root, "chrom_sizes.txt"),
        quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
    )
    dir.create(file.path(tmp_root, "tracks", "intervs", "global"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(tmp_root, "seq"), recursive = TRUE, showWarnings = FALSE)
    gdb.init(tmp_root)

    withr::defer(
        {
            if (dir.exists(old_root)) {
                try(gsetroot(old_root), silent = TRUE)
                try(gdb.reload(), silent = TRUE)
            }
        },
        envir = parent.frame()
    )
    gsetroot(tmp_root)
    gdb.reload()

    withr::defer(gintervals.rm("intervs.global.rmsk_tib", force = TRUE), envir = parent.frame())

    intervs <- tibble::tibble(
        chrom = c("chr1", "chr1", "chr2", "chr2", "chr3"),
        start = c(0L, 100L, 0L, 200L, 0L),
        end = c(50L, 150L, 100L, 300L, 50L),
        strand = c(1L, -1L, 1L, 1L, -1L),
        name = c("a", "b", "c", "d", "e")
    )

    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("intervs.global.rmsk_tib", intervs)
    })

    res <- gintervals.load("intervs.global.rmsk_tib")
    expect_equal(nrow(res), nrow(intervs))
    expect_setequal(as.character(res$chrom), as.character(intervs$chrom))
})
