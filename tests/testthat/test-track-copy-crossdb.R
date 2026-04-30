# tests/testthat/test-track-copy-crossdb.R
test_that(".gdb.is_indexed_at and .gdb.chrom_names_at probe a db without loading it", {
    withr::with_tempdir({
        create_test_db("perchrom_db")
        expect_false(misha:::.gdb.is_indexed_at(normalizePath("perchrom_db")))
        expect_equal(
            misha:::.gdb.chrom_names_at(normalizePath("perchrom_db")),
            c("chr1", "chr2")
        )
    })
})

test_that(".gdb.is_indexed_at returns TRUE for an indexed db", {
    withr::with_tempdir({
        create_test_db("idx_db")
        gdb.init("idx_db")
        gdb.convert_to_indexed(force = TRUE, verbose = FALSE)
        expect_true(misha:::.gdb.is_indexed_at(normalizePath("idx_db")))
    })
})
