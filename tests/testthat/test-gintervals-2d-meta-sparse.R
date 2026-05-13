# Phase 7a: verify that GIntervalsMeta2D uses a sparse in-memory
# representation. Before this change, a 2D interval set living in a database
# with many contigs allocated O(N^2) entries at load time. We exercise the
# load + query path here with a moderate-sized N (1000 contigs) and only a
# handful of populated chrom-pairs.

build_many_contig_db <- function(n_contigs, defer_envir) {
    skip_if_not_installed("withr")
    tmp_root <- withr::local_tempdir(.local_envir = defer_envir)
    chrom_df <- data.frame(
        chrom = as.character(seq_len(n_contigs)),
        size = rep(100000L, n_contigs)
    )
    utils::write.table(
        chrom_df,
        file.path(tmp_root, "chrom_sizes.txt"),
        quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
    )
    dir.create(file.path(tmp_root, "tracks", "test"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(tmp_root, "seq"), recursive = TRUE, showWarnings = FALSE)
    gdb.init(tmp_root)
    tmp_root
}

switch_to_db <- function(tmp_root, defer_envir) {
    old_gwd <- get("GWD", envir = .misha)
    old_root <- dirname(old_gwd)
    withr::defer(
        {
            if (dir.exists(old_root)) {
                try(gsetroot(old_root), silent = TRUE)
                try(gdb.reload(), silent = TRUE)
            }
        },
        envir = defer_envir
    )
    gsetroot(tmp_root)
    gdb.reload()
}

test_that("GIntervalsMeta2D loads OK on many-contig DB with few populated pairs", {
    skip_if_not_installed("withr")
    withr::local_options(list(gmulticontig.indexed_format = FALSE))

    n_contigs <- 1000L

    tmp_root <- build_many_contig_db(n_contigs, parent.frame())
    switch_to_db(tmp_root, parent.frame())
    withr::defer(
        try(gintervals.rm("test.sparse_2d", force = TRUE), silent = TRUE),
        envir = parent.frame()
    )

    # Only 3 populated chrom-pairs out of 1000*1000 = 1e6 possible pairs.
    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 500),
        starts1 = c(0, 1000, 0),
        ends1 = c(500, 1500, 500),
        chroms2 = c(1, 2, 999),
        starts2 = c(0, 0, 0),
        ends2 = c(500, 500, 500)
    )
    withr::with_options(list(gmax.data.size = 1), {
        misha::gintervals.save("test.sparse_2d", intervs)
    })

    # Load: with the prior dense N^2 layout this allocated ~24 MB per
    # interval set at N=1000 just for m_chroms2size + m_surfaces +
    # m_contains_overlaps. At N=1e6 that's 32 TB and OOMs. We exercise
    # the load path here as a sanity check; the real win is at larger N.
    withr::with_options(list(gmax.data.size = 1e9), {
        loaded <- misha::gintervals.load("test.sparse_2d")
        expect_equal(nrow(loaded), nrow(intervs))
        # Order may be canonicalized; compare sorted unique pair tuples.
        expect_setequal(
            paste(loaded$chrom1, loaded$start1, loaded$chrom2, loaded$start2),
            paste(intervs$chrom1, intervs$start1, intervs$chrom2, intervs$start2)
        )
    })
})

test_that("gintervals.2d.all enumerates contigs on many-contig DB", {
    skip_if_not_installed("withr")
    withr::local_options(list(gmulticontig.indexed_format = FALSE))

    n_contigs <- 500L

    tmp_root <- build_many_contig_db(n_contigs, parent.frame())
    switch_to_db(tmp_root, parent.frame())

    # gintervals.2d.all enumerates chrom-pair cells through the
    # GIntervalsFetcher2D get_next_chroms walk. Make sure that path is
    # still happy after the sparse rewrite.
    all2d <- gintervals.2d.all()
    expect_true(nrow(all2d) > 0)
})

test_that("masked copy preserves populated pair stats", {
    skip_if_not_installed("withr")
    withr::local_options(list(gmulticontig.indexed_format = FALSE))

    n_contigs <- 200L

    tmp_root <- build_many_contig_db(n_contigs, parent.frame())
    switch_to_db(tmp_root, parent.frame())
    withr::defer(
        try(gintervals.rm("test.masked_2d", force = TRUE), silent = TRUE),
        envir = parent.frame()
    )

    # 5 populated chrom-pairs across the contig range.
    chroms1 <- c(1, 1, 50, 100, 150)
    chroms2 <- c(1, 100, 1, 50, 199)
    intervs <- gintervals.2d(
        chroms1 = chroms1,
        starts1 = rep(0L, 5),
        ends1 = rep(1000L, 5),
        chroms2 = chroms2,
        starts2 = rep(0L, 5),
        ends2 = rep(1000L, 5)
    )
    withr::with_options(list(gmax.data.size = 1), {
        misha::gintervals.save("test.masked_2d", intervs)
    })

    # gintervals.intersect on a 2D bigset takes the meta load + walk path,
    # exercising the sparse iterator end-to-end. Intersect with a 2D mask
    # covering the first chrom-pair only.
    mask <- gintervals.2d(
        chroms1 = 1, starts1 = 0L, ends1 = 1000L,
        chroms2 = 1, starts2 = 0L, ends2 = 1000L
    )
    withr::with_options(list(gmax.data.size = 1e9), {
        hit <- gintervals.intersect("test.masked_2d", mask)
        # The first row is (1, 1) -> matched; the others on other pairs are excluded.
        expect_equal(nrow(hit), 1)
        expect_equal(as.character(hit$chrom1), "1")
        expect_equal(as.character(hit$chrom2), "1")
    })
})
