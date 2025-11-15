test_that("gtrack.liftover multi-target aggregation policies", {
    local_db_state()

    # Source genome with a single chromosome
    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 400), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Dense track with three bins (values 1, 2, 3; binsize 10)
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 3),
        start = c(0, 10, 20),
        end = c(10, 20, 30),
        stringsAsFactors = FALSE
    )
    src_values <- c(1, 2, 3)
    gtrack.create_dense("agg_dense_src", "Aggregation dense source", src_intervals, src_values, 10, NaN)
    src_track_dir <- file.path(source_db, "tracks", "agg_dense_src.track")

    # Variant with NA in the middle bin
    src_values_na <- c(1, NaN, 3)
    gtrack.create_dense("agg_dense_src_na", "Aggregation dense source NA", src_intervals, src_values_na, 10, NaN)
    src_track_dir_na <- file.path(source_db, "tracks", "agg_dense_src_na.track")

    # Target genome
    setup_db(list(paste0(">chrA\n", paste(rep("T", 400), collapse = ""), "\n")))

    # Chain mappings: all three bins map (partially) into chrA[0,10)
    chain_file <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 400, "+", 0, 10, "chrA", 400, "+", 0, 10, 1) # coverage len 10, value 1
    write_chain_entry(chain_file, "chrsource1", 400, "+", 10, 20, "chrA", 400, "+", 3, 13, 2) # coverage len 7, value 2
    write_chain_entry(chain_file, "chrsource1", 400, "+", 20, 30, "chrA", 400, "+", 7, 17, 3) # coverage len 3, value 3

    lifted_track <- "agg_dense_lifted"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    extract_bin <- function(track) {
        res <- gextract(track, gintervals("chrA", 0, 10))
        unique(as.numeric(res[[track]]))
    }

    liftover_with <- function(..., track_dir = src_track_dir, agg = "mean", params = NULL,
                              na.rm = TRUE, min_n = NULL, tgt_policy = "keep", desc = "lifted") {
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
        gtrack.liftover(
            lifted_track, desc, track_dir, chain_file,
            tgt_overlap_policy = tgt_policy,
            multi_target_agg = agg,
            params = params,
            na.rm = na.rm,
            min_n = min_n
        )
        extract_bin(lifted_track)
    }

    expect_equal(liftover_with(), 2) # mean of 1,2,3
    expect_equal(liftover_with(agg = "sum"), 6)
    expect_equal(liftover_with(agg = "min"), 1)
    expect_equal(liftover_with(agg = "max"), 3)
    expect_equal(liftover_with(agg = "median"), 2)
    expect_equal(liftover_with(agg = "count"), 3)

    expect_equal(liftover_with(agg = "first"), 1)
    expect_equal(liftover_with(agg = "last"), 3)
    expect_equal(liftover_with(agg = "nth", params = 2), 2)
    expect_equal(liftover_with(agg = "nth", params = list(n = 3)), 3)

    expect_equal(liftover_with(agg = "max.coverage_len"), 1)
    expect_equal(liftover_with(agg = "min.coverage_len"), 3)
    expect_equal(liftover_with(agg = "max.coverage_frac"), 1)
    expect_equal(liftover_with(agg = "min.coverage_frac"), 3)

    # NA handling
    expect_equal(liftover_with(track_dir = src_track_dir_na, agg = "mean", na.rm = TRUE), 2)
    expect_true(is.nan(liftover_with(track_dir = src_track_dir_na, agg = "mean", na.rm = FALSE)))

    # min_n gating
    expect_true(is.nan(liftover_with(track_dir = src_track_dir_na, agg = "mean", na.rm = TRUE, min_n = 3)))
    expect_equal(liftover_with(track_dir = src_track_dir_na, agg = "mean", na.rm = TRUE, min_n = 2), 2)
})

test_that("gtrack.liftover aggregation edge cases: all NAs", {
    local_db_state()

    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 100), collapse = ""), "\n", sep = "", file = source_fasta)
    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # All NA values
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 3),
        start = c(0, 10, 20),
        end = c(10, 20, 30),
        stringsAsFactors = FALSE
    )
    src_values <- c(NaN, NaN, NaN)
    gtrack.create_dense("all_na_src", "All NA source", src_intervals, src_values, 10, NaN)
    src_track_dir <- file.path(source_db, "tracks", "all_na_src.track")

    setup_db(list(paste0(">chrA\n", paste(rep("T", 100), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 10, "chrA", 100, "+", 0, 10, 1)
    write_chain_entry(chain_file, "chrsource1", 100, "+", 10, 20, "chrA", 100, "+", 3, 13, 2)
    write_chain_entry(chain_file, "chrsource1", 100, "+", 20, 30, "chrA", 100, "+", 7, 17, 3)

    lifted_track <- "all_na_lifted"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    # All aggregators should return NA when all values are NA (with na.rm=TRUE)
    for (agg in c("mean", "sum", "min", "max", "median", "first", "last", "max.coverage_len")) {
        gtrack.liftover(
            lifted_track, "test", src_track_dir, chain_file,
            tgt_overlap_policy = "keep",
            multi_target_agg = agg,
            na.rm = TRUE
        )
        res <- gextract(lifted_track, gintervals("chrA", 0, 10))
        expect_true(all(is.nan(res[[lifted_track]])), info = paste("aggregator:", agg))
        gtrack.rm(lifted_track, force = TRUE)
    }

    # count should return 0 for all NAs
    gtrack.liftover(
        lifted_track, "test", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "count",
        na.rm = TRUE
    )
    res <- gextract(lifted_track, gintervals("chrA", 0, 10))
    expect_equal(unique(res[[lifted_track]]), 0)
})

test_that("gtrack.liftover aggregation edge cases: ties and sorting with dense tracks", {
    local_db_state()

    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)
    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create dense track with values 5, 3, 5, 2 (testing ties)
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 4),
        start = c(0, 10, 20, 30),
        end = c(10, 20, 30, 40),
        stringsAsFactors = FALSE
    )
    src_values <- c(5, 3, 5, 2)
    gtrack.create_dense("tie_src", "Tie source", src_intervals, src_values, 10, NaN)
    src_track_dir <- file.path(source_db, "tracks", "tie_src.track")

    setup_db(list(paste0(">chrB\n", paste(rep("G", 200), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Map all to overlapping target locus [50,60)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 10, "chrB", 200, "+", 50, 60, 1)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 10, 20, "chrB", 200, "+", 50, 60, 2)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 20, 30, "chrB", 200, "+", 50, 60, 3)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 30, 40, "chrB", 200, "+", 50, 60, 4)

    lifted_track <- "tie_lifted"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    # first: all start at 50, end at 60, sorted by value descending: 5,5,3,2 -> pick first (5)
    gtrack.liftover(
        lifted_track, "first", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "first"
    )
    res <- gextract(lifted_track, gintervals("chrB", 50, 60))
    expect_equal(unique(res[[lifted_track]]), 5)
    gtrack.rm(lifted_track, force = TRUE)

    # last: pick last in sorted order (2)
    gtrack.liftover(
        lifted_track, "last", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "last"
    )
    res <- gextract(lifted_track, gintervals("chrB", 50, 60))
    expect_equal(unique(res[[lifted_track]]), 2)
    gtrack.rm(lifted_track, force = TRUE)

    # nth with n=2: sorted by value desc is 5,5,3,2, second is 5
    gtrack.liftover(
        lifted_track, "nth", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "nth",
        params = 2
    )
    res <- gextract(lifted_track, gintervals("chrB", 50, 60))
    expect_equal(unique(res[[lifted_track]]), 5)
    gtrack.rm(lifted_track, force = TRUE)

    # min: should be 2
    gtrack.liftover(
        lifted_track, "min", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "min"
    )
    res <- gextract(lifted_track, gintervals("chrB", 50, 60))
    expect_equal(unique(res[[lifted_track]]), 2)
    gtrack.rm(lifted_track, force = TRUE)

    # max: should be 5
    gtrack.liftover(
        lifted_track, "max", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "max"
    )
    res <- gextract(lifted_track, gintervals("chrB", 50, 60))
    expect_equal(unique(res[[lifted_track]]), 5)
})

test_that("gtrack.liftover median with even/odd numbers of contributors", {
    local_db_state()

    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)
    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create track with 4 values for even-count median test
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 4),
        start = c(0, 10, 20, 30),
        end = c(10, 20, 30, 40),
        stringsAsFactors = FALSE
    )
    src_values <- c(1, 2, 3, 4)
    gtrack.create_dense("median_even_src", "Median even", src_intervals, src_values, 10, NaN)
    src_track_dir_even <- file.path(source_db, "tracks", "median_even_src.track")

    # Create track with 3 values for odd-count median test
    src_intervals_odd <- data.frame(
        chrom = rep("chrsource1", 3),
        start = c(50, 60, 70),
        end = c(60, 70, 80),
        stringsAsFactors = FALSE
    )
    src_values_odd <- c(5, 7, 9)
    gtrack.create_dense("median_odd_src", "Median odd", src_intervals_odd, src_values_odd, 10, NaN)
    src_track_dir_odd <- file.path(source_db, "tracks", "median_odd_src.track")

    setup_db(list(paste0(">chrD\n", paste(rep("T", 200), collapse = ""), "\n")))

    # Even count test
    chain_file_even <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file_even, "chrsource1", 200, "+", 0, 10, "chrD", 200, "+", 10, 20, 1)
    write_chain_entry(chain_file_even, "chrsource1", 200, "+", 10, 20, "chrD", 200, "+", 10, 20, 2)
    write_chain_entry(chain_file_even, "chrsource1", 200, "+", 20, 30, "chrD", 200, "+", 10, 20, 3)
    write_chain_entry(chain_file_even, "chrsource1", 200, "+", 30, 40, "chrD", 200, "+", 10, 20, 4)

    lifted_track <- "median_lifted"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    # Median of 1,2,3,4 should be (2+3)/2 = 2.5
    gtrack.liftover(
        lifted_track, "median even", src_track_dir_even, chain_file_even,
        tgt_overlap_policy = "keep",
        multi_target_agg = "median"
    )
    res <- gextract(lifted_track, gintervals("chrD", 10, 20))
    expect_equal(unique(res[[lifted_track]]), 2.5)
    gtrack.rm(lifted_track, force = TRUE)

    # Odd count test
    chain_file_odd <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file_odd, "chrsource1", 200, "+", 50, 60, "chrD", 200, "+", 30, 40, 1)
    write_chain_entry(chain_file_odd, "chrsource1", 200, "+", 60, 70, "chrD", 200, "+", 30, 40, 2)
    write_chain_entry(chain_file_odd, "chrsource1", 200, "+", 70, 80, "chrD", 200, "+", 30, 40, 3)

    # Median of 5,7,9 should be 7
    gtrack.liftover(
        lifted_track, "median odd", src_track_dir_odd, chain_file_odd,
        tgt_overlap_policy = "keep",
        multi_target_agg = "median"
    )
    res <- gextract(lifted_track, gintervals("chrD", 30, 40))
    expect_equal(unique(res[[lifted_track]]), 7)
})

test_that("gtrack.liftover aggregation for sparse tracks", {
    local_db_state()

    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)
    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 3),
        start = c(0, 10, 20),
        end = c(10, 20, 30),
        stringsAsFactors = FALSE
    )
    src_vals <- c(5, 7, 11)
    gtrack.create_sparse("agg_sparse_src", "Aggregation sparse source", src_intervals, src_vals)
    src_track_dir <- file.path(source_db, "tracks", "agg_sparse_src.track")

    setup_db(list(paste0(">chrB\n", paste(rep("G", 200), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Map each interval onto the same target locus
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 10, "chrB", 200, "+", 100, 110, 1)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 10, 20, "chrB", 200, "+", 100, 110, 2)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 20, 30, "chrB", 200, "+", 100, 110, 3)

    lifted_track <- "agg_sparse_lifted"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(
        lifted_track, "sum sparse", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "sum"
    )
    vals <- gextract(lifted_track, gintervals("chrB", 100, 110))[[lifted_track]]
    expect_true(all(vals == 23))

    gtrack.rm(lifted_track, force = TRUE)
    gtrack.liftover(
        lifted_track, "count sparse", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "count"
    )
    vals <- gextract(lifted_track, gintervals("chrB", 100, 110))[[lifted_track]]
    expect_true(all(vals == 3))

    gtrack.rm(lifted_track, force = TRUE)
    gtrack.liftover(
        lifted_track, "first sparse", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "first"
    )
    vals <- gextract(lifted_track, gintervals("chrB", 100, 110))[[lifted_track]]
    expect_true(all(vals == 11))

    gtrack.rm(lifted_track, force = TRUE)
    gtrack.liftover(
        lifted_track, "nth sparse", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "nth",
        params = 2
    )
    vals <- gextract(lifted_track, gintervals("chrB", 100, 110))[[lifted_track]]
    expect_true(all(vals == 7))

    gtrack.rm(lifted_track, force = TRUE)
    gtrack.liftover(
        lifted_track, "min_n sparse", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "sum",
        min_n = 4
    )
    vals <- gextract(lifted_track, gintervals("chrB", 100, 110))[[lifted_track]]
    expect_true(all(is.nan(vals)))
})

test_that("gtrack.liftover nth aggregator validates params", {
    local_db_state()

    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 100), collapse = ""), "\n", sep = "", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 1),
        start = 0,
        end = 10,
        stringsAsFactors = FALSE
    )
    src_vals <- 1
    gtrack.create_sparse("agg_sparse_single", "single", src_intervals, src_vals)
    src_track_dir <- file.path(source_db, "tracks", "agg_sparse_single.track")

    setup_db(list(paste0(">chrC\n", paste(rep("C", 100), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 10, "chrC", 100, "+", 0, 10, 1)

    expect_error(
        gtrack.liftover(
            "nth_bad", "desc", src_track_dir, chain_file,
            multi_target_agg = "nth",
            tgt_overlap_policy = "keep"
        ),
        "params must be supplied for 'nth' aggregation"
    )

    expect_error(
        gtrack.liftover(
            "nth_bad2", "desc", src_track_dir, chain_file,
            multi_target_agg = "nth",
            params = list(),
            tgt_overlap_policy = "keep"
        ),
        "params list must contain an element 'n' for 'nth'",
        fixed = TRUE
    )
})

test_that("gtrack.liftover aggregation with non-consecutive overlapping chains", {
    local_db_state()

    # Test aggregation when overlapping chains are separated by non-overlapping ones
    # This tests that the aggregation finds ALL overlapping chains, not just consecutive ones

    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)
    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create track with interval that overlaps first two chains (but not the third)
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = "chrsource1",
        start = 20,
        end = 30,
        stringsAsFactors = FALSE
    )
    gtrack.create_sparse("nc_overlap_src", "Non-consecutive overlap", src_intervals, 100)
    src_track_dir <- file.path(source_db, "tracks", "nc_overlap_src.track")

    setup_db(list(
        ">chr1\n", paste(rep("A", 100), collapse = ""), "\n",
        ">chr2\n", paste(rep("C", 100), collapse = ""), "\n",
        ">chr3\n", paste(rep("G", 100), collapse = ""), "\n"
    ))

    chain_file <- new_chain_file()
    # Chains with same start (will be consecutive in sorted array)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 10, 40, "chr1", 100, "+", 0, 30, 1)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 10, 40, "chr2", 100, "+", 0, 30, 2)
    # Non-overlapping chain that separates the above in iteration order
    write_chain_entry(chain_file, "chrsource1", 200, "+", 50, 60, "chr3", 100, "+", 0, 10, 3)

    lifted_track <- "nc_overlap_lifted"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    # Should find both overlapping chains and aggregate (count=2)
    gtrack.liftover(
        lifted_track, "test", src_track_dir, chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "keep",
        multi_target_agg = "count"
    )

    result <- gextract(lifted_track, gintervals.all())

    # Both chr1 and chr2 should appear (the source interval overlaps both chains)
    expect_true("chr1" %in% result$chrom)
    expect_true("chr2" %in% result$chrom)
    expect_false("chr3" %in% result$chrom)
})

test_that("gtrack.liftover aggregation with chains of varying lengths starting at same position", {
    local_db_state()

    # Multiple chains starting at same source position but with different lengths
    # Tests deterministic ordering and correct aggregation

    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 300), collapse = ""), "\n", sep = "", file = source_fasta)
    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Interval that overlaps all chains
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = "chrsource1",
        start = 105,
        end = 115,
        stringsAsFactors = FALSE
    )
    gtrack.create_sparse("varylen_src", "Varying length", src_intervals, 50)
    src_track_dir <- file.path(source_db, "tracks", "varylen_src.track")

    setup_db(list(
        ">chr1\n", paste(rep("T", 200), collapse = ""), "\n",
        ">chr2\n", paste(rep("C", 200), collapse = ""), "\n",
        ">chr3\n", paste(rep("G", 200), collapse = ""), "\n"
    ))

    chain_file <- new_chain_file()
    # All start at 100 but have different lengths
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 300, "+", 100, 120, "chr1", 200, "+", 0, 20, 1) # len=20
    write_chain_entry(chain_file, "chrsource1", 300, "+", 100, 103, "chr2", 200, "+", 0, 3, 2) # len=3 (won't overlap interval at 105)
    write_chain_entry(chain_file, "chrsource1", 300, "+", 100, 130, "chr3", 200, "+", 0, 30, 3) # len=30

    lifted_track <- "varylen_lifted"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    # Our interval [105,115) overlaps chains 1 and 3, but not 2 (which ends at 103)
    gtrack.liftover(
        lifted_track, "test", src_track_dir, chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "keep",
        multi_target_agg = "count"
    )

    result <- gextract(lifted_track, gintervals.all())

    # Should get 2 results (chr1 and chr3)
    expect_equal(nrow(result), 2)
    expect_true("chr1" %in% result$chrom)
    expect_true("chr3" %in% result$chrom)
    expect_false("chr2" %in% result$chrom)
})

test_that("gtrack.liftover aggregation with reverse strand preserves correct values", {
    local_db_state()

    # Test that aggregation works correctly with reverse strand mappings

    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)
    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create track with three distinct values
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 3),
        start = c(10, 50, 100),
        end = c(20, 60, 110),
        stringsAsFactors = FALSE
    )
    src_vals <- c(10, 20, 30)
    gtrack.create_sparse("rev_strand_src", "Reverse strand", src_intervals, src_vals)
    src_track_dir <- file.path(source_db, "tracks", "rev_strand_src.track")

    setup_db(list(">chr1\n", paste(rep("T", 300), collapse = ""), "\n"))

    chain_file <- new_chain_file()
    # Map with reverse strand
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 150, "chr1", 300, "-", 150, 300, 1)

    lifted_track <- "rev_strand_lifted"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    # Test sum aggregation
    gtrack.liftover(
        lifted_track, "test", src_track_dir, chain_file,
        tgt_overlap_policy = "keep",
        multi_target_agg = "sum"
    )

    result <- gextract(lifted_track, gintervals.all())

    # All three intervals should be mapped
    expect_equal(nrow(result), 3)
    # Sum of all values in result should equal sum of source values
    expect_equal(sum(result[[lifted_track]]), sum(src_vals))
})

test_that("gtrack.liftover aggregation finds earlier long overlap when hint is to the right", {
    local_db_state()

    # Tests the prefix-max fallback logic for finding earlier overlapping chains
    # when the hint points past them

    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)
    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Two intervals: Q1 and Q2
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = "chrsource1",
        start = c(70, 90),
        end = c(71, 91),
        stringsAsFactors = FALSE
    )
    src_vals <- c(100, 200)
    gtrack.create_sparse("long_overlap_src", "Long overlap", src_intervals, src_vals)
    src_track_dir <- file.path(source_db, "tracks", "long_overlap_src.track")

    setup_db(list(
        ">chr1\n", paste(rep("A", 200), collapse = ""), "\n",
        ">chr2\n", paste(rep("C", 200), collapse = ""), "\n",
        ">chr3\n", paste(rep("G", 200), collapse = ""), "\n"
    ))

    chain_file <- new_chain_file()
    # A: [0,100) -> chr1 (overlaps both Q1 and Q2)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 100, "chr1", 200, "+", 0, 100, 1)
    # B: [15,16) -> chr2 (does not overlap Q1 or Q2)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 15, 16, "chr2", 200, "+", 0, 1, 2)
    # C: [80,110) -> chr3 (overlaps Q2 only)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 80, 110, "chr3", 200, "+", 0, 30, 3)

    lifted_track <- "long_overlap_lifted"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(
        lifted_track, "test", src_track_dir, chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "keep",
        multi_target_agg = "sum"
    )

    result <- gextract(lifted_track, gintervals.all())

    # Q1 (val=100) should map to chr1 only: 1 occurrence
    # Q2 (val=200) should map to chr1 AND chr3: 2 occurrences
    # Total: 3 rows
    expect_equal(nrow(result), 3)
    expect_equal(sum(result[[lifted_track]] == 100), 1)
    expect_equal(sum(result[[lifted_track]] == 200), 2)
    expect_true("chr1" %in% result$chrom[result[[lifted_track]] == 200])
    expect_true("chr3" %in% result$chrom[result[[lifted_track]] == 200])
})

test_that("gtrack.liftover aggregation with dense cluster of same-start chains", {
    local_db_state()

    # Many chains starting at same position with varying lengths
    # Tests that aggregation handles large numbers of contributors correctly

    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 1000), collapse = ""), "\n", sep = "", file = source_fasta)
    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Interval that overlaps many chains
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = "chrsource1",
        start = 60,
        end = 61,
        stringsAsFactors = FALSE
    )
    gtrack.create_sparse("dense_cluster_src", "Dense cluster", src_intervals, 77)
    src_track_dir <- file.path(source_db, "tracks", "dense_cluster_src.track")

    setup_db(list(
        ">chr1\n", paste(rep("A", 5000), collapse = ""), "\n",
        ">chr2\n", paste(rep("C", 5000), collapse = ""), "\n"
    ))

    chain_file <- new_chain_file()
    # Create 50 chains all starting at source1[50] with varying lengths
    for (i in 1:50) {
        len <- i + 10 # lengths 11 to 60
        target_chrom <- if (i %% 2 == 0) "chr1" else "chr2"
        # Chain file uses "chrsource1" to match database chromosome name
        write_chain_entry(
            chain_file, "chrsource1", 1000, "+", 50, 50 + len,
            target_chrom, 5000, "+", i * 20, i * 20 + len, i
        )
    }

    lifted_track <- "dense_cluster_lifted"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    # Query at 60 overlaps chains where 50+length > 60, i.e., length > 10
    # All 50 chains have length >= 11, so all should overlap
    gtrack.liftover(
        lifted_track, "test", src_track_dir, chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "keep",
        multi_target_agg = "mean"
    )

    result <- gextract(lifted_track, gintervals.all())

    # Should find many overlapping chains
    expect_true(nrow(result) >= 40) # At least most chains
    expect_true(nrow(result) <= 50) # At most all chains

    # All values should be 77 (mean of 77 is 77)
    expect_true(all(result[[lifted_track]] == 77))
})

test_that("gtrack.liftover aggregation with partial overlaps and varying coverage", {
    local_db_state()

    # Test coverage-based aggregators with varying overlap lengths
    # Create source intervals that will have different coverage lengths when lifted

    # Filename must start with "chr" for .gseq.import() to process it
    source_fasta <- file.path(tempdir(), "chrsource1.fasta")
    cat(">source1\n", paste(rep("A", 200), collapse = ""), "\n", sep = "", file = source_fasta)
    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta)
    gdb.init(source_db)

    # Create track with intervals that will map with different coverage
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 3),
        start = c(10, 40, 70),
        end = c(40, 70, 100), # all 30bp long
        stringsAsFactors = FALSE
    )
    src_vals <- c(100, 200, 300)
    gtrack.create_sparse("partial_overlap_src", "Partial overlap", src_intervals, src_vals)
    src_track_dir <- file.path(source_db, "tracks", "partial_overlap_src.track")

    setup_db(list(">chr1\n", paste(rep("C", 300), collapse = ""), "\n"))

    chain_file <- new_chain_file()
    # Create overlapping chains that will create different coverage lengths
    # Chain 1: source1[0,50) -> chr1[0,50) - will cover all of interval 1 (30bp) and part of interval 2
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 50, "chr1", 300, "+", 0, 50, 1)
    # Chain 2: chrsource1[30,80) -> chr1[30,80) - overlaps all three source intervals partially
    write_chain_entry(chain_file, "chrsource1", 200, "+", 30, 80, "chr1", 300, "+", 30, 80, 2)
    # Chain 3: chrsource1[60,120) -> chr1[60,120) - will cover all of interval 3 (30bp) and part of interval 2
    write_chain_entry(chain_file, "chrsource1", 200, "+", 60, 120, "chr1", 300, "+", 60, 120, 3)

    lifted_track <- "partial_overlap_lifted"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    # Use max.coverage_len - should pick the contributor with longest overlap
    # Use auto policy to handle target overlaps properly
    gtrack.liftover(
        lifted_track, "test", src_track_dir, chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "auto",
        multi_target_agg = "max.coverage_len"
    )

    result <- gextract(lifted_track, gintervals.all())

    # Should have at least one result
    expect_true(nrow(result) > 0)
    # All values should be from our source track
    expect_true(all(result[[lifted_track]] %in% c(100, 200, 300)))
})
