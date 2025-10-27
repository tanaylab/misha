test_that("near.count counts nearby intervals correctly", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- rbind(
        gintervals(1, 100, 110),
        gintervals(1, 300, 320),
        gintervals(1, 305, 320)
    )
    gvtrack.create("near10_simple", src, "near.count", 10)

    query <- rbind(
        gintervals(1, 90, 100),
        gintervals(1, 295, 305),
        gintervals(1, 600, 700)
    )
    res <- gextract("near10_simple", query, iterator = query)

    expect_equal(res$near10_simple, c(1, 2, 0))
})

test_that("near.count defaults to zero distance", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- gintervals(1, 100, 110)
    gvtrack.create("near0_default", src, "near.count")

    query <- rbind(
        gintervals(1, 95, 110),
        gintervals(1, 120, 140)
    )
    res <- gextract("near0_default", query, iterator = query)

    expect_equal(res$near0_default, c(1, 0))
})

test_that("near.count counts overlapping sources separately", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- rbind(
        gintervals(1, 100, 120),
        gintervals(1, 105, 125)
    )
    gvtrack.create("near_multi", src, "near.count", 5)

    res <- gextract("near_multi", gintervals(1, 110, 115), iterator = gintervals(1, 110, 115))

    expect_equal(res$near_multi, 2)
})

test_that("near.count honors iterator modifiers", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- gintervals(1, 100, 110)
    gvtrack.create("near_shift", src, "near.count", 0)

    query <- gintervals(1, 200, 210)
    res <- gextract("near_shift", query, iterator = query)
    expect_equal(res$near_shift, 0)

    gvtrack.iterator("near_shift", sshift = -100, eshift = -100)
    res_shifted <- gextract("near_shift", query, iterator = query)
    expect_equal(res_shifted$near_shift, 1)
})

test_that("near.count respects chromosome boundaries and separation", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    genome <- gintervals.all()
    skip_if(nrow(genome) < 2, "requires at least two chromosomes")
    chrom1 <- genome[1, ]
    chrom_label <- as.character(chrom1$chrom)
    chrom_start <- as.numeric(chrom1$start)
    chrom_end <- as.numeric(chrom1$end)
    other_chrom_label <- as.character(genome$chrom[2])

    src <- rbind(
        gintervals(chrom_label, chrom_start, chrom_start + 5),
        gintervals(chrom_label, chrom_end - 5, chrom_end),
        gintervals(other_chrom_label, 200, 220)
    )
    gvtrack.create("near_bounds", src, "near.count", 50)

    query <- rbind(
        gintervals(chrom_label, chrom_start, chrom_start + 10),
        gintervals(chrom_label, chrom_end - 10, chrom_end),
        gintervals(other_chrom_label, 0, 100)
    )
    res <- gextract("near_bounds", query, iterator = query)

    expect_equal(res$near_bounds, c(1, 1, 0))
})

test_that("near.count rejects invalid parameters", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- gintervals(1, 0, 10)

    expect_error(gvtrack.create("near_bad_neg", src, "near.count", -1))
    expect_error(gvtrack.create("near_bad_vec", src, "near.count", c(1, 2)))
    expect_error(gvtrack.create("near_bad_type", src, "near.count", "foo"))
})

test_that("near.count handles touching intervals correctly", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Source intervals that touch each other
    src <- rbind(
        gintervals(1, 100, 150),
        gintervals(1, 150, 200),
        gintervals(1, 300, 350)
    )
    gvtrack.create("near_touch", src, "near.count", 1)

    # With params=1: touching intervals (distance=0) count due to expansion overlap
    # In half-open coordinates, distance=1 doesn't count (expanded intervals touch but don't overlap)
    query <- rbind(
        gintervals(1, 150, 160), # overlaps [150,200), also [100,150) expanded to [99,151) overlaps → count=2
        gintervals(1, 200, 250), # distance=0 from [150,200), expansion makes them overlap → count=1
        gintervals(1, 350, 400), # distance=0 from [300,350), expansion makes them overlap → count=1
        gintervals(1, 351, 400) # distance=1 from [300,350), [299,351) touches [351,400) but no overlap → count=0
    )
    res <- gextract("near_touch", query, iterator = query)

    expect_equal(res$near_touch, c(2, 1, 1, 0))
})

test_that("near.count with varying distances", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- gintervals(1, 500, 600)

    # Test different distance thresholds (distance <= params, but half-open overlap semantics apply)
    gvtrack.create("near_dist5", src, "near.count", 5)
    gvtrack.create("near_dist20", src, "near.count", 20)
    gvtrack.create("near_dist101", src, "near.count", 101)

    query <- rbind(
        gintervals(1, 300, 400), # distance=100 from [500,600)
        gintervals(1, 480, 490), # distance=10 from [500,600)
        gintervals(1, 700, 800) # distance=100 from [500,600)
    )

    res5 <- gextract("near_dist5", query, iterator = query)
    res20 <- gextract("near_dist20", query, iterator = query)
    res101 <- gextract("near_dist101", query, iterator = query)

    # distance=100: needs params>=101 due to half-open interval overlap semantics
    # distance=10: counts with params>=11
    expect_equal(res5$near_dist5, c(0, 0, 0))
    expect_equal(res20$near_dist20, c(0, 1, 0))
    expect_equal(res101$near_dist101, c(1, 1, 1))
})

test_that("near.count handles completely overlapping query", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 150, 250),
        gintervals(1, 400, 500)
    )
    gvtrack.create("near_overlap", src, "near.count", 5)

    # Query completely overlaps first two sources
    query <- gintervals(1, 50, 300)
    res <- gextract("near_overlap", query, iterator = query)

    expect_equal(res$near_overlap, 2)
})

test_that("near.count with nested source intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Nested intervals
    src <- rbind(
        gintervals(1, 100, 500),
        gintervals(1, 200, 300),
        gintervals(1, 250, 275)
    )
    gvtrack.create("near_nest", src, "near.count", 0)

    query <- rbind(
        gintervals(1, 225, 260), # overlaps all three
        gintervals(1, 350, 400), # overlaps only first
        gintervals(1, 600, 700) # overlaps none
    )
    res <- gextract("near_nest", query, iterator = query)

    expect_equal(res$near_nest, c(3, 1, 0))
})

test_that("near.count with empty source", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Empty source - create empty data frame
    src <- data.frame(chrom = character(0), start = numeric(0), end = numeric(0))
    gvtrack.create("near_empty_src", src, "near.count", 10)

    query <- gintervals(1, 100, 200)
    res <- gextract("near_empty_src", query, iterator = query)

    expect_equal(res$near_empty_src, 0)
})

test_that("near.count across multiple chromosomes", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    genome <- gintervals.all()
    skip_if(nrow(genome) < 3, "requires at least three chromosomes")

    chrom1 <- as.character(genome$chrom[1])
    chrom2 <- as.character(genome$chrom[2])
    chrom3 <- as.character(genome$chrom[3])

    src <- rbind(
        gintervals(chrom1, 100, 200),
        gintervals(chrom1, 500, 600),
        gintervals(chrom2, 100, 200),
        gintervals(chrom3, 100, 200)
    )
    gvtrack.create("near_multi_chrom", src, "near.count", 51)

    query <- rbind(
        gintervals(chrom1, 150, 160), # overlaps 1 on chrom1
        gintervals(chrom1, 400, 450), # distance=50 from [500,600), expansion to [449,651) overlaps
        gintervals(chrom2, 50, 90), # distance=10 from [100,200), expansion to [49,251) overlaps
        gintervals(chrom3, 300, 400), # distance=100 from [100,200), expansion to [49,251) doesn't reach
        gintervals(chrom2, 150, 160) # overlaps 1 on chrom2
    )
    res <- gextract("near_multi_chrom", query, iterator = query)

    # Results in genomic order: chrom1[150,160], chrom1[400,450], chrom2[50,90], chrom2[150,160], chrom3[300,400]
    expect_equal(res$near_multi_chrom, c(1, 1, 1, 1, 0))
})

test_that("near.count with very large distance", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- rbind(
        gintervals(1, 1000, 2000),
        gintervals(1, 100000, 110000),
        gintervals(1, 500000, 510000)
    )
    gvtrack.create("near_large_dist", src, "near.count", 100000)

    query <- gintervals(1, 50000, 60000)
    res <- gextract("near_large_dist", query, iterator = query)

    # Query [50000, 60000): distance=48000 from [1000,2000), distance=40000 from [100000,110000)
    # Both distances <= 100000, so count=2
    expect_equal(res$near_large_dist, 2)
})

test_that("near.count with half-open coordinate edge cases", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- gintervals(1, 100, 200)
    gvtrack.create("near_exact", src, "near.count", 50)

    query <- rbind(
        gintervals(1, 250, 300), # distance=50, [100,200) expanded to [50,250) touches but doesn't overlap
        gintervals(1, 249, 300), # distance=49, [100,200) expanded to [50,250) overlaps with [249,300)
        gintervals(1, 251, 300), # distance=51 > 50, doesn't count
        gintervals(1, 0, 50), # distance=50, touches expanded interval but doesn't overlap
        gintervals(1, 0, 51) # distance=49, overlaps expanded interval [50,250)
    )
    res <- gextract("near_exact", query, iterator = query)

    # Results in genomic order: [0,50), [0,51), [249,300), [250,300), [251,300)
    # With half-open coordinates: distance=50 doesn't count (touching ≠ overlapping), distance=49 counts
    expect_equal(res$near_exact, c(0, 1, 1, 0, 0))
})

test_that("near.count with many nearby sources", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create many nearby intervals
    starts <- seq(100, 1000, by = 100)
    src <- gintervals(1, starts, starts + 50)
    gvtrack.create("near_many", src, "near.count", 60)

    query <- rbind(
        gintervals(1, 500, 510), # should overlap and be near several
        gintervals(1, 50, 60), # near first
        gintervals(1, 2000, 2100) # far from all
    )
    res <- gextract("near_many", query, iterator = query)

    # Results sorted by genomic position:
    # Query [50,60): distance=40 to [100,150) → count=1
    # Query [500,510): overlaps [500,550), distance=50 to [400,450) → count=2
    # Query [2000,2100): far from all → count=0
    expect_equal(res$near_many, c(1, 2, 0))
})

test_that("near.count with single base intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Single base source intervals
    src <- rbind(
        gintervals(1, 100, 101),
        gintervals(1, 200, 201),
        gintervals(1, 300, 301)
    )
    gvtrack.create("near_single", src, "near.count", 10)

    query <- rbind(
        gintervals(1, 90, 111), # overlaps first
        gintervals(1, 190, 211), # overlaps second
        gintervals(1, 150, 170) # between first and second, far from both
    )
    res <- gextract("near_single", query, iterator = query)

    # Results are sorted by genomic position
    expect_equal(res$near_single, c(1, 0, 1))
})

test_that("near.count with identical source and query intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    intervals <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 300, 400),
        gintervals(1, 500, 600)
    )

    gvtrack.create("near_identical", intervals, "near.count", 0)

    # Query with same intervals
    res <- gextract("near_identical", intervals, iterator = intervals)

    # Each interval should overlap itself
    expect_equal(res$near_identical, c(1, 1, 1))
})

test_that("near.count distance calculation with gaps", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Source with specific distances between intervals
    src <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 250, 350), # distance=50 from first
        gintervals(1, 500, 600) # distance=150 from second
    )

    gvtrack.create("near_gap25", src, "near.count", 25)
    gvtrack.create("near_gap51", src, "near.count", 51)
    gvtrack.create("near_gap91", src, "near.count", 91)

    # Query in the middle gaps
    query <- rbind(
        gintervals(1, 220, 230), # distance=20 from both [100,200) and [250,350)
        gintervals(1, 400, 410) # distance=50 from [250,350), distance=90 from [500,600)
    )

    res25 <- gextract("near_gap25", query, iterator = query)
    res51 <- gextract("near_gap51", query, iterator = query)
    res91 <- gextract("near_gap91", query, iterator = query)

    # Query 1: distances=20,20 → counts with all params (20<=25, 20<=51, 20<=91)
    # Query 2: distances=50,90 → needs params>=51 for first (half-open semantics), params>=91 for second
    expect_equal(res25$near_gap25, c(2, 0))
    expect_equal(res51$near_gap51, c(2, 1))
    expect_equal(res91$near_gap91, c(2, 2))
})
