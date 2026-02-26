create_isolated_test_db()

test_that("distance.center works with overlapping input intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Source intervals - these must CONTAIN the query midpoint for distance.center to return non-NaN
    # distance.center returns the distance from the query midpoint to the center of the containing source interval
    src <- rbind(
        gintervals(1, 100, 300), # center at 200
        gintervals(1, 500, 700), # center at 600
        gintervals(1, 1000, 1200) # center at 1100
    )
    gvtrack.create("dc_overlap", src, "distance.center")

    # Overlapping INPUT regions - these overlap with each other
    # This is the key: overlapping regions cause backward access in the iterator
    regions <- rbind(
        gintervals(1, 50, 350), # overlaps with next
        gintervals(1, 200, 750), # overlaps with prev and next
        gintervals(1, 600, 1250) # overlaps with prev
    )

    # Extract with small iterator so bins go backward between overlapping regions
    res_overlap <- gextract("dc_overlap", regions, iterator = 20)

    # Extract one region at a time (no backward access possible)
    res_individual <- do.call(rbind, lapply(1:nrow(regions), function(i) {
        gextract("dc_overlap", regions[i, , drop = FALSE], iterator = 20)
    }))

    # Results should be identical - merge by intervalID to properly compare
    merged <- merge(res_overlap, res_individual, by = c("chrom", "start", "end"), suffixes = c(".overlap", ".individual"))
    expect_equal(merged$dc_overlap.overlap, merged$dc_overlap.individual)
})

test_that("distance.center works with heavily overlapping input intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # More source intervals spread across the region
    src <- rbind(
        gintervals(1, 0, 200), # center at 100
        gintervals(1, 400, 600), # center at 500
        gintervals(1, 800, 1000), # center at 900
        gintervals(1, 1200, 1400), # center at 1300
        gintervals(1, 1800, 2000) # center at 1900
    )
    gvtrack.create("dc_heavy", src, "distance.center")

    # Many heavily overlapping input regions
    regions <- rbind(
        gintervals(1, 0, 500),
        gintervals(1, 100, 700),
        gintervals(1, 300, 900),
        gintervals(1, 500, 1100),
        gintervals(1, 700, 1500),
        gintervals(1, 1000, 2000)
    )

    res_overlap <- gextract("dc_heavy", regions, iterator = 20)

    res_individual <- do.call(rbind, lapply(1:nrow(regions), function(i) {
        gextract("dc_heavy", regions[i, , drop = FALSE], iterator = 20)
    }))

    merged <- merge(res_overlap, res_individual, by = c("chrom", "start", "end"), suffixes = c(".overlap", ".individual"))
    expect_equal(merged$dc_heavy.overlap, merged$dc_heavy.individual)
})

test_that("distance.center works with overlapping input intervals across multiple chromosomes", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    genome <- gintervals.all()
    skip_if(nrow(genome) < 2, "requires at least two chromosomes")

    chrom1 <- as.character(genome$chrom[1])
    chrom2 <- as.character(genome$chrom[2])

    # Source intervals on two chromosomes
    src <- rbind(
        gintervals(chrom1, 100, 300), # center at 200
        gintervals(chrom1, 500, 700), # center at 600
        gintervals(chrom1, 1000, 1200), # center at 1100
        gintervals(chrom2, 200, 400), # center at 300
        gintervals(chrom2, 800, 1000) # center at 900
    )
    gvtrack.create("dc_multi_chrom", src, "distance.center")

    # Overlapping input regions within each chromosome
    regions <- rbind(
        gintervals(chrom1, 50, 350),
        gintervals(chrom1, 200, 750),
        gintervals(chrom1, 600, 1250),
        gintervals(chrom2, 100, 500),
        gintervals(chrom2, 300, 1050)
    )

    res_overlap <- gextract("dc_multi_chrom", regions, iterator = 20)

    res_individual <- do.call(rbind, lapply(1:nrow(regions), function(i) {
        gextract("dc_multi_chrom", regions[i, , drop = FALSE], iterator = 20)
    }))

    merged <- merge(res_overlap, res_individual, by = c("chrom", "start", "end"), suffixes = c(".overlap", ".individual"))
    expect_equal(merged$dc_multi_chrom.overlap, merged$dc_multi_chrom.individual)
})
