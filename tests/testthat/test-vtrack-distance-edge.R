create_isolated_test_db()

test_that("distance.edge computes edge-to-edge distance correctly", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Source intervals
    src <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 500, 600)
    )
    gvtrack.create("dist_edge", src, "distance.edge")

    # Query intervals with known distances
    query <- rbind(
        gintervals(1, 150, 160), # overlaps first -> distance = 0
        gintervals(1, 250, 300), # between: dist to [100,200) = 50, dist to [500,600) = 200
        gintervals(1, 700, 800) # after second -> distance = 100
    )
    res <- gextract("dist_edge", query, iterator = query)

    expect_equal(res$dist_edge, c(0, 50, 100))
})

test_that("distance.edge respects strand for sign", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Source with positive strand
    src_pos <- gintervals(1, 500, 600)
    src_pos$strand <- 1
    gvtrack.create("dist_edge_pos", src_pos, "distance.edge")

    # Source with negative strand
    src_neg <- gintervals(1, 500, 600)
    src_neg$strand <- -1
    gvtrack.create("dist_edge_neg", src_neg, "distance.edge")

    query <- rbind(
        gintervals(1, 100, 200), # before source
        gintervals(1, 700, 800) # after source
    )

    res_pos <- gextract("dist_edge_pos", query, iterator = query)
    res_neg <- gextract("dist_edge_neg", query, iterator = query)

    # + strand: before source = negative, after source = positive
    expect_equal(res_pos$dist_edge_pos[1], -300) # 100-200 to 500-600: -300
    expect_equal(res_pos$dist_edge_pos[2], 100) # 700-800 to 500-600: +100

    # - strand: signs are flipped
    expect_equal(res_neg$dist_edge_neg[1], 300)
    expect_equal(res_neg$dist_edge_neg[2], -100)
})

test_that("distance.edge returns unsigned when no strand", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Source without strand column
    src <- gintervals(1, 500, 600)
    gvtrack.create("dist_edge_unsigned", src, "distance.edge")

    query <- rbind(
        gintervals(1, 100, 200), # before source
        gintervals(1, 700, 800) # after source
    )
    res <- gextract("dist_edge_unsigned", query, iterator = query)

    # Without strand, distances should be absolute values
    expect_equal(res$dist_edge_unsigned, c(300, 100))
})

test_that("distance.edge finds nearest among multiple intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 400, 500),
        gintervals(1, 900, 1000)
    )
    gvtrack.create("dist_edge_multi", src, "distance.edge")

    # Query between intervals
    query <- gintervals(1, 300, 350) # dist to [100,200) = 100, dist to [400,500) = 50
    res <- gextract("dist_edge_multi", query, iterator = query)

    expect_equal(res$dist_edge_multi, 50) # nearest is [400,500)
})

test_that("distance.edge returns NA for empty chromosome", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- gintervals(1, 100, 200)
    gvtrack.create("dist_edge_na", src, "distance.edge")

    query <- gintervals(2, 100, 200) # different chromosome
    res <- gextract("dist_edge_na", query, iterator = query)

    expect_true(is.na(res$dist_edge_na))
})

test_that("distance.edge honors iterator modifiers", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- gintervals(1, 500, 600)
    gvtrack.create("dist_edge_shift", src, "distance.edge")
    gvtrack.iterator("dist_edge_shift", sshift = -100, eshift = 100)

    # Query [300,400) with shift becomes [200,500) which touches [500,600)
    query <- gintervals(1, 300, 400)
    res <- gextract("dist_edge_shift", query, iterator = query)

    # After shift, [200,500) touches [500,600) => distance = 0 (like gintervals.neighbors)
    expect_equal(res$dist_edge_shift, 0)
})

test_that("distance.edge matches gintervals.neighbors", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 500, 600),
        gintervals(1, 900, 1000)
    )
    src$strand <- c(1, -1, 1)

    gvtrack.create("dist_edge_verify", src, "distance.edge")

    query <- rbind(
        gintervals(1, 50, 80),
        gintervals(1, 300, 400),
        gintervals(1, 1100, 1200)
    )

    res_vtrack <- gextract("dist_edge_verify", query, iterator = query)

    # Verify against gintervals.neighbors
    for (i in seq_len(nrow(query))) {
        neighbors <- gintervals.neighbors(query[i, ], src, maxneighbors = 1)
        if (!is.null(neighbors) && nrow(neighbors) > 0) {
            expect_equal(res_vtrack$dist_edge_verify[i], neighbors$dist[1])
        }
    }
})

test_that("distance.edge rejects parameters", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- gintervals(1, 100, 200)
    expect_error(gvtrack.create("dist_edge_param", src, "distance.edge", params = 100))
})

test_that("distance.edge works with overlapping source intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Overlapping sources are allowed
    src <- rbind(
        gintervals(1, 100, 300),
        gintervals(1, 200, 400)
    )
    gvtrack.create("dist_edge_overlap", src, "distance.edge")

    query <- rbind(
        gintervals(1, 150, 160), # overlaps both -> distance = 0
        gintervals(1, 500, 600) # distance to [200,400) = 100
    )
    res <- gextract("dist_edge_overlap", query, iterator = query)

    expect_equal(res$dist_edge_overlap, c(0, 100))
})

test_that("distance.edge handles touching intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- gintervals(1, 100, 200)
    gvtrack.create("dist_edge_touch", src, "distance.edge")

    # Query that exactly touches source
    query <- gintervals(1, 200, 300) # query starts where source ends
    res <- gextract("dist_edge_touch", query, iterator = query)

    # Touching intervals have distance 0 (like gintervals.neighbors)
    expect_equal(res$dist_edge_touch, 0)
})

test_that("distance.edge matches neighbors with maxneighbors=1", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create diverse set of source intervals
    src <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 500, 600),
        gintervals(1, 900, 1000),
        gintervals(2, 50, 150),
        gintervals(2, 300, 400)
    )

    gvtrack.create("dist_cmp", src, "distance.edge")

    # Test various query intervals
    query <- rbind(
        gintervals(1, 50, 80), # before first on chr1
        gintervals(1, 150, 160), # overlaps first
        gintervals(1, 250, 300), # between intervals
        gintervals(1, 700, 800), # between intervals
        gintervals(1, 1100, 1200), # after last on chr1
        gintervals(2, 10, 20), # before first on chr2
        gintervals(2, 200, 250), # between on chr2
        gintervals(2, 500, 600) # after last on chr2
    )

    res_vtrack <- gextract("dist_cmp", query, iterator = query)

    # Compare with gintervals.neighbors
    for (i in seq_len(nrow(query))) {
        neighbors <- gintervals.neighbors(query[i, ], src, maxneighbors = 1)
        vtrack_dist <- res_vtrack$dist_cmp[i]

        if (is.null(neighbors) || nrow(neighbors) == 0) {
            expect_true(is.na(vtrack_dist), info = paste("Row", i))
        } else {
            expect_equal(vtrack_dist, neighbors$dist[1], info = paste("Row", i))
        }
    }
})

test_that("distance.edge matches neighbors with strand variations", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Test with different strand values
    src <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 500, 600),
        gintervals(1, 900, 1000)
    )
    src$strand <- c(1, -1, 1) # Mix of strands

    gvtrack.create("dist_strand", src, "distance.edge")

    query <- rbind(
        gintervals(1, 50, 80), # -20 for +strand, +20 for -strand
        gintervals(1, 250, 300), # between
        gintervals(1, 650, 700), # between
        gintervals(1, 1050, 1100) # after
    )

    res_vtrack <- gextract("dist_strand", query, iterator = query)

    for (i in seq_len(nrow(query))) {
        neighbors <- gintervals.neighbors(query[i, ], src, maxneighbors = 1)
        expect_equal(res_vtrack$dist_strand[i], neighbors$dist[1], info = paste("Row", i))
    }
})

test_that("distance.edge matches neighbors with use_intervals1_strand=TRUE", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Source intervals with strand
    src <- rbind(
        gintervals(1, 500, 600),
        gintervals(1, 900, 1000)
    )

    # Query intervals with strand
    query <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 700, 800)
    )
    query$strand <- c(1, -1)

    gvtrack.create("dist_query", src, "distance.edge")

    res_vtrack <- gextract("dist_query", query, iterator = query)

    # Note: distance.edge uses SOURCE strand (like use_intervals1_strand=FALSE)
    # This is the expected behavior based on "like other distance vtracks"
    for (i in seq_len(nrow(query))) {
        neighbors <- gintervals.neighbors(query[i, ], src, maxneighbors = 1, use_intervals1_strand = FALSE, warn.ignored.strand = FALSE)
        expect_equal(res_vtrack$dist_query[i], neighbors$dist[1], info = paste("Row", i))
    }
})

test_that("distance.edge matches neighbors with distance ranges", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 500, 600),
        gintervals(1, 900, 1000)
    )

    gvtrack.create("dist_range", src, "distance.edge")

    # Test query that is equidistant from two intervals
    # Should return the one with lower coordinate (implementation detail)
    query <- gintervals(1, 350, 351) # Equidistant from [100,200) and [500,600)

    res_vtrack <- gextract("dist_range", query, iterator = query)
    neighbors <- gintervals.neighbors(query, src, maxneighbors = 1)

    expect_equal(res_vtrack$dist_range, neighbors$dist[1])
})

test_that("distance.edge matches neighbors with multiple chromosomes", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 500, 600),
        gintervals(2, 100, 200),
        gintervals(2, 500, 600),
        gintervals("X", 100, 200)
    )

    gvtrack.create("dist_multi_chr", src, "distance.edge")

    query <- rbind(
        gintervals(1, 300, 400),
        gintervals(2, 300, 400),
        gintervals("X", 50, 80),
        gintervals("Y", 100, 200) # No source intervals on Y
    )

    res_vtrack <- gextract("dist_multi_chr", query, iterator = query)

    for (i in seq_len(nrow(query))) {
        neighbors <- gintervals.neighbors(query[i, ], src, maxneighbors = 1)
        vtrack_dist <- res_vtrack$dist_multi_chr[i]

        if (is.null(neighbors) || nrow(neighbors) == 0) {
            expect_true(is.na(vtrack_dist), info = paste("Row", i, query$chrom[i]))
        } else {
            expect_equal(vtrack_dist, neighbors$dist[1], info = paste("Row", i, query$chrom[i]))
        }
    }
})

test_that("distance.edge matches neighbors with overlapping intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create overlapping source intervals
    src <- rbind(
        gintervals(1, 100, 300),
        gintervals(1, 200, 400),
        gintervals(1, 250, 350)
    )

    gvtrack.create("dist_overlap_src", src, "distance.edge")

    query <- rbind(
        gintervals(1, 50, 80), # Before all
        gintervals(1, 150, 180), # Overlaps first
        gintervals(1, 270, 280), # Overlaps all three
        gintervals(1, 320, 380), # Overlaps middle two
        gintervals(1, 500, 600) # After all
    )

    res_vtrack <- gextract("dist_overlap_src", query, iterator = query)

    for (i in seq_len(nrow(query))) {
        neighbors <- gintervals.neighbors(query[i, ], src, maxneighbors = 1)
        expect_equal(res_vtrack$dist_overlap_src[i], neighbors$dist[1], info = paste("Row", i))
    }
})

test_that("distance.edge matches neighbors with na.if.notfound behavior", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- gintervals(1, 100, 200)
    gvtrack.create("dist_na", src, "distance.edge")

    # Query on different chromosome
    query <- gintervals(2, 100, 200)

    res_vtrack <- gextract("dist_na", query, iterator = query)
    neighbors <- gintervals.neighbors(query, src, maxneighbors = 1, na.if.notfound = TRUE)

    # Both should return NA
    expect_true(is.na(res_vtrack$dist_na))
    expect_true(is.na(neighbors$dist))
})

test_that("distance.edge matches neighbors with very close intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Test edge cases with intervals 1bp apart
    src <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 201, 300), # 1bp gap
        gintervals(1, 500, 600)
    )

    gvtrack.create("dist_close", src, "distance.edge")

    query <- rbind(
        gintervals(1, 200, 201), # Exactly in the 1bp gap
        gintervals(1, 300, 301), # Right after second interval
        gintervals(1, 400, 450) # Between second and third
    )

    res_vtrack <- gextract("dist_close", query, iterator = query)

    for (i in seq_len(nrow(query))) {
        neighbors <- gintervals.neighbors(query[i, ], src, maxneighbors = 1)
        expect_equal(res_vtrack$dist_close[i], neighbors$dist[1], info = paste("Row", i))
    }
})
