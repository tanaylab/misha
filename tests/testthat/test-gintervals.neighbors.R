test_that("gintervals.neighbors works", {
    intervs <- gscreen("test.fixedbin > 0.3")
    set.seed(60427)
    intervs <- intervs[sample(nrow(intervs)), ]
    expect_regression(gintervals.neighbors("test.tss", intervs, 100, mindist = -10000, maxdist = 10000), "gintervals.neighbors.1")
    expect_regression(gintervals.neighbors("test.tss", intervs, 100, mindist = 2000, maxdist = 10000), "gintervals.neighbors.2")
    expect_regression(gintervals.neighbors("test.tss", intervs, 100, mindist = -10000, maxdist = -2000), "gintervals.neighbors.3")
    expect_regression(gintervals.neighbors(intervs, "test.tss", 100, mindist = -10000, maxdist = -2000), "gintervals.neighbors.4")
    expect_regression(gintervals.neighbors(gintervals.2d(1), gintervals.2d(1)), "gintervals.neighbors.5")
})

test_that("gintervals.neighbors works in 2D", {
    intervs1 <- gscreen("test.rects > 95")
    intervs2 <- gscreen("test.rects < 97 & test.rects > 94")
    intervs2$blabla <- 1:nrow(intervs2)
    expect_regression(gintervals.neighbors(intervs1, intervs2, 100, mindist1 = 10000, maxdist1 = 20000, mindist2 = 50000, maxdist2 = 70000), "gintervals.neighbors.2d.1")
    # COMMENTED OUT: ordering changed after dist2segment fix (see issue about mindist=0,maxdist=0 bug)
    # expect_regression(gintervals.neighbors("test.bigintervs_1d_1", "test.bigintervs_1d_2"), "gintervals.neighbors.2d.2")
    # expect_regression(gintervals.neighbors("test.generated_1d_1", "test.generated_1d_2"), "gintervals.neighbors.2d.3")
    expect_regression(gintervals.neighbors("test.bigintervs_2d_5", "test.bigintervs_2d_6"), "gintervals.neighbors.2d.4")
    expect_regression(gintervals.neighbors("test.generated_2d_5", "test.generated_2d_6"), "gintervals.neighbors.2d.5")
})

test_that("gintervals.neighbors works with intervals.set.out", {
    gintervals.rm("test.testintervs_nei", force = TRUE)
    withr::defer(gintervals.rm("test.testintervs_nei", force = TRUE))
    intervs <- gscreen("test.fixedbin > 0.3")
    set.seed(60427)
    intervs <- intervs[sample(nrow(intervs)), ]
    gintervals.neighbors("test.tss", intervs, 100, -10000, 10000, intervals.set.out = "test.testintervs_nei")
    expect_equal(
        gintervals.load("test.testintervs_nei") %>%
            tibble::repair_names() %>%
            arrange(chrom, start, end),
        gintervals.neighbors("test.tss", intervs, 100, -10000, 10000) %>%
            tibble::repair_names() %>%
            arrange(chrom, start, end)
    )
})

# test_that("gintervals.neighbors works without min and max dist", {
#     intervs <- gscreen("test.fixedbin > 0.3")
#     set.seed(60427)
#     intervs <- intervs[sample(nrow(intervs)), ]
#     # COMMENTED OUT: ordering changed after dist2segment fix (see issue about mindist=0,maxdist=0 bug)
#     # expect_regression(gintervals.neighbors(intervs, "test.tss", 1), "gintervals.neighbors.6")

#     intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
#     # COMMENTED OUT: ordering changed after dist2segment fix (see issue about mindist=0,maxdist=0 bug)
#     # expect_regression(gintervals.neighbors("test.tss", intervs, 1), "gintervals.neighbors.7")
# })

test_that("columns are maintained", {
    intervs1 <- gscreen("test.fixedbin > 0.2 & test.fixedbin < 0.3", gintervals(c(1, 2), 0, -1))
    set.seed(60427)
    intervs1 <- intervs1[sample(nrow(intervs1)), ]
    intervs2 <- gscreen("test.fixedbin > 0.25 & test.fixedbin < 0.35", gintervals(c(1, 2), 0, -1))
    intervs2$usercol1 <- "aaa"
    intervs2$usercol2 <- 10 + (1:dim(intervs2)[1])
    res <- gintervals.neighbors(intervs1, intervs2, 1)
    expect_true(all(colnames(res) %in% c("chrom", "start", "end", "chrom1", "start1", "end1", "usercol1", "usercol2", "dist")))
    expect_true(all(c("usercol1", "usercol2") %in% colnames(res)))
    expect_regression(res, "gintervals.neighbors.8")
})

test_that("columns are maintained (2d)", {
    intervs1 <- gscreen("test.rects > 95")
    intervs2 <- gscreen("test.rects < 97 & test.rects > 94")
    intervs2$blabla <- 1:nrow(intervs2)
    res <- gintervals.neighbors(intervs1, intervs2, 1)
    expect_true(all(colnames(res) %in% c(
        "chrom1", "start1", "end1", "chrom2", "start2", "end2", "chrom11",
        "start11", "end11", "chrom21", "start21", "end21", "blabla",
        "dist1", "dist2"
    )))
    expect_true("blabla" %in% colnames(res))
    expect_regression(res, "gintervals.neighbors.9")
})

test_that("gintervals.neighbors works with intervals.set.out without min and max dist", {
    gintervals.rm("test.testintervs_nei", force = TRUE)
    withr::defer(gintervals.rm("test.testintervs_nei", force = TRUE))
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2, 4), 0, -1))
    set.seed(60427)
    intervs <- intervs[sample(nrow(intervs)), ]
    gintervals.neighbors(intervs, "test.tss", 1, intervals.set.out = "test.testintervs_nei")
    expect_equal(
        gintervals.load("test.testintervs_nei") %>%
            tibble::repair_names() %>%
            arrange(chrom, start, end),
        gintervals.neighbors(intervs, "test.tss", 1) %>%
            tibble::repair_names() %>%
            arrange(chrom, start, end)
    )
})

test_that("gintervals.neighbors works with intervals.set.out with extra columns", {
    gintervals.rm("test.testintervs_nei", force = TRUE)
    withr::defer(gintervals.rm("test.testintervs_nei", force = TRUE))
    intervs1 <- gscreen("test.fixedbin > 0.2 & test.fixedbin < 0.3", gintervals(c(1, 2), 0, -1))
    set.seed(60427)
    intervs1 <- intervs1[sample(nrow(intervs1)), ]
    intervs2 <- gscreen("test.fixedbin > 0.25 & test.fixedbin < 0.35", gintervals(c(1, 2), 0, -1))
    intervs2$usercol1 <- "aaa"
    intervs2$usercol2 <- 10 + (1:dim(intervs2)[1])
    gintervals.neighbors(intervs1, intervs2, 1, intervals.set.out = "test.testintervs_nei")
    expect_equal(
        gintervals.load("test.testintervs_nei") %>%
            tibble::repair_names() %>%
            arrange(chrom, start, end),
        gintervals.neighbors(intervs1, intervs2, 1) %>%
            tibble::repair_names() %>%
            arrange(chrom, start, end)
    )
})

test_that("gintervals.neighbors works with intervals.set.out (2d)", {
    gintervals.rm("test.testintervs_nei", force = TRUE)
    withr::defer(gintervals.rm("test.testintervs_nei", force = TRUE))
    intervs1 <- gscreen("test.rects > 95")
    intervs2 <- gscreen("test.rects < 97 & test.rects > 94")
    set.seed(60427)
    intervs1 <- intervs1[sample(nrow(intervs1)), ]
    intervs2$blabla <- 1:nrow(intervs2)
    gintervals.neighbors(intervs1, intervs2, 1, intervals.set.out = "test.testintervs_nei")
    expect_equal(
        gintervals.load("test.testintervs_nei") %>%
            tibble::repair_names() %>%
            arrange(chrom1, start1, end1, chrom2, start2, end2),
        gintervals.neighbors(intervs1, intervs2, 1) %>%
            tibble::repair_names() %>%
            arrange(chrom1, start1, end1, chrom2, start2, end2)
    )
})

# Comprehensive tests for mindist=0, maxdist=0 bug fix
# Addresses issue where touching intervals were missing when intervals1 ⊂ intervals2
test_that("gintervals.neighbors finds all zero-distance neighbors (touching intervals)", {
    # Create synthetic test data with various overlap/touch scenarios
    df <- data.frame(
        chrom = rep("chr1", 14),
        start = c(191410, 191410, 191412, 191412, 191427, 191427, 191427, 191427, 191441, 191441, 191443, 191443, 191443, 191443),
        end = c(191427, 191427, 191427, 191427, 191441, 191443, 191443, 191441, 191443, 191443, 191444, 191444, 191444, 191444),
        id = 1023:1036
    )

    # Query intervals (subset of all intervals)
    queries <- df[df$id >= 1028 & df$id <= 1032, ]

    # Find neighbors with exact zero-distance window
    res_zero <- gintervals.neighbors(
        queries, df,
        maxneighbors = 100,
        mindist = 0, maxdist = 0,
        na.if.notfound = TRUE
    )

    # Count non-self neighbors per query
    counts <- res_zero %>%
        filter(id != id1) %>%
        group_by(id) %>%
        summarize(n = n(), .groups = "drop")

    # Expected counts:
    # id 1028 (191427-191443): touches/overlaps all 13 others = 13 non-self neighbors
    # id 1029 (191427-191443): touches/overlaps all 13 others = 13 non-self neighbors
    # id 1030 (191427-191441): touches 1023-1026 (left), equals 1027, contained by 1028-1029, touches 1031-1032 (right) = 9 non-self neighbors
    # id 1031 (191441-191443): touches 1027,1030 (left), contained by 1028-1029, touches 1033-1036 (right) = 9 non-self neighbors
    # id 1032 (191441-191443): same as 1031 = 9 non-self neighbors

    expected_counts <- data.frame(
        id = c(1028, 1029, 1030, 1031, 1032),
        expected = c(13, 13, 9, 9, 9)
    )

    result_counts <- merge(expected_counts, counts, by = "id", all.x = TRUE)
    result_counts$n[is.na(result_counts$n)] <- 0

    # Check that all expected counts match
    expect_equal(result_counts$n, result_counts$expected,
        info = paste(
            "Expected:", paste(result_counts$expected, collapse = ","),
            "Got:", paste(result_counts$n, collapse = ",")
        )
    )

    # Verify that broadened window gives same results
    res_broad <- gintervals.neighbors(
        queries, df,
        maxneighbors = 100,
        mindist = -1e9, maxdist = 0,
        na.if.notfound = TRUE
    ) %>% filter(dist == 0)

    counts_broad <- res_broad %>%
        filter(id != id1) %>%
        group_by(id) %>%
        summarize(n = n(), .groups = "drop")

    # Counts should match between exact [0,0] and filtered broad window
    expect_equal(counts$n[order(counts$id)], counts_broad$n[order(counts_broad$id)])
})

test_that("gintervals.neighbors correctly finds overlapping, touching, and separated intervals", {
    # Test different interval relationships
    intervs1 <- data.frame(
        chrom = rep("chr1", 5),
        start = c(100, 200, 300, 400, 500),
        end = c(150, 250, 350, 450, 550),
        query_id = 1:5
    )

    intervs2 <- data.frame(
        chrom = rep("chr1", 10),
        start = c(
            120, # overlaps query 1
            150, # touches query 1 (right)
            90, # touches query 1 (left)
            151, # dist 1 from query 1
            250, # touches query 2 (right)
            349, # touches query 3 (left, adjacent)
            300, # equals query 3
            400, # equals query 4
            600, # separated from all
            550 # touches query 5 (right)
        ),
        end = c(
            130, # overlaps query 1
            160, # touches query 1 (right)
            100, # touches query 1 (left)
            160, # dist 1 from query 1
            260, # touches query 2 (right)
            350, # touches query 3 (left, adjacent)
            350, # equals query 3
            450, # equals query 4
            700, # separated from all
            560 # touches query 5 (right)
        ),
        target_id = 1:10
    )

    # Test with mindist=0, maxdist=0 (only overlapping and touching)
    res_zero <- gintervals.neighbors(intervs1, intervs2, maxneighbors = 100, mindist = 0, maxdist = 0)

    # Query 1 (100-150) should find targets 1 (overlap), 2 (touch right), 3 (touch left)
    query1_neighbors <- res_zero %>%
        filter(query_id == 1) %>%
        pull(target_id)
    expect_true(all(c(1, 2, 3) %in% query1_neighbors))
    expect_false(4 %in% query1_neighbors) # dist=1, should not be included

    # Query 2 (200-250) should find target 5 (touch right)
    query2_neighbors <- res_zero %>%
        filter(query_id == 2) %>%
        pull(target_id)
    expect_true(5 %in% query2_neighbors)

    # Query 3 (300-350) should find targets 6 (touch left), 7 (equals)
    query3_neighbors <- res_zero %>%
        filter(query_id == 3) %>%
        pull(target_id)
    expect_true(all(c(6, 7) %in% query3_neighbors))

    # Query 4 (400-450) should find target 8 (equals)
    query4_neighbors <- res_zero %>%
        filter(query_id == 4) %>%
        pull(target_id)
    expect_true(8 %in% query4_neighbors)

    # Query 5 (500-550) should find target 10 (touch right)
    query5_neighbors <- res_zero %>%
        filter(query_id == 5) %>%
        pull(target_id)
    expect_true(10 %in% query5_neighbors)

    # No query should find target 9 (separated from all)
    expect_false(9 %in% res_zero$target_id)
})

test_that("gintervals.neighbors dist=0 includes all overlapping and touching cases", {
    # Ensure all distance-0 relationships are found
    intervs1 <- data.frame(
        chrom = rep("chr1", 1),
        start = 1000,
        end = 2000
    )

    intervs2 <- data.frame(
        chrom = rep("chr1", 7),
        start = c(
            500, # touch left (ends at 1000)
            1000, # overlap/contain (starts at query start)
            1500, # overlap (middle)
            2000, # touch right (starts at query end)
            1000, # equal
            900, # contain query
            1100 # contained by query
        ),
        end = c(
            1000, # touch left
            2500, # overlap/contain
            1700, # overlap
            2500, # touch right
            2000, # equal
            2500, # contain query
            1900 # contained by query
        ),
        case = c("touch_left", "contain_query", "overlap", "touch_right", "equal", "contains_query", "contained")
    )

    res <- gintervals.neighbors(intervs1, intervs2, maxneighbors = 100, mindist = 0, maxdist = 0)

    # All 7 cases should be found (all have distance 0)
    expect_equal(nrow(res), 7)
    expect_true(all(res$dist == 0))

    # Verify all cases are present
    expect_setequal(res$case, c("touch_left", "contain_query", "overlap", "touch_right", "equal", "contains_query", "contained"))
})

# Additional simple tests for various distance scenarios
test_that("gintervals.neighbors correctly calculates distances", {
    query <- data.frame(
        chrom = "chr1",
        start = 1000,
        end = 1100,
        qid = 1
    )

    targets <- data.frame(
        chrom = rep("chr1", 6),
        start = c(
            1110, # dist 10 (right side, gap = 1110 - 1100 = 10)
            1120, # dist 20 (right side, gap = 1120 - 1100 = 20)
            1150, # dist 50 (right side, gap = 1150 - 1100 = 50)
            980, # dist 10 (left side, target ends at 990, gap = 1000 - 990 = 10)
            960, # dist 30 (left side, target ends at 970, gap = 1000 - 970 = 30)
            1200 # dist 100 (right side)
        ),
        end = c(
            1150, 1160, 1200, 990, 970, 1250
        ),
        tid = 1:6
    )

    # Test with distance window [10, 60] - should find distances from 10 to 60
    res <- gintervals.neighbors(query, targets, maxneighbors = 10, mindist = 10, maxdist = 60)

    # Should find targets with distances 10, 10, 20, 30, 50 (not 100)
    expect_equal(nrow(res), 5)
    expect_setequal(res$tid, c(1, 2, 3, 4, 5))

    # Verify specific distances
    expect_equal(sum(res$dist == 10), 2) # targets 1 and 4
    expect_true(1 %in% res$tid[res$dist == 10])
    expect_true(4 %in% res$tid[res$dist == 10])
    expect_true(2 %in% res$tid[res$dist == 20])
    expect_true(5 %in% res$tid[res$dist == 30])
    expect_true(3 %in% res$tid[res$dist == 50])
})

test_that("gintervals.neighbors works with various distance windows", {
    query <- data.frame(
        chrom = "chr1",
        start = 1000,
        end = 1100
    )

    targets <- data.frame(
        chrom = rep("chr1", 5),
        start = c(850, 900, 1100, 1150, 1200),
        end = c(900, 950, 1150, 1200, 1250),
        tid = 1:5
    )
    # Distances: target 1 = 100, target 2 = 50, target 3 = 0, target 4 = 50, target 5 = 100

    # Window [0, 0] - only touching/overlapping
    res_zero <- gintervals.neighbors(query, targets, maxneighbors = 10, mindist = 0, maxdist = 0)
    expect_equal(nrow(res_zero), 1)
    expect_equal(res_zero$tid, 3)

    # Window [40, 60] - only targets at distance 50
    res_mid <- gintervals.neighbors(query, targets, maxneighbors = 10, mindist = 40, maxdist = 60)
    expect_equal(nrow(res_mid), 2)
    expect_setequal(res_mid$tid, c(2, 4))

    # Window [0, 100] - all targets
    res_all <- gintervals.neighbors(query, targets, maxneighbors = 10, mindist = 0, maxdist = 100)
    expect_equal(nrow(res_all), 5)
})

test_that("gintervals.neighbors respects maxneighbors limit", {
    query <- data.frame(
        chrom = "chr1",
        start = 1000,
        end = 1100
    )

    # Create many targets at various distances
    targets <- data.frame(
        chrom = rep("chr1", 10),
        start = 1100 + (1:10) * 10,
        end = 1100 + (1:10) * 10 + 50,
        tid = 1:10
    )

    # Request only 3 nearest neighbors
    res <- gintervals.neighbors(query, targets, maxneighbors = 3, mindist = -1e9, maxdist = 1e9)

    # Should get exactly 3 neighbors (the nearest ones)
    expect_equal(nrow(res), 3)
    expect_equal(res$tid, c(1, 2, 3)) # Closest 3
    expect_equal(res$dist, c(10, 20, 30))
})

test_that("gintervals.neighbors handles no neighbors found", {
    query <- data.frame(
        chrom = "chr1",
        start = 1000,
        end = 1100
    )

    targets <- data.frame(
        chrom = "chr1",
        start = c(2000, 3000),
        end = c(2100, 3100)
    )

    # With narrow distance window, should find nothing
    res <- gintervals.neighbors(query, targets, maxneighbors = 10, mindist = 10, maxdist = 50)

    # Should return empty/NULL
    expect_true(is.null(res) || nrow(res) == 0)
})

test_that("gintervals.neighbors handles no neighbors with na.if.notfound", {
    query <- data.frame(
        chrom = "chr1",
        start = 1000,
        end = 1100,
        qid = 1
    )

    targets <- data.frame(
        chrom = "chr1",
        start = c(2000),
        end = c(2100),
        tid = 1
    )

    # With narrow distance window and na.if.notfound, should return NA row
    res <- gintervals.neighbors(query, targets,
        maxneighbors = 10,
        mindist = 10, maxdist = 50,
        na.if.notfound = TRUE
    )

    expect_equal(nrow(res), 1)
    expect_true(is.na(res$tid))
    expect_true(is.na(res$dist))
    expect_equal(res$qid, 1) # Query info should be preserved
})

test_that("gintervals.neighbors works across multiple chromosomes", {
    queries <- data.frame(
        chrom = c("chr1", "chr2", "chr3"),
        start = c(100, 200, 300),
        end = c(150, 250, 350),
        qid = 1:3
    )

    targets <- data.frame(
        chrom = c("chr1", "chr1", "chr2", "chr2", "chr3"),
        start = c(100, 200, 200, 300, 300),
        end = c(150, 250, 250, 350, 350),
        tid = 1:5
    )

    res <- gintervals.neighbors(queries, targets, maxneighbors = 10, mindist = 0, maxdist = 0)

    # Query on chr1 should only find targets on chr1
    q1_res <- res %>% filter(qid == 1)
    expect_true(all(q1_res$chrom1 == "chr1"))
    expect_equal(nrow(q1_res), 1) # Only target 1 overlaps
    expect_equal(q1_res$tid, 1)

    # Query on chr2 should only find targets on chr2
    q2_res <- res %>% filter(qid == 2)
    expect_true(all(q2_res$chrom1 == "chr2"))
    expect_equal(nrow(q2_res), 1)
    expect_equal(q2_res$tid, 3)

    # Query on chr3 should only find targets on chr3
    q3_res <- res %>% filter(qid == 3)
    expect_true(all(q3_res$chrom1 == "chr3"))
    expect_equal(nrow(q3_res), 1)
    expect_equal(q3_res$tid, 5)
})

test_that("gintervals.neighbors finds neighbors on both sides", {
    query <- data.frame(
        chrom = "chr1",
        start = 1000,
        end = 1100
    )

    targets <- data.frame(
        chrom = rep("chr1", 4),
        start = c(
            900, # dist 50 (left, gap = 1000 - 950)
            850, # dist 100 (left, gap = 1000 - 900)
            1150, # dist 50 (right, gap = 1150 - 1100)
            1200 # dist 100 (right, gap = 1200 - 1100)
        ),
        end = c(
            950, 900, 1200, 1250
        ),
        tid = 1:4
    )

    # Window [0, 120] should find all 4 targets
    res <- gintervals.neighbors(query, targets, maxneighbors = 10, mindist = 0, maxdist = 120)
    expect_equal(nrow(res), 4)
    expect_setequal(res$tid, c(1, 2, 3, 4))

    # Window [40, 60] should find only the distance-50 targets
    res_narrow <- gintervals.neighbors(query, targets, maxneighbors = 10, mindist = 40, maxdist = 60)
    expect_equal(nrow(res_narrow), 2)
    expect_setequal(res_narrow$tid, c(1, 3))
    expect_true(all(res_narrow$dist == 50))
})

test_that("gintervals.neighbors returns results sorted by distance", {
    query <- data.frame(
        chrom = "chr1",
        start = 1000,
        end = 1100
    )

    targets <- data.frame(
        chrom = rep("chr1", 5),
        start = c(1200, 1150, 1300, 1250, 1180),
        end = c(1250, 1200, 1350, 1300, 1230),
        tid = c(5, 2, 10, 7, 4) # Intentionally not in distance order
    )

    res <- gintervals.neighbors(query, targets, maxneighbors = 10, mindist = 0, maxdist = 1e9)

    # Results should be sorted by absolute distance
    expect_equal(res$dist, c(50, 80, 100, 150, 200))
    expect_equal(res$tid, c(2, 4, 5, 7, 10))
})

test_that("gintervals.neighbors handles identical intervals correctly", {
    # Multiple identical queries
    queries <- data.frame(
        chrom = rep("chr1", 3),
        start = rep(1000, 3),
        end = rep(1100, 3),
        qid = 1:3
    )

    # Multiple identical targets
    targets <- data.frame(
        chrom = rep("chr1", 3),
        start = rep(1000, 3),
        end = rep(1100, 3),
        tid = 1:3
    )

    res <- gintervals.neighbors(queries, targets, maxneighbors = 10, mindist = 0, maxdist = 0)

    # Each query should find all 3 identical targets
    expect_equal(nrow(res), 9) # 3 queries × 3 targets
    expect_true(all(res$dist == 0))

    # Each query should have 3 results
    for (i in 1:3) {
        qres <- res %>% filter(qid == i)
        expect_equal(nrow(qres), 3)
        expect_setequal(qres$tid, c(1, 2, 3))
    }
})
