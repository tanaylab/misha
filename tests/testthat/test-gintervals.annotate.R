load_test_db()
test_that("gintervals.annotate basic annotation with distance", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")

    res <- gintervals.annotate(
        intervals = intervs,
        annotation_intervals = ann,
        annotation_columns = c("remark"),
        dist_column = "dist"
    )

    expect_equal(nrow(res), 2)
    expect_true(all(c("chrom", "start", "end", "remark", "dist") %in% colnames(res)))

    nei <- gintervals.neighbors(intervs, ann, 1, na.if.notfound = TRUE)
    expect_equal(as.character(res$remark), as.character(nei$remark))
    expect_equal(res$dist, nei$dist)
})

test_that("gintervals.annotate with max_dist and na_value", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")

    res <- gintervals.annotate(
        intervals = intervs,
        annotation_intervals = ann,
        annotation_columns = c("remark"),
        max_dist = 200,
        na_value = list(remark = "none")
    )

    expect_equal(nrow(res), 2)
    nei <- gintervals.neighbors(intervs, ann, 1, na.if.notfound = TRUE)
    expect_equal(res$dist, nei$dist)
    expect_equal(as.character(res$remark), ifelse(abs(nei$dist) > 200, "none", as.character(nei$remark)))
})

test_that("gintervals.annotate custom column and distance names", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")

    res <- gintervals.annotate(
        intervals = intervs,
        annotation_intervals = ann,
        annotation_columns = c("remark"),
        column_names = c("ann_remark"),
        dist_column = "ann_dist"
    )

    expect_true(all(c("ann_remark", "ann_dist") %in% colnames(res)))
    nei <- gintervals.neighbors(intervs, ann, 1, na.if.notfound = TRUE)
    expect_equal(as.character(res$ann_remark), as.character(nei$remark))
    expect_equal(res$ann_dist, nei$dist)
})

test_that("gintervals.annotate multiple neighbors and tie-breaking", {
    # One interval with two equidistant neighbors: left [800,900], right [1100,1200]
    intervs <- gintervals(1, 1000, 1100)
    # Make both neighbors 100bp away: left [800,900] (distance 100), right [1200,1300] (distance 100)
    ann <- gintervals(1, c(800, 1200), c(900, 1300))
    ann$label <- c("left", "right")

    # Ask for 2 neighbors, use tie method by min.start => left comes before right
    res_min_start <- gintervals.annotate(
        intervals = intervs,
        annotation_intervals = ann,
        annotation_columns = c("label"),
        maxneighbors = 2,
        tie_method = "min.start"
    )
    expect_equal(as.character(res_min_start$label), c("left", "right"))

    # With min.end, left interval ends earlier (900 vs 1200) => left then right as well
    res_min_end <- gintervals.annotate(
        intervals = intervs,
        annotation_intervals = ann,
        annotation_columns = c("label"),
        maxneighbors = 2,
        tie_method = "min.end"
    )
    expect_equal(as.character(res_min_end$label), c("left", "right"))
})

test_that("gintervals.annotate overwrite behavior", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    intervs$remark <- c("orig1", "orig2")
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")

    # By default, should error due to column conflict
    expect_error(gintervals.annotate(
        intervals = intervs,
        annotation_intervals = ann,
        annotation_columns = c("remark")
    ))

    # With overwrite=TRUE, proceeds and replaces existing column
    res <- gintervals.annotate(
        intervals = intervs,
        annotation_intervals = ann,
        annotation_columns = c("remark"),
        overwrite = TRUE
    )
    expect_equal(as.character(res$remark), c("a", "b"))
})

test_that("gintervals.annotate preserves original order when keep_order=TRUE", {
    intervs <- gintervals(1, c(1000, 5000, 2000), c(1100, 5050, 2100))
    ord <- c(3, 1, 2)
    intervs <- intervs[ord, ]
    ann <- gintervals(1, c(900, 5400, 1500), c(950, 5500, 1600))
    ann$remark <- c("a", "b", "c")

    res <- gintervals.annotate(
        intervals = intervs,
        annotation_intervals = ann,
        annotation_columns = c("remark"),
        keep_order = TRUE
    )

    # Map back to nearest annotations to confirm correspondence
    expect_equal(as.character(res$remark), c("c", "a", "b"))
})


test_that("gintervals.annotate selects all non-basic columns by default", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")
    ann$score <- c(10, 20)

    res <- gintervals.annotate(intervs, ann, annotation_columns = NULL)

    expect_true(all(c("remark", "score") %in% colnames(res)))
    nei <- gintervals.neighbors(intervs, ann, 1, na.if.notfound = TRUE)
    expect_equal(as.character(res$remark), as.character(nei$remark))
    expect_equal(res$score, nei$score)
})

test_that("gintervals.annotate handles empty intervals", {
    intervs <- data.frame(
        chrom = factor(character(0), levels = gintervals.all()$chrom),
        start = numeric(0),
        end = numeric(0)
    )
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")

    res <- gintervals.annotate(intervs, ann, annotation_columns = "remark")

    expect_equal(nrow(res), 0)
    expect_true(all(c("chrom", "start", "end", "remark", "dist") %in% colnames(res)))
})

test_that("gintervals.annotate handles empty annotation set", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- data.frame(
        chrom = factor(character(0), levels = gintervals.all()$chrom),
        start = numeric(0),
        end = numeric(0),
        remark = character(0)
    )

    res <- gintervals.annotate(intervs, ann, annotation_columns = "remark")

    expect_equal(nrow(res), nrow(intervs))
    expect_true(all(is.na(res$remark)))
    expect_true(all(is.na(res$dist)))
})

test_that("gintervals.annotate can suppress distance column", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")

    res <- gintervals.annotate(intervs, ann, annotation_columns = "remark", dist_column = NULL)

    expect_false("dist" %in% colnames(res))
})

test_that("gintervals.annotate errors when column_names length mismatches", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")
    ann$score <- c(1, 2)

    expect_error(gintervals.annotate(
        intervs, ann,
        annotation_columns = c("remark", "score"),
        column_names = c("only_one")
    ))
})

test_that("gintervals.annotate errors when distance column collides", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    intervs$dist <- c(1, 2)
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")

    expect_error(gintervals.annotate(intervs, ann, annotation_columns = "remark", dist_column = "dist"))
})

test_that("gintervals.annotate errors when annotation column missing", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")

    expect_error(gintervals.annotate(intervs, ann, annotation_columns = c("does_not_exist")))
})

test_that("gintervals.annotate limits neighbors count to available", {
    intervs <- gintervals(1, 1000, 1100)
    ann <- gintervals(1, 1200, 1300)
    ann$label <- "right"

    res <- gintervals.annotate(intervs, ann, annotation_columns = "label", maxneighbors = 3)
    expect_equal(nrow(res), 1)
    expect_equal(as.character(res$label), "right")
})

test_that("gintervals.annotate supports scalar na_value with threshold", {
    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")

    res <- gintervals.annotate(
        intervs, ann,
        annotation_columns = "remark",
        max_dist = 200, na_value = "no_ann"
    )

    nei <- gintervals.neighbors(intervs, ann, 1, na.if.notfound = TRUE)
    expect_equal(as.character(res$remark), ifelse(abs(nei$dist) > 200, "no_ann", as.character(nei$remark)))
})

test_that("gintervals.annotate works with intervals.set.out", {
    gintervals.rm("test.testintervs_ann", force = TRUE)
    withr::defer(gintervals.rm("test.testintervs_ann", force = TRUE))

    intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
    ann <- gintervals(1, c(900, 5400), c(950, 5500))
    ann$remark <- c("a", "b")

    res_mem <- gintervals.annotate(intervs, ann, annotation_columns = "remark") %>%
        tibble::repair_names() %>%
        arrange(chrom, start, end)

    gintervals.annotate(intervs, ann, annotation_columns = "remark", intervals.set.out = "test.testintervs_ann")

    res_file <- gintervals.load("test.testintervs_ann") %>%
        tibble::repair_names() %>%
        arrange(chrom, start, end)

    expect_equal(res_file, res_mem)
})
