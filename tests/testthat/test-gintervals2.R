create_isolated_test_db()
test_that("gintervals.is.bigset checks bigintervs1d", {
    withr::with_options(list(gmax.data.size = 100), {
        expect_true(gintervals.is.bigset("bigintervs1d"))
    })
})

test_that("gintervals.is.bigset checks bigintervs2d", {
    withr::with_options(list(gmax.data.size = 100), {
        expect_true(gintervals.is.bigset("bigintervs2d"))
    })
})

test_that("gintervals.is.bigset checks test.tss", {
    withr::with_options(list(gmax.data.size = 100), {
        expect_false(gintervals.is.bigset("test.tss"))
    })
})

test_that("gintervals.rbind with saved data works", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))
    intervs1 <- gextract("test.fixedbin", gintervals(c(1, 2), 1000, 4000))
    intervs2 <- gextract("test.fixedbin", gintervals(c(2, "X"), 2000, 5000))
    gintervals.save(temp_track_name, intervs2)
    r <- gintervals.rbind(intervs1, temp_track_name)
    expect_true(!is.null(r))
    expect_regression(r, "gintervals.rbind.1")
})

test_that("gintervals.ls updates after save and remove operations", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))
    gintervals.save(temp_track_name, gintervals(c(1, 2)))
    r1 <- as.character(gintervals.ls())
    r2 <- as.character(gintervals.ls())
    # Filter out temporary track names and legacy test interval sets for regression comparison
    r1_filtered <- r1[!grepl("^test\\.tmptrack_", r1) & r1 != "test.testintervs" & r1 != "testintervs"]
    r2_filtered <- r2[!grepl("^test\\.tmptrack_", r2) & r2 != "test.testintervs" & r2 != "testintervs"]
    # Sort to ensure consistent ordering for comparison
    r1_filtered <- sort(r1_filtered)
    r2_filtered <- sort(r2_filtered)
    expect_regression(list(r1_filtered, r2_filtered), "gintervals.rbind.2")
})

test_that("gextract with removed interval gives an error", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))
    gintervals.save(temp_track_name, gintervals(c(1, 2), 1000, 2000))
    gintervals.rm(temp_track_name, force = TRUE)
    expect_error(gextract("test.fixedbin", temp_track_name))
})

test_that("gintervals.rm handles non-existent data without error when using force", {
    expect_silent(gintervals.rm("test.aaaaaaaaaaaaaaaaaaa", force = TRUE))
    expect_error(gintervals.rm("test.aaaaaaaaaaaaaaaaaaa"))
})

test_that("gintervals.ls reflects changes after save", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))
    r1 <- as.character(gintervals.ls())
    gintervals.save(temp_track_name, gintervals(c(1, 2), 1000, 2000))
    r2 <- as.character(gintervals.ls())
    # Filter out temporary track names and legacy test interval sets for regression comparison
    r1_filtered <- r1[!grepl("^test\\.tmptrack_", r1) & r1 != "test.testintervs" & r1 != "testintervs"]
    r2_filtered <- r2[!grepl("^test\\.tmptrack_", r2) & r2 != "test.testintervs" & r2 != "testintervs"]
    # Sort to ensure consistent ordering for comparison
    r1_filtered <- sort(r1_filtered)
    r2_filtered <- sort(r2_filtered)
    expect_regression(list(r1_filtered, r2_filtered), "gintervals.ls.1")
})

test_that("gextract works with saved intervals", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))
    gintervals.save(temp_track_name, gintervals(c(1, 2), 1000, 2000))
    r <- gextract("test.fixedbin", temp_track_name)
    expect_true(!is.null(r))
    expect_regression(r, "gintervals.save.1")
})

test_that("gscreen and gintervals.union works correctly", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    intervs1 <- gscreen("test.fixedbin > 0.1 & test.fixedbin < 0.3", gintervals(c(1, 2, 4, 8, 9), 0, -1))
    intervs2 <- gscreen("test.fixedbin < 0.2", gintervals(c(1, 2, 4, 7, 9), 0, -1))

    withr::with_options(list(gmax.data.size = 1000000), {
        gintervals.union(intervs1, intervs2, intervals.set.out = temp_track_name)
    })

    r <- gintervals.load(temp_track_name)
    expect_true(!is.null(r))
    expect_regression(r, "gscreen_and_gintervals.union.1")
})

test_that("gintervals.update with saved data using chrom 1", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs1d")
    })

    expect_error(gintervals.update(temp_track_name, gintervals(1), chrom = 1))
})

test_that("gintervals.update with loaded data using chrom 1", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs1d")
    })

    r <- gintervals.load(temp_track_name, chrom = 2)
    expect_error(gintervals.update(temp_track_name, r[c(2, 3), ], chrom1 = 1))
})

test_that("gintervals.update with loaded data using chrom 2", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs1d")
    })

    r <- gintervals.load(temp_track_name, chrom = 2)
    gintervals.update(temp_track_name, r[c(2, 3), ], chrom = 2)
    result <- list(gintervals.load(temp_track_name, chrom = 2), gintervals.chrom_sizes(temp_track_name))
    expect_regression(result, "gintervals.update.3")
})

test_that("gintervals.update removes chrom 2 from saved data", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs1d")
    })

    gintervals.update(temp_track_name, NULL, chrom = 2)
    result <- list(gintervals.load(temp_track_name, chrom = 2), gintervals.chrom_sizes(temp_track_name))
    expect_regression(result, "gintervals.update.4")
})

test_that("gintervals.update with saved 2d data using chrom1 1", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs2d")
    })

    expect_error(gintervals.update(temp_track_name, gintervals.2d(1), chrom1 = 1))
})

test_that("gintervals.update with loaded 2d data using chrom 1", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs2d")
    })

    r <- gintervals.load(temp_track_name, chrom1 = 2, chrom2 = 2)
    expect_error(gintervals.update(temp_track_name, r[c(2, 3), ], chrom = 1))
})

test_that("gintervals.update with loaded 2d data using chrom1 2 and chrom2 2", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs2d")
    })

    r <- gintervals.load(temp_track_name, chrom1 = 2, chrom2 = 2)
    gintervals.update(temp_track_name, r[c(2, 3), ], chrom1 = 2, chrom2 = 2)
    result <- list(gintervals.load(temp_track_name, chrom1 = 2, chrom2 = 2), gintervals.chrom_sizes(temp_track_name))
    expect_regression(result, "gintervals.update.2d.3")
})

test_that("gintervals.update removes chrom1 2 and chrom2 2 from saved 2d data", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs2d")
    })

    gintervals.update(temp_track_name, NULL, chrom1 = 2, chrom2 = 2)
    result <- list(gintervals.load(temp_track_name, chrom1 = 2, chrom2 = 2), gintervals.chrom_sizes(temp_track_name))
    expect_regression(result, "gintervals.update.2d.4")
})

# gintervals.rbind() gives the columns base rbind() gives: values, types and factor levels
expect_rbind_like_base <- function(...) {
    expect_identical(as.list(gintervals.rbind(...)), as.list(rbind(...)))
}

test_that("gintervals.rbind binds factors with different levels by their values", {
    a <- data.frame(chrom = factor("chr2"), start = 0, end = 100)
    b <- data.frame(chrom = factor("chr1"), start = 500, end = 600)
    expect_equal(as.character(gintervals.rbind(a, b)$chrom), c("chr2", "chr1"))
    expect_rbind_like_base(a, b)

    # an extra factor column of saved sets: their chrom column is normalized on save, this one is not
    a$name <- factor("x")
    b$name <- factor("y")
    set_a <- paste0("test.tmprbind_", sample(1:1e9, 1))
    set_b <- paste0("test.tmprbind_", sample(1:1e9, 1))
    withr::defer({
        gintervals.rm(set_a, force = TRUE)
        gintervals.rm(set_b, force = TRUE)
    })
    gintervals.save(set_a, a)
    gintervals.save(set_b, b)
    expect_equal(as.character(gintervals.rbind(set_a, set_b)$name), c("x", "y"))
})

test_that("gintervals.rbind binds a factor with a character column and an integer with a double one", {
    fac <- gintervals("chr1", 0, 10)
    chr <- data.frame(chrom = "chr1", start = 5L, end = 20L)
    r <- gintervals.rbind(fac, chr)
    expect_equal(as.character(r$chrom), c("chr1", "chr1"))
    expect_identical(r$start, c(0, 5))
    expect_rbind_like_base(fac, chr)
    expect_rbind_like_base(chr, fac)
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
})

test_that("gintervals.rbind matches rbind on factor and character columns with NA", {
    set.seed(17)
    chroms <- c("chr1", "chr2", "chrX")
    random_frame <- function() {
        n <- sample(1:4, 1)
        d <- data.frame(
            chrom = sample(chroms, n, replace = TRUE), start = seq_len(n) * 10, end = seq_len(n) * 10 + 5,
            lab = sample(c("p", "q", "r", NA), n, replace = TRUE)
        )
        if (runif(1) < 0.7) {
            d$chrom <- factor(d$chrom, levels = sample(chroms))
        }
        kind <- sample(c("character", "factor", "factor with an NA level"), 1)
        if (kind == "factor") {
            d$lab <- factor(d$lab, levels = sample(c("p", "q", "r", "s")))
        } else if (kind == "factor with an NA level") {
            d$lab <- factor(d$lab, exclude = NULL)
        }
        if (runif(1) < 0.5) {
            d$start <- as.integer(d$start)
        }
        d
    }
    for (i in 1:200) {
        objs <- replicate(sample(2:4, 1), random_frame(), simplify = FALSE)
        expect_identical(as.list(do.call(gintervals.rbind, objs)), as.list(do.call(rbind, objs)), info = i)
    }
})

test_that("gintervals.rbind refuses columns of different types with an R error", {
    a <- data.frame(chrom = "chr1", start = 0, end = 10, x = 1)
    b <- data.frame(chrom = "chr1", start = 5, end = 20, x = "a")
    expect_error(gintervals.rbind(a, b), "column 'x' is numeric in object 1 and character in object 2")
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
})

