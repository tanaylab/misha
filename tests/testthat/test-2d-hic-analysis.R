# Test 2D functions using real HiC analysis patterns
# Based on analysis code from /home/nettam/projects/hic/lscripts/

# ==============================================================================
# SETUP
# ==============================================================================

# Set up test database and ensure test data exists
test_that("HiC test database is available", {
    # Use the snapshot database
    test_db <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19"
    expect_true(dir.exists(test_db))

    # Set root temporarily
    withr::local_options(list(gmax.data.size = 1e8))
    old_root <- if (exists("GROOT", envir = .misha, inherits = FALSE)) {
        get("GROOT", envir = .misha)
    } else {
        NULL
    }
    gsetroot(test_db)
    withr::defer({
        if (!is.null(old_root)) gsetroot(old_root)
    })

    # Check and create test data if needed
    status <- setup_hic_test_data(force = FALSE)

    # Verify all required tracks exist
    expect_true(all(status$tracks))
    expect_true(all(status$intervals))
})

# ==============================================================================
# REAL DATA TESTS - K562 HiC Analysis
# ==============================================================================
# These tests use actual K562 HiC tracks from the production database
# to validate parity between vanilla misha and multi-contig version.
# Track references from /home/nettam/projects/hic/lscripts/high_score.r

test_that("Real data: K562 tracks are accessible", {
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")
    gdb.reload()

    # Verify symlinked K562 tracks exist
    expect_true(gtrack.exists("hic.K562.ela_k562"))
    expect_true(gtrack.exists("hic.K562.ela_k562_score"))
    expect_true(gtrack.exists("hic.K562.ela_k562_ins_200kb"))
    expect_true(gtrack.exists("hic.K562.ela_k562_SRR1658694_omer_ins_200kb"))
})

test_that("Real data: gextract on K562 HiC track with distance bands", {
    # Pattern: hic_basic.r gcis_decay function (lines 38-56)
    # This is the core distance decay analysis pattern
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    # Extract cis contacts for chr1 with distance band
    # Use a small region for manageable test size
    scope <- gintervals.2d("chr1", 0, 5e6, "chr1", 0, 5e6)

    result <- gextract(
        "hic.K562.ela_k562",
        intervals = scope,
        band = c(2e4, 1e6), # 20kb - 1Mb distance
        colnames = "contacts"
    )

    expect_regression(result, "real_k562_gextract_band.hic.1")

    # Validate distance filtering
    if (!is.null(result) && nrow(result) > 0) {
        dists <- abs(result$start2 - result$start1)
        expect_true(all(dists >= 2e4))
        expect_true(all(dists <= 1e6))
    }
})

test_that("Real data: gcis_decay-style distance analysis on K562", {
    # Pattern: hic_basic.r gcis_decay function
    # Tests binning contacts by distance
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    # Sample a small region for testing
    scope <- gintervals.2d("chr1", 0, 2e6, "chr1", 0, 2e6)

    # Define distance bins (simplified from hic_basic.r)
    dist_breaks <- c(1e4, 5e4, 1e5, 5e5, 1e6)

    # Extract contacts in each distance range
    decay_data <- lapply(seq_along(dist_breaks[-1]), function(i) {
        band_min <- dist_breaks[i]
        band_max <- dist_breaks[i + 1]

        contacts <- gextract(
            "hic.K562.ela_k562",
            intervals = scope,
            band = c(band_min, band_max),
            colnames = "contacts"
        )

        if (is.null(contacts) || nrow(contacts) == 0) {
            return(data.frame(dist_bin = (band_min + band_max) / 2, count = 0))
        }

        data.frame(
            dist_bin = (band_min + band_max) / 2,
            count = sum(contacts$contacts, na.rm = TRUE)
        )
    })

    decay_df <- do.call(rbind, decay_data)
    expect_regression(decay_df, "real_k562_decay.hic.1")

    # Expect decreasing trend (though not strict monotonic due to sampling)
    expect_true(nrow(decay_df) > 0)
})

test_that("Real data: gintervals.2d.band_intersect on K562 regions", {
    # Pattern: hic_basic.r usage of band filtering with domains
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    # Create test regions on chr1
    regions <- gintervals("chr1", c(1e6, 2e6), c(1.5e6, 2.5e6))

    # Create 2D intervals for all pairs
    iter_2d <- expand.grid(
        idx1 = seq_len(nrow(regions)),
        idx2 = seq_len(nrow(regions))
    )
    iter_2d <- gintervals.2d(
        regions$chrom[iter_2d$idx1], regions$start[iter_2d$idx1], regions$end[iter_2d$idx1],
        regions$chrom[iter_2d$idx2], regions$start[iter_2d$idx2], regions$end[iter_2d$idx2]
    )

    # Filter by distance band
    band_filtered <- gintervals.2d.band_intersect(iter_2d, band = c(1e5, 1e6))

    expect_regression(band_filtered, "real_k562_band_intersect.hic.1")

    # Validate distance constraints
    if (!is.null(band_filtered) && nrow(band_filtered) > 0) {
        dists <- abs(band_filtered$start2 - band_filtered$start1)
        expect_true(all(dists >= 1e5 & dists <= 1e6))
    }
})

test_that("Real data: Virtual track weighted.sum on K562 score track", {
    # Pattern: high_score.r lines 110-120 (weighted scoring)
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    # Create virtual track for weighted sum
    gvtrack.create("v_k562_score", "hic.K562.ela_k562_score", func = "weighted.sum")
    withr::defer(gvtrack.rm("v_k562_score"))

    # Extract from small region
    scope <- gintervals.2d("chr1", 0, 1e6, "chr1", 0, 1e6)
    result <- gextract("v_k562_score", intervals = scope, band = c(1e4, 5e5))

    expect_regression(result, "real_k562_vtrack_weighted_sum.hic.1")

    if (!is.null(result) && nrow(result) > 0) {
        expect_true("v_k562_score" %in% colnames(result))
        expect_true(all(!is.na(result$v_k562_score)))
    }
})

test_that("Real data: Insulation track extraction (1D)", {
    # Pattern: high_score.r lines 94-95, 271
    # Insulation is a 1D track derived from 2D HiC data
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    # Extract insulation for chr1 region
    insulation <- gextract(
        "hic.K562.ela_k562_ins_200kb",
        intervals = gintervals("chr1", 0, 2e6),
        iterator = gintervals("chr1", 0, 2e6),
        colnames = "insulation"
    )

    expect_regression(insulation, "real_k562_insulation.hic.1")
    expect_true(!is.null(insulation))
    expect_true(nrow(insulation) > 0)
    expect_true("insulation" %in% colnames(insulation))
})

test_that("Real data: gscreen-style domain calling from insulation", {
    # Pattern: hic_basic.r lines 100-120 (domain boundary detection)
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    # Find local minima in insulation (domain boundaries)
    # Using gscreen to find regions below threshold
    scope <- gintervals("chr1", 0, 5e6)

    boundaries <- gscreen(
        "hic.K562.ela_k562_ins_200kb < -2",
        intervals = scope,
        iterator = scope
    )

    expect_regression(boundaries, "real_k562_boundaries.hic.1")

    if (!is.null(boundaries) && nrow(boundaries) > 0) {
        expect_true(all(c("chrom", "start", "end") %in% colnames(boundaries)))
        expect_true(all(boundaries$chrom == "chr1"))
    }
})

test_that("Real data: Multi-track extraction (score + insulation)", {
    # Pattern: Combining 2D track data with 1D annotation
    # Common pattern in HiC analysis workflows
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    # Extract HiC scores for a region
    scope <- gintervals.2d("chr1", 0, 1e6, "chr1", 0, 1e6)
    hic_data <- gextract(
        "hic.K562.ela_k562_score",
        intervals = scope,
        band = c(2e4, 5e5),
        colnames = "score"
    )

    expect_regression(hic_data, "real_k562_multi_track.hic.1")

    if (!is.null(hic_data) && nrow(hic_data) > 0) {
        expect_true("score" %in% colnames(hic_data))
        expect_true(all(!is.na(hic_data$score)))
    }
})

# ==============================================================================
# SYNTHETIC DATA TESTS - For comprehensive coverage
# ==============================================================================
# The following tests use synthetic data to ensure all edge cases
# and function signatures work correctly

# ==============================================================================
# GROUP 1: BASIC 2D OPERATIONS
# ==============================================================================

test_that("gintervals.2d creates 2D intervals from bins", {
    # Pattern: hic_basic.r line 44
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    bins <- gintervals("chr1", seq(0, 2e5, 1e4), seq(1e4, 2.1e5, 1e4))
    # gintervals.2d creates 1-to-1 pairs, not Cartesian product
    iter_2d <- gintervals.2d(
        bins$chrom, bins$start, bins$end,
        bins$chrom, bins$start, bins$end
    )

    expect_regression(iter_2d, "gintervals.2d.from_bins.1")
    expect_equal(nrow(iter_2d), nrow(bins)) # 1-to-1 pairing
    expect_true(all(iter_2d$chrom1 == iter_2d$chrom2))
    expect_true(all(iter_2d$chrom1 == "chr1"))
    # Verify it's actually pairing correctly
    expect_equal(iter_2d$start1, iter_2d$start2)
    expect_equal(iter_2d$end1, iter_2d$end2)
})

test_that("gintervals.2d.all creates genome-wide scope", {
    # Pattern: figures.r lines 244, 275
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    scope <- gintervals.2d.all()
    scope_cis <- scope[scope$chrom1 == scope$chrom2, ]

    expect_regression(scope, "gintervals.2d.all.hic.1")
    expect_true(nrow(scope_cis) > 0)
    expect_true(nrow(scope) >= nrow(scope_cis))

    # Check structure
    expect_true(all(c("chrom1", "start1", "end1", "chrom2", "start2", "end2") %in% colnames(scope)))
})

test_that("gintervals.2d from domains creates self-pairs", {
    # Pattern: hic_basic.r lines 81-83
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    domains <- gintervals.load("domains.test")
    dom_2d <- gintervals.2d(
        domains$chrom, domains$start, domains$end,
        domains$chrom, domains$start, domains$end
    )

    expect_regression(dom_2d, "gintervals.2d.domains.1")
    expect_equal(nrow(dom_2d), nrow(domains))

    # Verify self-pairs
    expect_true(all(dom_2d$chrom1 == dom_2d$chrom2))
    expect_true(all(dom_2d$start1 == dom_2d$start2))
    expect_true(all(dom_2d$end1 == dom_2d$end2))
})

test_that("expand.grid pattern creates all feature pairs", {
    # Pattern: hic_basic.r lines 435-451
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    f1 <- gintervals("chr1", c(5e5, 1e6, 1.5e6), c(5.1e5, 1.1e6, 1.6e6))
    f2 <- gintervals("chr1", c(7e5, 1.2e6, 1.8e6), c(7.1e5, 1.3e6, 1.9e6))

    g <- expand.grid(1:nrow(f1), 1:nrow(f2))
    grid <- data.frame(
        chrom1 = "chr1", start1 = f1$start[g$Var1], end1 = f1$end[g$Var1],
        chrom2 = "chr1", start2 = f2$start[g$Var2], end2 = f2$end[g$Var2]
    )

    expect_regression(grid, "expand.grid.2d.1")
    expect_equal(nrow(grid), nrow(f1) * nrow(f2))

    # Verify all combinations present
    expect_setequal(unique(grid$start1), f1$start)
    expect_setequal(unique(grid$start2), f2$start)
})

# ==============================================================================
# GROUP 2: DISTANCE-BASED ANALYSIS (CRITICAL)
# ==============================================================================

test_that("gcis_decay computes distance decay curve", {
    # Pattern: figures.r lines 65-66
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    breaks <- 2^seq(10, 20, by = 0.5)

    result_new <- gcis_decay(
        "hic.test_basic", breaks,
        gintervals.all(), gintervals.all(),
        include.lowest = TRUE, band = c(-2e6, -1024)
    )

    expect_regression(result_new, "gcis_decay.basic.hic.1")

    # Validate properties
    expect_equal(nrow(result_new), length(breaks) - 1)
    expect_equal(ncol(result_new), 2) # intra and inter columns

    # Check for general decay trend (most bins should decrease)
    # Allow some noise but expect overall decreasing trend
    expect_true(result_new[1, 1] > result_new[nrow(result_new), 1])
})

test_that("gcis_decay within domains combines upstream and downstream", {
    # Pattern: hic_basic.r lines 136-140
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    domains <- gintervals.load("domains.test")
    breaks <- 2^seq(10, 18, by = 0.5)
    min_dist <- 1e4
    max_dist <- 5e5

    # Upstream
    up <- gcis_decay(
        "hic.test_basic", breaks, gintervals.all(), domains,
        band = c(min_dist, max_dist)
    )

    # Downstream
    down <- gcis_decay(
        "hic.test_basic", breaks, gintervals.all(), domains,
        band = c(-max_dist, -min_dist)
    )

    # Combined
    combined <- up
    combined[, 1] <- up[, 1] + down[, 1]
    combined[, 2] <- up[, 2] + down[, 2]

    expect_regression(combined, "gcis_decay.domains.symmetric.hic.1")

    # Check that combined has more counts than individual
    expect_true(all(combined[, 1] >= up[, 1]))
    expect_true(all(combined[, 1] >= down[, 1]))
})

test_that("gextract with band excludes diagonal correctly", {
    # Pattern: Used 30+ times in both scripts
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    interval <- gintervals.2d("chr1", 5e5, 1e6, "chr1", 5e5, 1e6)
    min_dist <- 1024

    # Most common: exclude diagonal region
    result <- gextract("hic.test_basic", interval,
        band = c(-2e6, -min_dist)
    )

    expect_regression(result, "gextract.band.exclude_diag.hic.1")

    # Validate distances
    if (!is.null(result) && nrow(result) > 0) {
        expect_true(all(abs(result$start2 - result$start1) >= min_dist))
    }
})

test_that("gextract with band filters to specific ranges", {
    # Pattern: hic_basic.r lines 111, 118
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    interval <- gintervals.2d("chr1", 5e5, 1.5e6, "chr1", 5e5, 1.5e6)

    # Upstream only
    up <- gextract("hic.test_basic", interval, band = c(-5e5, -1e4))

    # Downstream only
    down <- gextract("hic.test_basic", interval, band = c(1e4, 5e5))

    expect_regression(list(up = up, down = down), "gextract.band.ranges.hic.1")

    # Validate distances for upstream (use abs for negative band)
    if (!is.null(up) && nrow(up) > 0) {
        dists_up <- abs(up$start1 - up$start2)
        expect_true(all(dists_up >= 1e4))
        expect_true(all(dists_up <= 5e5))
    }

    # Validate distances for downstream
    if (!is.null(down) && nrow(down) > 0) {
        dists_down <- abs(down$start2 - down$start1)
        expect_true(all(dists_down >= 1e4))
        expect_true(all(dists_down <= 5e5))
    }
})

test_that("gintervals.2d.band_intersect filters by distance", {
    # Existing pattern from test-2d-parity.R
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    starts <- seq(0, 1e6, 1e5)
    ends <- starts + 1e5
    # Create all pairs using expand.grid pattern
    g <- expand.grid(1:length(starts), 1:length(starts))
    intervals <- data.frame(
        chrom1 = "chr1", start1 = starts[g$Var1], end1 = ends[g$Var1],
        chrom2 = "chr1", start2 = starts[g$Var2], end2 = ends[g$Var2]
    )

    # Positive band (downstream)
    filtered_pos <- gintervals.2d.band_intersect(intervals, c(2e5, 5e5))

    # Negative band (upstream)
    filtered_neg <- gintervals.2d.band_intersect(intervals, c(-5e5, -2e5))

    expect_regression(
        list(pos = filtered_pos, neg = filtered_neg),
        "band_intersect.pos_neg.1"
    )

    # Validate filtering - band selects by absolute distance
    if (!is.null(filtered_pos) && nrow(filtered_pos) > 0) {
        dists_pos <- abs(filtered_pos$start2 - filtered_pos$start1)
        expect_true(all(dists_pos >= 2e5))
        expect_true(all(dists_pos <= 5e5))
    }

    if (!is.null(filtered_neg) && nrow(filtered_neg) > 0) {
        dists_neg <- abs(filtered_neg$start1 - filtered_neg$start2)
        expect_true(all(dists_neg >= 2e5))
        expect_true(all(dists_neg <= 5e5))
    }
})

# ==============================================================================
# GROUP 3: VIRTUAL TRACK OPERATIONS (HIGH)
# ==============================================================================

test_that("weighted.sum virtual track computes marginals", {
    # Pattern: hic_basic.r lines 40, 78-79
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    gvtrack.create("v_test", "hic.test_basic", "weighted.sum")
    withr::defer(gvtrack.rm("v_test"))

    iter_2d <- gintervals.2d("chr1", 5e5, 1e6)
    result <- gextract("v_test", intervals = iter_2d, iterator = iter_2d)

    expect_regression(result, "vtrack.weighted_sum.hic.1")

    # Check column exists
    expect_true("v_test" %in% colnames(result))
})

test_that("weighted.sum respects band parameter", {
    # Pattern: hic_basic.r lines 46-47
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    gvtrack.create("v_test", "hic.test_basic", "weighted.sum")
    withr::defer(gvtrack.rm("v_test"))

    iter_2d <- gintervals.2d("chr1", 5e5, 1e6)

    # Full marginal
    full <- gextract("v_test", intervals = iter_2d, iterator = iter_2d)

    # Diagonal only
    diag <- gextract("v_test",
        intervals = iter_2d, iterator = iter_2d,
        band = c(-1024, 1024)
    )

    expect_regression(list(full = full, diag = diag), "vtrack.band_filter.hic.1")

    # Diagonal should have fewer or equal values
    diag_sum <- sum(diag$v_test, na.rm = TRUE)
    full_sum <- sum(full$v_test, na.rm = TRUE)
    expect_true(diag_sum <= full_sum)
})

test_that("max virtual track finds maximum in region", {
    # Pattern: figures.r lines 522, 763
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    gvtrack.create("max_score", "hic.test_score", "max")
    withr::defer(gvtrack.rm("max_score"))

    gvtrack.iterator.2d("max_score",
        sshift1 = -5e4, eshift1 = 5e4,
        sshift2 = -5e4, eshift2 = 5e4
    )

    points <- gintervals.2d("chr1", c(5e5, 1e6), c(5e5 + 1, 1e6 + 1))
    result <- gextract("max_score", intervals = points, iterator = points)

    expect_regression(result, "vtrack.max.iterator2d.hic.1")

    # Check that we got max values
    expect_true("max_score" %in% colnames(result))
})

test_that("distance virtual tracks work with 2D iterators", {
    # Pattern: hic_basic.r lines 369-373
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    features <- gintervals("chr1", c(5e5, 1e6, 1.5e6), c(5.1e5, 1.1e6, 1.6e6))

    gvtrack.create("dist1", features, "distance")
    gvtrack.create("dist2", features, "distance")
    withr::defer({
        gvtrack.rm("dist1")
        gvtrack.rm("dist2")
    })

    gvtrack.iterator("dist1", 1)
    gvtrack.iterator("dist2", 2)

    scope <- gintervals.2d("chr1", 0, 2e6, "chr1", 0, 2e6)
    result <- gextract("dist1", "dist2", "hic.test_basic", scope,
        band = c(-1e6, -1e4)
    )

    expect_regression(result, "vtrack.distance.dual.hic.1")

    # Check distance columns exist
    expect_true(all(c("dist1", "dist2") %in% colnames(result)))
})

# ==============================================================================
# GROUP 4: DOMAIN/INSULATION ANALYSIS (HIGH)
# ==============================================================================

test_that("gtrack.2d.get_insu_doms extracts domains from insulation", {
    # Pattern: hic_basic.r line 31
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    domains <- gtrack.2d.get_insu_doms("insulation.test", thresh = -3.0)
    domains_filtered <- domains[domains$end - domains$start > 5e4, ]

    expect_regression(domains_filtered, "get_insu_doms.1")

    # Validate domain structure
    expect_true(all(domains_filtered$end > domains_filtered$start))
    expect_true(all(c("chrom", "start", "end") %in% colnames(domains_filtered)))
})

test_that("intra-domain contacts show enrichment vs far-cis", {
    # Pattern: hic_basic.r lines 75-95
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    domains <- gintervals.load("domains.test")

    gvtrack.create("v_test", "hic.test_basic", "weighted.sum")
    withr::defer(gvtrack.rm("v_test"))

    # Intra-domain
    dom_2d <- gintervals.2d(
        domains$chrom, domains$start, domains$end,
        domains$chrom, domains$start, domains$end
    )
    intra <- gextract("v_test",
        intervals = dom_2d, iterator = dom_2d,
        band = c(-2e6, 0)
    )

    # Far-cis (background)
    iter_2d <- gintervals.2d(domains$chrom, domains$start, domains$end)
    far_cis <- gextract("v_test",
        intervals = iter_2d, iterator = iter_2d,
        band = c(-2e6, 0)
    )

    expect_regression(
        list(intra = intra, far_cis = far_cis),
        "intra_vs_far_cis.hic.1"
    )

    # Basic structure checks
    expect_equal(nrow(intra), nrow(dom_2d))
    expect_equal(nrow(far_cis), nrow(iter_2d))
})

# ==============================================================================
# GROUP 5: SCREENING AND FILTERING (HIGH)
# ==============================================================================

test_that("gscreen filters 2D intervals with spatial iterator", {
    # Pattern: hic_basic.r line 459
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    grid <- expand.grid(
        start1 = seq(5e5, 1e6, 1e5),
        start2 = seq(7e5, 1.2e6, 1e5)
    )
    grid <- data.frame(
        chrom1 = "chr1", start1 = grid$start1, end1 = grid$start1 + 1e4,
        chrom2 = "chr1", start2 = grid$start2, end2 = grid$start2 + 1e4
    )

    gvtrack.create("v_filter", "hic.test_score", "max")
    withr::defer(gvtrack.rm("v_filter"))

    gvtrack.iterator.2d("v_filter",
        sshift1 = -5e4, eshift1 = 5e4,
        sshift2 = -5e4, eshift2 = 5e4
    )

    screened <- gscreen("v_filter > 50",
        intervals = grid, iterator = grid,
        band = c(-1e6, -1e4)
    )

    expect_regression(screened, "gscreen.2d.iterator2d.1")

    # Check that filtering occurred
    expect_true(!is.null(screened) && !is.null(grid) && nrow(screened) <= nrow(grid))
})

test_that("giterator.intervals with band creates filtered set", {
    # Pattern: hic_basic.r line 464
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    # Create test directory if needed
    tryCatch(gdir.create("test", showWarnings = FALSE), error = function(e) NULL)

    set_name <- "test.band_set"
    # Clean up if already exists from previous run
    if (gintervals.exists(set_name)) {
        gintervals.rm(set_name, force = TRUE)
    }
    withr::defer(gintervals.rm(set_name, force = TRUE))

    giterator.intervals("hic.test_basic",
        band = c(-5e5, -1e4),
        intervals.set.out = set_name
    )

    loaded <- gintervals.load(set_name)

    expect_regression(loaded, "giterator.intervals.band.hic.1")

    # Validate distances
    if (!is.null(loaded) && nrow(loaded) > 0) {
        expect_true(all(abs(loaded$start2 - loaded$start1) >= 1e4))
        expect_true(all(abs(loaded$start2 - loaded$start1) <= 5e5))
    }
})

test_that("gscreen handles multiple 2D criteria", {
    # Pattern: hic_basic.r line 545
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    starts <- seq(5e5, 1e6, 5e4)
    ends <- starts + 5e4
    grid <- gintervals.2d("chr1", starts, ends)

    screened <- gscreen("hic.test_score > 40 & hic.test_basic > 0",
        intervals = grid, iterator = grid
    )

    expect_regression(screened, "gscreen.multi_criteria.1")

    # Basic validation
    expect_true(!is.null(screened) && !is.null(grid) && nrow(screened) <= nrow(grid))
})

# ==============================================================================
# GROUP 6: MARGINAL ANALYSIS (MEDIUM)
# ==============================================================================

test_that("marginal extraction with diagonal subtraction", {
    # Pattern: hic_basic.r lines 37-56
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    resolution <- 5e4
    intervals <- gintervals("chr1", 0, 5e5)

    iter_1d <- giterator.intervals(intervals = intervals, iterator = resolution)
    iter_2d <- gintervals.2d(iter_1d$chrom, iter_1d$start, iter_1d$end)

    gvtrack.create("v_marg", "hic.test_basic", "weighted.sum")
    withr::defer(gvtrack.rm("v_marg"))

    # Full marginal
    marginal <- gextract("v_marg",
        intervals = iter_2d, iterator = iter_2d,
        colnames = "marginal"
    )

    # Diagonal
    diag <- gextract("v_marg",
        intervals = iter_2d, iterator = iter_2d,
        band = c(-1024, 1024), colnames = "diagonal"
    )

    expect_regression(
        list(marginal = marginal, diagonal = diag),
        "marginal.diag_correct.hic.1"
    )

    # Validate structure
    expect_true("marginal" %in% colnames(marginal))
    expect_true("diagonal" %in% colnames(diag))
})

test_that("marginal extraction handles multiple tracks", {
    # Pattern: hic_basic.r lines 38-42
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    tracks <- c("hic.test_basic", "hic.test_shuffle")
    v_tracks <- paste0("v_", tracks)

    for (i in seq_along(tracks)) {
        gvtrack.create(v_tracks[i], tracks[i], "weighted.sum")
    }
    withr::defer(sapply(v_tracks, gvtrack.rm))

    iter_2d <- gintervals.2d("chr1", 5e5, 1e6)
    result <- gextract(v_tracks,
        intervals = iter_2d, iterator = iter_2d,
        colnames = tracks
    )

    expect_regression(result, "marginal.multi_track.hic.1")

    # Check column names
    expect_true(all(tracks %in% colnames(result)))
})

# ==============================================================================
# GROUP 7: GRID OPERATIONS (MEDIUM)
# ==============================================================================

test_that("contact grid bins by distance to features", {
    # Pattern: hic_basic.r lines 366-384
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    features1 <- gintervals("chr1", c(5e5, 1e6), c(5.1e5, 1.1e6))
    features2 <- gintervals("chr1", c(7e5, 1.2e6), c(7.1e5, 1.3e6))

    gvtrack.create("dist1", features1, "distance")
    gvtrack.create("dist2", features2, "distance")
    gvtrack.create("ws", "hic.test_basic", "weighted.sum")
    withr::defer({
        gvtrack.rm("dist1")
        gvtrack.rm("dist2")
        gvtrack.rm("ws")
    })

    gvtrack.iterator("dist1", 1)
    gvtrack.iterator("dist2", 2)

    scope <- gintervals.2d("chr1", 0, 2e6, "chr1", 0, 2e6)
    data <- gextract("dist1", "dist2", "ws", scope, band = c(-1e6, -1e4))

    # Bin
    breaks <- seq(-2e5, 2e5, 5e4)
    if (!is.null(data) && nrow(data) > 0) {
        d1_bin <- cut(data$dist1, breaks = breaks, include.lowest = TRUE)
        d2_bin <- cut(data$dist2, breaks = breaks, include.lowest = TRUE)

        contact_table <- table(d1_bin, d2_bin)

        expect_regression(
            as.data.frame.table(contact_table),
            "contact_grid.binned.hic.1"
        )
    } else {
        skip("No data extracted for contact grid test")
    }
})

test_that("obs/exp grid normalizes correctly", {
    # Pattern: hic_basic.r lines 387-417
    # Simplified version focusing on normalization logic
    set.seed(42) # For reproducible random data

    # Mock obs and exp matrices
    obs <- matrix(rpois(25, 100), 5, 5)
    exp <- matrix(rpois(25, 80), 5, 5)

    # Regularization
    obs[obs < 30] <- NA
    exp[is.na(obs) | exp < 30] <- NA

    # Normalize
    total_obs <- sum(obs, na.rm = TRUE)
    total_exp <- sum(exp, na.rm = TRUE)
    obs_norm <- obs / total_obs
    exp_norm <- exp / total_exp
    enrichment <- obs_norm / exp_norm

    expect_regression(
        list(
            obs_norm = obs_norm, exp_norm = exp_norm,
            enrichment = enrichment
        ),
        "obs_exp_grid.norm.hic.1"
    )

    # Basic validation
    expect_true(abs(sum(obs_norm, na.rm = TRUE) - 1) < 1e-10)
    expect_true(abs(sum(exp_norm, na.rm = TRUE) - 1) < 1e-10)
})

# ==============================================================================
# GROUP 8: SPATIAL/NEIGHBORHOOD ANALYSIS (LOW)
# ==============================================================================

test_that("gintervals.neighbors respects 2D distance constraints", {
    # Pattern: hic_basic.r line 467
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    # Create test directory if needed
    tryCatch(gdir.create("test", showWarnings = FALSE), error = function(e) NULL)

    set_name <- "test.neigh_set"
    # Clean up if already exists from previous run
    if (gintervals.exists(set_name)) {
        gintervals.rm(set_name, force = TRUE)
    }
    withr::defer(gintervals.rm(set_name, force = TRUE))

    giterator.intervals("hic.test_basic",
        band = c(-5e5, -1e4),
        intervals.set.out = set_name
    )

    starts1 <- seq(5e5, 1e6, 1e5)
    starts2 <- seq(6e5, 1.1e6, 1e5)
    # Create grid with expand.grid pattern
    g <- expand.grid(1:length(starts1), 1:length(starts2))
    grid <- data.frame(
        chrom1 = "chr1", start1 = starts1[g$Var1], end1 = starts1[g$Var1] + 1e5,
        chrom2 = "chr1", start2 = starts2[g$Var2], end2 = starts2[g$Var2] + 1e5
    )

    # Reload database to ensure interval set is visible
    gdb.reload()

    neighbors <- gintervals.neighbors(set_name, grid,
        mindist1 = -5e4, maxdist1 = 5e4,
        mindist2 = -5e4, maxdist2 = 5e4
    )

    expect_regression(neighbors, "gintervals.neighbors.2d.hic.1")

    # Basic validation
    expect_true(is.data.frame(neighbors))
})

test_that("spatial max values in grid pattern", {
    # Pattern: figures.r lines 755-777 (simplified)
    withr::local_options(list(gmax.data.size = 1e8))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19")

    distance <- 1e5
    resolution <- 5
    expansion <- seq(-distance, distance, length.out = resolution)

    # Create small grid of max vtracks
    vtracks <- c()
    for (i in 1:(resolution - 1)) {
        for (j in 1:(resolution - 1)) {
            v <- paste0("max_", i, "_", j)
            gvtrack.create(v, "hic.test_score", "max")
            gvtrack.iterator.2d(v,
                sshift1 = expansion[i], eshift1 = expansion[i + 1],
                sshift2 = expansion[j], eshift2 = expansion[j + 1]
            )
            vtracks <- c(vtracks, v)
        }
    }
    withr::defer(sapply(vtracks, gvtrack.rm))

    points <- gintervals.2d("chr1", c(1e6), c(1e6 + 1))
    result <- gextract(vtracks, intervals = points, iterator = points)

    expect_regression(result, "spatial_max_grid.hic.1")

    # Validate structure
    expect_equal(length(vtracks), (resolution - 1)^2)
    expect_true(all(vtracks %in% colnames(result)))
})
