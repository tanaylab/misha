# Stress tests for auto-configuration
# These tests verify edge cases and stress scenarios

test_that("low RAM systems respect total memory budget", {
    # Edge case: 2GB RAM, 4 cores
    # Per-process: (2GB * 0.7) / 2 = 0.7GB
    # Minimum: 1GB
    # Problem: 2 * 1GB = 2GB > 1.4GB (70% of 2GB)

    with_mocked_bindings(
        detectCores = function(...) 4,
        .package = "parallel",
        {
            size <- misha:::calculate_optimal_gmax_data_size(2e9)
            procs <- as.integer(4 * 0.7) # 2 processes

            total_usage <- procs * size
            budget <- 2e9 * 0.7

            expect_true(total_usage <= budget,
                info = sprintf(
                    "Total usage (%.2f GB) exceeds budget (%.2f GB)",
                    total_usage / 1e9, budget / 1e9
                )
            )
        }
    )
})

test_that("very low RAM systems stay within budget", {
    # Extreme edge case: 1GB RAM, 8 cores
    with_mocked_bindings(
        detectCores = function(...) 8,
        .package = "parallel",
        {
            size <- misha:::calculate_optimal_gmax_data_size(1e9)
            procs <- as.integer(8 * 0.7) # 5 processes

            total_usage <- procs * size
            budget <- 1e9 * 0.7

            expect_true(total_usage <= budget,
                info = sprintf(
                    "Total usage (%.2f GB) exceeds budget (%.2f GB)",
                    total_usage / 1e9, budget / 1e9
                )
            )
        }
    )
})

test_that("budget constraint holds for various system configurations", {
    # Test many different RAM/core combinations
    test_configs <- expand.grid(
        ram_gb = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024),
        cores = c(1, 2, 4, 8, 16, 32, 64, 128)
    )

    for (i in 1:nrow(test_configs)) {
        ram <- test_configs$ram_gb[i] * 1e9
        cores <- test_configs$cores[i]

        with_mocked_bindings(
            detectCores = function(...) cores,
            .package = "parallel",
            {
                size <- misha:::calculate_optimal_gmax_data_size(ram)
                procs <- as.integer(cores * 0.7)
                if (procs < 1) procs <- 1

                total_usage <- procs * size
                budget <- ram * 0.7

                expect_true(total_usage <= budget * 1.01, # Allow 1% tolerance for rounding
                    info = sprintf(
                        "Config: %.0f GB RAM, %d cores -> %d procs * %.2f GB = %.2f GB > %.2f GB budget",
                        ram / 1e9, cores, procs, size / 1e9,
                        total_usage / 1e9, budget / 1e9
                    )
                )
            }
        )
    }
})

test_that("parallel execution respects buffer size on real database", {
    skip_if_not(
        dir.exists("/net/mraid20/ifs/wisdom/tanay_lab/tgdata/db/tgdb/misha_test_db_indexed/"),
        "Real test database not available"
    )

    # Load real database
    gdb.init("/net/mraid20/ifs/wisdom/tanay_lab/tgdata/db/tgdb/misha_test_db_indexed/")

    tracks <- gtrack.ls()
    skip_if(length(tracks) == 0, "No tracks in database")

    test_track <- tracks[1]

    # Try a moderately large query with multitasking enabled
    options(gmultitasking = TRUE)

    # This should succeed without exceeding buffer size
    result <- tryCatch(
        {
            gextract(test_track, gintervals.all(), iterator = 10000)
        },
        error = function(e) {
            if (grepl("Result size.*exceeded", e$message)) {
                fail(sprintf("Buffer size exceeded with auto-configured settings: %s", e$message))
            }
            stop(e)
        }
    )

    expect_true(nrow(result) > 0)
})

test_that("single-core systems get reasonable configuration", {
    with_mocked_bindings(
        detectCores = function(...) 1,
        .package = "parallel",
        {
            # 8GB RAM, 1 core
            size <- misha:::calculate_optimal_gmax_data_size(8e9)
            procs <- as.integer(1 * 0.7) # 0 -> should be 1
            if (procs < 1) procs <- 1

            # Should get 70% of RAM but capped at 10GB
            expect_equal(size, min(8e9 * 0.7, 10e9))

            # Budget check
            total_usage <- procs * size
            budget <- 8e9 * 0.7
            expect_true(total_usage <= budget)
        }
    )
})

test_that("high-core low-RAM systems are handled correctly", {
    # Pathological case: 128 cores, 16GB RAM
    with_mocked_bindings(
        detectCores = function(...) 128,
        .package = "parallel",
        {
            size <- misha:::calculate_optimal_gmax_data_size(16e9)
            procs <- as.integer(128 * 0.7) # 89 processes

            # Per-process: (16GB * 0.7) / 89 = 0.126 GB
            # Should be enforced to minimum 1GB
            # But total would be 89 * 1GB = 89GB >> 16GB!

            total_usage <- procs * size
            budget <- 16e9 * 0.7

            expect_true(total_usage <= budget * 1.01,
                info = sprintf(
                    "High-core low-RAM: %d procs * %.2f GB = %.2f GB > %.2f GB budget",
                    procs, size / 1e9, total_usage / 1e9, budget / 1e9
                )
            )
        }
    )
})

test_that("formula handles integer truncation correctly", {
    # Test that we don't lose too much precision with integer truncation
    test_cases <- list(
        list(cores = 3, expected_procs = 2), # 3 * 0.7 = 2.1 -> 2
        list(cores = 7, expected_procs = 4), # 7 * 0.7 = 4.9 -> 4
        list(cores = 11, expected_procs = 7), # 11 * 0.7 = 7.7 -> 7
        list(cores = 13, expected_procs = 9) # 13 * 0.7 = 9.1 -> 9
    )

    for (tc in test_cases) {
        procs <- as.integer(tc$cores * 0.7)
        expect_equal(procs, tc$expected_procs,
            info = sprintf("cores=%d", tc$cores)
        )
    }
})
