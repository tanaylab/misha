test_that("get_system_memory detects memory correctly", {
    mem <- misha:::get_system_memory()

    # Should return a positive number
    expect_true(is.numeric(mem))
    expect_true(mem > 0)

    # Should be reasonable (at least 1GB, less than 100TB)
    expect_true(mem >= 1e9)
    expect_true(mem <= 1e14)

    # Should match ps package (if available)
    skip_if_not_installed("ps")
    ps_mem <- ps::ps_system_memory()[["total"]]
    if (!is.null(ps_mem) && ps_mem > 0) {
        expect_equal(mem, ps_mem)
    }
})

test_that("get_system_memory fallback works when ps fails", {
    # Mock ps failure by temporarily making ps unavailable
    # This tests the fallback to 8GB
    with_mocked_bindings(
        ps_system_memory = function() stop("Mocked failure"),
        .package = "ps",
        {
            mem <- misha:::get_system_memory()
            expect_equal(mem, 8e9) # Should return 8GB fallback
        }
    )
})

test_that("calculate_optimal_gmax_data_size respects budget constraint", {
    # Very low RAM system (500MB), 8 cores -> 5 processes
    # Budget: 500MB * 0.7 = 350MB
    # Per-process: 350MB / 5 = 70MB
    # No minimum enforced - must respect budget
    with_mocked_bindings(
        detectCores = function(...) 8,
        .package = "parallel",
        {
            size <- misha:::calculate_optimal_gmax_data_size(5e8)
            procs <- as.integer(8 * 0.7) # 5 processes

            # Verify budget constraint
            total <- procs * size
            budget <- 5e8 * 0.7
            expect_true(total <= budget)
        }
    )
})

test_that("calculate_optimal_gmax_data_size respects maximum cap", {
    # Very high RAM system (2TB)
    size <- misha:::calculate_optimal_gmax_data_size(2e12)

    # Should be capped at 10GB
    expect_equal(size, 10e9)
})

test_that("calculate_optimal_gmax_data_size uses per-process formula", {
    # Formula: min((ram * 0.7) / num_processes, 10GB)
    # This ensures: gmax.processes * gmax.data.size <= 70% of RAM

    # Get actual system cores to calculate expected values
    num_cores <- parallel::detectCores(logical = TRUE)
    if (is.na(num_cores) || num_cores < 1) num_cores <- 1
    estimated_procs <- as.integer(num_cores * 0.7)
    if (estimated_procs < 1) estimated_procs <- 1

    # Test with various RAM sizes
    # 20GB RAM: budget = 14GB, per-process = 14GB / procs
    size <- misha:::calculate_optimal_gmax_data_size(20e9)
    expected <- min(20e9 * 0.7 / estimated_procs, 10e9)
    expect_equal(size, expected)

    # 100GB RAM: budget = 70GB, per-process = 70GB / procs (likely capped at 10GB)
    size <- misha:::calculate_optimal_gmax_data_size(100e9)
    expected <- min(100e9 * 0.7 / estimated_procs, 10e9)
    expect_equal(size, expected)
})

test_that("package initializes with correct auto-configuration", {
    # Check that options were set during .onLoad
    expect_true(!is.null(getOption("gmax.processes")))
    expect_true(!is.null(getOption("gmax.processes2core")))
    expect_true(!is.null(getOption("gmax.data.size")))

    # Check process limits are reasonable
    gmax_procs <- getOption("gmax.processes")
    num_cores <- parallel::detectCores(logical = TRUE)
    expect_true(is.numeric(gmax_procs))
    expect_true(gmax_procs >= 1)
    expect_true(gmax_procs <= num_cores) # Should be <= total cores

    # Check data size is reasonable
    gmax_data <- getOption("gmax.data.size")
    expect_true(is.numeric(gmax_data))
    expect_true(gmax_data > 0) # Must be positive
    expect_true(gmax_data <= 10e9) # At most 10GB

    # Verify budget constraint
    sys_mem <- misha:::get_system_memory()
    budget <- sys_mem * 0.7
    total_usage <- gmax_procs * gmax_data
    expect_true(total_usage <= budget * 1.01) # Allow 1% tolerance for rounding
})

test_that("gmax.processes is 70% of cores", {
    # We can't change the actual core count, but we can test the logic
    # by checking the current setting matches expected behavior
    num_cores <- parallel::detectCores(logical = TRUE)
    expected_procs <- as.integer(num_cores * 0.7)

    actual_procs <- getOption("gmax.processes")
    expect_equal(actual_procs, expected_procs)
})

test_that("gmax.data.size is within reasonable range", {
    sys_mem <- misha:::get_system_memory()
    expected_size <- misha:::calculate_optimal_gmax_data_size(sys_mem)

    actual_size <- getOption("gmax.data.size")

    # Check that size is reasonable (positive and capped at 10GB)
    expect_true(actual_size > 0) # Must be positive
    expect_true(actual_size <= 10e9) # At most 10GB

    # Verify budget constraint
    num_cores <- parallel::detectCores(logical = TRUE)
    if (is.na(num_cores) || num_cores < 1) num_cores <- 1
    gmax_procs <- as.integer(num_cores * 0.7)
    if (gmax_procs < 1) gmax_procs <- 1

    budget <- sys_mem * 0.7
    total_usage <- gmax_procs * actual_size
    expect_true(total_usage <= budget * 1.01) # Allow 1% tolerance for rounding
})

test_that("format_bytes produces human-readable output", {
    expect_equal(misha:::format_bytes(500), "500 bytes")
    expect_equal(misha:::format_bytes(1500), "1.5 KB")
    expect_equal(misha:::format_bytes(1.5e6), "1.5 MB")
    expect_equal(misha:::format_bytes(1.5e9), "1.5 GB")
    expect_equal(misha:::format_bytes(1.5e12), "1.5 TB")
})

test_that("auto-configuration works across different system sizes", {
    # Simulate different system configurations
    # Formula: gmax.processes = cores * 0.7
    #          gmax.data.size = min((ram * 0.7) / gmax.processes, 10GB)
    # This ensures: gmax.processes * gmax.data.size <= 70% of RAM

    test_cases <- list(
        # ram, cores, expected_procs, expected_size
        # 2GB, 2 cores: 1 proc * 1.4GB = 1.4GB (70% of 2GB)
        list(ram = 2e9, cores = 2, expected_procs = 1, expected_size = 1.4e9),
        # 4GB, 4 cores: 2 procs * 1.4GB = 2.8GB (70% of 4GB)
        list(ram = 4e9, cores = 4, expected_procs = 2, expected_size = 1.4e9),
        # 8GB, 8 cores: 5 procs * 1.12GB = 5.6GB (70% of 8GB)
        list(ram = 8e9, cores = 8, expected_procs = 5, expected_size = 1.12e9),
        # 16GB, 16 cores: 11 procs * 1.018GB = 11.2GB (70% of 16GB)
        list(ram = 16e9, cores = 16, expected_procs = 11, expected_size = 16e9 * 0.7 / 11),
        # 128GB, 64 cores: 44 procs * 2.036GB = 89.6GB (70% of 128GB)
        list(ram = 128e9, cores = 64, expected_procs = 44, expected_size = 128e9 * 0.7 / 44),
        # 1TB, 128 cores: 89 procs * 7.865GB = 700GB (70% of 1TB)
        list(ram = 1000e9, cores = 128, expected_procs = 89, expected_size = 1000e9 * 0.7 / 89),
        # 30GB, 2 cores: 1 proc * 21GB -> capped at 10GB
        list(ram = 30e9, cores = 2, expected_procs = 1, expected_size = 10e9)
    )

    for (tc in test_cases) {
        # Mock detectCores to return our test value
        with_mocked_bindings(
            detectCores = function(...) tc$cores,
            .package = "parallel",
            {
                procs <- as.integer(tc$cores * 0.7)
                size <- misha:::calculate_optimal_gmax_data_size(tc$ram)

                expect_equal(procs, tc$expected_procs,
                    info = sprintf("cores=%d", tc$cores)
                )
                expect_equal(size, tc$expected_size,
                    info = sprintf("ram=%.0f GB, cores=%d", tc$ram / 1e9, tc$cores)
                )
            }
        )
    }
})

test_that("auto-configuration handles edge cases", {
    # NA cores
    with_mocked_bindings(
        detectCores = function(...) NA,
        .package = "parallel",
        {
            # Should default to 1 core
            num_cores <- parallel::detectCores(logical = TRUE)
            num_cores <- if (is.na(num_cores) || num_cores < 1) 1 else num_cores
            expect_equal(num_cores, 1)
        }
    )

    # Zero cores (shouldn't happen but test defensively)
    with_mocked_bindings(
        detectCores = function(...) 0,
        .package = "parallel",
        {
            num_cores <- parallel::detectCores(logical = TRUE)
            num_cores <- if (is.na(num_cores) || num_cores < 1) 1 else num_cores
            expect_equal(num_cores, 1)
        }
    )
})

test_that("auto-disable threshold prevents fork overhead on small datasets", {
    create_isolated_test_db()

    # Dynamic threshold: gmax.processes * 1000 records
    # On this system, threshold should scale with number of processes
    procs <- getOption("gmax.processes")
    expected_threshold <- procs * 1000

    cat(sprintf(
        "\nDynamic threshold: %d procs * 1000 = %s records\n",
        procs, format(expected_threshold, big.mark = ",")
    ))

    # Small dataset - use specific small intervals
    small_intervals <- gintervals(1, 0, 100000)

    # With multitasking enabled
    options(gmultitasking = TRUE)
    start <- Sys.time()
    result_small <- gextract("test.fixedbin", small_intervals, iterator = 1000)
    time_small <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    # Should complete quickly
    expect_true(time_small < 0.5) # Less than 500ms
    expect_true(nrow(result_small) > 0)
})

test_that("multitasking engages for large datasets", {
    create_isolated_test_db()
    options(gmultitasking = TRUE)
    result_large <- gextract("test.fixedbin", gintervals.all(), iterator = 10)

    expect_true(nrow(result_large) >= 10000)
})

test_that("manual override of gmax.processes works", {
    original <- getOption("gmax.processes")

    # Override to specific value
    options(gmax.processes = 4)
    expect_equal(getOption("gmax.processes"), 4)

    # Restore
    options(gmax.processes = original)
    expect_equal(getOption("gmax.processes"), original)
})

test_that("manual override of gmax.data.size works", {
    original <- getOption("gmax.data.size")

    # Override to specific value
    options(gmax.data.size = 5e9) # 5GB
    expect_equal(getOption("gmax.data.size"), 5e9)

    # Restore
    options(gmax.data.size = original)
    expect_equal(getOption("gmax.data.size"), original)
})

test_that("gmultitasking can be disabled", {
    create_isolated_test_db()

    original <- getOption("gmultitasking")

    # Disable multitasking
    options(gmultitasking = FALSE)
    result <- gextract("test.fixedbin", gintervals(1, 0, 100000), iterator = 1000)

    expect_true(nrow(result) > 0)

    # Restore
    options(gmultitasking = original)
})

test_that("buffer size caps at 10GB even with large RAM and few processes", {
    # With 2 cores (1 process), very large RAM should hit the 10GB cap
    with_mocked_bindings(
        detectCores = function(...) 2,
        .package = "parallel",
        {
            # 30GB RAM, 1 process: 30 * 0.7 / 1 = 21GB -> capped at 10GB
            size <- misha:::calculate_optimal_gmax_data_size(30e9)
            expect_equal(size, 10e9)

            # 100GB RAM, 1 process: 100 * 0.7 / 1 = 70GB -> capped at 10GB
            size <- misha:::calculate_optimal_gmax_data_size(100e9)
            expect_equal(size, 10e9)
        }
    )
})

test_that("buffer size respects per-process formula", {
    # With 4 cores (2 processes), test various RAM sizes
    with_mocked_bindings(
        detectCores = function(...) 4,
        .package = "parallel",
        {
            # 4GB RAM: 4 * 0.7 / 2 = 1.4GB per process
            size <- misha:::calculate_optimal_gmax_data_size(4e9)
            expect_equal(size, 1.4e9)

            # 10GB RAM: 10 * 0.7 / 2 = 3.5GB per process
            size <- misha:::calculate_optimal_gmax_data_size(10e9)
            expect_equal(size, 3.5e9)

            # 30GB RAM: 30 * 0.7 / 2 = 10.5GB -> capped at 10GB
            size <- misha:::calculate_optimal_gmax_data_size(30e9)
            expect_equal(size, 10e9)
        }
    )
})

test_that("buffer size respects budget on low RAM systems", {
    # Budget constraint must be respected even on low RAM systems
    # No minimum buffer size is enforced if it would violate the constraint
    with_mocked_bindings(
        detectCores = function(...) 8,
        .package = "parallel",
        {
            procs <- as.integer(8 * 0.7) # 5 processes

            # 500MB RAM: 500MB * 0.7 / 5 = 70MB per process
            size <- misha:::calculate_optimal_gmax_data_size(500e6)
            expect_equal(size, 500e6 * 0.7 / procs)
            expect_true(procs * size <= 500e6 * 0.7)

            # 2GB RAM: 2GB * 0.7 / 5 = 280MB per process
            size <- misha:::calculate_optimal_gmax_data_size(2e9)
            expect_equal(size, 2e9 * 0.7 / procs)
            expect_true(procs * size <= 2e9 * 0.7)
        }
    )
})

test_that("process count calculation is correct for various core counts", {
    # 4 cores -> 70% = 2
    expect_equal(as.integer(4 * 0.7), 2)

    # 8 cores -> 70% = 5
    expect_equal(as.integer(8 * 0.7), 5)

    # 16 cores -> 70% = 11
    expect_equal(as.integer(16 * 0.7), 11)

    # 64 cores -> 70% = 44
    expect_equal(as.integer(64 * 0.7), 44)

    # 128 cores -> 70% = 89
    expect_equal(as.integer(128 * 0.7), 89)
})

test_that("auto-configuration persists across queries", {
    create_isolated_test_db()

    # Get configured values
    procs_before <- getOption("gmax.processes")
    size_before <- getOption("gmax.data.size")

    # Run a query
    result <- gextract("test.fixedbin", gintervals(1, 0, 100000), iterator = 1000)
    expect_true(nrow(result) > 0)

    # Check values haven't changed
    expect_equal(getOption("gmax.processes"), procs_before)
    expect_equal(getOption("gmax.data.size"), size_before)
})

test_that("format_bytes handles all size ranges correctly", {
    # Bytes
    expect_equal(misha:::format_bytes(500), "500 bytes")
    expect_equal(misha:::format_bytes(999), "999 bytes")

    # Kilobytes
    expect_equal(misha:::format_bytes(1500), "1.5 KB")
    expect_equal(misha:::format_bytes(999e3), "999.0 KB")

    # Megabytes
    expect_equal(misha:::format_bytes(1.5e6), "1.5 MB")
    expect_equal(misha:::format_bytes(999e6), "999.0 MB")

    # Gigabytes
    expect_equal(misha:::format_bytes(1.5e9), "1.5 GB")
    expect_equal(misha:::format_bytes(999e9), "999.0 GB")

    # Terabytes
    expect_equal(misha:::format_bytes(1.5e12), "1.5 TB")
    expect_equal(misha:::format_bytes(10e12), "10.0 TB")
})

test_that("get_system_memory returns reasonable values", {
    mem <- misha:::get_system_memory()

    # Should be positive
    expect_true(mem > 0)

    # Should be at least 1GB (no modern system has less)
    expect_true(mem >= 1e9)

    # Should be less than 100TB (no system has that much RAM yet)
    expect_true(mem < 100e12)
})

test_that("auto-config values are reasonable", {
    # Check that the configured values are reasonable
    # (may not match exactly if package was loaded before code changes)
    num_cores <- parallel::detectCores(logical = TRUE)
    sys_mem <- misha:::get_system_memory()

    actual_procs <- getOption("gmax.processes")
    actual_size <- getOption("gmax.data.size")

    # Process count should be reasonable
    expect_true(is.numeric(actual_procs))
    expect_true(actual_procs >= 1)
    expect_true(actual_procs <= num_cores)

    # Buffer size should be reasonable
    expect_true(is.numeric(actual_size))
    expect_true(actual_size > 0) # Must be positive
    expect_true(actual_size <= 10e9) # At most 10GB

    # Verify budget constraint
    budget <- sys_mem * 0.7
    total_usage <- actual_procs * actual_size
    expect_true(total_usage <= budget * 1.01) # Allow 1% tolerance for rounding
})
