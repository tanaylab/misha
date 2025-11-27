# Memory utilities for auto-sizing gmax.data.size

#' Get system memory in bytes
#'
#' @return Total system memory in bytes, or 8GB fallback if detection fails
#' @keywords internal
#' @noRd
get_system_memory <- function() {
    tryCatch(
        {
            # Use ps package for cross-platform system memory detection
            mem_info <- ps::ps_system_memory()

            # ps_system_memory() returns a named list with 'total' in bytes
            if (!is.null(mem_info$total) && mem_info$total > 0) {
                return(mem_info$total)
            }
        },
        error = function(e) {
            # Silently fall through to fallback
        }
    )

    # Fallback: conservative 8GB default if ps fails
    return(8e9)
}

#' Calculate optimal gmax.data.size based on system memory and cores
#'
#' @param sys_mem Total system memory in bytes (from get_system_memory())
#' @return Recommended gmax.data.size value
#' @keywords internal
#' @noRd
calculate_optimal_gmax_data_size <- function(sys_mem) {
    # Get number of cores for process estimation
    num_cores <- parallel::detectCores(logical = TRUE)
    if (is.na(num_cores) || num_cores < 1) {
        num_cores <- 1
    }

    # Estimate max processes that will be used
    # Auto-configured to 70% of cores
    estimated_processes <- as.integer(num_cores * 0.7)
    if (estimated_processes < 1) {
        estimated_processes <- 1
    }

    # Calculate buffer size:
    # - Total memory budget: 70% of RAM
    # - Divide by number of processes (since each can use gmax.data.size)
    # - Cap at practical mmap limit (10GB proven to work)
    #
    # This ensures: gmax.processes * gmax.data.size <= 70% of RAM
    #
    # Examples:
    # - 4 cores, 8GB RAM: 2 procs * 2.8GB = 5.6GB (70% of 8GB)
    # - 16 cores, 32GB RAM: 11 procs * 2.0GB = 22.4GB (70% of 32GB)
    # - 128 cores, 256GB RAM: 89 procs * 2.0GB = 179.2GB (70% of 256GB)
    #
    # The cap of 10GB is empirically determined - larger allocations fail
    # even though system limits permit them
    total_ram_budget <- sys_mem * 0.7
    per_process_buffer <- total_ram_budget / estimated_processes
    optimal_size <- min(per_process_buffer, 10e9)

    # Note: We do NOT enforce a minimum buffer size (e.g., 1GB) because that
    # would violate the budget constraint on high-core/low-RAM systems.
    # Example: 128 cores, 16GB RAM would give 89 * 1GB = 89GB > 11.2GB budget
    # It's better to have a small buffer that respects the constraint than to
    # violate the memory budget and risk OOM.

    return(optimal_size)
}

#' Format bytes as human-readable string
#'
#' @param bytes Number of bytes
#' @return Formatted string like "1.5 GB"
#' @keywords internal
#' @noRd
format_bytes <- function(bytes) {
    if (bytes >= 1e12) {
        sprintf("%.1f TB", bytes / 1e12)
    } else if (bytes >= 1e9) {
        sprintf("%.1f GB", bytes / 1e9)
    } else if (bytes >= 1e6) {
        sprintf("%.1f MB", bytes / 1e6)
    } else if (bytes >= 1e3) {
        sprintf("%.1f KB", bytes / 1e3)
    } else {
        sprintf("%d bytes", as.integer(bytes))
    }
}
