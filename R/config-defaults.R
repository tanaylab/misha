# Centralized configuration defaults for misha package
#
# This file defines the single source of truth for all default configuration values.
# Both R code and C++ code should use consistent defaults matching these values.
#
# NOTE: If you modify a default here, ensure the corresponding constant in
# src/ConfigurationDefaults.h is also updated to maintain consistency.

#' Centralized configuration defaults
#'
#' This list contains all default values for misha configuration options.
#' Use `.ggetOption()` to retrieve options with these defaults as fallbacks.
#'
#' @keywords internal
#' @noRd
.misha_defaults <- list(
    # Memory and buffer sizes
    gmax.data.size = NULL, # Auto-calculated at load time based on system RAM
    gbig.intervals.size = 1000000, # Threshold for converting to disk-based "big" format
    gmax.mem.usage = NULL, # Auto-calculated: 80% of system RAM (min 10 GB)
    gbuf.size = 1000, # Evaluation buffer size for vectorized R expression evaluation
    gmultitask.max.records.factor = 64, # Inflate multitask max_records estimates to avoid under-allocation

    # Quantile computation
    gquantile.edge.data.size = 100000, # Buffer size for high-precision edge values
    gpv.middle.size = 0.96, # P-value middle region size
    gpv.middle.precision = 1e-4, # P-value middle region precision
    gpv.edge.precision = 1e-9, # P-value edge precision

    # Track chunk sizes (for 2D tracks)
    gtrack.chunk.size = 100000, # Chunk size in base pairs
    gtrack.num.chunks = 0, # Number of chunks (0 = auto)

    # Parallelization
    gmax.processes = NULL, # Auto-calculated: 70% of cores
    gmax.processes2core = 2, # Max processes per core
    gmin.scope4process = 10000, # Min genomic scope per parallel process
    gmultitasking = TRUE # Enable/disable multitasking
)
