# Genome Synthesis: Generate synthetic genomes from stratified Markov models

# Helper function to compute flat bin indices from per-dimension indices
# Uses vectorized matrix multiplication instead of loop-based stride computation
.compute_flat_indices <- function(per_dim_indices, dim_sizes) {
    n_dims <- length(dim_sizes)
    n_positions <- nrow(per_dim_indices)

    flat_indices <- rep(NA_integer_, n_positions)

    # Check which positions have valid indices in all dimensions
    valid_positions <- rowSums(is.na(per_dim_indices)) == 0L

    if (any(valid_positions)) {
        # Precompute strides: c(1, s1, s1*s2, ...)
        strides <- c(1L, cumprod(dim_sizes[-n_dims]))

        # Get valid indices matrix (already 1-based)
        valid_per_dim <- per_dim_indices[valid_positions, , drop = FALSE]

        # Convert to 0-based for flat index computation: (idx - 1) * stride
        # Then sum across dimensions and convert back to 1-based
        # flat_idx = sum((idx_d - 1) * stride_d) + 1
        #          = sum(idx_d * stride_d) - sum(stride_d) + 1
        flat_idx <- as.integer(valid_per_dim %*% strides - sum(strides) + 1L)

        flat_indices[valid_positions] <- flat_idx
    }

    flat_indices
}

# Default chunk size threshold for parallel processing (1 billion bases)
.GSYNTH_MAX_CHUNK_SIZE <- 1e9

#' Process large genomes in parallel chunks
#'
#' Internal helper function that extracts the common parallel processing pattern
#' from gsynth functions. Handles both file and vector output modes.
#'
#' @param intervals Genomic intervals to process
#' @param output_format Output format: "misha" (binary), "fasta" (text), or "vector"
#' @param output_path Path to output file (ignored for vector mode)
#' @param process_chunk_fn Function to process each chunk. Should accept:
#'   \itemize{
#'     \item chunk_interval: Single-row intervals data frame
#'     \item chunk_index: Integer index of the chunk (1-based)
#'     \item chunk_seed: Optional seed for this chunk (if seed provided)
#'     \item ...: Additional arguments passed through
#'   }
#'   For file output, should return path to temp file.
#'   For vector output, should return character vector.
#' @param max_chunk_size Maximum total bases before triggering parallel processing
#' @param seed Optional seed for reproducible chunk seeds
#' @param ... Additional arguments passed to process_chunk_fn
#'
#' @return NULL if parallel processing not needed (genome too small).
#'         For file output: invisible(NULL) after writing combined file.
#'         For vector output: Character vector of combined results.
#'
#' @noRd
.gsynth_process_parallel <- function(intervals,
                                     output_format,
                                     output_path,
                                     process_chunk_fn,
                                     max_chunk_size = .GSYNTH_MAX_CHUNK_SIZE,
                                     seed = NULL,
                                     ...) {
    # Calculate total bases
    total_bases <- sum(intervals$end - intervals$start)

    # Return NULL if genome is small enough to process normally
    if (total_bases <= max_chunk_size) {
        return(NULL)
    }

    # Number of chunks = number of intervals (one per row)
    n_chunks <- nrow(intervals)

    message(sprintf(
        "Large genome detected (%s bases). Processing %d chunks in parallel...",
        format(total_bases, big.mark = ","), n_chunks
    ))

    # Determine number of cores to use
    n_cores <- min(getOption("gmax.processes", 1L), n_chunks)
    message(sprintf("Using %d cores for parallel processing", n_cores))

    # Generate reproducible seeds for each chunk if seed is provided
    chunk_seeds <- if (!is.null(seed)) {
        set.seed(seed)
        sample.int(.Machine$integer.max, n_chunks)
    } else {
        rep(NA, n_chunks)
    }

    # Capture extra arguments
    extra_args <- list(...)

    if (output_format == "vector") {
        # Vector output mode: process chunks and combine results
        all_results <- parallel::mclapply(seq_len(n_chunks), function(i) {
            chunk_interval <- intervals[i, , drop = FALSE]
            chunk_seed <- if (!is.na(chunk_seeds[i])) chunk_seeds[i] else NULL

            # Call process function with all arguments
            do.call(process_chunk_fn, c(
                list(
                    chunk_interval = chunk_interval,
                    chunk_index = i,
                    chunk_seed = chunk_seed
                ),
                extra_args
            ))
        }, mc.cores = n_cores, mc.preschedule = FALSE)

        # Combine vector results
        result <- unlist(all_results)
        return(result)
    }

    # File output mode: process to temp files and combine
    temp_files <- parallel::mclapply(seq_len(n_chunks), function(i) {
        chunk_interval <- intervals[i, , drop = FALSE]
        chunk_seed <- if (!is.na(chunk_seeds[i])) chunk_seeds[i] else NULL

        # Call process function with all arguments
        do.call(process_chunk_fn, c(
            list(
                chunk_interval = chunk_interval,
                chunk_index = i,
                chunk_seed = chunk_seed
            ),
            extra_args
        ))
    }, mc.cores = n_cores, mc.preschedule = FALSE)

    # Combine results in order
    message("Combining results...")

    for (i in seq_along(temp_files)) {
        temp_file <- temp_files[[i]]

        if (!file.exists(temp_file)) {
            stop(sprintf("Failed to generate chunk %d", i), call. = FALSE)
        }

        if (i == 1) {
            # First chunk: copy to output
            file.copy(temp_file, output_path, overwrite = TRUE)
        } else {
            # Subsequent chunks: append
            if (output_format == "fasta") {
                # Text mode: read lines and append
                temp_content <- readLines(temp_file)
                cat(temp_content, file = output_path, sep = "\n", append = TRUE)
            } else {
                # Binary mode (misha): read raw and append
                temp_content <- readBin(temp_file, "raw", n = file.info(temp_file)$size)
                con <- file(output_path, "ab")
                writeBin(temp_content, con)
                close(con)
            }
        }

        # Clean up temp file
        unlink(temp_file)
    }

    invisible(NULL)
}

#' Create a bin mapping from value-based merge specifications
#'
#' Converts value-based bin merge specifications into a bin_map named vector
#' that can be used with \code{\link{gsynth.train}}. This allows you to
#' specify merges using actual track values rather than bin indices.
#'
#' @param breaks Numeric vector of bin boundaries (same as used in
#'        \code{\link{gsynth.train}})
#' @param merge_ranges List of merge specifications. Each specification is a
#'        named list with:
#'   \describe{
#'     \item{from}{Numeric vector of length 2 \code{c(min, max)} defining the
#'           source value range to merge. Use \code{-Inf} or \code{Inf} for
#'           open-ended ranges. Can also be a single number (shorthand for
#'           \code{c(value, Inf)}).}
#'     \item{to}{Numeric vector of length 2 \code{c(min, max)} defining the
#'           target bin that source bins should map to. Must match an existing
#'           bin defined by \code{breaks}.}
#'   }
#'
#' @return A named vector (bin_map) compatible with \code{bin_map} parameter
#'         in \code{\link{gsynth.train}}. The names are source bin indices
#'         (1-based), and values are target bin indices (1-based).
#'
#' @examples
#' # Define breaks for GC content [0, 1] in 0.025 increments
#' breaks <- seq(0, 1, 0.025)
#'
#' # Merge all GC content above 70% (0.7) into the bin (0.675, 0.7]
#' bin_map <- gsynth.bin_map(
#'     breaks = breaks,
#'     merge_ranges = list(
#'         list(from = 0.7, to = c(0.675, 0.7))
#'     )
#' )
#'
#' # Multiple merges: merge low GC (< 0.3) and high GC (> 0.7) into middle bins
#' bin_map2 <- gsynth.bin_map(
#'     breaks = breaks,
#'     merge_ranges = list(
#'         list(from = c(-Inf, 0.3), to = c(0.4, 0.425)), # low GC -> (0.4, 0.425]
#'         list(from = 0.7, to = c(0.675, 0.7)) # high GC -> (0.675, 0.7]
#'     )
#' )
#'
#' @seealso \code{\link{gsynth.train}}, \code{\link{gsynth.cell_merge}},
#'          \code{\link{gsynth.sample}}
#' @export
gsynth.bin_map <- function(breaks, merge_ranges = NULL) {
    if (!is.numeric(breaks) || length(breaks) < 2) {
        stop("breaks must be a numeric vector with at least 2 elements", call. = FALSE)
    }
    breaks <- sort(breaks)
    num_bins <- length(breaks) - 1

    # Start with identity mapping
    bin_map_vec <- seq_len(num_bins)
    names(bin_map_vec) <- as.character(bin_map_vec)

    if (is.null(merge_ranges) || length(merge_ranges) == 0) {
        return(bin_map_vec)
    }

    # Find target bin indices for each target range
    find_bin_for_range <- function(range) {
        # range should be c(min, max)
        if (length(range) != 2) {
            stop("'to' range must be a vector of length 2: c(min, max)", call. = FALSE)
        }
        min_val <- range[1]
        max_val <- range[2]

        # Find which bin this range corresponds to
        # Bins are [breaks[i], breaks[i+1]) for i=1..num_bins, except last is closed
        for (i in seq_len(num_bins)) {
            bin_min <- breaks[i]
            bin_max <- breaks[i + 1]
            # Check if the range matches this bin (within tolerance)
            if (abs(min_val - bin_min) < 1e-10 && abs(max_val - bin_max) < 1e-10) {
                return(i)
            }
        }
        stop(sprintf(
            "Target range [%.6f, %.6f] does not match any bin defined by breaks",
            min_val, max_val
        ), call. = FALSE)
    }

    # Process each merge specification
    for (spec in merge_ranges) {
        if (!is.list(spec) || !("from" %in% names(spec)) || !("to" %in% names(spec))) {
            stop("Each merge specification must be a list with 'from' and 'to' elements", call. = FALSE)
        }

        from_range <- spec$from
        to_range <- spec$to

        # Handle shorthand: single number means [value, Inf)
        if (length(from_range) == 1 && is.numeric(from_range)) {
            from_range <- c(from_range, Inf)
        }

        if (length(from_range) != 2) {
            stop("'from' must be a vector of length 1 or 2", call. = FALSE)
        }

        from_min <- from_range[1]
        from_max <- from_range[2]

        # Find target bin index
        target_bin <- find_bin_for_range(to_range)

        # Find all source bins whose center is in from_range and map them to target_bin
        for (i in seq_len(num_bins)) {
            bin_min <- breaks[i]
            bin_max <- breaks[i + 1]
            bin_center <- (bin_min + bin_max) / 2

            # Check if bin center is within from_range
            # Handle Inf/-Inf by treating them as always true for that bound
            in_range <- TRUE
            if (is.finite(from_min)) {
                in_range <- in_range && (bin_center >= from_min)
            }
            if (is.finite(from_max)) {
                in_range <- in_range && (bin_center <= from_max)
            }

            if (in_range) {
                bin_map_vec[i] <- target_bin
            }
        }
    }

    # Return as named vector (source -> target)
    result <- as.integer(bin_map_vec)
    names(result) <- as.character(seq_len(num_bins))
    result
}

#' Resolve a cell-level merge specification into flat bin indices
#'
#' Unlike \code{\link{gsynth.bin_map}}, which merges bins independently along
#' each dimension, \code{gsynth.cell_merge} resolves per-\emph{joint-cell}
#' redirects: each entry redirects one specific training cell (identified by
#' its per-dimension values) to another specific training cell. This lets
#' callers redirect arbitrary cells whose Cartesian position cannot be
#' expressed via per-axis merges (e.g., \dQuote{cell (GC=0.725, CG=0.05) -->
#' cell (GC=0.70, CG=0.08)}).
#'
#' This function is primarily a utility for inspecting / debugging what
#' \code{\link{gsynth.sample}} will do when invoked with the \code{cell_merge}
#' argument. It returns a data frame describing, for every entry, the resolved
#' source and target cells in both per-dimension-bin and flat-bin space.
#'
#' @param model A \code{gsynth.model} object from \code{\link{gsynth.train}}.
#' @param cell_merge A list of redirect specifications. Each entry is a named
#'        list with:
#'   \describe{
#'     \item{from}{Numeric vector of length \code{n_dims}; one representative
#'           value per dimension that identifies the \emph{source} cell.}
#'     \item{to}{Numeric vector of length \code{n_dims}; one representative
#'           value per dimension that identifies the \emph{target} cell.}
#'   }
#'   A single data frame with columns \code{from_1, from_2, ..., to_1, to_2,
#'   ...} is also accepted and converted internally.
#' @param bin_merge Optional sampling-time bin merge specification (same format
#'        as in \code{\link{gsynth.sample}}). When supplied, source and target
#'        per-dimension bin indices are remapped through the resulting
#'        per-axis maps before being combined into flat indices, so that cell
#'        values reference post-bin_merge cells.
#'
#' @return A data frame with one row per \code{cell_merge} entry and columns:
#'   \itemize{
#'     \item \code{from_<d>}, \code{to_<d>} for each dimension \code{d}:
#'           1-based bin index after any bin_merge remapping.
#'     \item \code{source_flat}, \code{target_flat}: 1-based flat bin index
#'           into \code{model$model_data$cdf}.
#'   }
#'
#' @examples
#' \dontrun{
#' # Resolve a redirect table before handing it to gsynth.sample:
#' redirects <- list(
#'     list(from = c(0.725, 0.05), to = c(0.70, 0.08)),
#'     list(from = c(0.75, 0.06), to = c(0.70, 0.08))
#' )
#' resolved <- gsynth.cell_merge(model, redirects)
#' print(resolved)
#'
#' gsynth.sample(model, "out.fa",
#'     output_format = "fasta",
#'     cell_merge = redirects
#' )
#' }
#'
#' @seealso \code{\link{gsynth.bin_map}}, \code{\link{gsynth.sample}}
#' @export
gsynth.cell_merge <- function(model, cell_merge, bin_merge = NULL) {
    if (!inherits(model, "gsynth.model")) {
        stop("model must be a gsynth.model object", call. = FALSE)
    }
    n_dims <- model$n_dims
    if (is.null(n_dims) || n_dims < 1) {
        stop("cell_merge requires a stratified model (n_dims >= 1)", call. = FALSE)
    }

    cell_merge <- .cell_merge_normalize(cell_merge, n_dims)
    n_entries <- length(cell_merge)
    if (n_entries == 0) {
        return(.cell_merge_empty_result(n_dims))
    }

    sample_bin_maps <- .cell_merge_sample_bin_maps(model, bin_merge)
    dim_sizes <- model$dim_sizes

    source_dims <- matrix(NA_integer_, nrow = n_entries, ncol = n_dims)
    target_dims <- matrix(NA_integer_, nrow = n_entries, ncol = n_dims)

    for (d in seq_len(n_dims)) {
        spec <- model$dim_specs[[d]]
        breaks <- spec$breaks
        num_bins <- spec$num_bins

        from_vals <- vapply(cell_merge, function(e) as.numeric(e$from[d]), numeric(1))
        to_vals <- vapply(cell_merge, function(e) as.numeric(e$to[d]), numeric(1))

        from_bins <- findInterval(from_vals, breaks, rightmost.closed = TRUE)
        to_bins <- findInterval(to_vals, breaks, rightmost.closed = TRUE)

        bad_from <- from_bins < 1 | from_bins > num_bins | is.na(from_bins)
        bad_to <- to_bins < 1 | to_bins > num_bins | is.na(to_bins)
        if (any(bad_from)) {
            stop(sprintf(
                "cell_merge: 'from' values out of range in dimension %d (entries %s)",
                d, paste(which(bad_from), collapse = ", ")
            ), call. = FALSE)
        }
        if (any(bad_to)) {
            stop(sprintf(
                "cell_merge: 'to' values out of range in dimension %d (entries %s)",
                d, paste(which(bad_to), collapse = ", ")
            ), call. = FALSE)
        }

        # Remap through sample-time bin_map so cell coordinates reference
        # post-bin_merge training cells (the cells whose CDFs actually exist).
        source_dims[, d] <- sample_bin_maps[[d]][from_bins]
        target_dims[, d] <- sample_bin_maps[[d]][to_bins]
    }

    source_flat <- .compute_flat_indices(source_dims, dim_sizes)
    target_flat <- .compute_flat_indices(target_dims, dim_sizes)

    resolved <- data.frame(
        source_flat = as.integer(source_flat),
        target_flat = as.integer(target_flat),
        stringsAsFactors = FALSE
    )
    for (d in seq_len(n_dims)) {
        resolved[[paste0("from_", d)]] <- as.integer(source_dims[, d])
        resolved[[paste0("to_", d)]] <- as.integer(target_dims[, d])
    }
    attr(resolved, "n_dims") <- n_dims
    resolved
}

# --- cell_merge helpers ------------------------------------------------------

# Normalize user input into a list-of-lists with $from and $to numeric vectors
# of length n_dims.
.cell_merge_normalize <- function(cell_merge, n_dims) {
    if (is.null(cell_merge) || (is.list(cell_merge) && length(cell_merge) == 0)) {
        return(list())
    }

    if (is.data.frame(cell_merge)) {
        from_cols <- paste0("from_", seq_len(n_dims))
        to_cols <- paste0("to_", seq_len(n_dims))
        missing_cols <- setdiff(c(from_cols, to_cols), colnames(cell_merge))
        if (length(missing_cols) > 0) {
            stop(sprintf(
                "cell_merge data frame missing required columns: %s",
                paste(missing_cols, collapse = ", ")
            ), call. = FALSE)
        }
        cell_merge <- lapply(seq_len(nrow(cell_merge)), function(i) {
            list(
                from = as.numeric(unlist(cell_merge[i, from_cols])),
                to   = as.numeric(unlist(cell_merge[i, to_cols]))
            )
        })
    }

    if (!is.list(cell_merge)) {
        stop("cell_merge must be a list of entries or a data frame", call. = FALSE)
    }

    for (i in seq_along(cell_merge)) {
        entry <- cell_merge[[i]]
        if (!is.list(entry) || !all(c("from", "to") %in% names(entry))) {
            stop(sprintf(
                "cell_merge entry %d must be a list with 'from' and 'to' elements",
                i
            ), call. = FALSE)
        }
        if (length(entry$from) != n_dims || length(entry$to) != n_dims) {
            stop(sprintf(
                "cell_merge entry %d: 'from' and 'to' must each have length %d (n_dims)",
                i, n_dims
            ), call. = FALSE)
        }
    }
    cell_merge
}

# Empty resolved result with the right columns (for zero-entry input).
.cell_merge_empty_result <- function(n_dims) {
    resolved <- data.frame(
        source_flat = integer(0),
        target_flat = integer(0),
        stringsAsFactors = FALSE
    )
    for (d in seq_len(n_dims)) {
        resolved[[paste0("from_", d)]] <- integer(0)
        resolved[[paste0("to_", d)]] <- integer(0)
    }
    attr(resolved, "n_dims") <- n_dims
    resolved
}

# Compute sample-time per-dim bin_maps from a (possibly NULL) bin_merge spec.
# Falls back to training-time bin_map for dimensions not overridden. This
# duplicates the logic inside gsynth.sample() so gsynth.cell_merge() can be
# called standalone and still see the same remapping.
.cell_merge_sample_bin_maps <- function(model, bin_merge) {
    n_dims <- model$n_dims
    out <- vector("list", n_dims)
    if (!is.null(bin_merge)) {
        if (!is.list(bin_merge) || length(bin_merge) != n_dims) {
            stop(sprintf(
                "bin_merge must be a list with %d elements (one per dimension)",
                n_dims
            ), call. = FALSE)
        }
    }
    for (d in seq_len(n_dims)) {
        spec <- model$dim_specs[[d]]
        dim_merge <- if (is.null(bin_merge)) NULL else bin_merge[[d]]
        if (is.null(dim_merge)) {
            out[[d]] <- spec$bin_map
        } else {
            bin_map_result <- gsynth.bin_map(spec$breaks, dim_merge)
            new_bin_map <- seq_len(spec$num_bins)
            for (j in seq_along(bin_map_result)) {
                src_bin <- as.integer(names(bin_map_result)[j])
                tgt_bin <- as.integer(bin_map_result[j])
                if (!is.na(src_bin) && !is.na(tgt_bin) &&
                    src_bin >= 1 && src_bin <= spec$num_bins &&
                    tgt_bin >= 1 && tgt_bin <= spec$num_bins) {
                    new_bin_map[src_bin] <- tgt_bin
                }
            }
            out[[d]] <- as.integer(new_bin_map)
        }
    }
    out
}

#' Forbid a k-mer pattern in a trained gsynth model
#'
#' Returns a new \code{gsynth.model} whose samples are guaranteed not to contain
#' \code{pattern} as a substring (subject to the seeding caveat below).
#' Analytically equivalent to rejection sampling the output, implemented by
#' zeroing every transition that would produce the pattern and renormalizing
#' per state-row.
#'
#' Useful for building CpG-null, motif-null, or repeat-class-null synthetic
#' backgrounds from a standard \code{gsynth.train()} model without retraining.
#'
#' @param model A \code{gsynth.model} from \code{\link{gsynth.train}}.
#' @param pattern Character scalar, uppercase DNA (\code{ACGT} only), with
#'        \code{nchar(pattern) <= model$k + 1}. Patterns longer than one
#'        transition cannot be forbidden locally and error.
#' @param check Logical. If \code{TRUE} (default), print a short summary of how
#'        many transitions and how many bins were affected.
#'
#' @details
#' \strong{Seeding caveat.} \code{\link{gsynth.sample}} initializes the first
#' \code{k} bases of each sampling interval by uniform random draw, so those
#' seed bases may themselves contain \code{pattern}. If the seed lands on a
#' state k-mer that already contains \code{pattern} as a substring, every
#' possible next base would extend that occurrence and thus be forbidden; such
#' "trapped" states fall back to uniform sampling (not the forbid'd CDF) until
#' the pattern slides out of the state window. The guarantee applies to the
#' Markov-sampled bases downstream of the trap-escape window, not to the first
#' few bases of the interval. Expected residual per interval is small but
#' nonzero; for strict pattern-free output, pass \code{mask_copy} to
#' \code{\link{gsynth.sample}} to seed from a known pattern-free reference, or
#' scrub residuals after sampling.
#'
#' @return A new \code{gsynth.model} with modified \code{model_data$counts} and
#'         \code{model_data$cdf}. The original model is not mutated.
#'
#' @examples
#' \dontrun{
#' # CpG-null synthetic background: train on the genome, then forbid CG.
#' model <- gsynth.train(
#'     list(expr = "gc_vt", breaks = seq(0, 1, 0.05)),
#'     intervals = gintervals.all(),
#'     iterator = 200
#' )
#' model_no_cg <- gsynth.forbid_kmer(model, "CG")
#' seqs <- gsynth.sample(model_no_cg,
#'     output_format = "vector",
#'     intervals = some_regions, seed = 42
#' )
#'
#' # Motif-null background: forbid a 4-mer TF consensus substring.
#' model_no_ebox <- gsynth.forbid_kmer(model, "CACG")
#' }
#'
#' @seealso \code{\link{gsynth.sample}}, \code{\link{gsynth.train}}
#' @export
gsynth.forbid_kmer <- function(model, pattern, check = TRUE) {
    if (!inherits(model, "gsynth.model")) {
        stop("model must be a gsynth.model object", call. = FALSE)
    }
    if (!is.character(pattern) || length(pattern) != 1L || is.na(pattern)) {
        stop("pattern must be a single non-NA character string", call. = FALSE)
    }
    pattern <- toupper(pattern)
    if (!grepl("^[ACGT]+$", pattern)) {
        stop("pattern must be non-empty DNA over ACGT (got: '", pattern, "')",
            call. = FALSE
        )
    }
    L <- nchar(pattern)
    k <- if (is.null(model$k)) 5L else as.integer(model$k)
    if (L > k + 1L) {
        stop(sprintf(
            "pattern length %d exceeds model$k + 1 = %d; cannot forbid a pattern longer than a single transition",
            L, k + 1L
        ), call. = FALSE)
    }

    base_map <- c(A = 0L, C = 1L, G = 2L, T = 3L)
    pat_vec <- unname(base_map[strsplit(pattern, "")[[1]]])
    n_states <- 4L^k

    # flag_mat[state + 1, next + 1] = TRUE iff the (k+1)-mer formed by
    # concatenating the k-mer state with `next` contains `pattern` as any
    # L-length substring. Layout matches the counts/cdf matrices (row=state,
    # col=next).
    flag_mat <- matrix(FALSE, nrow = n_states, ncol = 4L)
    n_windows <- k + 2L - L # number of L-length windows inside a (k+1)-mer
    for (state in 0L:(n_states - 1L)) {
        state_bases <- integer(k)
        tmp <- state
        for (i in k:1L) {
            state_bases[i] <- tmp %% 4L
            tmp <- tmp %/% 4L
        }
        for (next_base in 0L:3L) {
            bases <- c(state_bases, next_base)
            hit <- FALSE
            for (start in seq_len(n_windows)) {
                if (all(bases[start:(start + L - 1L)] == pat_vec)) {
                    hit <- TRUE
                    break
                }
            }
            flag_mat[state + 1L, next_base + 1L] <- hit
        }
    }

    # Precompute per-state legal-next masks for empty-row fallback.
    legal_next <- !flag_mat

    # Any state whose k-mer already contains `pattern` as a substring has no
    # legal next base (every extension still contains the pattern). Such a
    # state is unreachable from a pattern-free Markov walk but IS reachable
    # via the uniform-random seeding step in gsynth.sample (see seeding caveat
    # above). Fall back to uniform over all 4 bases so the chain can eventually
    # slide the pattern out of the k-mer state window.
    trapped_state <- rowSums(legal_next) == 0L

    out <- model
    out_counts <- model$model_data$counts
    out_cdf <- model$model_data$cdf

    n_zeroed_cells <- 0L
    n_affected_bins <- 0L
    n_empty_rows <- 0L

    for (b in seq_along(out_counts)) {
        co <- out_counts[[b]]
        if (is.null(co) || all(is.na(co))) next
        zeroed_here <- sum(flag_mat & co != 0)
        if (zeroed_here > 0L) {
            n_affected_bins <- n_affected_bins + 1L
            n_zeroed_cells <- n_zeroed_cells + zeroed_here
        }
        co[flag_mat] <- 0

        rs <- rowSums(co)
        probs <- co
        nz <- rs > 0
        probs[nz, ] <- co[nz, ] / rs[nz]
        zr <- which(!nz)
        if (length(zr) > 0L) {
            n_empty_rows <- n_empty_rows + length(zr)
            for (r in zr) {
                if (trapped_state[r]) {
                    # State already contains pattern (seeding residual). Use
                    # uniform to let the chain escape the pattern window.
                    probs[r, ] <- 0.25
                } else {
                    legal <- legal_next[r, ]
                    probs[r, legal] <- 1 / sum(legal)
                    probs[r, !legal] <- 0
                }
            }
        }

        new_cdf <- t(apply(probs, 1L, cumsum))
        new_cdf[, 4L] <- 1

        out_counts[[b]] <- co
        out_cdf[[b]] <- new_cdf
    }

    out$model_data$counts <- out_counts
    out$model_data$cdf <- out_cdf

    if (isTRUE(check)) {
        message(sprintf(
            "gsynth.forbid_kmer('%s'): zeroed %s transitions across %d bins; %d empty rows used uniform fallback.",
            pattern, format(n_zeroed_cells, big.mark = ","),
            n_affected_bins, n_empty_rows
        ))
    }
    out
}

#' Train a stratified Markov model from genome sequences
#'
#' Computes a Markov model of order \code{k} (default 5) optionally stratified by
#' bins of one or more track expressions (e.g., GC content and CG dinucleotide
#' frequency). This model can be used to generate synthetic genomes that preserve
#' the k-mer statistics of the original genome within each stratification bin.
#' When called with no dimension specifications, trains a single unstratified model.
#'
#' @param ... Zero or more dimension specifications. Each specification is a list
#'        containing:
#'   \describe{
#'     \item{expr}{Track expression for this dimension (required)}
#'     \item{breaks}{Numeric vector of bin boundaries for this dimension (required)}
#'     \item{bin_merge}{Optional list of merge specifications for merging sparse bins.
#'           Each specification is a named list with 'from' and 'to' elements.}
#'   }
#'   If no dimensions are provided, trains an unstratified model with a single bin.
#' @param mask Optional intervals to exclude from training. Regions in the mask
#'        will not contribute to k-mer counts. Can be computed using \code{gscreen()}.
#' @param intervals Genomic intervals to process. If NULL, uses all chromosomes.
#' @param iterator Iterator for track evaluation, determines the resolution at which
#'        track values are computed.
#' @param pseudocount Pseudocount added to all k-mer counts to avoid zero probabilities.
#'        Default is 1.
#' @param min_obs Minimum number of observations ((k+1)-mers) required per bin. Bins
#'        with fewer observations will be marked as NA (not learned) and a warning will
#'        be issued. Default is 0 (no minimum). During sampling, NA bins will fall back
#'        to uniform sampling unless merged via \code{bin_merge}.
#' @param k Integer Markov order (1--10). Default is 5, which models 6-mer
#'        (context of length 5 plus the emitted base) transition probabilities.
#'        Higher values capture longer-range sequence dependencies but require
#'        exponentially more memory (\eqn{4^k} context states).
#'
#' @details
#' \strong{Strand symmetry:} The training process counts both the forward strand
#' (k+1)-mer and its reverse complement for each position, ensuring strand-symmetric
#' transition probabilities. This means the reported total_kmers is approximately
#' double the number of genomic positions processed.
#'
#' \strong{N bases:} Positions where the (k+1)-mer contains any N (unknown) bases
#' are skipped during training and counted in \code{total_n}. The model only learns
#' from valid A/C/G/T sequences.
#'
#' @return A \code{gsynth.model} object containing:
#'   \describe{
#'     \item{k}{Markov order used for training}
#'     \item{num_kmers}{Number of context states (\eqn{4^k})}
#'     \item{n_dims}{Number of stratification dimensions}
#'     \item{dim_specs}{List of dimension specifications (expr, breaks, num_bins, bin_map)}
#'     \item{dim_sizes}{Vector of bin counts per dimension}
#'     \item{total_bins}{Total number of bins (product of dim_sizes)}
#'     \item{total_kmers}{Total number of valid (k+1)-mers counted}
#'     \item{per_bin_kmers}{Number of (k+1)-mers counted per bin}
#'     \item{total_masked}{Number of positions skipped due to mask}
#'     \item{total_n}{Number of positions skipped due to N bases}
#'     \item{model_data}{Internal model data (counts and CDFs)}
#'   }
#'
#' @examples
#' gdb.init_examples()
#'
#' # Create virtual tracks for stratification
#' gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
#' gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
#' gvtrack.create("cg_frac", NULL, "kmer.frac", kmer = "CG")
#' gvtrack.create("masked_frac", NULL, "masked.frac")
#'
#' # Define repeat mask
#' repeats <- gscreen("masked_frac > 0.5",
#'     intervals = gintervals.all(),
#'     iterator = 100
#' )
#'
#' # Train unstratified model (no stratification)
#' model_0d <- gsynth.train(
#'     mask = repeats,
#'     intervals = gintervals.all(),
#'     iterator = 200
#' )
#'
#' # Train model with 2D stratification (GC content and CG dinucleotide)
#' model <- gsynth.train(
#'     list(
#'         expr = "g_frac + c_frac",
#'         breaks = seq(0, 1, 0.025),
#'         bin_merge = list(list(from = 0.7, to = c(0.675, 0.7)))
#'     ),
#'     list(
#'         expr = "cg_frac",
#'         breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.2),
#'         bin_merge = list(list(from = 0.04, to = c(0.03, 0.04)))
#'     ),
#'     mask = repeats,
#'     intervals = gintervals.all(),
#'     iterator = 200
#' )
#'
#' @seealso \code{\link{gsynth.sample}}, \code{\link{gsynth.save}},
#'          \code{\link{gsynth.load}}, \code{\link{gsynth.bin_map}}
#' @export
gsynth.train <- function(...,
                         mask = NULL,
                         intervals = NULL,
                         iterator = NULL,
                         pseudocount = 1,
                         min_obs = 0,
                         k = 5L) {
    .gcheckroot()

    # Validate k (Markov order)
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k != as.integer(k) || k < 1L || k > 10L) {
        stop("k must be a single integer between 1 and 10", call. = FALSE)
    }
    k <- as.integer(k)
    num_kmers <- as.integer(4L^k)

    # Capture all dimension specs from ...
    args <- list(...)

    # Allow zero dimensions for unstratified model
    zero_dim_model <- (length(args) == 0)

    # Validate and process each dimension spec
    dim_specs <- list()
    for (i in seq_along(args)) {
        spec <- args[[i]]

        if (!is.list(spec)) {
            stop(sprintf("Dimension %d must be a list", i), call. = FALSE)
        }

        if (!("expr" %in% names(spec))) {
            stop(sprintf("Dimension %d must have an 'expr' element", i), call. = FALSE)
        }

        if (!("breaks" %in% names(spec))) {
            stop(sprintf("Dimension %d must have a 'breaks' element", i), call. = FALSE)
        }

        # Validate breaks
        breaks <- sort(spec$breaks)
        if (!is.numeric(breaks) || length(breaks) < 2) {
            stop(sprintf("Dimension %d breaks must be numeric with at least 2 elements", i),
                call. = FALSE
            )
        }

        num_bins <- length(breaks) - 1

        # Process bin_merge -> bin_map for this dimension
        bin_map <- seq_len(num_bins) # Identity mapping
        if (!is.null(spec$bin_merge)) {
            bin_map_result <- gsynth.bin_map(breaks, spec$bin_merge)
            # Convert from named vector format to integer vector
            for (j in seq_along(bin_map_result)) {
                src_bin <- as.integer(names(bin_map_result)[j])
                tgt_bin <- as.integer(bin_map_result[j])
                if (!is.na(src_bin) && !is.na(tgt_bin) &&
                    src_bin >= 1 && src_bin <= num_bins &&
                    tgt_bin >= 1 && tgt_bin <= num_bins) {
                    bin_map[src_bin] <- tgt_bin
                }
            }
        }

        dim_specs[[i]] <- list(
            expr = spec$expr,
            breaks = breaks,
            num_bins = num_bins,
            bin_map = as.integer(bin_map) # 1-based for R
        )
    }

    # Set up iterator (needed for both 0D and multi-D models)
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    # Get intervals (default to all chromosomes)
    if (is.null(intervals)) {
        intervals <- gintervals.all()
    }

    # Branch based on dimensionality
    if (zero_dim_model) {
        # Zero-dimensional model: single bin, no stratification
        n_dims <- 0L
        dim_sizes <- integer(0)
        total_bins <- 1L

        # Get number of positions for creating flat indices
        # Use a minimal gextract to get chromosome/position info
        message("Setting up iterator positions...")
        iter_info <- gextract("1", intervals = intervals, iterator = .iterator)

        if (is.null(iter_info) || nrow(iter_info) == 0) {
            stop("No positions extracted. Check that intervals and iterator are valid.", call. = FALSE)
        }

        n_positions <- nrow(iter_info)

        # All positions map to bin 1 in R (will be converted to 0-based for C++)
        flat_indices <- rep(1L, n_positions)

        # Get chromosome info for processing. Use levels() (full chromkey order
        # from gintervals_chrom_sizes) rather than chrom_sizes$chrom (only chroms
        # present in the input subset) so chromids match misha's internal chromkey
        # even when the input is missing chroms that sort earlier in the chromkey.
        chrom_sizes <- gintervals.chrom_sizes(intervals)
        chrom_key <- levels(chrom_sizes$chrom)
        chrom_ids <- match(as.character(intervals$chrom), chrom_key) - 1L
        chrom_starts <- intervals$start
        chrom_ends <- intervals$end

        # Prepare iterator position data
        iter_chroms <- match(as.character(iter_info$chrom), chrom_key) - 1L
        iter_starts <- as.integer(iter_info$start)
    } else {
        # Multi-dimensional model: existing logic
        n_dims <- length(dim_specs)
        dim_sizes <- sapply(dim_specs, function(d) d$num_bins)
        total_bins <- prod(dim_sizes)

        # Extract track values for ALL dimensions
        message("Extracting track values...")
        exprs <- sapply(dim_specs, function(d) d$expr)

        # For multiple expressions, we need to extract each one
        # gextract can handle multiple expressions
        track_data <- gextract(exprs, intervals = intervals, iterator = .iterator)

        if (is.null(track_data) || nrow(track_data) == 0) {
            stop("No track data extracted. Check that intervals and iterator are valid.", call. = FALSE)
        }

        n_positions <- nrow(track_data)

        # Compute per-dimension bin indices (1-based)
        per_dim_indices <- matrix(NA_integer_, nrow = n_positions, ncol = n_dims)

        for (d in seq_len(n_dims)) {
            spec <- dim_specs[[d]]
            track_values <- track_data[[spec$expr]]
            bin_idx <- findInterval(track_values, spec$breaks, rightmost.closed = TRUE)
            bin_idx[bin_idx == 0] <- NA # Values below first break
            bin_idx[bin_idx > spec$num_bins] <- NA # Values above last break

            # Apply bin mapping for this dimension
            valid_mask <- !is.na(bin_idx) & bin_idx >= 1 & bin_idx <= spec$num_bins
            if (any(valid_mask)) {
                bin_idx[valid_mask] <- spec$bin_map[bin_idx[valid_mask]]
            }

            per_dim_indices[, d] <- bin_idx
        }

        # Compute flat bin indices using vectorized helper
        flat_indices <- .compute_flat_indices(per_dim_indices, dim_sizes)

        # Get chromosome info for processing. See note above re: levels() vs
        # chrom_sizes$chrom.
        chrom_sizes <- gintervals.chrom_sizes(intervals)
        chrom_key <- levels(chrom_sizes$chrom)
        chrom_ids <- match(as.character(intervals$chrom), chrom_key) - 1L
        chrom_starts <- intervals$start
        chrom_ends <- intervals$end

        # Prepare iterator position data
        iter_chroms <- match(as.character(track_data$chrom), chrom_key) - 1L
        iter_starts <- as.integer(track_data$start)
    }

    # Handle NA flat indices (convert to -1 for C++)
    flat_indices[is.na(flat_indices)] <- 0L
    flat_indices <- as.integer(flat_indices) - 1L # 0-based for C++

    # Prepare bin_map for C++ (0-based, for applying during training)
    # We've already applied bin_map in R, so just pass identity mapping
    bin_map_vec <- seq_len(total_bins) - 1L # 0-based identity

    message("Training Markov model...")

    # Call C++ training function
    # We pass flat indices directly and create dummy breaks that give total_bins bins
    # The C++ code computes num_bins = length(breaks) - 1
    # So we need breaks with total_bins + 1 elements
    dummy_breaks <- seq(0, total_bins, length.out = total_bins + 1)

    result <- .gcall(
        "C_gsynth_train",
        as.integer(chrom_ids),
        as.integer(chrom_starts),
        as.integer(chrom_ends),
        flat_indices,
        iter_starts,
        iter_chroms,
        as.numeric(dummy_breaks),
        bin_map_vec,
        mask,
        as.numeric(pseudocount),
        as.integer(k),
        .misha_env()
    )

    # Override C++ result with our multi-dimensional metadata
    result$k <- k
    result$num_kmers <- num_kmers
    result$n_dims <- n_dims
    result$dim_specs <- dim_specs
    result$dim_sizes <- dim_sizes
    result$total_bins <- total_bins
    result$num_bins <- total_bins # For compatibility
    result$iterator <- .iterator

    # Store all breaks for reference
    result$breaks <- NULL # Remove single breaks

    # Check for sparse bins (fewer than min_obs observations)
    result$min_obs <- min_obs
    sparse_bins <- which(result$per_bin_kmers < min_obs)
    result$sparse_bins <- sparse_bins

    if (length(sparse_bins) > 0) {
        # Mark CDFs for sparse bins as NA
        for (bin_idx in sparse_bins) {
            result$model_data$cdf[[bin_idx]][] <- NA_real_
        }

        # Notify user about sparse bins
        n_sparse <- length(sparse_bins)
        warning(sprintf(
            "%d out of %d bins have fewer than %d observations and are marked as NA.\n%s",
            n_sparse,
            total_bins,
            min_obs,
            "Use bin_merge to merge sparse bins before sampling, or they will use uniform sampling."
        ), call. = FALSE)

        # If we have dimension info, try to identify which dimensions have issues
        if (n_dims > 1 && n_sparse <= 50) {
            # Convert flat indices to per-dimension indices for user reference
            sparse_coords <- matrix(NA_integer_, nrow = n_sparse, ncol = n_dims)
            for (i in seq_along(sparse_bins)) {
                flat_idx <- sparse_bins[i] - 1L # 0-based
                for (d in seq_len(n_dims)) {
                    sparse_coords[i, d] <- (flat_idx %% dim_sizes[d]) + 1L
                    flat_idx <- flat_idx %/% dim_sizes[d]
                }
            }
            colnames(sparse_coords) <- sapply(dim_specs, function(s) s$expr)
            result$sparse_coords <- sparse_coords
        }
    }

    class(result) <- "gsynth.model"

    if (n_dims == 0) {
        message(sprintf(
            "Trained unstratified Markov-%d model: %s %d-mers (no stratification)",
            k, format(result$total_kmers, big.mark = ","), k + 1L
        ))
    } else {
        message(sprintf(
            "Trained Markov-%d model: %s %d-mers across %d bins (%d dimensions)",
            k, format(result$total_kmers, big.mark = ","), k + 1L,
            result$total_bins,
            result$n_dims
        ))
    }

    if (length(sparse_bins) > 0) {
        message(sprintf(
            "  Note: %d bins have < %d observations (marked as NA)",
            length(sparse_bins), min_obs
        ))
    }

    result
}

#' Print summary of a gsynth.model
#'
#' @param x A gsynth.model object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.gsynth.model <- function(x, ...) {
    k <- if (is.null(x$k)) 5L else x$k
    cat(sprintf("Synthetic Genome Markov-%d Model\n", k))
    cat("----------------------------\n")
    cat(sprintf("Markov order (k): %d  [context states: %d]\n", k, 4L^k))

    if (!is.null(x$n_dims) && x$n_dims == 0) {
        # Zero-dimensional (unstratified) model
        cat("Stratification: None (single global model)\n")
        cat(sprintf("Total bins: 1\n"))
    } else if (!is.null(x$n_dims)) {
        # Multi-dimensional model
        cat(sprintf("Dimensions: %d\n", x$n_dims))

        for (d in seq_len(x$n_dims)) {
            spec <- x$dim_specs[[d]]
            cat(sprintf("\n  Dimension %d:\n", d))
            cat(sprintf("    Expression: %s\n", spec$expr))
            cat(sprintf("    Bins: %d\n", spec$num_bins))
            cat(sprintf("    Range: [%.4f, %.4f]\n", min(spec$breaks), max(spec$breaks)))
        }

        cat(sprintf("\nTotal bins: %d (", x$total_bins))
        cat(paste(x$dim_sizes, collapse = " x "))
        cat(")\n")
    } else {
        # Legacy single-dimension model
        cat(sprintf("Expression: %s\n", x$expr))
        cat(sprintf("Number of bins: %d\n", x$num_bins))
    }

    cat(sprintf("Total k-mers: %s\n", format(x$total_kmers, big.mark = ",")))
    cat(sprintf("Masked positions: %s\n", format(x$total_masked, big.mark = ",")))
    cat(sprintf("N positions: %s\n", format(x$total_n, big.mark = ",")))

    # Show per-bin k-mer counts (abbreviated for multi-dimensional)
    if (x$total_bins <= 20) {
        cat("\nPer-bin k-mer counts:\n")
        for (i in seq_along(x$per_bin_kmers)) {
            if (x$per_bin_kmers[i] > 0) {
                sparse_marker <- if (i %in% x$sparse_bins) " (NA - sparse)" else ""
                cat(sprintf(
                    "  Bin %d: %s%s\n",
                    i,
                    format(x$per_bin_kmers[i], big.mark = ","),
                    sparse_marker
                ))
            }
        }
    } else {
        non_empty <- sum(x$per_bin_kmers > 0)
        cat(sprintf("\nNon-empty bins: %d / %d\n", non_empty, x$total_bins))
    }

    # Show sparse bins info
    if (!is.null(x$sparse_bins) && length(x$sparse_bins) > 0) {
        cat(sprintf(
            "\nSparse bins (< %d obs): %d bins marked as NA\n",
            x$min_obs,
            length(x$sparse_bins)
        ))
        if (length(x$sparse_bins) <= 10) {
            cat(sprintf("  Bin indices: %s\n", paste(x$sparse_bins, collapse = ", ")))
        }
    }

    invisible(x)
}

#' Save a gsynth.model to disk in .gsm format
#'
#' Saves a trained Markov model in the cross-platform .gsm format, which
#' consists of a metadata YAML file and raw binary arrays for counts and CDFs.
#' The .gsm format can be stored as a directory (default) or a ZIP archive.
#'
#' @param model A gsynth.model object from \code{\link{gsynth.train}}
#' @param file Path to save the model (directory or .zip file)
#' @param compress Logical. If \code{TRUE}, save as a ZIP archive. If
#'        \code{FALSE} (default), save as a directory.
#'
#' @return Invisibly returns the file path.
#'
#' @seealso \code{\link{gsynth.load}}, \code{\link{gsynth.train}},
#'          \code{\link{gsynth.convert}}
#' @export
gsynth.save <- function(model, file, compress = FALSE) {
    if (!inherits(model, "gsynth.model")) {
        stop("model must be a gsynth.model object", call. = FALSE)
    }

    total_bins <- model$total_bins
    n_dims <- if (is.null(model$n_dims)) 0L else model$n_dims
    k <- if (is.null(model$k)) 5L else model$k
    num_kmers <- if (is.null(model$num_kmers)) as.integer(4L^k) else model$num_kmers

    # Build dim_specs for YAML
    yaml_dim_specs <- list()
    if (n_dims > 0) {
        for (d in seq_len(n_dims)) {
            spec <- model$dim_specs[[d]]
            bin_map_val <- spec$bin_map
            # Check if bin_map is identity (1:num_bins) — store as null
            if (!is.null(bin_map_val) && length(bin_map_val) == spec$num_bins &&
                all(bin_map_val == seq_len(spec$num_bins))) {
                bin_map_val <- NULL
            }
            yaml_dim_specs[[d]] <- list(
                expr = spec$expr,
                breaks = as.numeric(spec$breaks),
                num_bins = as.integer(spec$num_bins),
                bin_map = if (is.null(bin_map_val)) NULL else as.integer(bin_map_val)
            )
        }
    }

    # Build metadata
    metadata <- list(
        format = "gsynth_model",
        version = if (is.null(model$k) || model$k == 5L) 1L else 2L,
        markov_order = as.integer(if (is.null(model$k)) 5L else model$k),
        n_dims = as.integer(n_dims),
        dim_sizes = if (n_dims > 0) as.integer(model$dim_sizes) else list(),
        total_bins = as.integer(total_bins),
        pseudocount = if (!is.null(model$pseudocount)) as.numeric(model$pseudocount) else 1.0,
        min_obs = as.integer(if (!is.null(model$min_obs)) model$min_obs else 0L),
        total_kmers = as.numeric(model$total_kmers),
        total_masked = as.numeric(model$total_masked),
        total_n = as.numeric(model$total_n),
        per_bin_kmers = as.numeric(model$per_bin_kmers),
        dim_specs = if (n_dims > 0) yaml_dim_specs else list(),
        data = list(
            counts = list(
                dtype = "float64",
                shape = list(as.integer(total_bins), as.integer(num_kmers), 4L),
                order = "C",
                file = "counts.bin"
            ),
            cdf = list(
                dtype = "float64",
                shape = list(as.integer(total_bins), as.integer(num_kmers), 4L),
                order = "C",
                file = "cdf.bin"
            )
        )
    )

    # Helper to write binary data for a list of num_kmers x 4 matrices in row-major order
    .write_bin_matrices <- function(mat_list, filepath) {
        con <- file(filepath, "wb")
        on.exit(close(con))
        for (i in seq_along(mat_list)) {
            mat <- mat_list[[i]]
            # mat is num_kmers x 4 column-major in R
            # Writing t(mat) column-major = writing mat row-major (C order)
            writeBin(as.double(t(mat)), con, size = 8L, endian = "little")
        }
    }

    if (compress) {
        # Write to temp dir, then zip
        tmp_dir <- tempfile("gsm_")
        dir.create(tmp_dir, recursive = TRUE)
        on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

        yaml::write_yaml(metadata, file.path(tmp_dir, "metadata.yaml"), precision = 15)
        .write_bin_matrices(model$model_data$counts, file.path(tmp_dir, "counts.bin"))
        .write_bin_matrices(model$model_data$cdf, file.path(tmp_dir, "cdf.bin"))

        # Create zip — use setwd to avoid directory hierarchy in archive
        zip_path <- normalizePath(file, mustWork = FALSE)
        old_wd <- setwd(tmp_dir)
        on.exit(setwd(old_wd), add = TRUE)
        utils::zip(zip_path,
            files = c("metadata.yaml", "counts.bin", "cdf.bin"),
            flags = "-q"
        )
    } else {
        # Write directly to directory
        if (file.exists(file) && !dir.exists(file)) {
            stop(sprintf("'%s' exists and is not a directory", file), call. = FALSE)
        }
        dir.create(file, recursive = TRUE, showWarnings = FALSE)

        yaml::write_yaml(metadata, file.path(file, "metadata.yaml"), precision = 15)
        .write_bin_matrices(model$model_data$counts, file.path(file, "counts.bin"))
        .write_bin_matrices(model$model_data$cdf, file.path(file, "cdf.bin"))
    }

    invisible(file)
}

#' Load a gsynth.model from disk
#'
#' Loads a previously saved Markov model. Auto-detects the format:
#' \itemize{
#'   \item If \code{file} is a directory, reads the .gsm directory format
#'   \item If \code{file} is a file, tries ZIP .gsm format first, then
#'         falls back to legacy RDS format
#' }
#'
#' @param file Path to the saved model (directory, .gsm zip, or legacy .rds)
#'
#' @return A gsynth.model object
#'
#' @seealso \code{\link{gsynth.save}}, \code{\link{gsynth.train}},
#'          \code{\link{gsynth.convert}}
#' @export
gsynth.load <- function(file) {
    if (!file.exists(file)) {
        stop(sprintf("File not found: %s", file), call. = FALSE)
    }

    if (dir.exists(file)) {
        # Directory .gsm format
        return(.load_gsm_dir(file))
    }

    # File: try ZIP first, then legacy RDS
    result <- tryCatch(.load_gsm_zip(file), error = function(e) NULL)
    if (!is.null(result)) {
        return(result)
    }

    # Legacy RDS fallback
    model <- readRDS(file)
    if (!inherits(model, "gsynth.model")) {
        stop("File does not contain a valid gsynth.model", call. = FALSE)
    }
    model
}

#' Convert a legacy RDS gsynth model to .gsm format
#'
#' Reads a gsynth.model from a legacy RDS file and saves it in the
#' cross-platform .gsm format.
#'
#' @param input_file Path to the legacy RDS model file
#' @param output_file Path for the output .gsm model (directory or zip)
#' @param compress Logical. If \code{TRUE}, save as a ZIP archive. If
#'        \code{FALSE} (default), save as a directory.
#'
#' @return Invisibly returns the output file path.
#'
#' @seealso \code{\link{gsynth.save}}, \code{\link{gsynth.load}}
#' @export
gsynth.convert <- function(input_file, output_file, compress = FALSE) {
    model <- readRDS(input_file)
    if (!inherits(model, "gsynth.model")) {
        stop("Input file does not contain a valid gsynth.model", call. = FALSE)
    }
    gsynth.save(model, output_file, compress = compress)
}

# Internal: load .gsm from a directory
.load_gsm_dir <- function(dir_path) {
    meta_path <- file.path(dir_path, "metadata.yaml")
    if (!file.exists(meta_path)) {
        stop(sprintf("metadata.yaml not found in '%s'", dir_path), call. = FALSE)
    }
    metadata <- yaml::read_yaml(meta_path)
    .build_model_from_gsm(metadata, dir_path)
}

# Internal: load .gsm from a ZIP file
.load_gsm_zip <- function(zip_path) {
    # Quick check: ZIP files start with PK (0x50 0x4B)
    sig <- readBin(zip_path, "raw", n = 2L)
    if (length(sig) < 2L || sig[1] != as.raw(0x50) || sig[2] != as.raw(0x4B)) {
        stop("Not a ZIP file", call. = FALSE)
    }

    tmp_dir <- tempfile("gsm_unzip_")
    dir.create(tmp_dir, recursive = TRUE)
    on.exit(unlink(tmp_dir, recursive = TRUE))

    utils::unzip(zip_path, exdir = tmp_dir)

    meta_path <- file.path(tmp_dir, "metadata.yaml")
    if (!file.exists(meta_path)) {
        stop("Not a valid .gsm ZIP file (no metadata.yaml)", call. = FALSE)
    }
    metadata <- yaml::read_yaml(meta_path)
    .build_model_from_gsm(metadata, tmp_dir)
}

# Internal: reconstruct gsynth.model from metadata + binary files
.build_model_from_gsm <- function(metadata, dir_path) {
    if (is.null(metadata$format) || metadata$format != "gsynth_model") {
        stop("Invalid .gsm metadata: missing or wrong format field", call. = FALSE)
    }

    total_bins <- as.integer(metadata$total_bins)
    n_dims <- as.integer(metadata$n_dims)

    # Read Markov order (default to 5 for backward compatibility)
    k <- as.integer(if (is.null(metadata$markov_order)) 5L else metadata$markov_order)
    num_kmers <- as.integer(4L^k)

    # Read binary data
    counts_path <- file.path(dir_path, "counts.bin")
    cdf_path <- file.path(dir_path, "cdf.bin")

    if (!file.exists(counts_path) || !file.exists(cdf_path)) {
        stop("Binary data files (counts.bin, cdf.bin) not found", call. = FALSE)
    }

    expected_n <- as.double(total_bins) * as.double(num_kmers) * 4

    counts_flat <- readBin(counts_path, "double", n = expected_n, size = 8L, endian = "little")
    cdf_flat <- readBin(cdf_path, "double", n = expected_n, size = 8L, endian = "little")

    if (length(counts_flat) != expected_n || length(cdf_flat) != expected_n) {
        stop(sprintf(
            "Binary data size mismatch: expected %d doubles, got counts=%d, cdf=%d",
            expected_n, length(counts_flat), length(cdf_flat)
        ), call. = FALSE)
    }

    # Reshape: row-major flat -> list of num_kmers x 4 column-major R matrices
    chunk <- num_kmers * 4L
    counts_list <- vector("list", total_bins)
    cdf_list <- vector("list", total_bins)

    for (i in seq_len(total_bins)) {
        offset <- (i - 1L) * chunk
        # Data was written as t(mat) in column-major = row-major
        # Read as 4 x num_kmers matrix (column-major), then transpose to num_kmers x 4
        counts_list[[i]] <- t(matrix(counts_flat[(offset + 1L):(offset + chunk)],
            nrow = 4L, ncol = num_kmers
        ))
        cdf_list[[i]] <- t(matrix(cdf_flat[(offset + 1L):(offset + chunk)],
            nrow = 4L, ncol = num_kmers
        ))
    }

    model_data <- list(counts = counts_list, cdf = cdf_list)

    # Reconstruct dim_specs
    dim_specs <- list()
    if (n_dims > 0 && length(metadata$dim_specs) > 0) {
        for (d in seq_len(n_dims)) {
            raw_spec <- metadata$dim_specs[[d]]
            bm <- raw_spec$bin_map
            if (is.null(bm)) {
                bm <- seq_len(as.integer(raw_spec$num_bins))
            } else {
                bm <- as.integer(unlist(bm))
            }
            dim_specs[[d]] <- list(
                expr = raw_spec$expr,
                breaks = as.numeric(unlist(raw_spec$breaks)),
                num_bins = as.integer(raw_spec$num_bins),
                bin_map = bm
            )
        }
    }

    dim_sizes <- if (n_dims > 0) as.integer(unlist(metadata$dim_sizes)) else integer(0)

    # Reconstruct sparse_bins from CDF data
    min_obs <- as.integer(if (!is.null(metadata$min_obs)) metadata$min_obs else 0L)
    per_bin_kmers <- as.numeric(unlist(metadata$per_bin_kmers))
    sparse_bins <- which(per_bin_kmers < min_obs)

    # Build the model object
    pseudocount <- if (!is.null(metadata$pseudocount)) as.numeric(metadata$pseudocount) else 1.0
    model <- list(
        k = k,
        num_kmers = num_kmers,
        num_bins = total_bins,
        breaks = NULL,
        total_kmers = as.numeric(metadata$total_kmers),
        per_bin_kmers = per_bin_kmers,
        total_masked = as.numeric(metadata$total_masked),
        total_n = as.numeric(metadata$total_n),
        model_data = model_data,
        n_dims = n_dims,
        dim_specs = dim_specs,
        dim_sizes = dim_sizes,
        total_bins = total_bins,
        pseudocount = pseudocount,
        iterator = NULL,
        min_obs = min_obs,
        sparse_bins = sparse_bins
    )

    class(model) <- "gsynth.model"
    model
}

#' Internal: build a per-bin log-probability table from a gsynth.model.
#'
#' Differences each bin's cdf to recover P(base|context), then takes log.
#' Sparse bins (cdf rows are NA from training-time min_obs filter) yield
#' NaN log_p rows, which the C++ kernel uses as the "is sparse" signal.
#'
#' @param model A gsynth.model object.
#' @return list of length total_bins; each element is a num_kmers x 4
#'         numeric matrix of natural-log conditional probabilities.
#' @noRd
.gsynth_build_log_p <- function(model) {
    stopifnot(inherits(model, "gsynth.model"))
    cdf_list <- model$model_data$cdf
    lapply(cdf_list, function(cdf) {
        # cdf is num_kmers x 4 (cumulative across the 4 bases).
        # p[, 1] = cdf[, 1]; p[, j] = cdf[, j] - cdf[, j-1] for j>1.
        p <- cbind(
            cdf[, 1],
            cdf[, 2] - cdf[, 1],
            cdf[, 3] - cdf[, 2],
            cdf[, 4] - cdf[, 3]
        )
        log(p) # NA in cdf -> NaN here
    })
}

#' Sample a synthetic genome from a trained Markov model
#'
#' Generates a synthetic genome by sampling from a trained stratified Markov
#' model. The generated genome preserves the k-mer statistics of the original
#' genome within each stratification bin.
#'
#' @param model A gsynth.model object from \code{\link{gsynth.train}}
#' @param output_path Path to the output file (ignored when output_format = "vector")
#' @param output_format Output format:
#'   \itemize{
#'     \item "misha": .seq binary format (default)
#'     \item "fasta": FASTA text format
#'     \item "vector": Return sequences as a character vector (does not write to file)
#'   }
#' @param mask_copy Optional intervals to copy from the original genome instead of
#'        sampling. Use this to preserve specific regions (e.g., repeats, regulatory
#'        elements) exactly as they appear in the reference. Regions not in mask_copy
#'        will be sampled using the Markov model.
#'        Note: mask_copy intervals should be non-overlapping and sorted by start
#'        position within each chromosome. Overlapping intervals may result in
#'        only the first overlapping region being copied, with subsequent overlaps
#'        skipped due to cursor advancement during sequential processing.
#' @param seed Random seed for reproducibility. If NULL, uses current random state.
#' @param intervals Genomic intervals to sample. If NULL, uses all chromosomes.
#' @param n_samples Number of samples to generate per interval. Default is 1.
#'        When n_samples > 1 and output_format = "fasta", headers include "_sampleN".
#'        When output_format = "vector", returns n_samples * n_intervals sequences.
#' @param bin_merge Optional list of bin merge specifications to apply during sampling,
#'        one per dimension (length must equal model$n_dims). Each element should be:
#'        \itemize{
#'          \item A list of merge specifications (same format as in
#'                \code{\link{gsynth.train}}: each spec is \code{list(from = ..., to = ...)})
#'          \item Or NULL to use the bin mapping from training for that dimension
#'        }
#'        This allows merging sparse bins at sampling time without re-training.
#'        Example for a 2D model:
#'        \preformatted{
#' bin_merge = list(
#'   # Dimension 1: merge bins with values >= 0.8 to bin [0.7, 0.8)
#'   list(list(from = c(0.8, Inf), to = c(0.7, 0.8))),
#'   # Dimension 2: use training-time bin_map (no override)
#'   NULL
#' )
#'        }
#' @param cell_merge Optional per-\emph{joint-cell} redirect table, applied
#'        after \code{bin_merge} and after any sparse-bin uniform fallback.
#'        Each entry redirects one specific training cell to another specific
#'        training cell whose Cartesian position cannot be expressed as a
#'        per-axis merge. Format: a list where each element is
#'        \code{list(from = c(v1, v2, ...), to = c(v1, v2, ...))}, with one
#'        representative value per dimension; or a data frame with columns
#'        \code{from_1, from_2, ..., to_1, to_2, ...}. Internally resolved via
#'        \code{\link{gsynth.cell_merge}}; at sampling time each source cell's
#'        CDF is replaced with the target cell's CDF (a pointer-level copy
#'        inside \code{cdf_list} --- no matrix duplication and no change to
#'        the C++ hot path).
#'
#' @details
#' \strong{FASTA index (.fai):} When \code{output_format = "fasta"}, the
#' function also writes a samtools-compatible \code{.fai} file alongside the
#' FASTA (byte offsets are tracked during the write loop, so this is
#' effectively free). The index is suitable for any downstream tool that
#' expects a samtools-indexed reference.
#'
#' \strong{N bases during sampling:} When the sampler needs to initialize the
#' first k-mer context and encounters regions with only N bases, it falls back
#' to uniform random base selection until a valid context is established.
#' Similarly, if a bin has no learned statistics (sparse bin with NA CDF),
#' uniform sampling is used for that position.
#'
#' \strong{Sparse bins:} If the model has sparse bins (from \code{min_obs} during
#' training), a warning is issued when sampling regions that fall into these bins.
#' Consider using \code{bin_merge} to redirect sparse bins to well-populated ones.
#'
#' @return When output_format is "misha" or "fasta", returns invisible NULL and
#'         writes the synthetic genome to output_path. When output_format is "vector",
#'         returns a character vector of sequences (length = n_intervals * n_samples).
#'
#' @examples
#' gdb.init_examples()
#'
#' # Create virtual tracks
#' gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
#' gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
#' gvtrack.create("cg_frac", NULL, "kmer.frac", kmer = "CG")
#' gvtrack.create("masked_frac", NULL, "masked.frac")
#'
#' # Define repeat mask (regions to preserve from original)
#' repeats <- gscreen("masked_frac > 0.5",
#'     intervals = gintervals.all(),
#'     iterator = 100
#' )
#'
#' # Train model (excluding repeats from training)
#' model <- gsynth.train(
#'     list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.025)),
#'     list(expr = "cg_frac", breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.2)),
#'     mask = repeats,
#'     iterator = 200,
#'     min_obs = 1000
#' )
#'
#' # Sample with mask_copy to preserve repeats from original genome
#' temp_dir <- tempdir()
#' synthetic_genome_file <- file.path(temp_dir, "synthetic_genome.fa")
#' gsynth.sample(model, synthetic_genome_file,
#'     output_format = "fasta",
#'     mask_copy = repeats,
#'     seed = 60427,
#'     bin_merge = list(
#'         list(list(from = 0.7, to = c(0.675, 0.7))),
#'         list(list(from = 0.04, to = c(0.03, 0.04)))
#'     )
#' )
#'
#' @seealso \code{\link{gsynth.train}}, \code{\link{gsynth.save}}
#' @export
gsynth.sample <- function(model,
                          output_path = NULL,
                          output_format = c("misha", "fasta", "vector"),
                          mask_copy = NULL,
                          seed = NULL,
                          intervals = NULL,
                          n_samples = 1,
                          bin_merge = NULL,
                          cell_merge = NULL) {
    .gcheckroot()

    if (!inherits(model, "gsynth.model")) {
        stop("model must be a gsynth.model object", call. = FALSE)
    }

    # Extract Markov order from model (default k=5 for backward compatibility)
    k <- if (is.null(model$k)) 5L else model$k
    num_kmers <- if (is.null(model$num_kmers)) as.integer(4L^k) else model$num_kmers

    output_format <- match.arg(output_format)

    # Validate output_path requirement
    if (output_format != "vector" && is.null(output_path)) {
        stop("output_path is required when output_format is not 'vector'", call. = FALSE)
    }

    # Validate n_samples
    n_samples <- as.integer(n_samples)
    if (is.na(n_samples) || n_samples < 1) {
        stop("n_samples must be a positive integer", call. = FALSE)
    }

    # Set random seed if provided
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Get intervals (default to all chromosomes)
    if (is.null(intervals)) {
        intervals <- gintervals.all()
    }

    # Get iterator from model (with validation)
    .iterator <- model$iterator
    if (is.null(.iterator)) {
        stop("Model iterator is NULL. The model may be corrupted or from an incompatible version.", call. = FALSE)
    }

    n_dims <- model$n_dims
    dim_sizes <- model$dim_sizes

    message("Extracting track values for sampling...")

    if (n_dims == 0) {
        # Zero-dimensional model: all positions use bin 1 (R 1-based)

        # Still need iterator info to know positions
        iter_info <- gextract("1", intervals = intervals, iterator = .iterator)

        if (is.null(iter_info) || nrow(iter_info) == 0) {
            stop("No positions extracted. Check that intervals are valid.", call. = FALSE)
        }

        n_positions <- nrow(iter_info)

        # All positions map to bin 1 in R (will be converted to 0-based for C++)
        flat_indices <- rep(1L, n_positions)

        # Get chromosome info. Use levels() (full chromkey order) rather than
        # chrom_sizes$chrom (only chroms present in the input) so chromids match
        # misha's internal chromkey.
        chrom_sizes <- gintervals.chrom_sizes(intervals)
        chrom_key <- levels(chrom_sizes$chrom)
        iter_chroms <- match(as.character(iter_info$chrom), chrom_key) - 1L
        iter_starts <- as.integer(iter_info$start)
    } else {
        # Multi-dimensional model: extract track values and compute bins

        # Extract track values for all dimensions
        exprs <- sapply(model$dim_specs, function(d) d$expr)
        track_data <- gextract(exprs, intervals = intervals, iterator = .iterator)

        if (is.null(track_data) || nrow(track_data) == 0) {
            stop("No track data extracted. Check that intervals are valid.", call. = FALSE)
        }

        n_positions <- nrow(track_data)

        # Process sampling-time bin_merge if provided
        sample_bin_maps <- vector("list", n_dims)
        if (!is.null(bin_merge)) {
            if (!is.list(bin_merge) || length(bin_merge) != n_dims) {
                stop(sprintf(
                    "bin_merge must be a list with %d elements (one per dimension)",
                    n_dims
                ), call. = FALSE)
            }

            for (d in seq_len(n_dims)) {
                spec <- model$dim_specs[[d]]
                dim_merge <- bin_merge[[d]]

                if (is.null(dim_merge)) {
                    # Use training-time bin_map for this dimension
                    sample_bin_maps[[d]] <- spec$bin_map
                } else {
                    # Compute new bin_map from sampling-time bin_merge
                    bin_map_result <- gsynth.bin_map(spec$breaks, dim_merge)
                    new_bin_map <- seq_len(spec$num_bins)
                    for (j in seq_along(bin_map_result)) {
                        src_bin <- as.integer(names(bin_map_result)[j])
                        tgt_bin <- as.integer(bin_map_result[j])
                        if (!is.na(src_bin) && !is.na(tgt_bin) &&
                            src_bin >= 1 && src_bin <= spec$num_bins &&
                            tgt_bin >= 1 && tgt_bin <= spec$num_bins) {
                            new_bin_map[src_bin] <- tgt_bin
                        }
                    }
                    sample_bin_maps[[d]] <- as.integer(new_bin_map)
                }
            }
        } else {
            # Use training-time bin_maps for all dimensions
            for (d in seq_len(n_dims)) {
                sample_bin_maps[[d]] <- model$dim_specs[[d]]$bin_map
            }
        }

        # Compute per-dimension bin indices (1-based)
        per_dim_indices <- matrix(NA_integer_, nrow = n_positions, ncol = n_dims)

        for (d in seq_len(n_dims)) {
            spec <- model$dim_specs[[d]]
            track_values <- track_data[[spec$expr]]
            bin_idx <- findInterval(track_values, spec$breaks, rightmost.closed = TRUE)
            bin_idx[bin_idx == 0] <- NA
            bin_idx[bin_idx > spec$num_bins] <- NA

            # Apply bin mapping for this dimension (use sample-time bin_map)
            valid_mask <- !is.na(bin_idx) & bin_idx >= 1 & bin_idx <= spec$num_bins
            if (any(valid_mask)) {
                bin_idx[valid_mask] <- sample_bin_maps[[d]][bin_idx[valid_mask]]
            }

            per_dim_indices[, d] <- bin_idx
        }

        # Compute flat bin indices using vectorized helper
        flat_indices <- .compute_flat_indices(per_dim_indices, dim_sizes)

        # Get chromosome info. See note above re: levels() vs chrom_sizes$chrom.
        chrom_sizes <- gintervals.chrom_sizes(intervals)
        chrom_key <- levels(chrom_sizes$chrom)
        iter_chroms <- match(as.character(track_data$chrom), chrom_key) - 1L
        iter_starts <- as.integer(track_data$start)
    }

    # Handle NA flat indices (convert to -1 for C++)
    flat_indices[is.na(flat_indices)] <- 0L
    flat_indices <- as.integer(flat_indices) - 1L # 0-based

    # Get CDF data from model (make a copy to avoid modifying the original)
    cdf_list <- model$model_data$cdf

    # Check for sparse bins that will be used during sampling
    sparse_bins <- model$sparse_bins
    if (!is.null(sparse_bins) && length(sparse_bins) > 0) {
        # Find which sparse bins will be hit
        unique_flat <- unique(flat_indices[flat_indices >= 0]) + 1L # 1-based
        sparse_used <- intersect(unique_flat, sparse_bins)

        if (length(sparse_used) > 0) {
            warning(sprintf(
                "%d sparse bins (with < %d observations) will be used during sampling.\n%s",
                length(sparse_used),
                model$min_obs,
                "These bins will use uniform base distribution. Consider using bin_merge to merge them."
            ), call. = FALSE)

            # Replace NA CDFs with uniform distribution for sampling
            # Uniform CDF: [0.25, 0.5, 0.75, 1.0]
            uniform_cdf <- matrix(
                rep(c(0.25, 0.5, 0.75, 1.0), each = num_kmers),
                nrow = num_kmers, ncol = 4
            )
            for (bin_idx in sparse_used) {
                cdf_list[[bin_idx]] <- uniform_cdf
            }
        }
    }

    # Apply cell_merge: per-joint-cell CDF redirects. Runs AFTER bin_merge and
    # the sparse-bin uniform fallback, so callers can layer the mechanisms
    # (cell_merge wins). Each redirect is a pointer-level list-element copy;
    # the underlying CDF matrix is shared, no duplication, and the C++ hot
    # path sees identical SEXPs for source and target flat bins.
    if (!is.null(cell_merge) && model$n_dims >= 1) {
        resolved <- gsynth.cell_merge(model, cell_merge, bin_merge = bin_merge)

        if (nrow(resolved) > 0) {
            self <- resolved$source_flat == resolved$target_flat
            if (any(self)) {
                warning(sprintf(
                    "cell_merge: %d %s a self-redirect and %s ignored.",
                    sum(self),
                    if (sum(self) == 1L) "entry is" else "entries are",
                    if (sum(self) == 1L) "was" else "were"
                ), call. = FALSE)
                resolved <- resolved[!self, , drop = FALSE]
            }
        }

        if (nrow(resolved) > 0) {
            dup <- duplicated(resolved$source_flat)
            if (any(dup)) {
                warning(sprintf(
                    "cell_merge: %d duplicate source cell(s); later entries override earlier.",
                    sum(dup)
                ), call. = FALSE)
            }

            if (!is.null(sparse_bins) && length(sparse_bins) > 0 &&
                any(resolved$target_flat %in% sparse_bins)) {
                warning("cell_merge: at least one target cell is itself a sparse bin (uniform fallback).",
                    call. = FALSE
                )
            }

            message(sprintf(
                "cell_merge: redirecting %d source cell(s).",
                nrow(resolved)
            ))

            for (i in seq_len(nrow(resolved))) {
                cdf_list[[resolved$source_flat[i]]] <-
                    cdf_list[[resolved$target_flat[i]]]
            }
        }
    } else if (!is.null(cell_merge) && model$n_dims < 1) {
        stop("cell_merge requires a stratified model (n_dims >= 1)", call. = FALSE)
    }

    # Output format: 0 = misha, 1 = fasta, 2 = vector
    output_format_int <- switch(output_format,
        misha = 0L,
        fasta = 1L,
        vector = 2L
    )

    if (n_samples > 1 || output_format == "vector") {
        message(sprintf("Sampling synthetic genome (%d samples per interval)...", n_samples))
    } else {
        message("Sampling synthetic genome...")
    }

    # Call C++ sampling function
    # Create dummy breaks that give total_bins bins (length = total_bins + 1)
    dummy_breaks <- seq(0, model$total_bins, length.out = model$total_bins + 1)

    # Use empty string for output_path if vector mode
    output_path_str <- if (is.null(output_path)) "" else as.character(output_path)

    result <- .gcall(
        "C_gsynth_sample",
        cdf_list,
        as.numeric(dummy_breaks),
        flat_indices,
        iter_starts,
        iter_chroms,
        intervals,
        mask_copy,
        output_path_str,
        output_format_int,
        n_samples,
        as.integer(k),
        as.integer(model$iterator),
        .misha_env()
    )

    if (output_format == "vector") {
        message(sprintf("Generated %d sequence(s)", length(result)))
        return(result)
    }

    message(sprintf("Synthetic genome written to: %s", output_path))

    invisible(NULL)
}

#' Score the genome under a trained gsynth model
#'
#' Writes a misha fixed-bin dense track whose value at each output bin
#' is the summed natural-log conditional probability of the reference
#' sequence under the trained Markov model:
#' \deqn{T(a) = \sum_{i=a}^{a+r-1} \log P_M(\mathrm{seq}[i] \mid
#'             \mathrm{seq}[i-k..i-1], b(i))}
#' where \eqn{b(i)} is the model's stratum bin (constant within each
#' \code{model$iterator}-bp window). The first \eqn{k} bp of every
#' chromosome are NA (no upstream context available); under default
#' policies, any output bin containing an NA per-bp contribution is NA.
#'
#' Two models scored against the same reference give a windowed log
#' Bayes factor as a track expression:
#' \preformatted{
#'     gsynth.score(genome_model, "genome_score")
#'     gsynth.score(cre_model,    "cre_score")
#'     gextract("cre_score - genome_score", iterator=1000, ...)
#' }
#'
#' @param model       A \code{gsynth.model} object (from
#'                    \code{\link{gsynth.train}}).
#' @param track       Name of the misha track to create.
#' @param description Optional track description.
#' @param intervals   Intervals to score. Defaults to
#'                    \code{gintervals.all()}. Best results when interval
#'                    starts are aligned to multiples of
#'                    \code{model$iterator}; otherwise the first stratum
#'                    window is shorter than \code{model$iterator} and
#'                    its bin label may differ from training.
#' @param mask        Optional intervals to NA-out in the output (e.g.
#'                    repeats). Intersects with \code{intervals} per bp;
#'                    every output bin containing a masked bp becomes
#'                    \code{NaN}.
#' @param resolution  Output bin size in bp. Defaults to
#'                    \code{model$iterator}. \code{1} produces a per-bp
#'                    track; any positive integer is allowed.
#' @param sparse_policy How to score positions whose stratum bin is
#'                    sparse in the model:
#'                    \code{"NA"} (default) propagates NA;
#'                    \code{"uniform"} contributes \code{log(1/4)} per
#'                    bp.
#' @param n_policy    How to score positions whose k-mer \emph{context}
#'                    contains an N: \code{"NA"} (default) or
#'                    \code{"uniform"} (\code{log(1/4)} per bp). The
#'                    predicted base itself is always NA when N — the
#'                    model has no \eqn{\log P} for non-ACGT bases.
#' @param overwrite   If \code{TRUE}, replace an existing track of the
#'                    same name.
#'
#' @return Invisibly \code{NULL}. Side effect: creates the named misha
#'         track. Output bins are written as \code{NaN} where the sum
#'         is undefined (out of intervals, or any per-bp NA).
#'
#' @seealso \code{\link{gsynth.train}}, \code{\link{gsynth.sample}}
#' @export
gsynth.score <- function(model,
                         track,
                         description = NULL,
                         intervals = NULL,
                         mask = NULL,
                         resolution = NULL,
                         sparse_policy = c("NA", "uniform"),
                         n_policy = c("NA", "uniform"),
                         overwrite = FALSE) {
    .gcheckroot()

    if (!inherits(model, "gsynth.model")) {
        stop("model must be a gsynth.model object", call. = FALSE)
    }
    if (missing(track) || !is.character(track) || length(track) != 1L ||
        is.na(track) || nchar(track) == 0L) {
        stop("track is required (single non-empty string)", call. = FALSE)
    }
    if (is.null(resolution)) {
        resolution <- as.integer(model$iterator)
    }
    if (!is.numeric(resolution) || length(resolution) != 1L ||
        is.na(resolution) || resolution != as.integer(resolution) ||
        resolution < 1L) {
        stop("resolution must be a positive integer (>= 1)", call. = FALSE)
    }
    resolution <- as.integer(resolution)

    sparse_policy <- match.arg(sparse_policy)
    n_policy <- match.arg(n_policy)

    if (is.null(intervals)) {
        intervals <- gintervals.all()
    }

    # ---- stratum extraction (mirrors gsynth.sample) ----
    .iterator <- model$iterator
    if (is.null(.iterator)) {
        stop("Model iterator is NULL. Model may be corrupted.", call. = FALSE)
    }
    n_dims <- model$n_dims
    dim_sizes <- model$dim_sizes

    # Bin lookup queries pos - k (the leftmost base of the (k+1)-mer
    # context) to match training. Extend gextract upstream by one
    # iter_size so the first k bp of every interval get bin info from
    # the prior iter window — matching what training saw on the same
    # genome. Clamped to chromosome start.
    .k <- as.integer(model$k)
    strata_intervals <- intervals
    if (is.numeric(.iterator) && length(.iterator) == 1L && .iterator > 0L) {
        strata_intervals$start <- as.integer(pmax(
            0L, as.integer(strata_intervals$start) - as.integer(.iterator)
        ))
    }

    if (n_dims == 0) {
        iter_info <- gextract("1",
            intervals = strata_intervals, iterator = .iterator
        )
        if (is.null(iter_info) || nrow(iter_info) == 0) {
            stop("No positions extracted. Check intervals.", call. = FALSE)
        }
        n_positions <- nrow(iter_info)
        flat_indices <- rep(1L, n_positions)
        chrom_sizes <- gintervals.chrom_sizes(intervals)
        chrom_key <- levels(chrom_sizes$chrom)
        iter_chroms <- match(as.character(iter_info$chrom), chrom_key) - 1L
        iter_starts <- as.integer(iter_info$start)
    } else {
        exprs <- sapply(model$dim_specs, function(d) d$expr)
        track_data <- gextract(exprs,
            intervals = strata_intervals, iterator = .iterator
        )
        if (is.null(track_data) || nrow(track_data) == 0) {
            stop("No track data extracted. Check intervals.", call. = FALSE)
        }
        n_positions <- nrow(track_data)
        per_dim_indices <- matrix(NA_integer_, nrow = n_positions, ncol = n_dims)
        for (d in seq_len(n_dims)) {
            spec <- model$dim_specs[[d]]
            tv <- track_data[[spec$expr]]
            bi <- findInterval(tv, spec$breaks, rightmost.closed = TRUE)
            bi[bi == 0] <- NA
            bi[bi > spec$num_bins] <- NA
            valid <- !is.na(bi) & bi >= 1 & bi <= spec$num_bins
            if (any(valid)) bi[valid] <- spec$bin_map[bi[valid]]
            per_dim_indices[, d] <- bi
        }
        flat_indices <- .compute_flat_indices(per_dim_indices, dim_sizes)
        chrom_sizes <- gintervals.chrom_sizes(intervals)
        chrom_key <- levels(chrom_sizes$chrom)
        iter_chroms <- match(as.character(track_data$chrom), chrom_key) - 1L
        iter_starts <- as.integer(track_data$start)
    }

    flat_indices[is.na(flat_indices)] <- 0L
    flat_indices <- as.integer(flat_indices) - 1L # 0-based for C++

    log_p_list <- .gsynth_build_log_p(model)

    # Encode policies for C++: 0 = NA, 1 = uniform.
    n_policy_int <- if (n_policy == "uniform") 1L else 0L
    sparse_policy_int <- if (sparse_policy == "uniform") 1L else 0L

    # ---- Track-creation flow: parent creates the dir once (idempotent
    # across forked workers in parallel mode), kernel writes per-chrom
    # files into it, parent registers the track on success.
    if (gtrack.exists(track)) {
        if (!overwrite) {
            stop(sprintf(
                "Track '%s' already exists; pass overwrite=TRUE to replace.",
                track
            ), call. = FALSE)
        }
        gtrack.rm(track, force = TRUE)
    }

    track_path <- .track_dir(track)
    if (!dir.exists(track_path)) {
        dir.create(track_path, recursive = TRUE, mode = "0777")
    }

    # ---- Decide serial vs parallel chrom dispatch.
    # Mirrors misha's general convention: gmultitasking + gmax.processes.
    # Each chrom is independent (own state, own output file), so chrom-
    # parallelism scales near-linearly until n_workers ~ #chroms.
    all_chromids <- seq_along(chrom_key) - 1L # 0-based
    use_mt <- isTRUE(.ggetOption("gmultitasking", TRUE)) &&
        length(all_chromids) > 1L
    n_workers <- if (use_mt) {
        as.integer(min(.ggetOption("gmax.processes", 1L), length(all_chromids)))
    } else {
        1L
    }

    success <- FALSE
    tryCatch(
        {
            if (n_workers <= 1L) {
                # Serial path — single .gcall covering all chroms.
                .gcall(
                    "C_gsynth_score",
                    log_p_list,
                    flat_indices,
                    iter_starts,
                    iter_chroms,
                    intervals,
                    track,
                    all_chromids,
                    as.integer(resolution),
                    as.integer(model$k),
                    as.integer(.iterator),
                    n_policy_int,
                    sparse_policy_int,
                    mask,
                    .misha_env()
                )
            } else {
                # Parallel path — split chroms across workers, each worker
                # writes its assigned chrom files into the shared track dir.
                # Round-robin split keeps the largest chroms spread across
                # workers (chromkey is roughly size-sorted).
                groups <- split(all_chromids, seq_along(all_chromids) %% n_workers)
                results <- parallel::mclapply(groups, function(cids) {
                    tryCatch(
                        {
                            .gcall(
                                "C_gsynth_score",
                                log_p_list,
                                flat_indices,
                                iter_starts,
                                iter_chroms,
                                intervals,
                                track,
                                as.integer(cids),
                                as.integer(resolution),
                                as.integer(model$k),
                                as.integer(.iterator),
                                n_policy_int,
                                sparse_policy_int,
                                mask,
                                .misha_env()
                            )
                            list(ok = TRUE)
                        },
                        error = function(e) list(ok = FALSE, msg = conditionMessage(e))
                    )
                }, mc.cores = n_workers, mc.preschedule = FALSE)

                fails <- vapply(results, function(r) {
                    if (inherits(r, "try-error")) {
                        return(TRUE)
                    }
                    if (is.null(r) || !isTRUE(r$ok)) {
                        return(TRUE)
                    }
                    FALSE
                }, logical(1))
                if (any(fails)) {
                    msgs <- vapply(results[fails], function(r) {
                        if (is.list(r) && !is.null(r$msg)) r$msg else "<unknown>"
                    }, character(1))
                    stop(sprintf(
                        "gsynth.score parallel kernel failed in %d/%d worker(s): %s",
                        sum(fails), length(results),
                        paste(unique(msgs), collapse = "; ")
                    ), call. = FALSE)
                }
            }

            .gdb.add_track(track)

            .gtrack.attr.set(
                track, "created.by",
                sprintf(
                    "gsynth.score(model, '%s', resolution=%d)",
                    track, resolution
                ),
                TRUE
            )
            .gtrack.attr.set(track, "created.date", date(), TRUE)
            .gtrack.attr.set(track, "created.user", Sys.getenv("USER"), TRUE)
            .gtrack.attr.set(
                track, "description",
                if (is.null(description)) "" else as.character(description),
                TRUE
            )
            .gtrack.attr.set(track, "type", "dense", TRUE)
            .gtrack.attr.set(track, "binsize", as.integer(resolution), TRUE)

            success <- TRUE
        },
        finally = {
            if (!success) {
                # Best-effort cleanup if kernel or registration failed.
                track_path <- tryCatch(.track_dir(track),
                    error = function(e) NULL
                )
                if (!is.null(track_path) && dir.exists(track_path)) {
                    unlink(track_path, recursive = TRUE, force = TRUE)
                }
            }
        }
    )

    invisible(NULL)
}

#' Generate random genome sequences
#'
#' Generates random DNA sequences based on nucleotide probabilities without
#' using a trained Markov model. Each nucleotide is sampled independently
#' according to the specified probabilities.
#'
#' @param intervals Genomic intervals to sample. If NULL, uses all chromosomes.
#' @param output_path Path to the output file (ignored when output_format = "vector")
#' @param output_format Output format:
#'   \itemize{
#'     \item "misha": .seq binary format (default)
#'     \item "fasta": FASTA text format
#'     \item "vector": Return sequences as a character vector (does not write to file)
#'   }
#' @param nuc_probs Nucleotide probabilities. Can be specified as:
#'   \itemize{
#'     \item A named vector: \code{c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)}
#'     \item An unnamed vector in A, C, G, T order: \code{c(0.3, 0.2, 0.2, 0.3)}
#'   }
#'   Probabilities are automatically normalized to sum to 1. Default is uniform
#'   (0.25 each).
#' @param mask_copy Optional intervals to copy from the original genome instead of
#'        random sampling. Use this to preserve specific regions exactly as they
#'        appear in the reference.
#' @param seed Random seed for reproducibility. If NULL, uses current random state.
#' @param n_samples Number of samples to generate per interval. Default is 1.
#' @param iterator Iterator for position resolution. Default is 1 (base-pair resolution).
#'        Larger values may speed up processing but are typically not needed for
#'        random sampling.
#'
#' @details
#' Unlike \code{\link{gsynth.sample}} which uses a trained Markov model to generate
#' sequences that preserve k-mer statistics, \code{gsynth.random} generates purely
#' random sequences where each nucleotide is sampled independently. This is useful
#' for generating baseline random sequences or sequences with specific GC content.
#'
#' \strong{Nucleotide ordering:} When using an unnamed vector for \code{nuc_probs},
#' the order is A, C, G, T. Named vectors can be in any order.
#'
#' @return When output_format is "misha" or "fasta", returns invisible NULL and
#'         writes the random sequences to output_path. When output_format is "vector",
#'         returns a character vector of sequences (length = n_intervals * n_samples).
#'
#' @examples
#' gdb.init_examples()
#'
#' # Generate random sequences with uniform nucleotide probabilities
#' seqs <- gsynth.random(
#'     intervals = gintervals(1, 0, 1000),
#'     output_format = "vector",
#'     seed = 42
#' )
#'
#' # Generate GC-rich sequences (60% GC)
#' gc_rich <- gsynth.random(
#'     intervals = gintervals(1, 0, 1000),
#'     output_format = "vector",
#'     nuc_probs = c(A = 0.2, C = 0.3, G = 0.3, T = 0.2),
#'     seed = 42
#' )
#'
#' # Generate AT-rich sequences
#' at_rich <- gsynth.random(
#'     intervals = gintervals(1, 0, 1000),
#'     output_format = "vector",
#'     nuc_probs = c(A = 0.35, C = 0.15, G = 0.15, T = 0.35),
#'     seed = 42
#' )
#'
#' @seealso \code{\link{gsynth.sample}}, \code{\link{gsynth.train}}
#' @export
gsynth.random <- function(intervals = NULL,
                          output_path = NULL,
                          output_format = c("misha", "fasta", "vector"),
                          nuc_probs = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                          mask_copy = NULL,
                          seed = NULL,
                          n_samples = 1,
                          iterator = 1) {
    .gcheckroot()

    output_format <- match.arg(output_format)

    # Validate output_path requirement
    if (output_format != "vector" && is.null(output_path)) {
        stop("output_path is required when output_format is not 'vector'", call. = FALSE)
    }

    # Validate n_samples
    n_samples <- as.integer(n_samples)
    if (is.na(n_samples) || n_samples < 1) {
        stop("n_samples must be a positive integer", call. = FALSE)
    }

    # Validate and normalize nuc_probs
    if (!is.numeric(nuc_probs) || length(nuc_probs) != 4) {
        stop("nuc_probs must be a numeric vector of length 4 (A, C, G, T)", call. = FALSE)
    }
    if (any(nuc_probs < 0)) {
        stop("nuc_probs values must be non-negative", call. = FALSE)
    }
    if (sum(nuc_probs) == 0) {
        stop("nuc_probs must have at least one positive value", call. = FALSE)
    }

    # Handle named vector - reorder to A, C, G, T
    if (!is.null(names(nuc_probs))) {
        expected_names <- c("A", "C", "G", "T")
        if (!all(toupper(names(nuc_probs)) %in% expected_names)) {
            stop("nuc_probs names must be A, C, G, T", call. = FALSE)
        }
        # Reorder to A, C, G, T
        nuc_probs <- nuc_probs[match(expected_names, toupper(names(nuc_probs)))]
    }

    # Normalize to sum to 1
    nuc_probs <- nuc_probs / sum(nuc_probs)

    # Set random seed if provided
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Get intervals (default to all chromosomes)
    if (is.null(intervals)) {
        intervals <- gintervals.all()
    }

    # Check if we need to process in chunks for large genomes
    # Define chunk processor for parallel execution
    # Captures output_format, nuc_probs, mask_copy, n_samples, iterator from enclosing scope
    .process_random_chunk <- function(chunk_interval, chunk_index, chunk_seed = NULL, ...) {
        # Filter mask to this chunk's chromosome
        chunk_mask <- if (!is.null(mask_copy)) {
            mask_copy[mask_copy$chrom == chunk_interval$chrom[1], ]
        } else {
            NULL
        }

        if (output_format == "vector") {
            # Vector mode: return sequences directly
            gsynth.random(
                intervals = chunk_interval,
                output_path = NULL,
                output_format = "vector",
                nuc_probs = nuc_probs,
                mask_copy = chunk_mask,
                seed = chunk_seed,
                n_samples = n_samples,
                iterator = iterator
            )
        } else {
            # File mode: write to temp file and return path
            temp_file <- tempfile(fileext = if (output_format == "fasta") ".fa" else ".seq")
            gsynth.random(
                intervals = chunk_interval,
                output_path = temp_file,
                output_format = output_format,
                nuc_probs = nuc_probs,
                mask_copy = chunk_mask,
                seed = chunk_seed,
                n_samples = n_samples,
                iterator = iterator
            )
            temp_file
        }
    }

    # Try parallel processing for large genomes
    parallel_result <- .gsynth_process_parallel(
        intervals = intervals,
        output_format = output_format,
        output_path = output_path,
        process_chunk_fn = .process_random_chunk,
        max_chunk_size = .GSYNTH_MAX_CHUNK_SIZE,
        seed = seed
    )

    # If parallel processing handled it, return the result
    if (!is.null(parallel_result)) {
        message(sprintf("Generated %d random sequence(s)", length(parallel_result)))
        return(parallel_result)
    }
    if (is.null(parallel_result) && sum(intervals$end - intervals$start) > .GSYNTH_MAX_CHUNK_SIZE) {
        # File output was written by parallel processing
        message(sprintf("Random sequences written to: %s", output_path))
        return(invisible(NULL))
    }

    # Set up iterator
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    # Get position info using gextract
    message("Setting up random sampling positions...")
    iter_info <- gextract("1", intervals = intervals, iterator = .iterator)

    if (is.null(iter_info) || nrow(iter_info) == 0) {
        stop("No positions extracted. Check that intervals are valid.", call. = FALSE)
    }

    n_positions <- nrow(iter_info)

    # All positions map to bin 0 (single bin)
    flat_indices <- rep(0L, n_positions)

    # Get chromosome info. Use levels() (full chromkey order) rather than
    # chrom_sizes$chrom (only chroms present in the input) so chromids match
    # misha's internal chromkey.
    chrom_sizes <- gintervals.chrom_sizes(intervals)
    chrom_key <- levels(chrom_sizes$chrom)
    iter_chroms <- match(as.character(iter_info$chrom), chrom_key) - 1L
    iter_starts <- as.integer(iter_info$start)

    # Create CDF from probabilities: [P(A), P(A)+P(C), P(A)+P(C)+P(G), 1]
    # Order is A=0, C=1, G=2, T=3 (same as C++ code)
    cdf_row <- cumsum(as.numeric(nuc_probs))

    # Use k=5 for context-independent random sampling (all rows identical)
    random_k <- 5L
    random_num_kmers <- as.integer(4L^random_k)

    # Create CDF matrix: same for all contexts (no context dependency)
    # Each row is the same CDF - pure random, context-independent sampling
    cdf_matrix <- matrix(rep(cdf_row, random_num_kmers), nrow = random_num_kmers, ncol = 4, byrow = TRUE)
    cdf_list <- list(cdf_matrix)

    # Output format: 0 = misha, 1 = fasta, 2 = vector
    output_format_int <- switch(output_format,
        misha = 0L,
        fasta = 1L,
        vector = 2L
    )

    if (n_samples > 1 || output_format == "vector") {
        message(sprintf("Generating random sequences (%d samples per interval)...", n_samples))
    } else {
        message("Generating random sequences...")
    }

    # Use empty string for output_path if vector mode
    output_path_str <- if (is.null(output_path)) "" else as.character(output_path)

    # Create dummy breaks for single bin
    dummy_breaks <- c(0, 1)

    # Derive iter_size for C_gsynth_sample. For gsynth.random all positions
    # map to the single bin (index 0), so iter_size only controls whether a
    # position inside an interval counts as bin-covered. For numeric iterators
    # this is just the iterator itself; for track/intervals iterators the bin
    # extent isn't fixed, so pass .Machine$integer.max to mark every position
    # inside any interval as bin-covered (preserves pre-fix behavior).
    iter_size_int <- if (is.numeric(.iterator) && length(.iterator) == 1 && .iterator > 0) {
        as.integer(.iterator)
    } else {
        .Machine$integer.max
    }

    result <- .gcall(
        "C_gsynth_sample",
        cdf_list,
        as.numeric(dummy_breaks),
        flat_indices,
        iter_starts,
        iter_chroms,
        intervals,
        mask_copy,
        output_path_str,
        output_format_int,
        n_samples,
        as.integer(random_k),
        iter_size_int,
        .misha_env()
    )

    if (output_format == "vector") {
        message(sprintf("Generated %d random sequence(s)", length(result)))
        return(result)
    }

    message(sprintf("Random sequences written to: %s", output_path))

    invisible(NULL)
}

#' Iteratively replace a k-mer in the genome
#'
#' Performs an iterative replacement of a \code{target} k-mer with a
#' \code{replacement} sequence. This is useful for creating synthetic genomes
#' with specific motifs removed (e.g., creating a CpG-null genome by iteratively
#' swapping CG to GC).
#'
#' @param target The k-mer sequence to remove (e.g., "CG").
#' @param replacement The replacement sequence (e.g., "GC").
#' @param output_path Path to the output file (ignored when output_format = "vector").
#' @param output_format Output format:
#'   \itemize{
#'     \item "misha": .seq binary format (default)
#'     \item "fasta": FASTA text format
#'     \item "vector": Return sequences as a character vector (does not write to file)
#'   }
#' @param intervals Genomic intervals to process. If NULL, uses all chromosomes.
#' @param check_composition Logical. If TRUE (default), ensures target and
#'        replacement have the same nucleotide composition (preserving exact
#'        base counts).
#'
#' @details
#' \strong{Bubble Sort / Iterative Logic:} The function scans the sequence and
#' replaces occurrences of \code{target} with \code{replacement}. If a replacement
#' creates a new instance of \code{target} (e.g., removing "CG" with "GC" in
#' the sequence "CCG" -> "CGC"), the new instance is also replaced. This continues
#' until the sequence is free of the \code{target} k-mer.
#'
#' When \code{target} and \code{replacement} are permutations of each other
#' (e.g., "CG" and "GC"), this acts as a "bubble sort" of nucleotides, moving
#' bases locally without altering the total GC content or base counts of the genome.
#'
#' @return When output_format is "misha" or "fasta", returns invisible NULL and
#'         writes to output_path. When output_format is "vector", returns a
#'         character vector of modified sequences.
#'
#' @examples
#' \dontrun{
#' # Robust removal of all CpG dinucleotides (preserving GC%)
#' gsynth.replace_kmer(
#'     target = "CG",
#'     replacement = "GC",
#'     output_path = "genome_no_cpg.seq",
#'     output_format = "misha"
#' )
#' }
#' @export
gsynth.replace_kmer <- function(target,
                                replacement,
                                output_path = NULL,
                                output_format = c("misha", "fasta", "vector"),
                                intervals = NULL,
                                check_composition = TRUE) {
    .gcheckroot()
    output_format <- match.arg(output_format)

    # Validate inputs
    if (!is.character(target) || length(target) != 1 || nchar(target) == 0) {
        stop("target must be a non-empty string", call. = FALSE)
    }
    if (!is.character(replacement) || length(replacement) != 1 || nchar(replacement) == 0) {
        stop("replacement must be a non-empty string", call. = FALSE)
    }
    if (target == replacement) {
        stop("target and replacement cannot be identical", call. = FALSE)
    }

    # Composition check for robustness
    if (check_composition) {
        t_counts <- table(strsplit(target, "")[[1]])
        r_counts <- table(strsplit(replacement, "")[[1]])
        if (!identical(t_counts, r_counts)) {
            stop(sprintf(
                "Composition mismatch between '%s' and '%s'. Set check_composition=FALSE to override.",
                target, replacement
            ), call. = FALSE)
        }
    }

    # Validate output path
    if (output_format != "vector" && is.null(output_path)) {
        stop("output_path is required when output_format is not 'vector'", call. = FALSE)
    }

    if (is.null(intervals)) {
        intervals <- gintervals.all()
    }

    # Check if we need to process in chunks for large genomes
    # Define chunk processor for parallel execution
    # Captures target, replacement, output_format from enclosing scope
    .process_replace_chunk <- function(chunk_interval, chunk_index, chunk_seed = NULL, ...) {
        if (output_format == "vector") {
            # Vector mode: return sequences directly
            gsynth.replace_kmer(
                target = target,
                replacement = replacement,
                output_path = NULL,
                output_format = "vector",
                intervals = chunk_interval,
                check_composition = FALSE # Already checked above
            )
        } else {
            # File mode: write to temp file and return path
            temp_file <- tempfile(fileext = if (output_format == "fasta") ".fa" else ".seq")
            gsynth.replace_kmer(
                target = target,
                replacement = replacement,
                output_path = temp_file,
                output_format = output_format,
                intervals = chunk_interval,
                check_composition = FALSE
            )
            temp_file
        }
    }

    # Try parallel processing for large genomes
    parallel_result <- .gsynth_process_parallel(
        intervals = intervals,
        output_format = output_format,
        output_path = output_path,
        process_chunk_fn = .process_replace_chunk,
        max_chunk_size = .GSYNTH_MAX_CHUNK_SIZE
    )

    # If parallel processing handled it, return the result
    if (!is.null(parallel_result)) {
        return(parallel_result)
    }
    if (is.null(parallel_result) && sum(intervals$end - intervals$start) > .GSYNTH_MAX_CHUNK_SIZE) {
        # File output was written by parallel processing
        message(sprintf("Modified genome written to: %s", output_path))
        return(invisible(NULL))
    }

    # Map format to integer for C++
    fmt_int <- switch(output_format,
        misha = 0L,
        fasta = 1L,
        vector = 2L
    )

    # Handle NULL path for vector output
    path_str <- if (is.null(output_path)) "" else as.character(output_path)

    message(sprintf("Iteratively replacing '%s' with '%s'...", target, replacement))

    res <- .gcall(
        "C_gsynth_replace_kmer",
        as.character(target),
        as.character(replacement),
        intervals,
        path_str,
        fmt_int,
        .misha_env()
    )

    if (output_format == "vector") {
        return(res)
    }

    message(sprintf("Modified genome written to: %s", output_path))
    invisible(NULL)
}
