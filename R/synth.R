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
#' @seealso \code{\link{gsynth.train}}
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

#' Train a stratified Markov-5 model from genome sequences
#'
#' Computes a 5th-order Markov model optionally stratified by bins of one or more
#' track expressions (e.g., GC content and CG dinucleotide frequency). This model
#' can be used to generate synthetic genomes that preserve the k-mer statistics of
#' the original genome within each stratification bin. When called with no dimension
#' specifications, trains a single unstratified model.
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
#' @param min_obs Minimum number of observations (6-mers) required per bin. Bins with
#'        fewer observations will be marked as NA (not learned) and a warning will be
#'        issued. Default is 0 (no minimum). During sampling, NA bins will fall back
#'        to uniform sampling unless merged via \code{bin_merge}.
#'
#' @details
#' \strong{Strand symmetry:} The training process counts both the forward strand
#' 6-mer and its reverse complement for each position, ensuring strand-symmetric
#' transition probabilities. This means the reported total_kmers is approximately
#' double the number of genomic positions processed.
#'
#' \strong{N bases:} Positions where the 6-mer contains any N (unknown) bases are
#' skipped during training and counted in \code{total_n}. The model only learns
#' from valid A/C/G/T sequences.
#'
#' @return A \code{gsynth.model} object containing:
#'   \describe{
#'     \item{n_dims}{Number of stratification dimensions}
#'     \item{dim_specs}{List of dimension specifications (expr, breaks, num_bins, bin_map)}
#'     \item{dim_sizes}{Vector of bin counts per dimension}
#'     \item{total_bins}{Total number of bins (product of dim_sizes)}
#'     \item{total_kmers}{Total number of valid 6-mers counted}
#'     \item{per_bin_kmers}{Number of 6-mers counted per bin}
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
                         min_obs = 0) {
    .gcheckroot()

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

        # Get chromosome info for processing
        chrom_sizes <- gintervals.chrom_sizes(intervals)
        chrom_ids <- match(as.character(intervals$chrom), chrom_sizes$chrom) - 1L
        chrom_starts <- intervals$start
        chrom_ends <- intervals$end

        # Prepare iterator position data
        iter_chroms <- match(as.character(iter_info$chrom), chrom_sizes$chrom) - 1L
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

        # Get chromosome info for processing
        chrom_sizes <- gintervals.chrom_sizes(intervals)
        chrom_ids <- match(as.character(intervals$chrom), chrom_sizes$chrom) - 1L
        chrom_starts <- intervals$start
        chrom_ends <- intervals$end

        # Prepare iterator position data
        iter_chroms <- match(as.character(track_data$chrom), chrom_sizes$chrom) - 1L
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
        .misha_env()
    )

    # Override C++ result with our multi-dimensional metadata
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
            "Trained unstratified model: %s 6-mers (no stratification)",
            format(result$total_kmers, big.mark = ",")
        ))
    } else {
        message(sprintf(
            "Trained model: %s 6-mers across %d bins (%d dimensions)",
            format(result$total_kmers, big.mark = ","),
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
    cat("Synthetic Genome Markov-5 Model\n")
    cat("----------------------------\n")

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

#' Save a gsynth.model to disk
#'
#' Saves a trained Markov model to an RDS file for later use.
#'
#' @param model A gsynth.model object from \code{\link{gsynth.train}}
#' @param file Path to save the model
#'
#' @seealso \code{\link{gsynth.load}}, \code{\link{gsynth.train}}
#' @export
gsynth.save <- function(model, file) {
    if (!inherits(model, "gsynth.model")) {
        stop("model must be a gsynth.model object", call. = FALSE)
    }
    saveRDS(model, file)
    invisible(file)
}

#' Load a gsynth.model from disk
#'
#' Loads a previously saved Markov model from an RDS file.
#'
#' @param file Path to the saved model file
#'
#' @return A gsynth.model object
#'
#' @seealso \code{\link{gsynth.save}}, \code{\link{gsynth.train}}
#' @export
gsynth.load <- function(file) {
    if (!file.exists(file)) {
        stop(sprintf("File not found: %s", file), call. = FALSE)
    }
    model <- readRDS(file)
    if (!inherits(model, "gsynth.model")) {
        stop("File does not contain a valid gsynth.model", call. = FALSE)
    }
    model
}

#' Sample a synthetic genome from a trained Markov model
#'
#' Generates a synthetic genome by sampling from a trained stratified Markov-5
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
#'
#' @details
#' \strong{N bases during sampling:} When the sampler needs to initialize the
#' first 5-mer context and encounters regions with only N bases, it falls back
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
#' gsynth.sample(model, "synthetic_genome.fa",
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
                          bin_merge = NULL) {
    .gcheckroot()

    if (!inherits(model, "gsynth.model")) {
        stop("model must be a gsynth.model object", call. = FALSE)
    }

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

        # Get chromosome info
        chrom_sizes <- gintervals.chrom_sizes(intervals)
        iter_chroms <- match(as.character(iter_info$chrom), chrom_sizes$chrom) - 1L
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

        # Get chromosome info
        chrom_sizes <- gintervals.chrom_sizes(intervals)
        iter_chroms <- match(as.character(track_data$chrom), chrom_sizes$chrom) - 1L
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
                rep(c(0.25, 0.5, 0.75, 1.0), each = 1024),
                nrow = 1024, ncol = 4
            )
            for (bin_idx in sparse_used) {
                cdf_list[[bin_idx]] <- uniform_cdf
            }
        }
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
        .misha_env()
    )

    if (output_format == "vector") {
        message(sprintf("Generated %d sequence(s)", length(result)))
        return(result)
    }

    message(sprintf("Synthetic genome written to: %s", output_path))

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

    # Get chromosome info
    chrom_sizes <- gintervals.chrom_sizes(intervals)
    iter_chroms <- match(as.character(iter_info$chrom), chrom_sizes$chrom) - 1L
    iter_starts <- as.integer(iter_info$start)

    # Create CDF from probabilities: [P(A), P(A)+P(C), P(A)+P(C)+P(G), 1]
    # Order is A=0, C=1, G=2, T=3 (same as C++ code)
    cdf_row <- cumsum(as.numeric(nuc_probs))

    # Create CDF matrix: same for all 1024 contexts (no context dependency)
    # Each row is the same CDF - pure random, context-independent sampling
    cdf_matrix <- matrix(rep(cdf_row, 1024), nrow = 1024, ncol = 4, byrow = TRUE)
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
        .misha_env()
    )

    if (output_format == "vector") {
        message(sprintf("Generated %d random sequence(s)", length(result)))
        return(result)
    }

    message(sprintf("Random sequences written to: %s", output_path))

    invisible(NULL)
}
