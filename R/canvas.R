# Genome Canvas: Generate synthetic genomes from stratified Markov models

#' Create a bin mapping from value-based merge specifications
#'
#' Converts value-based bin merge specifications into a bin_map named vector
#' that can be used with \code{\link{gcanvas.train}}. This allows you to
#' specify merges using actual track values rather than bin indices.
#'
#' @param breaks Numeric vector of bin boundaries (same as used in
#'        \code{\link{gcanvas.train}})
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
#'         in \code{\link{gcanvas.train}}. The names are source bin indices
#'         (1-based), and values are target bin indices (1-based).
#'
#' @examples
#' # Define breaks for GC content [0, 1] in 0.025 increments
#' breaks <- seq(0, 1, 0.025)
#'
#' # Merge all GC content above 70% (0.7) into the bin (0.675, 0.7]
#' bin_map <- gcanvas.bin_map(
#'     breaks = breaks,
#'     merge_ranges = list(
#'         list(from = 0.7, to = c(0.675, 0.7))
#'     )
#' )
#'
#' # Multiple merges: merge low GC (< 0.3) and high GC (> 0.7) into middle bins
#' bin_map2 <- gcanvas.bin_map(
#'     breaks = breaks,
#'     merge_ranges = list(
#'         list(from = c(-Inf, 0.3), to = c(0.4, 0.425)), # low GC -> (0.4, 0.425]
#'         list(from = 0.7, to = c(0.675, 0.7)) # high GC -> (0.675, 0.7]
#'     )
#' )
#'
#' @seealso \code{\link{gcanvas.train}}
#' @export
gcanvas.bin_map <- function(breaks, merge_ranges = NULL) {
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
#' Computes a 5th-order Markov model stratified by bins of a track expression
#' (e.g., GC content). This model can be used to generate synthetic genomes
#' that preserve the k-mer statistics of the original genome within each bin.
#'
#' @param expr Track expression for stratification (e.g., a GC virtual track)
#' @param breaks Numeric vector of bin boundaries. For example, \code{seq(0, 1, 0.025)}
#'        creates 40 bins of 2.5\% width.
#' @param bin_merge Optional list of merge specifications for merging sparse bins.
#'        Each specification is a named list with:
#'   \describe{
#'     \item{from}{Numeric vector of length 2 \code{c(min, max)} defining the
#'           source value range to merge. Use \code{-Inf} or \code{Inf} for
#'           open-ended ranges. Can also be a single number (shorthand for
#'           \code{c(value, Inf)}).}
#'     \item{to}{Numeric vector of length 2 \code{c(min, max)} defining the
#'           target bin that source bins should map to. Must exactly match a bin
#'           defined by \code{breaks}.}
#'   }
#'        For example, to map all GC content above 70\% to the (0.675, 0.7] bin:
#'        \code{bin_merge = list(list(from = 0.7, to = c(0.675, 0.7)))}
#' @param mask Optional intervals to exclude from training. Regions in the mask
#'        will not contribute to k-mer counts. Can be computed using \code{gscreen()}.
#' @param intervals Genomic intervals to process. If NULL, uses all chromosomes.
#' @param iterator Iterator for track evaluation, determines the resolution at which
#'        track values are computed.
#' @param pseudocount Pseudocount added to all k-mer counts to avoid zero probabilities.
#'        Default is 1.
#'
#' @return A \code{gcanvas.model} object containing:
#'   \describe{
#'     \item{num_bins}{Number of stratification bins}
#'     \item{breaks}{Bin boundaries}
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
#' # Create GC virtual track at 200bp resolution
#' gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
#' gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
#' gvtrack.create("masked_frac", NULL, "masked.frac")
#'
#' # Define repeat mask
#' repeats <- gscreen("masked_frac > 0.5",
#'     intervals = gintervals.all(),
#'     iterator = 100
#' )
#'
#' # Train model with GC stratification
#' model <- gcanvas.train(
#'     expr = "g_frac + c_frac",
#'     breaks = seq(0, 1, 0.025),
#'     bin_merge = list(
#'         list(from = 0.7, to = c(0.675, 0.7)) # Map high GC (>70%) to (0.675, 0.7]
#'     ),
#'     mask = repeats,
#'     iterator = 200
#' )
#'
#' @seealso \code{\link{gcanvas.sample}}, \code{\link{gcanvas.save}},
#'          \code{\link{gcanvas.load}}, \code{\link{gcanvas.bin_map}}
#' @export
gcanvas.train <- function(expr,
                          breaks,
                          bin_merge = NULL,
                          mask = NULL,
                          intervals = NULL,
                          iterator = NULL,
                          pseudocount = 1) {
    .gcheckroot()

    # Validate breaks
    if (!is.numeric(breaks) || length(breaks) < 2) {
        stop("breaks must be a numeric vector with at least 2 elements", call. = FALSE)
    }
    breaks <- sort(breaks)
    num_bins <- length(breaks) - 1

    # Convert bin_merge to bin_map_vec if provided
    bin_map_vec <- seq_len(num_bins) # Identity mapping
    if (!is.null(bin_merge)) {
        bin_map <- gcanvas.bin_map(breaks, bin_merge)
        # Convert from named vector format to integer vector
        for (i in seq_along(bin_map)) {
            src_bin <- as.integer(names(bin_map)[i])
            tgt_bin <- as.integer(bin_map[i])
            if (!is.na(src_bin) && !is.na(tgt_bin) &&
                src_bin >= 1 && src_bin <= num_bins &&
                tgt_bin >= 1 && tgt_bin <= num_bins) {
                bin_map_vec[src_bin] <- tgt_bin
            }
        }
    }
    # Convert to 0-based for C++
    bin_map_vec <- bin_map_vec - 1L

    # Convert expression to string
    expr_str <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())

    # Set up iterator
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    # Get intervals (default to all chromosomes)
    if (is.null(intervals)) {
        intervals <- gintervals.all()
    }

    # Extract track values for each iterator position
    message("Extracting track values...")
    track_data <- gextract(expr_str, intervals = intervals, iterator = .iterator)

    if (is.null(track_data) || nrow(track_data) == 0) {
        stop("No track data extracted. Check that intervals and iterator are valid.", call. = FALSE)
    }

    # Assign bin indices to each position
    # Use findInterval to map track values to bins
    track_values <- track_data[[expr_str]]
    bin_indices <- findInterval(track_values, breaks, rightmost.closed = TRUE)
    bin_indices[bin_indices == 0] <- NA # Values below first break
    bin_indices[bin_indices > num_bins] <- NA # Values above last break

    # Get chromosome info for processing
    chrom_sizes <- gintervals.chrom_sizes(intervals)
    chrom_ids <- match(as.character(intervals$chrom), chrom_sizes$chrom) - 1L
    chrom_starts <- intervals$start
    chrom_ends <- intervals$end

    # Prepare iterator position data
    iter_chroms <- match(as.character(track_data$chrom), chrom_sizes$chrom) - 1L
    iter_starts <- as.integer(track_data$start)

    # Handle NA bin indices (convert to -1 for C++)
    bin_indices[is.na(bin_indices)] <- 0
    bin_indices <- as.integer(bin_indices) - 1L # 0-based

    message("Training Markov model...")

    # Call C++ training function
    result <- .gcall(
        "C_gcanvas_train",
        as.integer(chrom_ids),
        as.integer(chrom_starts),
        as.integer(chrom_ends),
        bin_indices,
        iter_starts,
        iter_chroms,
        as.numeric(breaks),
        bin_map_vec,
        mask,
        as.numeric(pseudocount),
        .misha_env()
    )

    # Add class and metadata
    result$expr <- expr_str
    result$iterator <- .iterator

    # Store bin mapping (convert back to 1-based for R)
    # bin_map_vec is 0-based for C++, convert to 1-based for storage
    result$bin_map <- bin_map_vec + 1L

    class(result) <- "gcanvas.model"

    message(sprintf(
        "Trained model: %s 6-mers across %d bins",
        format(result$total_kmers, big.mark = ","),
        result$num_bins
    ))

    result
}

#' Print summary of a gcanvas.model
#'
#' @param x A gcanvas.model object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.gcanvas.model <- function(x, ...) {
    cat("Genome Canvas Markov-5 Model\n")
    cat("----------------------------\n")
    cat(sprintf("Expression: %s\n", x$expr))
    cat(sprintf("Number of bins: %d\n", x$num_bins))
    cat(sprintf("Total k-mers: %s\n", format(x$total_kmers, big.mark = ",")))
    cat(sprintf("Masked positions: %s\n", format(x$total_masked, big.mark = ",")))
    cat(sprintf("N positions: %s\n", format(x$total_n, big.mark = ",")))
    cat("\nPer-bin k-mer counts:\n")
    for (i in seq_along(x$per_bin_kmers)) {
        if (x$per_bin_kmers[i] > 0) {
            cat(sprintf(
                "  Bin %d [%.3f-%.3f]: %s\n",
                i, x$breaks[i], x$breaks[i + 1],
                format(x$per_bin_kmers[i], big.mark = ",")
            ))
        }
    }
    invisible(x)
}

#' Save a gcanvas.model to disk
#'
#' Saves a trained Markov model to an RDS file for later use.
#'
#' @param model A gcanvas.model object from \code{\link{gcanvas.train}}
#' @param file Path to save the model
#'
#' @seealso \code{\link{gcanvas.load}}, \code{\link{gcanvas.train}}
#' @export
gcanvas.save <- function(model, file) {
    if (!inherits(model, "gcanvas.model")) {
        stop("model must be a gcanvas.model object", call. = FALSE)
    }
    saveRDS(model, file)
    invisible(file)
}

#' Load a gcanvas.model from disk
#'
#' Loads a previously saved Markov model from an RDS file.
#'
#' @param file Path to the saved model file
#'
#' @return A gcanvas.model object
#'
#' @seealso \code{\link{gcanvas.save}}, \code{\link{gcanvas.train}}
#' @export
gcanvas.load <- function(file) {
    if (!file.exists(file)) {
        stop(sprintf("File not found: %s", file), call. = FALSE)
    }
    model <- readRDS(file)
    if (!inherits(model, "gcanvas.model")) {
        stop("File does not contain a valid gcanvas.model", call. = FALSE)
    }
    model
}

#' Sample a synthetic genome from a trained Markov model
#'
#' Generates a synthetic genome by sampling from a trained stratified Markov-5
#' model. The generated genome preserves the k-mer statistics of the original
#' genome within each stratification bin.
#'
#' @param model A gcanvas.model object from \code{\link{gcanvas.train}}
#' @param output_path Path to the output file
#' @param output_format Output format: "misha" for .seq binary format,
#'        "fasta" for FASTA text format
#' @param mask Optional intervals for special handling during sampling
#' @param mask_mode How to handle masked regions:
#'   \itemize{
#'     \item "sample": Sample masked regions like any other (ignore mask)
#'     \item "copy": Copy masked regions from the original genome
#'   }
#' @param seed Random seed for reproducibility. If NULL, uses current random state.
#' @param intervals Genomic intervals to sample. If NULL, uses all chromosomes.
#'
#' @return Invisible NULL. The synthetic genome is written to output_path.
#'
#' @examples
#' gdb.init_examples()
#'
#' # Train a model first
#' gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
#' gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
#' gvtrack.create("masked_frac", NULL, "masked.frac")
#'
#' # Define repeat mask
#' repeats <- gscreen("masked_frac > 0.5",
#'     intervals = gintervals.all(),
#'     iterator = 100
#' )
#' model <- gcanvas.train(
#'     expr = "g_frac + c_frac",
#'     breaks = seq(0, 1, 0.025),
#'     mask = repeats,
#'     iterator = 200
#' )
#'
#' # Sample a synthetic genome
#' gcanvas.sample(model, "synthetic_genome.fa",
#'     output_format = "fasta",
#'     mask = repeats,
#'     mask_mode = "copy",
#'     seed = 60427
#' )
#'
#' @seealso \code{\link{gcanvas.train}}, \code{\link{gcanvas.save}}
#' @export
gcanvas.sample <- function(model,
                           output_path,
                           output_format = c("misha", "fasta"),
                           mask = NULL,
                           mask_mode = c("sample", "copy"),
                           seed = NULL,
                           intervals = NULL) {
    .gcheckroot()

    if (!inherits(model, "gcanvas.model")) {
        stop("model must be a gcanvas.model object", call. = FALSE)
    }

    output_format <- match.arg(output_format)
    mask_mode <- match.arg(mask_mode)

    # Set random seed if provided
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Get intervals (default to all chromosomes)
    if (is.null(intervals)) {
        intervals <- gintervals.all()
    }

    # Re-extract track values for sampling (using the same expression and iterator)
    .iterator <- model$iterator
    expr_str <- model$expr

    message("Extracting track values for sampling...")
    track_data <- gextract(expr_str, intervals = intervals, iterator = .iterator)

    if (is.null(track_data) || nrow(track_data) == 0) {
        stop("No track data extracted. Check that intervals are valid.", call. = FALSE)
    }

    # Assign bin indices to each position
    track_values <- track_data[[expr_str]]
    bin_indices <- findInterval(track_values, model$breaks, rightmost.closed = TRUE)
    bin_indices[bin_indices == 0] <- NA
    bin_indices[bin_indices > model$num_bins] <- NA

    # Apply bin mapping if present (from bin_merge during training)
    # This remaps bins that were merged during training to their target bins
    if (!is.null(model$bin_map)) {
        # bin_indices are 1-based here, bin_map is 1-based
        valid_mask <- !is.na(bin_indices) &
            bin_indices >= 1 &
            bin_indices <= length(model$bin_map)
        if (any(valid_mask)) {
            # Apply mapping: bin_indices[i] -> bin_map[bin_indices[i]]
            mapped_bins <- model$bin_map[bin_indices[valid_mask]]
            # Ensure mapped bins are valid (should be, but check for safety)
            mapped_bins <- pmin(pmax(mapped_bins, 1L), model$num_bins)
            bin_indices[valid_mask] <- mapped_bins
        }
    }

    # Get chromosome info
    chrom_sizes <- gintervals.chrom_sizes(intervals)
    iter_chroms <- match(as.character(track_data$chrom), chrom_sizes$chrom) - 1L
    iter_starts <- as.integer(track_data$start)

    # Handle NA bin indices (convert to -1 for C++)
    bin_indices[is.na(bin_indices)] <- 0
    bin_indices <- as.integer(bin_indices) - 1L # 0-based

    # Get CDF data from model
    cdf_list <- model$model_data$cdf

    # Output format: 0 = misha, 1 = fasta
    output_format_int <- if (output_format == "fasta") 1L else 0L

    # Mask mode: 0 = sample, 1 = copy
    mask_mode_int <- if (mask_mode == "copy") 1L else 0L

    message("Sampling synthetic genome...")

    # Call C++ sampling function
    .gcall(
        "C_gcanvas_sample",
        cdf_list,
        as.numeric(model$breaks),
        bin_indices,
        iter_starts,
        iter_chroms,
        intervals,
        mask,
        mask_mode_int,
        as.character(output_path),
        output_format_int,
        .misha_env()
    )

    message(sprintf("Synthetic genome written to: %s", output_path))

    invisible(NULL)
}
