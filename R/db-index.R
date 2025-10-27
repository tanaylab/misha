#' Change Database to Indexed Genome Format
#'
#' Converts a per-chromosome database to indexed genome format
#' with a single consolidated genome.seq file and genome.idx index.
#' Optionally also converts tracks and interval sets to indexed format.
#'
#' @param groot Root directory of the database to change to indexed format. If NULL, uses the currently active database.
#' @param remove_old_files Logical. If TRUE, removes old per-chromosome files after successful conversion. Default: FALSE.
#' @param force Logical. If TRUE, forces the conversion without confirmation. Default: FALSE.
#' @param validate Logical. If TRUE, validates the conversion by comparing sequences. Default: TRUE.
#' @param convert_tracks Logical. If TRUE, also converts all eligible tracks to indexed format. Default: FALSE.
#' @param convert_intervals Logical. If TRUE, also converts all eligible interval sets to indexed format. Default: FALSE.
#' @param verbose Logical. If TRUE, prints verbose messages. Default: FALSE.
#'
#' @return Invisible NULL
#'
#' @details
#' This function converts a per-chromosome database (with separate .seq files per contig) to
#' indexed format (single genome.seq + genome.idx). The indexed format
#' provides better performance and scalability, especially for genomes with many contigs.
#'
#' The conversion process:
#' \enumerate{
#'   \item Checks if database is already in indexed format
#'   \item Reads existing chromosome information from chrom_sizes.txt
#'   \item Consolidates all per-chromosome .seq files into genome.seq
#'   \item Creates genome.idx with CRC64 checksum
#'   \item Optionally validates the conversion
#'   \item Optionally removes old .seq files
#'   \item If convert_tracks=TRUE, converts all eligible 1D tracks (dense, sparse, array)
#'   \item If convert_intervals=TRUE, converts all eligible interval sets (1D and 2D)
#' }
#'
#' Tracks and intervals that cannot be converted (and are skipped):
#' \itemize{
#'   \item Tracks: 2D tracks, virtual tracks, single-file tracks, already converted tracks
#'   \item Intervals: Single-file interval sets, already converted interval sets
#' }
#'
#' @examples
#' \dontrun{
#' # Convert current database to indexed format (genome only)
#' gdb.convert_to_indexed()
#'
#' # Convert specific database to indexed format
#' gdb.convert_to_indexed(groot = "/path/to/database")
#'
#' # Convert genome and all tracks to indexed format
#' gdb.convert_to_indexed(convert_tracks = TRUE)
#'
#' # Convert genome, tracks, and intervals to indexed format
#' gdb.convert_to_indexed(convert_tracks = TRUE, convert_intervals = TRUE)
#'
#' # Full conversion with cleanup
#' gdb.convert_to_indexed(convert_tracks = TRUE, convert_intervals = TRUE, remove_old_files = TRUE)
#' }
#'
#' @seealso \code{\link{gdb.create}}, \code{\link{gdb.init}}, \code{\link{gtrack.convert_to_indexed}}, \code{\link{gintervals.convert_to_indexed}}, \code{\link{gintervals.2d.to_indexed_format}}
#' @export
gdb.convert_to_indexed <- function(groot = NULL, remove_old_files = FALSE, force = FALSE, validate = TRUE, convert_tracks = FALSE, convert_intervals = FALSE, verbose = FALSE) {
    # Validate database and get setup information
    setup_info <- .gdb.convert_to_indexed.validate_and_setup(groot, verbose)

    # Return early if already indexed
    if (setup_info$already_indexed) {
        return(invisible(NULL))
    }

    # Get user confirmation
    if (!.gdb.convert_to_indexed.get_confirmation(setup_info$groot, setup_info$chrom_sizes, remove_old_files, force)) {
        return(invisible(NULL))
    }

    # Convert genome sequences
    .gdb.convert_to_indexed.genome(setup_info, validate, remove_old_files, verbose)

    # Convert tracks if requested
    if (convert_tracks) {
        .gdb.convert_to_indexed.tracks(setup_info$groot, verbose)
    }

    # Convert intervals if requested
    if (convert_intervals) {
        .gdb.convert_to_indexed.intervals(setup_info$groot, remove_old_files, verbose)
    }

    if (verbose) message("\n=== Conversion Complete ===")

    invisible(NULL)
}

# Helper function to validate database and get chromosome information
.gdb.convert_to_indexed.validate_and_setup <- function(groot, verbose = FALSE) {
    # Use current database if not specified
    if (is.null(groot)) {
        groot <- get("GROOT", envir = .misha)
        if (is.null(groot) || groot == "") {
            stop("No database is currently active. Please call gdb.init() or specify groot parameter.", call. = FALSE)
        }
    }

    # Check if database exists
    if (!dir.exists(groot)) {
        stop(sprintf("Database directory does not exist: %s", groot), call. = FALSE)
    }

    seq_dir <- file.path(groot, "seq")
    if (!dir.exists(seq_dir)) {
        stop(sprintf("seq directory does not exist: %s", seq_dir), call. = FALSE)
    }

    # Check if already in indexed format
    index_path <- file.path(seq_dir, "genome.idx")
    genome_seq_path <- file.path(seq_dir, "genome.seq")

    if (file.exists(index_path) && file.exists(genome_seq_path)) {
        if (verbose) message("Database is already in indexed format.")
        return(list(already_indexed = TRUE, groot = groot))
    }

    # Get canonical chromosome names from ALLGENOME (if database is initialized)
    # This preserves the "chr" prefix for per-chromosome databases
    canonical_names <- NULL
    if (exists("ALLGENOME", envir = .misha, inherits = FALSE) &&
        !is.null(get("ALLGENOME", envir = .misha)) &&
        !is.null(get("GROOT", envir = .misha)) &&
        get("GROOT", envir = .misha) == groot) {
        allgenome <- get("ALLGENOME", envir = .misha)
        if (!is.null(allgenome[[1]]) && "chrom" %in% colnames(allgenome[[1]])) {
            canonical_names <- as.character(allgenome[[1]]$chrom)
        }
    }

    # Read chromosome information
    chrom_sizes_path <- file.path(groot, "chrom_sizes.txt")
    if (!file.exists(chrom_sizes_path)) {
        stop(sprintf("chrom_sizes.txt not found: %s", chrom_sizes_path), call. = FALSE)
    }

    chrom_sizes <- read.table(chrom_sizes_path, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    colnames(chrom_sizes) <- c("chrom", "size")

    # If we have canonical names from ALLGENOME, use them
    # This preserves the "chr" prefix for per-chromosome databases
    if (!is.null(canonical_names) && length(canonical_names) == nrow(chrom_sizes)) {
        # Sort to match ALLGENOME order
        chrom_sizes <- chrom_sizes[order(chrom_sizes$chrom), ]
        # Replace with canonical names (which are already sorted in ALLGENOME)
        chrom_sizes$chrom <- canonical_names
    } else {
        # Fallback: just sort by chromosome name
        chrom_sizes <- chrom_sizes[order(chrom_sizes$chrom), ]
    }

    # Check that per-chromosome .seq files exist
    # Handle chr prefix mismatch between chrom_sizes.txt and .seq files
    seq_files <- character(nrow(chrom_sizes))
    for (i in seq_len(nrow(chrom_sizes))) {
        chrom <- chrom_sizes$chrom[i]

        # Try the chromosome name as-is first
        seq_file <- file.path(seq_dir, paste0(chrom, ".seq"))
        if (file.exists(seq_file)) {
            seq_files[i] <- seq_file
            next
        }

        # If not found, try with chr prefix
        if (!startsWith(chrom, "chr")) {
            chr_seq_file <- file.path(seq_dir, paste0("chr", chrom, ".seq"))
            if (file.exists(chr_seq_file)) {
                seq_files[i] <- chr_seq_file
                next
            }
        }

        # If not found, try without chr prefix
        if (startsWith(chrom, "chr")) {
            no_chr_seq_file <- file.path(seq_dir, paste0(sub("^chr", "", chrom), ".seq"))
            if (file.exists(no_chr_seq_file)) {
                seq_files[i] <- no_chr_seq_file
                next
            }
        }

        # If still not found, this is a missing file
        seq_files[i] <- seq_file # Use original name for error reporting
    }

    missing_files <- seq_files[!file.exists(seq_files)]
    if (length(missing_files) > 0) {
        stop(sprintf("Missing sequence files: %s", paste(basename(missing_files), collapse = ", ")), call. = FALSE)
    }

    return(list(
        already_indexed = FALSE,
        groot = groot,
        seq_dir = seq_dir,
        chrom_sizes = chrom_sizes,
        seq_files = seq_files,
        index_path = index_path,
        genome_seq_path = genome_seq_path,
        chrom_sizes_path = chrom_sizes_path
    ))
}

# Helper function to get user confirmation for conversion
.gdb.convert_to_indexed.get_confirmation <- function(groot, chrom_sizes, remove_old_files, force) {
    if (interactive() && !force) {
        cat(sprintf("About to convert database to indexed format: %s\n", groot))
        cat(sprintf("  Chromosomes: %d\n", nrow(chrom_sizes)))
        cat(sprintf("  Total size: %.2f MB\n", sum(chrom_sizes$size) / 1024^2))
        if (remove_old_files) {
            cat("  Old .seq files will be REMOVED after conversion\n")
        }
        response <- readline("Proceed with conversion? (yes/no): ")
        if (!(tolower(response) %in% c("yes", "y"))) {
            message("Conversion cancelled.")
            return(FALSE)
        }
    }
    return(TRUE)
}

# Helper function to convert genome sequences to indexed format
.gdb.convert_to_indexed.genome <- function(setup_info, validate = TRUE, remove_old_files = FALSE, verbose = FALSE) {
    groot <- setup_info$groot
    chrom_sizes <- setup_info$chrom_sizes
    seq_files <- setup_info$seq_files
    index_path <- setup_info$index_path
    genome_seq_path <- setup_info$genome_seq_path
    chrom_sizes_path <- setup_info$chrom_sizes_path

    if (verbose) message("Converting database to indexed format...")

    # Create temporary FASTA file from .seq files
    temp_fasta <- tempfile(fileext = ".fasta")
    on.exit(unlink(temp_fasta), add = TRUE)

    tryCatch(
        {
            # Write multi-FASTA file
            if (verbose) message("Creating temporary multi-FASTA file...")
            fasta_con <- file(temp_fasta, "w")

            for (i in seq_len(nrow(chrom_sizes))) {
                chrom <- chrom_sizes$chrom[i]
                seq_file <- seq_files[i]

                # Write header
                cat(sprintf(">%s\n", chrom), file = fasta_con)

                # Read and write sequence (chunk by chunk to handle large chromosomes)
                seq_con <- file(seq_file, "rb")
                chunk_size <- 1048576 # 1 MB chunks

                repeat {
                    chunk <- readBin(seq_con, "raw", n = chunk_size)
                    if (length(chunk) == 0) break

                    # Convert to character and write
                    seq_str <- rawToChar(chunk)
                    # Write in lines of 80 characters (standard FASTA format)
                    seq_chars <- strsplit(seq_str, "")[[1]]
                    for (j in seq(1, length(seq_chars), 80)) {
                        end_idx <- min(j + 79, length(seq_chars))
                        cat(paste(seq_chars[j:end_idx], collapse = ""), "\n", file = fasta_con)
                    }
                }

                close(seq_con)

                if ((i %% 10) == 0 || i == nrow(chrom_sizes)) {
                    if (verbose) message(sprintf("  Processed %d/%d chromosomes", i, nrow(chrom_sizes)))
                }
            }

            close(fasta_con)

            # Call C++ import function
            if (verbose) message("Creating indexed format...")
            contig_info <- .gcall(
                "gseq_multifasta_import",
                temp_fasta,
                genome_seq_path,
                index_path,
                .misha_env()
            )

            if (verbose) message("Index created successfully")

            # Update chrom_sizes.txt with the canonical chromosome names from the index
            # This ensures consistency between chrom_sizes.txt and genome.idx
            if (verbose) message("Updating chrom_sizes.txt with canonical chromosome names...")
            updated_chrom_sizes <- data.frame(
                chrom = contig_info$name,
                size = contig_info$size
            )
            write.table(updated_chrom_sizes, chrom_sizes_path,
                quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
            )

            # Validate if requested
            if (validate) {
                if (verbose) message("Validating conversion...")

                # Sample validation: check first 100 bases of each chromosome
                validation_failed <- FALSE
                for (i in seq_len(min(10, nrow(chrom_sizes)))) { # Check first 10 chroms
                    chrom <- chrom_sizes$chrom[i]
                    seq_file <- seq_files[i]

                    # Read from old file
                    old_con <- file(seq_file, "rb")
                    old_seq <- rawToChar(readBin(old_con, "raw", n = 100))
                    close(old_con)

                    # Read from indexed format (need to reload database first)
                    # Save current state
                    old_groot <- NULL
                    if (exists("GROOT", envir = .misha, inherits = FALSE)) {
                        old_groot <- get("GROOT", envir = .misha)
                    }

                    tryCatch({
                        # Temporarily init the converted database
                        suppressMessages(gdb.init(groot))

                        # Extract from indexed format
                        new_seq <- gseq.extract(gintervals(chrom, 0, min(100, chrom_sizes$size[i])))

                        if (old_seq != new_seq) {
                            warning(sprintf("Validation failed for chromosome %s", chrom))
                            validation_failed <- TRUE
                        }
                    }, finally = {
                        # Restore old state
                        if (!is.null(old_groot) && old_groot != "") {
                            suppressMessages(gdb.init(old_groot))
                        }
                    })
                }

                if (validation_failed) {
                    stop("Validation failed! Conversion may be corrupted. Old files have NOT been removed.", call. = FALSE)
                } else {
                    if (verbose) message("Validation passed")
                }
            }

            # Remove old files if requested
            if (remove_old_files) {
                if (verbose) message("Removing old .seq files...")
                for (seq_file in seq_files) {
                    unlink(seq_file)
                }
                if (verbose) message(sprintf("Removed %d old .seq files", length(seq_files)))
            }

            if (verbose) message(sprintf("Database sequence conversion complete: %s", groot))
        },
        error = function(e) {
            # Clean up partial files on error
            if (file.exists(genome_seq_path)) {
                unlink(genome_seq_path)
            }
            if (file.exists(index_path)) {
                unlink(index_path)
            }
            stop(sprintf("Conversion failed: %s", conditionMessage(e)), call. = FALSE)
        }
    )
}

# Helper function to convert tracks to indexed format
.gdb.convert_to_indexed.tracks <- function(groot, verbose = FALSE) {
    if (verbose) message("\n=== Converting Tracks ===")

    # Temporarily init the database to get track list
    old_groot <- NULL
    if (exists("GROOT", envir = .misha, inherits = FALSE)) {
        old_groot <- get("GROOT", envir = .misha)
    }

    tryCatch({
        suppressMessages(gdb.init(groot))

        all_tracks <- gtrack.ls()

        if (length(all_tracks) == 0) {
            if (verbose) message("No tracks found in database")
        } else {
            if (verbose) message(sprintf("Found %d tracks in database", length(all_tracks)))

            # Filter to only 1D tracks that can be converted
            convertible_tracks <- c()
            skipped_tracks <- list()

            for (track in all_tracks) {
                # Get track info to determine type
                info <- tryCatch(gtrack.info(track), error = function(e) NULL)

                if (is.null(info)) {
                    skipped_tracks[[track]] <- "failed to get info"
                    next
                }

                # Only 1D tracks (dense, sparse, array) can be converted
                if (!info$type %in% c("dense", "sparse", "array")) {
                    skipped_tracks[[track]] <- sprintf("unsupported type (%s)", info$type)
                    next
                }

                # Check if it's a Big Set track (directory format)
                trackstr <- gsub("\\.", "/", track)
                trackdir <- sprintf("%s.track", paste(groot, "tracks", trackstr, sep = "/"))

                if (!dir.exists(trackdir)) {
                    skipped_tracks[[track]] <- "single-file format"
                    next
                }

                # Check if already converted
                idx_path <- file.path(trackdir, "track.idx")
                if (file.exists(idx_path)) {
                    skipped_tracks[[track]] <- "already converted"
                    next
                }

                convertible_tracks <- c(convertible_tracks, track)
            }

            # Report what we found
            if (length(convertible_tracks) > 0) {
                if (verbose) message(sprintf("  Convertible: %d tracks", length(convertible_tracks)))
            } else {
                if (verbose) message("  No tracks need conversion")
            }

            if (length(skipped_tracks) > 0) {
                if (verbose) message(sprintf("  Skipped: %d tracks", length(skipped_tracks)))
                for (track in names(skipped_tracks)) {
                    if (verbose) message(sprintf("    - %s: %s", track, skipped_tracks[[track]]))
                }
            }

            # Convert tracks
            if (length(convertible_tracks) > 0) {
                converted_count <- 0
                failed_tracks <- c()

                for (track in convertible_tracks) {
                    if (verbose) message(sprintf("  Converting track: %s", track))

                    tryCatch(
                        {
                            gtrack.convert_to_indexed(track)
                            converted_count <- converted_count + 1
                        },
                        error = function(e) {
                            warning(sprintf("Failed to convert track %s: %s", track, conditionMessage(e)))
                            failed_tracks <<- c(failed_tracks, track)
                        }
                    )
                }

                if (verbose) message(sprintf("Successfully converted %d/%d tracks", converted_count, length(convertible_tracks)))

                if (length(failed_tracks) > 0) {
                    warning(sprintf(
                        "Failed to convert %d tracks: %s",
                        length(failed_tracks),
                        paste(failed_tracks, collapse = ", ")
                    ))
                }
            }
        }
    }, finally = {
        # Restore old state
        if (!is.null(old_groot) && old_groot != "") {
            suppressMessages(gdb.init(old_groot))
        }
    })
}

# Helper function to convert interval sets to indexed format
.gdb.convert_to_indexed.intervals <- function(groot, remove_old_files = FALSE, verbose = FALSE) {
    if (verbose) message("\n=== Converting Interval Sets ===")

    # Temporarily init the database to get interval list
    old_groot <- NULL
    if (exists("GROOT", envir = .misha, inherits = FALSE)) {
        old_groot <- get("GROOT", envir = .misha)
    }

    tryCatch({
        suppressMessages(gdb.init(groot))

        all_intervals <- gintervals.ls()

        if (length(all_intervals) == 0) {
            if (verbose) message("No interval sets found in database")
        } else {
            if (verbose) message(sprintf("Found %d interval sets in database", length(all_intervals)))

            # Filter to only Big Set intervals that can be converted
            convertible_1d <- c()
            convertible_2d <- c()
            skipped_intervals <- list()

            for (intervset in all_intervals) {
                # Construct path
                path <- gsub("\\.", "/", intervset)
                intervset_path <- paste0(groot, "/tracks/", path, ".interv")

                # Check if it's a Big Set (directory format)
                if (!file.exists(intervset_path)) {
                    skipped_intervals[[intervset]] <- "does not exist"
                    next
                }

                if (!dir.exists(intervset_path)) {
                    skipped_intervals[[intervset]] <- "single-file format"
                    next
                }

                # Determine if 1D or 2D
                idx_path_1d <- file.path(intervset_path, "intervals.idx")
                idx_path_2d <- file.path(intervset_path, "intervals2d.idx")

                # Check for 2D files
                pair_files <- list.files(intervset_path, pattern = "-")
                has_2d <- length(pair_files) > 0

                if (has_2d) {
                    # 2D interval set
                    if (file.exists(idx_path_2d)) {
                        skipped_intervals[[intervset]] <- "already converted (2D)"
                    } else {
                        convertible_2d <- c(convertible_2d, intervset)
                    }
                } else {
                    # 1D interval set
                    if (file.exists(idx_path_1d)) {
                        skipped_intervals[[intervset]] <- "already converted (1D)"
                    } else {
                        convertible_1d <- c(convertible_1d, intervset)
                    }
                }
            }

            # Report what we found
            if (length(convertible_1d) > 0) {
                if (verbose) message(sprintf("  Convertible 1D: %d interval sets", length(convertible_1d)))
            }
            if (length(convertible_2d) > 0) {
                if (verbose) message(sprintf("  Convertible 2D: %d interval sets", length(convertible_2d)))
            }
            if (length(convertible_1d) == 0 && length(convertible_2d) == 0) {
                if (verbose) message("  No interval sets need conversion")
            }

            if (length(skipped_intervals) > 0) {
                if (verbose) message(sprintf("  Skipped: %d interval sets", length(skipped_intervals)))
                for (intervset in names(skipped_intervals)) {
                    if (verbose) message(sprintf("    - %s: %s", intervset, skipped_intervals[[intervset]]))
                }
            }

            # Convert 1D intervals
            converted_1d_count <- 0
            failed_1d <- c()

            if (length(convertible_1d) > 0) {
                for (intervset in convertible_1d) {
                    if (verbose) message(sprintf("  Converting 1D interval set: %s", intervset))

                    tryCatch(
                        {
                            gintervals.convert_to_indexed(intervset, remove.old = remove_old_files)
                            converted_1d_count <- converted_1d_count + 1
                        },
                        error = function(e) {
                            warning(sprintf("Failed to convert 1D interval set %s: %s", intervset, conditionMessage(e)))
                            failed_1d <<- c(failed_1d, intervset)
                        }
                    )
                }
            }

            # Convert 2D intervals
            converted_2d_count <- 0
            failed_2d <- c()

            if (length(convertible_2d) > 0) {
                for (intervset in convertible_2d) {
                    if (verbose) message(sprintf("  Converting 2D interval set: %s", intervset))

                    tryCatch(
                        {
                            gintervals.2d.to_indexed_format(intervset, remove.old = remove_old_files)
                            converted_2d_count <- converted_2d_count + 1
                        },
                        error = function(e) {
                            warning(sprintf("Failed to convert 2D interval set %s: %s", intervset, conditionMessage(e)))
                            failed_2d <<- c(failed_2d, intervset)
                        }
                    )
                }
            }

            # Report results
            total_converted <- converted_1d_count + converted_2d_count
            total_convertible <- length(convertible_1d) + length(convertible_2d)

            if (total_convertible > 0) {
                if (verbose) {
                    message(sprintf(
                        "Successfully converted %d/%d interval sets (%d 1D, %d 2D)",
                        total_converted,
                        total_convertible,
                        converted_1d_count,
                        converted_2d_count
                    ))
                }

                if (length(failed_1d) > 0 || length(failed_2d) > 0) {
                    all_failed <- c(failed_1d, failed_2d)
                    warning(sprintf(
                        "Failed to convert %d interval sets: %s",
                        length(all_failed),
                        paste(all_failed, collapse = ", ")
                    ))
                }
            }
        }
    }, finally = {
        # Restore old state
        if (!is.null(old_groot) && old_groot != "") {
            suppressMessages(gdb.init(old_groot))
        }
    })
}


#' Get Database Information
#'
#' Returns information about a misha genome database including format, number of chromosomes,
#' total genome size, and whether it uses the indexed format.
#'
#' @param groot Root directory of the database. If NULL, uses the currently active database.
#' @return A list with database information:
#' \itemize{
#'   \item \code{path} - Full path to the database
#'   \item \code{is_db} - TRUE if this is a valid misha database
#'   \item \code{format} - "indexed" or "per-chromosome"
#'   \item \code{num_chromosomes} - Number of chromosomes/contigs
#'   \item \code{genome_size} - Total length of genome in bases
#'   \item \code{chromosomes} - Data frame with chromosome names and sizes
#' }
#'
#' @examples
#' \dontrun{
#' # Get info about currently active database
#' info <- gdb.info()
#' cat("Database format:", info$format, "\n")
#' cat("Genome size:", info$genome_size / 1e6, "Mb\n")
#'
#' # Get info about specific database
#' info <- gdb.info("/path/to/database")
#' }
#'
#' @export gdb.info
gdb.info <- function(groot = NULL) {
    # Use current database if not specified
    if (is.null(groot)) {
        if (!exists("GROOT", envir = .misha) || is.null(get("GROOT", envir = .misha))) {
            stop("No database is currently active. Please call gdb.init() or specify groot parameter.", call. = FALSE)
        }
        groot <- get("GROOT", envir = .misha)
    }

    # Normalize path
    groot <- normalizePath(groot, mustWork = FALSE)

    # Check if directory exists
    if (!dir.exists(groot)) {
        return(list(
            path = groot,
            is_db = FALSE,
            error = "Directory does not exist"
        ))
    }

    # Check for chrom_sizes.txt
    chrom_sizes_path <- file.path(groot, "chrom_sizes.txt")
    if (!file.exists(chrom_sizes_path)) {
        return(list(
            path = groot,
            is_db = FALSE,
            error = "Not a misha database (chrom_sizes.txt not found)"
        ))
    }

    # Read chromosome information
    chrom_sizes <- tryCatch(
        read.csv(chrom_sizes_path,
            sep = "\t", header = FALSE,
            col.names = c("chrom", "size"), colClasses = c("character", "numeric")
        ),
        error = function(e) NULL
    )

    if (is.null(chrom_sizes)) {
        return(list(
            path = groot,
            is_db = FALSE,
            error = "Invalid chrom_sizes.txt format"
        ))
    }

    # Detect format
    idx_path <- file.path(groot, "seq", "genome.idx")
    genome_seq_path <- file.path(groot, "seq", "genome.seq")

    if (file.exists(idx_path) && file.exists(genome_seq_path)) {
        format <- "indexed"
    } else {
        format <- "per-chromosome"
    }

    # Calculate total genome size
    genome_size <- sum(chrom_sizes$size)

    list(
        path = groot,
        is_db = TRUE,
        format = format,
        num_chromosomes = nrow(chrom_sizes),
        genome_size = genome_size,
        chromosomes = chrom_sizes
    )
}

#' Convert a track to indexed format
#'
#' Converts a per-chromosome track to indexed format (track.dat + track.idx).
#'
#' This function converts a track from the per-chromosome file format to
#' single-file indexed format. The indexed format dramatically reduces file descriptor
#' usage for genomes with many contigs and provides better performance for parallel access.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates that all per-chromosome files have consistent metadata
#'   \item Creates track.dat by concatenating all per-chromosome files
#'   \item Creates track.idx with offset/length information for each chromosome
#'   \item Uses atomic operations (fsync + rename) to ensure data integrity
#'   \item Removes the old per-chromosome files after successful conversion
#' }
#'
#' @param track track name to convert
#' @return None
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.create_dense}}
#' @examples
#' \dontrun{
#' # Convert a track to indexed format
#' gtrack.convert_to_indexed("my_track")
#' }
#' @export gtrack.convert_to_indexed
gtrack.convert_to_indexed <- function(track = NULL) {
    if (is.null(substitute(track))) {
        stop("Usage: gtrack.convert_to_indexed(track)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    if (is.na(match(trackstr, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Track %s does not exist", trackstr), call. = FALSE)
    }

    trackdir <- sprintf("%s.track", paste(get("GWD", envir = .misha), gsub("\\.", "/", trackstr), sep = "/"))
    idx_path <- file.path(trackdir, "track.idx")
    dat_path <- file.path(trackdir, "track.dat")

    # Check if already converted
    if (file.exists(idx_path)) {
        message(sprintf("Track %s is already in indexed format.", trackstr))
        return(invisible(0))
    }

    # Get track info to determine type
    info <- gtrack.info(track)
    track_type <- info$type

    # Only 1D tracks can be converted
    if (!track_type %in% c("dense", "sparse", "array")) {
        stop(sprintf("Cannot convert track %s: only 1D tracks (dense, sparse, array) can be converted", trackstr), call. = FALSE)
    }

    # Call C++ function to perform the conversion (always remove old files)
    success <- FALSE
    tryCatch(
        {
            .gcall("gtrack_convert_to_indexed_format", trackstr, TRUE, .misha_env())
            success <- TRUE
        },
        error = function(e) {
            # Clean up temporary files on error
            if (file.exists(paste0(dat_path, ".tmp"))) {
                unlink(paste0(dat_path, ".tmp"))
            }
            if (file.exists(paste0(idx_path, ".tmp"))) {
                unlink(paste0(idx_path, ".tmp"))
            }
            stop(sprintf("Failed to convert track %s: %s", trackstr, e$message), call. = FALSE)
        }
    )

    invisible(0)
}

#' Convert 1D interval set to indexed format
#'
#' Converts a per-chromosome interval set to indexed format
#' (intervals.dat + intervals.idx) which reduces file descriptor usage.
#'
#' @param set.name name of interval set to convert
#' @param remove.old if TRUE, removes old per-chromosome files after successful conversion
#' @param force if TRUE, re-converts even if already in indexed format
#' @return invisible NULL
#' @details
#' The indexed format stores all chromosomes in a single intervals.dat file
#' with an intervals.idx index file. This reduces file descriptor usage from
#' N files (one per chromosome) to just 2 files.
#'
#' The conversion process:
#' \enumerate{
#'   \item Creates temporary intervals.dat.tmp and intervals.idx.tmp files
#'   \item Concatenates all per-chromosome files into intervals.dat.tmp
#'   \item Builds index with offsets and checksums
#'   \item Atomically renames temporary files to final names
#'   \item Optionally removes old per-chromosome files
#' }
#'
#' The indexed format is 100% backward compatible with all existing misha functions.
#'
#' @examples
#' \dontrun{
#' # Convert an interval set
#' gintervals.convert_to_indexed("my_intervals")
#'
#' # Convert and remove old files
#' gintervals.convert_to_indexed("my_intervals", remove.old = TRUE)
#'
#' # Force re-conversion
#' gintervals.convert_to_indexed("my_intervals", force = TRUE)
#' }
#' @seealso \code{\link{gintervals.save}}, \code{\link{gintervals.load}}
#' @export
gintervals.convert_to_indexed <- function(set.name = NULL, remove.old = FALSE, force = FALSE) {
    if (is.null(set.name) || !is.character(set.name) || length(set.name) != 1) {
        stop("Usage: gintervals.convert_to_indexed(set.name, remove.old = FALSE, force = FALSE)", call. = FALSE)
    }
    .gcheckroot()

    # Get interval set path - mimic C++ interv2path logic
    path <- gsub("\\.", "/", set.name)
    intervset_path <- paste0(get("GWD", envir = .misha), "/", path, ".interv")

    # Check if it's a Big Set (directory) or single-file format
    if (!file.exists(intervset_path)) {
        stop(sprintf("Interval set %s does not exist", set.name), call. = FALSE)
    }

    is_bigset <- dir.exists(intervset_path)

    if (!is_bigset) {
        message(sprintf("Interval set %s is in single-file format and does not need conversion.", set.name))
        return(invisible(NULL))
    }

    # Check if already converted (check index file instead of directory for robustness)
    idx_path <- file.path(intervset_path, "intervals.idx")
    dat_path <- file.path(intervset_path, "intervals.dat")

    if (file.exists(idx_path) && !force) {
        message(sprintf("Interval set %s is already in indexed format. Use force=TRUE to re-convert.", set.name))
        return(invisible(NULL))
    }

    # Call C++ function to perform the conversion
    tryCatch(
        {
            .gcall("ginterv_convert", set.name, remove.old, .misha_env())
        },
        error = function(e) {
            # Clean up temporary files on error
            tmp_dat <- paste0(dat_path, ".tmp")
            tmp_idx <- paste0(idx_path, ".tmp")
            if (file.exists(tmp_dat)) unlink(tmp_dat)
            if (file.exists(tmp_idx)) unlink(tmp_idx)
            stop(sprintf("Failed to convert interval set %s: %s", set.name, e$message), call. = FALSE)
        }
    )

    invisible(NULL)
}

#' Convert 2D interval set to indexed format
#'
#' Converts a per-chromosome interval set to indexed format
#' (intervals2d.dat + intervals2d.idx) which reduces file descriptor usage.
#'
#' @param set.name name of 2D interval set to convert
#' @param remove.old if TRUE, removes old per-chromosome files after successful conversion
#' @param force if TRUE, re-converts even if already in indexed format
#' @return invisible NULL
#' @details
#' The indexed format stores all chromosome pairs in a single intervals2d.dat file
#' with an intervals2d.idx index file. This dramatically reduces file descriptor
#' usage, especially for genomes with many chromosomes (N*(N-1)/2 files to just 2).
#'
#' Only non-empty pairs are stored in the index, avoiding O(N²) space overhead.
#'
#' The conversion process:
#' \enumerate{
#'   \item Scans directory for existing per-pair files
#'   \item Creates temporary intervals2d.dat.tmp and intervals2d.idx.tmp files
#'   \item Concatenates all per-pair files into intervals2d.dat.tmp
#'   \item Builds index with pair offsets and checksums
#'   \item Atomically renames temporary files to final names
#'   \item Optionally removes old per-pair files
#' }
#'
#' The indexed format is 100% backward compatible with all existing misha functions.
#'
#' @examples
#' \dontrun{
#' # Convert a 2D interval set
#' gintervals.2d.to_indexed_format("my_2d_intervals")
#'
#' # Convert and remove old files
#' gintervals.2d.to_indexed_format("my_2d_intervals", remove.old = TRUE)
#'
#' # Force re-conversion
#' gintervals.2d.to_indexed_format("my_2d_intervals", force = TRUE)
#' }
#'
#' @export
gintervals.2d.to_indexed_format <- function(set.name = NULL, remove.old = FALSE, force = FALSE) {
    if (is.null(set.name) || !is.character(set.name) || length(set.name) != 1) {
        stop("Usage: gintervals.2d.to_indexed_format(set.name, remove.old = FALSE, force = FALSE)", call. = FALSE)
    }
    .gcheckroot()

    # Get interval set path - mimic C++ interv2path logic
    path <- gsub("\\.", "/", set.name)
    intervset_path <- paste0(get("GWD", envir = .misha), "/", path, ".interv")

    # Check if it's a Big Set (directory) or single-file format
    if (!file.exists(intervset_path)) {
        stop(sprintf("2D interval set %s does not exist", set.name), call. = FALSE)
    }

    is_bigset <- dir.exists(intervset_path)

    if (!is_bigset) {
        message(sprintf("2D interval set %s is in single-file format and does not need conversion.", set.name))
        return(invisible(NULL))
    }

    # Check if already converted (check index file instead of directory for robustness)
    idx_path <- file.path(intervset_path, "intervals2d.idx")
    dat_path <- file.path(intervset_path, "intervals2d.dat")

    if (file.exists(idx_path) && !force) {
        message(sprintf("2D interval set %s is already in indexed format. Use force=TRUE to re-convert.", set.name))
        return(invisible(NULL))
    }

    # Call C++ function to perform the conversion
    tryCatch(
        {
            .gcall("ginterv2d_convert", set.name, remove.old, .misha_env())
        },
        error = function(e) {
            # Clean up temporary files on error
            tmp_dat <- paste0(dat_path, ".tmp")
            tmp_idx <- paste0(idx_path, ".tmp")
            if (file.exists(tmp_dat)) unlink(tmp_dat)
            if (file.exists(tmp_idx)) unlink(tmp_idx)
            stop(sprintf("Failed to convert 2D interval set %s: %s", set.name, e$message), call. = FALSE)
        }
    )

    invisible(NULL)
}
