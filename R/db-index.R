# Internal helper function to check if the current database is in indexed format
# @return Logical. TRUE if the database is indexed, FALSE otherwise
.gdb.is_indexed <- function() {
    groot <- get("GROOT", envir = .misha)
    if (is.null(groot) || groot == "") {
        return(FALSE)
    }

    seq_dir <- file.path(groot, "seq")
    if (!dir.exists(seq_dir)) {
        return(FALSE)
    }

    index_path <- file.path(seq_dir, "genome.idx")
    genome_seq_path <- file.path(seq_dir, "genome.seq")

    return(file.exists(index_path) && file.exists(genome_seq_path))
}

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
#' @param chunk_size Integer. The size of the chunk to read from the sequence files. Default: 104857600 (100MB). Reduce if
#' you are running into memory issues.
#' @param threads Integer or NULL. Number of parallel processes used when converting
#' tracks and interval sets (each worker handles one track/interval set via
#' \code{parallel::mclapply}). If \code{NULL} (default), uses
#' \code{min(parallel::detectCores(), 8)}. Set to \code{1} for serial execution.
#' Falls back to serial on non-Unix platforms (mclapply requires fork).
#'
#' @return Invisible NULL
#'
#' @details
#' This function converts a per-chromosome database (with separate .seq files per contig) to
#' indexed format (single genome.seq + genome.idx). The indexed format
#' provides better performance and scalability, especially for genomes with many contigs.
#'
#' The converted database keeps the chromosome names and order that \code{gsetroot()} gives it,
#' whether or not it is the loaded database, so its indexed tracks and interval sets, iteration
#' order and interval IDs stay the same. In a per-chromosome database whose chrom_sizes.txt names
#' lack the "chr" prefix of its .seq files, that order is the sorted chromosome names, not the
#' chrom_sizes.txt order, and chrom_sizes.txt is rewritten in it.
#' The conversion stops, leaving the database as it was, unless the indexed sequence holds exactly
#' the chromosomes of chrom_sizes.txt: a chromosome name with a character other than a letter,
#' digit, '_', '-' or '.' (which the index would replace with '_'), or a .seq file whose length
#' differs from chrom_sizes.txt, stops it, and so does a chromosome named 'genome' (its genome.seq
#' would be overwritten). A database whose seq/ or chrom_sizes.txt is in another database, such as a
#' dataset saved with \code{copy_seq = FALSE} or a \code{gdb.create_linked()} database, is not
#' converted: convert that database. A seq/ or chrom_sizes.txt linked to storage of its own elsewhere
#' is converted there. The session loaded before the conversion is put back as it was after it,
#' except the converted database if it was loaded as a dataset and the conversion changed it: it
#' is unloaded, with a warning. If
#' the conversion changed the sequence it reads (its seq/ is the converted database's, with a
#' chrom_sizes.txt of its own), a warning says so, and nothing is loaded when its chromosomes would
#' read other chromosomes' sequence. chrom_sizes.txt keeps its group when the converting user
#' belongs to it; otherwise it takes that user's group.
#'
#' The conversion process:
#' \enumerate{
#'   \item Checks if database is already in indexed format
#'   \item Gets the chromosome names and order that \code{gsetroot()} gives the database
#'   \item Consolidates all per-chromosome .seq files into genome.seq
#'   \item Creates genome.idx with CRC64 checksum
#'   \item Optionally validates the conversion
#'   \item Optionally removes old .seq files
#'   \item If convert_tracks=TRUE, converts all eligible tracks (1D: dense, sparse, array; 2D: rectangles, points)
#'   \item If convert_intervals=TRUE, converts all eligible interval sets (1D and 2D)
#' }
#'
#' Tracks and intervals that cannot be converted (and are skipped):
#' \itemize{
#'   \item Tracks: virtual tracks, single-file tracks, already converted tracks
#'   \item Intervals: Single-file interval sets, already converted interval sets
#' }
#'
#' @examples
#' \dontrun{
#' # Convert the loaded database with its tracks and interval sets
#' gsetroot("/path/to/database")
#' gdb.convert_to_indexed(
#'     convert_tracks = TRUE,
#'     convert_intervals = TRUE,
#'     remove_old_files = TRUE,
#'     verbose = TRUE
#' )
#'
#' # Convert current database to indexed format (genome only)
#' gdb.convert_to_indexed()
#'
#' # Convert specific database without loading it first
#' gdb.convert_to_indexed(groot = "/path/to/database")
#'
#' # Convert genome and all tracks to indexed format
#' gdb.convert_to_indexed(convert_tracks = TRUE)
#'
#' # Full conversion with validation and cleanup
#' gsetroot("/path/to/database")
#' gdb.convert_to_indexed(
#'     convert_tracks = TRUE,
#'     convert_intervals = TRUE,
#'     remove_old_files = TRUE,
#'     validate = TRUE,
#'     verbose = TRUE
#' )
#' }
#'
#' @seealso \code{\link{gdb.create}}, \code{\link{gdb.init}}, \code{\link{gtrack.convert_to_indexed}}, \code{\link{gintervals.convert_to_indexed}}, \code{\link{gintervals.2d.convert_to_indexed}}
#' @export
gdb.convert_to_indexed <- function(groot = NULL, remove_old_files = FALSE, force = FALSE, validate = TRUE, convert_tracks = FALSE, convert_intervals = FALSE, verbose = FALSE, chunk_size = 104857600, threads = NULL) {
    # Resolve thread count: NULL -> min(detectCores, 8). Non-Unix -> 1 (mclapply forks).
    threads <- .gdb.convert_to_indexed.resolve_threads(threads)

    # Validate database and get setup information
    setup_info <- .gdb.convert_to_indexed.validate_and_setup(groot, verbose)

    # Return early if already indexed
    if (setup_info$already_indexed && !force && !convert_tracks && !convert_intervals) {
        return(invisible(NULL))
    }

    # Get user confirmation
    if (!setup_info$already_indexed &&
        !.gdb.convert_to_indexed.get_confirmation(setup_info$groot, setup_info$chrom_sizes, remove_old_files, force)) {
        return(invisible(NULL))
    }

    # The steps below load the database being converted (to validate it, and to list its tracks and
    # interval sets). The session loaded before is put back as it was at the end, its datasets and
    # working directory included, with the directory caches dropped, as tracks may have changed
    # format under the same paths. The converted database, if it was a loaded dataset, is left out
    # once the conversion changed it (its tracks may be keyed by its own chrom ids now). If the
    # conversion changed the sequence the session reads (its seq/ is the converted one, with a
    # chrom_sizes.txt of its own), it is unloaded when its chrom ids no longer match the index, and
    # a warning says so; one that still reads right but that gsetroot() would not load now gets a
    # warning too.
    old_session <- as.list(.misha, all.names = TRUE)
    old_groot <- old_session$GROOT
    datasets <- as.character(old_session$GDATASETS)
    converted <- datasets[normalizePath(datasets, mustWork = FALSE) == normalizePath(setup_info$groot)]
    # what the conversion got to, for a converted loaded dataset: the genome step finished, the
    # track or interval step rewrote something (files keyed by chrom id), or chrom_sizes.txt was
    # replaced and not put back; and whether all of it finished
    changed <- FALSE
    finished <- FALSE
    if (length(converted)) {
        chrom_sizes_md5 <- unname(tools::md5sum(file.path(setup_info$groot, "chrom_sizes.txt")))
    }
    if (!is.null(old_groot) && nzchar(old_groot)) {
        on.exit(
            {
                rm(list = ls(.misha, all.names = TRUE), envir = .misha)
                list2env(old_session, envir = .misha)
                .gdb.clear_all_dir_caches()
                unload <- length(converted) &&
                    (changed || !identical(unname(tools::md5sum(file.path(setup_info$groot, "chrom_sizes.txt"))), chrom_sizes_md5))
                if (unload) {
                    gdataset.unload(converted[1])
                }
                genome <- get("ALLGENOME", envir = .misha)[[1]]
                problem <- tryCatch(
                    .gdb.check_genome_idx(old_groot, as.character(genome$chrom), genome$end),
                    warning = function(w) NULL,
                    error = conditionMessage
                )
                reads_own <- is.null(problem)
                if (!reads_own) {
                    gdb.unload()
                } else {
                    # as gsetroot() checks it; and its format as the files now have it (the loaded
                    # database may be the converted one)
                    problem <- tryCatch(
                        {
                            chromsizes <- utils::read.csv(file.path(old_groot, "chrom_sizes.txt"),
                                sep = "\t", header = FALSE, col.names = c("chrom", "size"), colClasses = c("character", "numeric")
                            )
                            chrom_order <- .gdb.chrom_order(old_groot, chromsizes)
                            assign("DB_IS_PER_CHROMOSOME", chrom_order$per_chromosome, envir = .misha)
                            .gdb.check_genome_idx(old_groot, chrom_order$names[chrom_order$id_order], chromsizes$size[chrom_order$id_order])
                        },
                        warning = function(w) NULL,
                        error = conditionMessage
                    )
                }
                # warnings last, once the session is in place: a handler that stops at one leaves it so
                if (unload) {
                    warning(sprintf(
                        if (finished) {
                            "%s was converted, so it is no longer loaded as a dataset."
                        } else {
                            "The conversion of %s failed after changing it, so it is no longer loaded as a dataset."
                        },
                        converted[1]
                    ), call. = FALSE)
                }
                if (!reads_own) {
                    warning(sprintf(
                        "%s, loaded before the conversion, no longer reads its own sequence: %s No database is loaded.",
                        old_groot, problem
                    ), call. = FALSE)
                } else if (!is.null(problem)) {
                    warning(sprintf(
                        "%s is loaded as before the conversion, but gsetroot() would not load it now: %s",
                        old_groot, problem
                    ), call. = FALSE)
                }
            },
            add = TRUE
        )
    }

    if (!setup_info$already_indexed) {
        # Convert genome sequences
        .gdb.convert_to_indexed.genome(setup_info, validate, remove_old_files, verbose = verbose, chunk_size = chunk_size)
        changed <- TRUE
    }


    # Convert tracks if requested
    if (convert_tracks && .gdb.convert_to_indexed.tracks(setup_info$groot, verbose, threads = threads) > 0) {
        changed <- TRUE
    }

    # Convert intervals if requested
    if (convert_intervals && .gdb.convert_to_indexed.intervals(setup_info$groot, remove_old_files, verbose, threads = threads) > 0) {
        changed <- TRUE
    }

    finished <- TRUE
    if (verbose) message("\n=== Conversion Complete ===")

    invisible(NULL)
}

# Resolve thread count for bulk conversion. NULL -> capped detectCores.
# Non-Unix forces serial since parallel::mclapply uses fork.
.gdb.convert_to_indexed.resolve_threads <- function(threads) {
    if (.Platform$OS.type != "unix") {
        if (!is.null(threads) && is.numeric(threads) && threads > 1L) {
            warning("parallel::mclapply requires fork (Unix); falling back to serial.", call. = FALSE)
        }
        return(1L)
    }
    if (is.null(threads)) {
        n <- tryCatch(parallel::detectCores(), error = function(e) 1L)
        if (is.na(n) || n < 1L) n <- 1L
        return(as.integer(min(n, 8L)))
    }
    if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) || threads < 1L) {
        stop("threads must be NULL or a positive integer", call. = FALSE)
    }
    as.integer(threads)
}

# Apply fn to xs serially or via parallel::mclapply, capturing per-element
# failures as ("error" + message) list entries rather than aborting. Returns
# a list of per-element result records of the form
#   list(item = <xs[[i]]>, status = "ok"|"error", error = NULL|<chr>, value = <result>).
.gdb.convert_to_indexed.parallel_apply <- function(xs, fn, threads) {
    wrap <- function(x) {
        tryCatch(
            {
                value <- fn(x)
                list(item = x, status = "ok", error = NULL, value = value)
            },
            error = function(e) {
                list(item = x, status = "error", error = conditionMessage(e), value = NULL)
            }
        )
    }

    if (threads <= 1L || length(xs) <= 1L) {
        return(lapply(xs, wrap))
    }

    n_workers <- min(threads, length(xs))
    results <- parallel::mclapply(xs, wrap,
        mc.cores = n_workers, mc.preschedule = FALSE
    )

    # mclapply may surface fork-level failures as "try-error" instead of our wrap()
    # record. Normalize so callers always see the same shape.
    for (i in seq_along(results)) {
        r <- results[[i]]
        if (inherits(r, "try-error") || !is.list(r) || is.null(r$status)) {
            msg <- if (inherits(r, "try-error")) {
                cond <- attr(r, "condition")
                if (!is.null(cond)) conditionMessage(cond) else as.character(r)
            } else {
                "worker returned an unexpected value"
            }
            results[[i]] <- list(
                item = xs[[i]], status = "error",
                error = msg, value = NULL
            )
        }
    }
    results
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

    # Read chromosome information as gsetroot() reads it
    chrom_sizes_path <- file.path(groot, "chrom_sizes.txt")
    if (!file.exists(chrom_sizes_path)) {
        stop(sprintf("chrom_sizes.txt not found: %s", chrom_sizes_path), call. = FALSE)
    }

    # The conversion writes seq/ and rewrites chrom_sizes.txt. In a dataset saved without copy_seq,
    # or a database made by gdb.create_linked(), they are links into another database (a directory
    # with tracks/, seq/ or chrom_sizes.txt besides them), which the conversion would change for
    # every database that shares it. A seq/ or chrom_sizes.txt linked to storage of its own
    # elsewhere is converted there.
    groot_real <- normalizePath(groot, mustWork = TRUE)
    for (path in c(seq_dir, chrom_sizes_path)) {
        real <- normalizePath(path, mustWork = TRUE)
        owner <- dirname(real)
        if (owner != groot_real && any(file.exists(setdiff(file.path(owner, c("tracks", "seq", "chrom_sizes.txt")), real)))) {
            stop(sprintf(
                "%s is in the database %s: convert that database instead.",
                path, owner
            ), call. = FALSE)
        }
    }

    chrom_sizes <- utils::read.csv(
        chrom_sizes_path,
        sep = "\t", header = FALSE, col.names = c("chrom", "size"), colClasses = c("character", "numeric")
    )

    # The converted database keeps the chromosome names and chrom ids that gsetroot() gives it
    # (.gdb.chrom_order), whether or not it is the loaded database, so that its indexed tracks
    # and interval sets, which are keyed by these ids, stay valid. In a per-chromosome database
    # the ids follow the sorted names, not the chrom_sizes.txt order.
    chrom_order <- .gdb.chrom_order(groot, chrom_sizes)
    chrom_sizes$chrom <- chrom_order$names
    chrom_sizes <- chrom_sizes[chrom_order$id_order, ]
    rownames(chrom_sizes) <- NULL

    # Check that per-chromosome .seq files exist
    # Handle chr prefix mismatch between chromosome names and .seq files
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

    # The indexed format's sequence file is seq/genome.seq, so a chromosome whose own file has that
    # name (in any case: macOS file systems ignore it) would lose its sequence to the conversion
    genome_chrom <- chrom_sizes$chrom[tolower(basename(seq_files)) == "genome.seq"]
    if (length(genome_chrom)) {
        stop(sprintf(
            "The sequence file of chromosome %s is seq/genome.seq, the name of the indexed format's sequence file; rename the chromosome before converting.",
            genome_chrom[1]
        ), call. = FALSE)
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
.gdb.convert_to_indexed.genome <- function(setup_info, validate = TRUE, remove_old_files = FALSE, verbose = FALSE, chunk_size = 104857600) {
    groot <- setup_info$groot
    chrom_sizes <- setup_info$chrom_sizes
    seq_files <- setup_info$seq_files
    index_path <- setup_info$index_path
    genome_seq_path <- setup_info$genome_seq_path
    chrom_sizes_path <- setup_info$chrom_sizes_path

    # genome.idx is written last, under a temporary name, and renamed into place only after
    # chrom_sizes.txt holds the new order: a database with genome.idx counts as converted (and the
    # sequence readers use the index), so a conversion stopped before that point leaves a
    # consistent per-chromosome database that the next run converts again.
    index_path_tmp <- paste0(index_path, ".tmp")
    # chrom_sizes.txt is replaced through a temporary file next to it, and the original is kept, to
    # put back if the conversion fails after replacing it: a hard link, or where that is refused (a
    # file of another user under fs.protected_hardlinks), a copy with its group, mode and bytes. Both
    # get names no other file has. A new file gets the original's group when the converting user is
    # in it (else chgrp fails and it keeps the user's group), then its mode.
    chrom_sizes_target <- normalizePath(chrom_sizes_path, mustWork = TRUE)
    chrom_sizes_tmp <- tempfile("chrom_sizes.txt.", tmpdir = dirname(chrom_sizes_target))
    chrom_sizes_orig <- tempfile("chrom_sizes.txt.", tmpdir = dirname(chrom_sizes_target))
    original <- file.info(chrom_sizes_target, extra_cols = TRUE)
    chrom_sizes_replaced <- FALSE
    # The temporary file goes on exit, an interrupt included. The backup goes once the conversion is
    # done (genome.idx in place) or if chrom_sizes.txt was not replaced; it stays when putting the
    # original back failed (the error names it) and when stopped between the two renames.
    on.exit(
        {
            unlink(chrom_sizes_tmp)
            if (!chrom_sizes_replaced || file.exists(index_path)) unlink(chrom_sizes_orig)
        },
        add = TRUE
    )
    if (!suppressWarnings(file.link(chrom_sizes_target, chrom_sizes_orig))) {
        copied <- file.copy(chrom_sizes_target, chrom_sizes_orig)
        if (copied && !is.na(original$gid)) {
            suppressWarnings(system2("chgrp", c(original$gid, shQuote(chrom_sizes_orig)), stdout = FALSE, stderr = FALSE))
        }
        if (!copied ||
            !Sys.chmod(chrom_sizes_orig, original$mode, use_umask = FALSE) ||
            !identical(
                readBin(chrom_sizes_orig, "raw", file.size(chrom_sizes_orig)),
                readBin(chrom_sizes_target, "raw", file.size(chrom_sizes_target))
            )) {
            stop(sprintf("Failed to keep a copy of %s", chrom_sizes_target), call. = FALSE)
        }
    }
    # R reports a failed write or close (a full disk, say) as a warning; here it stops the conversion
    as_error <- function(expr) {
        withCallingHandlers(expr, warning = function(w) stop(conditionMessage(w), call. = FALSE))
    }

    if (verbose) message("Converting database to indexed format...")

    # Create temporary FASTA file from .seq files
    temp_fasta <- tempfile(fileext = ".fasta")
    on.exit(unlink(temp_fasta), add = TRUE)
    # closed here if writing it fails midway (it is closed and set back to NULL when complete)
    fasta_con <- NULL
    on.exit(if (!is.null(fasta_con)) close(fasta_con), add = TRUE)

    tryCatch(
        {
            # Files of an interrupted conversion: a genome.idx left next to the genome.seq written
            # below would be read with it, so it goes first (an index alone counts as indexed)
            unlink(c(index_path, index_path_tmp))
            unlink(genome_seq_path)

            # Write multi-FASTA file
            if (verbose) message("Creating temporary multi-FASTA file...")
            fasta_con <- file(temp_fasta, "wb") # Use binary mode for faster writes

            for (i in seq_len(nrow(chrom_sizes))) {
                chrom <- chrom_sizes$chrom[i]
                seq_file <- seq_files[i]

                # Write header
                as_error(writeBin(charToRaw(sprintf(">%s\n", chrom)), fasta_con))

                # Copy sequence file directly - no line breaking needed
                # The C++ parser handles long lines just fine
                seq_con <- file(seq_file, "rb")
                tryCatch(
                    repeat {
                        chunk <- readBin(seq_con, "raw", n = chunk_size)
                        if (length(chunk) == 0) break
                        as_error(writeBin(chunk, fasta_con))
                    },
                    finally = close(seq_con)
                )

                # Write newline after sequence
                as_error(writeBin(charToRaw("\n"), fasta_con))

                if ((i %% 10) == 0 || i == nrow(chrom_sizes)) {
                    if (verbose) message(sprintf("  Processed %d/%d chromosomes", i, nrow(chrom_sizes)))
                }
            }

            con <- fasta_con
            fasta_con <- NULL
            as_error(close(con))

            # Call C++ import function
            # Use sort=FALSE to keep the chrom id order set up by validate_and_setup
            if (verbose) message("Creating indexed format...")
            contig_info <- .gcall(
                "gseq_multifasta_import",
                temp_fasta,
                genome_seq_path,
                index_path_tmp,
                FALSE, # sort=FALSE: keep the chrom id order
                .misha_env()
            )

            # The index has to hold the chromosomes of chrom_sizes.txt, names, sizes and chrom id
            # order, before anything is replaced or removed. A short temporary FASTA (a full disk),
            # a .seq file of another length, or a name the import changes (it replaces characters
            # other than letters, digits, '_', '-' and '.') would renumber or rename chromosomes.
            expected <- sprintf("%s (%.0f bp)", chrom_sizes$chrom, as.numeric(chrom_sizes$size))
            imported <- sprintf("%s (%.0f bp)", contig_info$name, as.numeric(contig_info$size))
            if (!identical(imported, expected)) {
                k <- seq_len(max(length(expected), length(imported)))
                i <- which(is.na(expected[k]) | is.na(imported[k]) | expected[k] != imported[k])[1]
                if (identical(chrom_sizes$chrom[i], contig_info$name[i]) &&
                    file.size(seq_files[i]) != as.numeric(chrom_sizes$size[i])) {
                    stop(sprintf(
                        "seq/%s has %.0f bytes and chrom_sizes.txt says %.0f; they must agree before converting",
                        basename(seq_files[i]), file.size(seq_files[i]), as.numeric(chrom_sizes$size[i])
                    ), call. = FALSE)
                }
                stop(sprintf(
                    "chromosome %d of chrom_sizes.txt, %s, was imported as %s",
                    i, if (is.na(expected[i])) "none" else expected[i], if (is.na(imported[i])) "nothing" else imported[i]
                ), call. = FALSE)
            }
            seq_bytes <- as.numeric(file.info(genome_seq_path)$size)
            if (!identical(seq_bytes, sum(as.numeric(chrom_sizes$size)))) {
                stop(sprintf(
                    "genome.seq has %.0f bytes, not the %.0f of chrom_sizes.txt",
                    seq_bytes, sum(as.numeric(chrom_sizes$size))
                ), call. = FALSE)
            }

            if (verbose) message("Index created successfully")

            # chrom_sizes.txt in chrom id order, as the index has it
            if (verbose) message("Updating chrom_sizes.txt with canonical chromosome names...")

            updated_chrom_sizes <- data.frame(
                chrom = chrom_sizes$chrom,
                size = as.numeric(chrom_sizes$size),
                stringsAsFactors = FALSE
            )

            .gwith_umask(as_error(write.table(updated_chrom_sizes, chrom_sizes_tmp,
                quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
            )))
            # the original's group and mode, before it replaces it
            if (!is.na(original$gid)) {
                suppressWarnings(system2("chgrp", c(original$gid, shQuote(chrom_sizes_tmp)), stdout = FALSE, stderr = FALSE))
            }
            suppressWarnings(Sys.chmod(chrom_sizes_tmp, original$mode, use_umask = FALSE))
            if (!file.rename(chrom_sizes_tmp, chrom_sizes_target)) {
                stop(sprintf("Failed to replace %s", chrom_sizes_target), call. = FALSE)
            }
            chrom_sizes_replaced <- TRUE
            if (!file.rename(index_path_tmp, index_path)) {
                stop(sprintf("Failed to move %s into place", index_path), call. = FALSE)
            }

            # Validate if requested
            if (validate) {
                if (verbose) message("Validating conversion...")

                # Sample validation: check first 100 bases of each chromosome
                validation_failed <- FALSE

                # The converted database is loaded for the validation (gdb.convert_to_indexed()
                # puts the session loaded before back)
                suppressMessages(gdb.init(groot))

                for (i in seq_len(min(10, nrow(chrom_sizes)))) { # Check first 10 chroms
                    chrom_name <- updated_chrom_sizes$chrom[i]
                    seq_file <- seq_files[i]

                    # Read from old file
                    old_con <- file(seq_file, "rb")
                    old_seq <- toupper(rawToChar(readBin(old_con, "raw", n = 100)))
                    close(old_con)

                    # Extract from indexed format
                    new_seq <- toupper(gseq.extract(gintervals(chrom_name, 0, min(100, updated_chrom_sizes$size[i]))))

                    if (old_seq != new_seq) {
                        if (verbose) {
                            message(sprintf("\nValidation mismatch for %s:", chrom_name))
                            message(sprintf("  Old (first 60): %s", substr(old_seq, 1, 60)))
                            message(sprintf("  New (first 60): %s", substr(new_seq, 1, 60)))
                            message(sprintf("  Seq file: %s", seq_file))
                        }
                        warning(sprintf("Validation failed for chromosome %s", chrom_name))
                        validation_failed <- TRUE
                    } else if (verbose && i <= 3) {
                        message(sprintf("  [OK] %s validated successfully", chrom_name))
                    }
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
            # Clean up partial files on error: genome.idx before genome.seq, since an index alone
            # counts as indexed, and the original chrom_sizes.txt back if it was replaced
            unlink(c(index_path, index_path_tmp))
            unlink(genome_seq_path)
            if (chrom_sizes_replaced && !file.rename(chrom_sizes_orig, chrom_sizes_target)) {
                stop(sprintf(
                    "Conversion failed: %s; and putting back chrom_sizes.txt failed too: the original is %s",
                    conditionMessage(e), chrom_sizes_orig
                ), call. = FALSE)
            }
            stop(sprintf("Conversion failed: %s", conditionMessage(e)), call. = FALSE)
        }
    )
}

# Helper function to convert tracks to indexed format
.gdb.convert_to_indexed.tracks <- function(groot, verbose = FALSE, threads = 1L) {
    if (verbose) message("\n=== Converting Tracks ===")

    # Init the database to get track list (gdb.convert_to_indexed() puts the session loaded before
    # back)
    suppressMessages(gdb.init(groot))

    all_tracks <- gtrack.ls()

    if (length(all_tracks) == 0) {
        if (verbose) message("No tracks found in database")
        return(invisible(0L))
    }

    if (verbose) message(sprintf("Found %d tracks in database", length(all_tracks)))

    # Classification (parallel): gtrack.info per track is independent.
    # Each worker returns the per-track verdict; main process aggregates.
    classify_one <- function(track) {
        info <- tryCatch(gtrack.info(track), error = function(e) NULL)
        if (is.null(info)) {
            return(list(track = track, verdict = "skip", reason = "failed to get info"))
        }

        is_1d <- info$type %in% c("dense", "sparse", "array")
        is_2d <- info$type %in% c("rectangles", "points")
        if (!is_1d && !is_2d) {
            return(list(
                track = track, verdict = "skip",
                reason = sprintf("unsupported type (%s)", info$type)
            ))
        }

        trackstr <- gsub("\\.", "/", track)
        trackdir <- sprintf("%s.track", paste(groot, "tracks", trackstr, sep = "/"))
        if (!dir.exists(trackdir)) {
            return(list(track = track, verdict = "skip", reason = "single-file format"))
        }
        if (file.exists(file.path(trackdir, "track.idx"))) {
            return(list(track = track, verdict = "skip", reason = "already converted"))
        }

        list(track = track, verdict = if (is_1d) "1d" else "2d")
    }

    classify_results <- .gdb.convert_to_indexed.parallel_apply(
        all_tracks, classify_one, threads
    )

    convertible_1d_tracks <- character(0)
    convertible_2d_tracks <- character(0)
    skipped_tracks <- list()

    for (r in classify_results) {
        track <- r$item
        if (r$status == "error") {
            skipped_tracks[[track]] <- sprintf("classification error: %s", r$error)
            next
        }
        v <- r$value
        if (v$verdict == "skip") {
            skipped_tracks[[track]] <- v$reason
        } else if (v$verdict == "1d") {
            convertible_1d_tracks <- c(convertible_1d_tracks, track)
        } else if (v$verdict == "2d") {
            convertible_2d_tracks <- c(convertible_2d_tracks, track)
        }
    }

    total_convertible <- length(convertible_1d_tracks) + length(convertible_2d_tracks)

    if (total_convertible > 0) {
        if (verbose) {
            message(sprintf(
                "  Convertible: %d tracks (%d 1D, %d 2D)",
                total_convertible, length(convertible_1d_tracks), length(convertible_2d_tracks)
            ))
            if (threads > 1L) {
                message(sprintf("  Running with %d parallel workers", threads))
            }
        }
    } else {
        if (verbose) message("  No tracks need conversion")
    }

    if (length(skipped_tracks) > 0) {
        if (verbose) message(sprintf("  Skipped: %d tracks", length(skipped_tracks)))
        for (track in names(skipped_tracks)) {
            if (verbose) message(sprintf("    - %s: %s", track, skipped_tracks[[track]]))
        }
    }

    # Convert tracks in parallel; per-track failures are captured as
    # warnings without aborting the batch. mclapply children fork from
    # this process so they inherit GROOT/GTRACKS/GWD - no need to gdb.init
    # inside each worker.
    convert_1d <- function(track) {
        if (verbose) message(sprintf("  Converting 1D track: %s", track))
        gtrack.convert_to_indexed(track)
        track
    }
    convert_2d <- function(track) {
        if (verbose) message(sprintf("  Converting 2D track: %s", track))
        gtrack.2d.convert_to_indexed(track, remove.old = TRUE)
        track
    }

    results_1d <- .gdb.convert_to_indexed.parallel_apply(
        convertible_1d_tracks, convert_1d, threads
    )
    results_2d <- .gdb.convert_to_indexed.parallel_apply(
        convertible_2d_tracks, convert_2d, threads
    )

    all_results <- c(results_1d, results_2d)
    converted_count <- sum(vapply(all_results, function(r) r$status == "ok", logical(1)))
    failed <- Filter(function(r) r$status == "error", all_results)

    for (r in failed) {
        warning(sprintf("Failed to convert track %s: %s", r$item, r$error),
            call. = FALSE
        )
    }

    if (total_convertible > 0) {
        if (verbose) message(sprintf("Successfully converted %d/%d tracks", converted_count, total_convertible))
        if (length(failed) > 0) {
            warning(sprintf(
                "Failed to convert %d tracks: %s",
                length(failed),
                paste(vapply(failed, function(r) r$item, character(1)), collapse = ", ")
            ), call. = FALSE)
        }
        if (verbose) {
            message(sprintf(
                "Track conversion summary: %d succeeded, %d failed",
                converted_count, length(failed)
            ))
        }
    }
    # the number of tracks rewritten, for gdb.convert_to_indexed()
    invisible(converted_count)
}

# Helper function to convert interval sets to indexed format
.gdb.convert_to_indexed.intervals <- function(groot, remove_old_files = FALSE, verbose = FALSE, threads = 1L) {
    if (verbose) message("\n=== Converting Interval Sets ===")

    # Init the database to get interval list (gdb.convert_to_indexed() puts the session loaded before
    # back)
    suppressMessages(gdb.init(groot))

    all_intervals <- gintervals.ls()

    if (length(all_intervals) == 0) {
        if (verbose) message("No interval sets found in database")
        return(invisible(0L))
    }

    if (verbose) message(sprintf("Found %d interval sets in database", length(all_intervals)))

    # Classify each interval set (parallel-safe: each call only reads
    # filesystem metadata for its own directory).
    classify_one <- function(intervset) {
        path <- gsub("\\.", "/", intervset)
        intervset_path <- paste0(groot, "/tracks/", path, ".interv")

        if (!file.exists(intervset_path)) {
            return(list(verdict = "skip", reason = "does not exist"))
        }
        if (!dir.exists(intervset_path)) {
            return(list(verdict = "skip", reason = "single-file format"))
        }

        idx_path_1d <- file.path(intervset_path, "intervals.idx")
        idx_path_2d <- file.path(intervset_path, "intervals2d.idx")
        pair_files <- list.files(intervset_path, pattern = "-")
        has_2d <- length(pair_files) > 0

        if (has_2d) {
            if (file.exists(idx_path_2d)) {
                list(verdict = "skip", reason = "already converted (2D)")
            } else {
                list(verdict = "2d")
            }
        } else {
            if (file.exists(idx_path_1d)) {
                list(verdict = "skip", reason = "already converted (1D)")
            } else {
                list(verdict = "1d")
            }
        }
    }

    classify_results <- .gdb.convert_to_indexed.parallel_apply(
        all_intervals, classify_one, threads
    )

    convertible_1d <- character(0)
    convertible_2d <- character(0)
    skipped_intervals <- list()

    for (r in classify_results) {
        intervset <- r$item
        if (r$status == "error") {
            skipped_intervals[[intervset]] <- sprintf("classification error: %s", r$error)
            next
        }
        v <- r$value
        if (v$verdict == "skip") {
            skipped_intervals[[intervset]] <- v$reason
        } else if (v$verdict == "1d") {
            convertible_1d <- c(convertible_1d, intervset)
        } else if (v$verdict == "2d") {
            convertible_2d <- c(convertible_2d, intervset)
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
    } else if (verbose && threads > 1L) {
        message(sprintf("  Running with %d parallel workers", threads))
    }

    if (length(skipped_intervals) > 0) {
        if (verbose) message(sprintf("  Skipped: %d interval sets", length(skipped_intervals)))
        for (intervset in names(skipped_intervals)) {
            if (verbose) message(sprintf("    - %s: %s", intervset, skipped_intervals[[intervset]]))
        }
    }

    # Convert in parallel; per-set failures captured as warnings.
    convert_1d <- function(intervset) {
        if (verbose) message(sprintf("  Converting 1D interval set: %s", intervset))
        gintervals.convert_to_indexed(intervset, remove.old = remove_old_files)
        intervset
    }
    convert_2d <- function(intervset) {
        if (verbose) message(sprintf("  Converting 2D interval set: %s", intervset))
        gintervals.2d.convert_to_indexed(intervset, remove.old = remove_old_files)
        intervset
    }

    results_1d <- .gdb.convert_to_indexed.parallel_apply(
        convertible_1d, convert_1d, threads
    )
    results_2d <- .gdb.convert_to_indexed.parallel_apply(
        convertible_2d, convert_2d, threads
    )

    converted_1d_count <- sum(vapply(results_1d, function(r) r$status == "ok", logical(1)))
    converted_2d_count <- sum(vapply(results_2d, function(r) r$status == "ok", logical(1)))
    failed_1d <- Filter(function(r) r$status == "error", results_1d)
    failed_2d <- Filter(function(r) r$status == "error", results_2d)

    for (r in c(failed_1d, failed_2d)) {
        warning(sprintf("Failed to convert interval set %s: %s", r$item, r$error),
            call. = FALSE
        )
    }

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
            all_failed_names <- c(
                vapply(failed_1d, function(r) r$item, character(1)),
                vapply(failed_2d, function(r) r$item, character(1))
            )
            warning(sprintf(
                "Failed to convert %d interval sets: %s",
                length(all_failed_names),
                paste(all_failed_names, collapse = ", ")
            ), call. = FALSE)
        }
    }
    # the number of interval sets rewritten, for gdb.convert_to_indexed()
    invisible(total_converted)
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

    # the session's chrom ids are not those of a dataset loaded with an order of its own
    if (any(.gtrack_db_path(trackstr) %in% names(get0("GDATASET_CHROMS", envir = .misha, ifnotfound = NULL)))) {
        stop(sprintf("Track %s is in a dataset that numbers the chromosomes differently from the working database, so it cannot be written in the indexed format in this session; load that database with gsetroot() to do it there.", trackstr), call. = FALSE)
    }

    trackdir <- .track_dir(trackstr)
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

    # Dispatch based on track type
    if (track_type %in% c("rectangles", "points")) {
        # 2D tracks: delegate to gtrack.2d.convert_to_indexed
        gtrack.2d.convert_to_indexed(trackstr, remove.old = TRUE)
        return(invisible(0))
    }

    if (!track_type %in% c("dense", "sparse", "array")) {
        stop(sprintf("Cannot convert track %s: unsupported type '%s'", trackstr, track_type), call. = FALSE)
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

    # Track layout flipped per-chrom -> indexed; any prior cache entry
    # (nullptr from the per-chrom era) is now wrong.
    .gdb.invalidate_dir_cache(trackdir)

    invisible(0)
}

.gtrack.split_indexed_to_per_chrom <- function(track_dir, chrom_names, remove_indexed = TRUE) {
    if (!dir.exists(track_dir)) {
        stop(sprintf("Track directory does not exist: %s", track_dir), call. = FALSE)
    }
    .gcall(
        "gtrack_split_indexed_to_per_chrom",
        track_dir, as.character(chrom_names), isTRUE(remove_indexed)
    )
    # Track layout flipped indexed -> per-chrom; the cached TrackIndex is
    # now invalid.
    .gdb.invalidate_dir_cache(track_dir)
    invisible()
}

.gtrack.pack_per_chrom_to_indexed <- function(track_dir, chrom_names, track_type) {
    if (!dir.exists(track_dir)) {
        stop(sprintf("Track directory does not exist: %s", track_dir), call. = FALSE)
    }
    .gcall(
        "gtrack_pack_per_chrom_to_indexed",
        track_dir, as.character(chrom_names), as.character(track_type)
    )
    .gdb.invalidate_dir_cache(track_dir)
    invisible()
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

    # the session's chrom ids are not those of a dataset loaded with an order of its own
    if (any(.gintervals_db_path(set.name) %in% names(get0("GDATASET_CHROMS", envir = .misha, ifnotfound = NULL)))) {
        stop(sprintf("Interval set %s is in a dataset that numbers the chromosomes differently from the working database, so it cannot be written in the indexed format in this session; load that database with gsetroot() to do it there.", set.name), call. = FALSE)
    }

    # Get interval set path using database-aware resolution
    intervset_path <- .intervals_dir(set.name)

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

    # Layout flipped per-chrom -> indexed; drop any cached IntervalsIndex1D
    # entry for this dir.
    .gdb.invalidate_dir_cache(intervset_path)

    invisible(NULL)
}

#' Convert 2D track to indexed format
#'
#' Converts a per-chromosome-pair 2D track (rectangles or points) to indexed format
#' (track.dat + track.idx). This reduces file descriptor usage from O(N^2) to O(1),
#' which is especially beneficial for genomes with many contigs.
#'
#' @param track track name to convert
#' @param remove.old Logical. If TRUE, removes old per-chromosome-pair files after
#' successful conversion. Default: FALSE.
#' @param force Logical. If TRUE, re-converts even if already in indexed format.
#' Default: FALSE.
#' @return None.
#' @seealso \code{\link{gtrack.2d.create}}, \code{\link{gtrack.2d.import}},
#' \code{\link{gtrack.2d.import_contacts}}, \code{\link{gtrack.convert_to_indexed}},
#' \code{\link{gdb.convert_to_indexed}}
#' @examples
#' \dontrun{
#' # Convert a 2D track to indexed format
#' gtrack.2d.convert_to_indexed("my_2d_track")
#'
#' # Convert and remove old per-pair files
#' gtrack.2d.convert_to_indexed("my_2d_track", remove.old = TRUE)
#'
#' # Force re-conversion
#' gtrack.2d.convert_to_indexed("my_2d_track", force = TRUE)
#' }
#' @export gtrack.2d.convert_to_indexed
gtrack.2d.convert_to_indexed <- function(track = NULL, remove.old = FALSE, force = FALSE) {
    if (is.null(substitute(track))) {
        stop("Usage: gtrack.2d.convert_to_indexed(track, remove.old = FALSE, force = FALSE)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
    if (is.na(match(trackstr, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Track %s does not exist", trackstr), call. = FALSE)
    }

    # the session's chrom ids are not those of a dataset loaded with an order of its own
    if (any(.gtrack_db_path(trackstr) %in% names(get0("GDATASET_CHROMS", envir = .misha, ifnotfound = NULL)))) {
        stop(sprintf("Track %s is in a dataset that numbers the chromosomes differently from the working database, so it cannot be written in the indexed format in this session; load that database with gsetroot() to do it there.", trackstr), call. = FALSE)
    }

    trackdir <- .track_dir(trackstr)
    idx_path <- file.path(trackdir, "track.idx")
    dat_path <- file.path(trackdir, "track.dat")

    # Check if already converted
    if (file.exists(idx_path) && !force) {
        message(sprintf("Track %s is already in indexed format.", trackstr))
        return(invisible(0))
    }

    # Get track info to determine type
    info <- gtrack.info(track)
    track_type <- info$type

    # Only 2D tracks can be converted with this function
    if (!track_type %in% c("rectangles", "points")) {
        stop(sprintf("Cannot convert track %s: only 2D tracks (rectangles, points) can be converted with this function", trackstr), call. = FALSE)
    }

    # Call C++ function to perform the conversion
    tryCatch(
        {
            .gcall("gtrack2d_convert_to_indexed", trackstr, remove.old, .misha_env())
        },
        error = function(e) {
            # Clean up temporary files on error
            if (file.exists(paste0(dat_path, ".tmp"))) {
                unlink(paste0(dat_path, ".tmp"))
            }
            if (file.exists(paste0(idx_path, ".tmp"))) {
                unlink(paste0(idx_path, ".tmp"))
            }
            stop(sprintf("Failed to convert 2D track %s: %s", trackstr, e$message), call. = FALSE)
        }
    )

    # 2D layout flipped per-pair -> indexed; drop any cached TrackIndex2D
    # for this dir. (The C++ already calls TrackIndex2D::clear_cache() as
    # a blanket reset; the targeted call here is harmless and keeps the
    # invariant uniform across paths.)
    .gdb.invalidate_dir_cache(trackdir)

    invisible(0)
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
#' Only non-empty pairs are stored in the index, avoiding O(N^2) space overhead.
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
#' gintervals.2d.convert_to_indexed("my_2d_intervals")
#'
#' # Convert and remove old files
#' gintervals.2d.convert_to_indexed("my_2d_intervals", remove.old = TRUE)
#'
#' # Force re-conversion
#' gintervals.2d.convert_to_indexed("my_2d_intervals", force = TRUE)
#' }
#'
#' @export
gintervals.2d.convert_to_indexed <- function(set.name = NULL, remove.old = FALSE, force = FALSE) {
    if (is.null(set.name) || !is.character(set.name) || length(set.name) != 1) {
        stop("Usage: gintervals.2d.convert_to_indexed(set.name, remove.old = FALSE, force = FALSE)", call. = FALSE)
    }
    .gcheckroot()

    # the session's chrom ids are not those of a dataset loaded with an order of its own
    if (any(.gintervals_db_path(set.name) %in% names(get0("GDATASET_CHROMS", envir = .misha, ifnotfound = NULL)))) {
        stop(sprintf("Interval set %s is in a dataset that numbers the chromosomes differently from the working database, so it cannot be written in the indexed format in this session; load that database with gsetroot() to do it there.", set.name), call. = FALSE)
    }

    # Get interval set path using database-aware resolution
    intervset_path <- .intervals_dir(set.name)

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

    # 2D bigset layout flipped per-pair -> indexed; drop the cached
    # IntervalsIndex2D entry for this dir.
    .gdb.invalidate_dir_cache(intervset_path)

    invisible(NULL)
}
