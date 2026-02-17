# Sequence functions

.gseq.import <- function(groot = NULL, path = NULL) {
    chroms <- c()
    files <- c()
    tmp.dirname <- tempfile(pattern = "", tmpdir = paste(groot, "/downloads", sep = ""))
    if (!dir.create(tmp.dirname, recursive = TRUE, mode = "0777")) {
        stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = FALSE)
    }

    tryCatch(
        {
            files <- c()

            ftp_files <- grep("^ftp://", path, perl = TRUE, value = TRUE)
            other_files <- setdiff(path, ftp_files)
            if (length(other_files) + length(ftp_files) != length(path)) {
                stop("Some paths are not supported", call. = FALSE)
            }
            if (length(ftp_files) > 0) {
                files <- c(files, gwget(ftp_files, tmp.dirname))
            }
            if (any(!file.exists(other_files))) {
                stop("Some files do not exist", call. = FALSE)
            }
            files <- c(files, other_files)

            message("Building Seq files...")
            # Accept FASTA files with any naming convention, not just chr* prefix
            # Look for common FASTA extensions
            fasta_pattern <- "\\.(fa|fasta|fna|seq)(\\.gz)?$"
            fastas <- files[grep(fasta_pattern, basename(files), perl = TRUE, ignore.case = TRUE)]

            for (fasta in fastas) {
                # Extract chromosome name from filename
                # Remove common FASTA extensions and .gz suffix
                # Use the filename (without extension) as the chromosome name directly
                base_name <- basename(fasta)
                chrom <- gsub("\\.(fa|fasta|fna|seq)(\\.gz)?$", "", base_name, perl = TRUE, ignore.case = TRUE)

                if (!is.na(match(chrom, chroms))) {
                    next
                }

                fasta.original <- fasta

                if (length(grep("^.+\\.gz$", fasta, perl = TRUE))) {
                    .chrom <- basename(gsub("^(.+)\\.gz$", "\\1", fasta, perl = TRUE))
                    fasta.unzipped <- paste(tmp.dirname, "/", .chrom, sep = "")
                    cmd <- paste("/bin/sh -c \"gunzip -q -c", fasta, ">", fasta.unzipped, "\"")
                    if (system(cmd)) {
                        stop(sprintf("Command failed: %s", cmd), call. = FALSE)
                    }
                    fasta <- fasta.unzipped
                }

                # Output without forcing chr prefix - use chromosome name as-is
                seq <- sprintf("%s/seq/%s.seq", groot, chrom)

                message(sprintf("%s", chrom))
                .gcall("gseqimport", fasta, seq, .misha_env())

                chroms <- c(chroms, chrom)
            }
        },
        finally = {
            unlink(tmp.dirname, recursive = TRUE)
        }
    )
    chroms
}

# Internal function: Import multi-FASTA file to indexed genome format
# Creates genome.seq and genome.idx files
.gseq.import_multifasta <- function(groot = NULL, fasta = NULL, verbose = TRUE) {
    if (is.null(groot) || is.null(fasta)) {
        stop("Usage: .gseq.import_multifasta(groot, fasta, verbose = TRUE)", call. = FALSE)
    }

    if (length(fasta) != 1) {
        stop("Multi-FASTA import requires exactly one input file", call. = FALSE)
    }

    if (!file.exists(fasta)) {
        stop(sprintf("FASTA file %s does not exist", fasta), call. = FALSE)
    }

    # Handle gzipped files
    fasta.original <- fasta
    tmp.dirname <- NULL
    cleanup_files <- c()

    tryCatch(
        {
            if (grepl("\\.gz$", fasta, perl = TRUE)) {
                tmp.dirname <- tempfile(pattern = "", tmpdir = paste(groot, "/downloads", sep = ""))
                if (!dir.create(tmp.dirname, recursive = TRUE, mode = "0777")) {
                    stop(sprintf("Failed to create temp directory %s", tmp.dirname), call. = FALSE)
                }

                fasta.unzipped <- file.path(tmp.dirname, gsub("\\.gz$", "", basename(fasta)))
                cmd <- sprintf("/bin/sh -c \"gunzip -q -c '%s' > '%s'\"", fasta, fasta.unzipped)
                if (system(cmd) != 0) {
                    stop(sprintf("Failed to decompress file: %s", cmd), call. = FALSE)
                }
                fasta <- fasta.unzipped
                cleanup_files <- c(cleanup_files, fasta.unzipped)
            }

            # Output paths
            seq_path <- file.path(groot, "seq", "genome.seq")
            index_path <- file.path(groot, "seq", "genome.idx")

            # Call C++ import function
            # It will: parse FASTA, sanitize headers, sort by name (if sort=TRUE),
            # assign chromids 0..N-1 in sorted order, write genome.seq and genome.idx,
            # and return a data frame with contig names and sizes
            if (verbose) message("Importing multi-FASTA file...")
            contig_info <- .gcall(
                "gseq_multifasta_import",
                fasta,
                seq_path,
                index_path,
                TRUE, # sort=TRUE (default): sort chromosomes alphabetically
                .misha_env()
            )

            if (nrow(contig_info) == 0) {
                stop("No contigs were imported from FASTA file", call. = FALSE)
            }

            if (verbose) message(sprintf("Successfully imported %d contigs", nrow(contig_info)))

            # Return contig info data frame with columns: name, size
            contig_info
        },
        finally = {
            # Cleanup temp files
            for (f in cleanup_files) {
                if (file.exists(f)) {
                    unlink(f)
                }
            }
            if (!is.null(tmp.dirname) && dir.exists(tmp.dirname)) {
                unlink(tmp.dirname, recursive = TRUE)
            }
        }
    )
}


#' Returns DNA sequences
#'
#' Returns DNA sequences for given intervals
#'
#' This function returns an array of sequence strings for each interval from
#' 'intervals'. If intervals contain an additional 'strand' column and its
#' value is '-1', the reverse-complementary sequence is returned.
#'
#' @param intervals intervals for which DNA sequence is returned
#' @return An array of character strings representing DNA sequence.
#' @seealso \code{\link{gextract}}
#' @keywords ~extract ~DNA ~sequence
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' intervs <- gintervals(c(1, 2), 10000, 10020)
#' gseq.extract(intervs)
#'
#' @export gseq.extract
gseq.extract <- function(intervals = NULL) {
    if (is.null(intervals)) {
        stop("Usage: gseq.extract(intervals)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    res <- .gcall("gseqread", intervals, .misha_env())
    res
}

#' Export a database genome as FASTA
#'
#' Writes all contigs from a misha database to a multi-FASTA file.
#'
#' By default, the currently active database is used. You can also provide
#' \code{groot} to export another database without changing the caller's active
#' database.
#'
#' @param file Output FASTA file path
#' @param groot Optional database root path. If NULL, uses current database.
#' @param line_width Number of bases per FASTA line. Default: 80.
#' @param chunk_size Number of bases to extract per chunk while writing.
#'   Default: 1000000.
#' @param overwrite Logical. If TRUE, overwrite existing output file.
#'   Default: FALSE.
#' @param verbose Logical. If TRUE, prints progress messages. Default: FALSE.
#'
#' @return Invisibly returns \code{file}.
#' @seealso \code{\link{gdb.init}}, \code{\link{gseq.extract}}
#' @examples
#' \dontrun{
#' gdb.init_examples()
#' out <- tempfile(fileext = ".fa")
#' gdb.export_fasta(out)
#' head(readLines(out))
#' }
#' @export
gdb.export_fasta <- function(file = NULL,
                             groot = NULL,
                             line_width = 80L,
                             chunk_size = 1000000L,
                             overwrite = FALSE,
                             verbose = FALSE) {
    if (is.null(file) || !is.character(file) || length(file) != 1 || nchar(file) == 0) {
        stop("Usage: gdb.export_fasta(file, groot = NULL, line_width = 80, chunk_size = 1000000, overwrite = FALSE, verbose = FALSE)", call. = FALSE)
    }

    line_width <- suppressWarnings(as.integer(line_width))
    if (is.na(line_width) || line_width < 1) {
        stop("line_width must be a positive integer", call. = FALSE)
    }

    chunk_size <- suppressWarnings(as.integer(chunk_size))
    if (is.na(chunk_size) || chunk_size < 1) {
        stop("chunk_size must be a positive integer", call. = FALSE)
    }

    out_dir <- dirname(file)
    if (!dir.exists(out_dir)) {
        stop(sprintf("Output directory does not exist: %s", out_dir), call. = FALSE)
    }

    if (file.exists(file) && !overwrite) {
        stop(sprintf("Output file already exists: %s. Use overwrite = TRUE to replace it.", file), call. = FALSE)
    }

    old_groot <- NULL
    if (exists("GROOT", envir = .misha, inherits = FALSE)) {
        old_groot <- get("GROOT", envir = .misha)
    }

    switched_root <- FALSE
    if (is.null(groot)) {
        .gcheckroot()
    } else {
        if (!dir.exists(groot)) {
            stop(sprintf("Database directory does not exist: %s", groot), call. = FALSE)
        }
        groot <- normalizePath(groot, mustWork = TRUE)

        needs_switch <- is.null(old_groot) || old_groot == ""
        if (!needs_switch) {
            needs_switch <- normalizePath(old_groot, mustWork = TRUE) != groot
        }

        if (needs_switch) {
            suppressMessages(gdb.init(groot))
            switched_root <- TRUE
        }
    }

    temp_file <- tempfile(pattern = "misha_export_", fileext = ".fasta", tmpdir = out_dir)
    con <- NULL

    on.exit(
        {
            if (!is.null(con) && isOpen(con)) {
                close(con)
            }
            if (file.exists(temp_file)) {
                unlink(temp_file)
            }
            if (switched_root && !is.null(old_groot) && old_groot != "") {
                suppressMessages(gdb.init(old_groot))
            }
        },
        add = TRUE
    )

    all_genome <- gintervals.all()
    if (is.null(all_genome) || nrow(all_genome) == 0) {
        stop("No chromosomes found in the database", call. = FALSE)
    }

    chroms <- as.character(all_genome$chrom)
    starts <- as.numeric(all_genome$start)
    ends <- as.numeric(all_genome$end)

    con <- file(temp_file, open = "wt")

    for (i in seq_along(chroms)) {
        chrom <- chroms[i]
        start_pos <- starts[i]
        end_pos <- ends[i]

        if (verbose) {
            message(sprintf("Exporting %s (%d bp)", chrom, as.integer(end_pos - start_pos)))
        }

        writeLines(sprintf(">%s", chrom), con)

        line_buffer <- ""
        while (start_pos < end_pos) {
            next_end <- min(start_pos + chunk_size, end_pos)
            seq_chunk <- gseq.extract(gintervals(chrom, start_pos, next_end))

            if (length(seq_chunk) != 1 || is.na(seq_chunk[1])) {
                stop(sprintf("Failed reading sequence for chromosome %s [%d,%d)", chrom, as.integer(start_pos), as.integer(next_end)), call. = FALSE)
            }

            chunk <- paste0(line_buffer, seq_chunk[1])
            chunk_len <- nchar(chunk, type = "chars")

            if (chunk_len >= line_width) {
                line_starts <- seq.int(1L, chunk_len - line_width + 1L, by = line_width)
                line_ends <- line_starts + line_width - 1L
                writeLines(substring(chunk, line_starts, line_ends), con)

                next_pos <- max(line_ends) + 1L
                if (next_pos <= chunk_len) {
                    line_buffer <- substr(chunk, next_pos, chunk_len)
                } else {
                    line_buffer <- ""
                }
            } else {
                line_buffer <- chunk
            }

            start_pos <- next_end
        }

        if (nchar(line_buffer, type = "chars") > 0) {
            writeLines(line_buffer, con)
        }
    }

    close(con)
    con <- NULL

    if (file.exists(file)) {
        unlink(file)
    }
    if (!file.rename(temp_file, file)) {
        stop(sprintf("Failed to move temporary FASTA to output path: %s", file), call. = FALSE)
    }

    invisible(file)
}


#' Computes auto-correlation between the strands for a file of mapped sequences
#'
#' Calculates auto-correlation between plus and minus strands for the given
#' chromosome in a file of mapped sequences.
#'
#' This function calculates auto-correlation between plus and minus strands for
#' the given chromosome in a file of mapped sequences. Each line in the file
#' describes one read. Each column is separated by a TAB character.
#'
#' The following columns must be presented in the file: sequence, chromosome,
#' coordinate and strand. The position of these columns are controlled by
#' 'cols.order' argument accordingly. The default value of 'cols.order' is a
#' vector (9,11,13,14) meaning that sequence is expected to be found at column
#' number 9, chromosome - at column 11, coordinate - at column 13 and strand -
#' at column 14. The first column should be referenced by 1 and not by 0.
#'
#' Coordinates that are not in [min.coord, max.coord] range are ignored.
#'
#' gcompute_strands_autocorr outputs the total statistics and the
#' auto-correlation given by bins. The size of the bin is indicated by
#' 'binsize' parameter. Statistics is calculated for bins in the range of
#' [-maxread, maxread].
#'
#' @param file the name of the file containing mapped sequences
#' @param chrom chromosome for which the auto-correlation is computed
#' @param binsize calculate the auto-correlation for bins in the range of
#' [-maxread, maxread]
#' @param maxread maximal length of the sequence used for statistics
#' @param cols.order order of sequence, chromosome, coordinate and strand
#' columns in file
#' @param min.coord minimal coordinate used for statistics
#' @param max.coord maximal coordinate used for statistics
#' @return Statistics for each strand and auto-correlation by given bins.
#' @keywords ~gcompute_strands_autocorr ~auto-correlation ~autocorrelation
#' ~correlation
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gcompute_strands_autocorr(paste(.misha$GROOT, "reads", sep = "/"),
#'     "chr1", 50,
#'     maxread = 300
#' )
#'
#' @export gcompute_strands_autocorr
gcompute_strands_autocorr <- function(file = NULL, chrom = NULL, binsize = NULL, maxread = 400, cols.order = c(9, 11, 13, 14), min.coord = 0, max.coord = 3e+8) {
    if (is.null(file) || is.null(chrom) || is.null(binsize)) {
        stop("Usage: gcompute_strands_autocorr(file, chrom, binsize, maxread = 400, cols.order = c(9, 11, 13, 14), min.coord = 0, max.coord = 3e+8)", call. = FALSE)
    }
    .gcheckroot()

    # Normalize chromosome name using aliases (handles per-chromosome databases)
    chrom <- as.character(.gchroms(as.character(chrom)))

    res <- .gcall("C_gcompute_strands_autocorr", file, chrom, binsize, maxread, cols.order, min.coord, max.coord, .misha_env())
    res
}
