# Genome Editing: Replace intervals in a reference genome with donor sequences

#' Read a FASTA file into a named list of sequences
#'
#' @param fasta_path Path to FASTA file
#' @return Named list with chromosome names as keys and full sequences as values
#' @noRd
.read_fasta <- function(fasta_path) {
    lines <- readLines(fasta_path)
    if (length(lines) == 0L) {
        stop("FASTA file is empty: ", fasta_path, call. = FALSE)
    }

    header_idx <- which(startsWith(lines, ">"))
    if (length(header_idx) == 0L) {
        stop("No FASTA headers found in: ", fasta_path, call. = FALSE)
    }

    # Pre-allocate result list
    chrom_names <- sub("^>\\s*", "", sub("\\s.*", "", lines[header_idx]))
    result <- vector("list", length(chrom_names))
    names(result) <- chrom_names

    # For each header, concatenate all sequence lines until the next header
    for (i in seq_along(header_idx)) {
        start <- header_idx[i] + 1L
        end <- if (i < length(header_idx)) header_idx[i + 1L] - 1L else length(lines)
        if (start > end) {
            result[[i]] <- ""
        } else {
            result[[i]] <- paste0(lines[start:end], collapse = "")
        }
    }

    result
}

#' Write sequences as a FASTA file with fixed-width lines
#'
#' @param sequences Named list of sequences (name -> sequence string)
#' @param output Path to output file
#' @param line_width Bases per line
#' @return Invisibly returns the output path
#' @noRd
.write_fasta <- function(sequences, output, line_width = 80L) {
    con <- file(output, open = "wt")
    on.exit(close(con), add = TRUE)

    for (chrom_name in names(sequences)) {
        writeLines(paste0(">", chrom_name), con)
        seq_str <- sequences[[chrom_name]]
        seq_len <- nchar(seq_str)
        if (seq_len == 0L) next

        line_starts <- seq.int(1L, seq_len, by = line_width)
        line_ends <- pmin(line_starts + line_width - 1L, seq_len)
        writeLines(substring(seq_str, line_starts, line_ends), con)
    }

    invisible(output)
}

#' Write a FASTA index (.fai) file
#'
#' Creates a samtools-compatible .fai index alongside a FASTA file.
#' The .fai format is: name\\tlength\\toffset\\tlinebases\\tlinewidth
#' where linewidth includes the newline character.
#'
#' @param fasta_path Path to the FASTA file to index
#' @param line_width Number of bases per line used when writing
#' @return Invisibly returns the .fai path
#' @noRd
.write_fai <- function(fasta_path, line_width = 80L) {
    fai_path <- paste0(fasta_path, ".fai")
    lines <- readLines(fasta_path)

    header_idx <- which(startsWith(lines, ">"))
    if (length(header_idx) == 0L) {
        stop("No FASTA headers found in: ", fasta_path, call. = FALSE)
    }

    fai_entries <- character(length(header_idx))

    for (i in seq_along(header_idx)) {
        chrom_name <- sub("^>\\s*", "", sub("\\s.*", "", lines[header_idx[i]]))

        seq_start <- header_idx[i] + 1L
        seq_end <- if (i < length(header_idx)) header_idx[i + 1L] - 1L else length(lines)

        if (seq_start > seq_end) {
            seq_length <- 0L
        } else {
            seq_lines <- lines[seq_start:seq_end]
            seq_length <- sum(nchar(seq_lines))
        }

        # Calculate byte offset to the first sequence character
        # Each previous line has its content + newline (\n)
        offset <- 0L
        for (j in seq_len(header_idx[i])) {
            offset <- offset + nchar(lines[j]) + 1L # +1 for newline
        }

        # linebases = actual bases per line, linewidth = bytes per line including newline
        linebases <- line_width
        linewidth <- line_width + 1L # +1 for newline character

        fai_entries[i] <- paste(chrom_name, seq_length, offset, linebases, linewidth, sep = "\t")
    }

    writeLines(fai_entries, fai_path)
    invisible(fai_path)
}


#' Implant donor sequences into a reference genome
#'
#' Replaces specified intervals in a reference genome with donor DNA sequences
#' and writes the result as a new FASTA file. Optionally creates a misha trackdb
#' from the output.
#'
#' This function fills the gap between manual sequence extraction and genome
#' perturbation by providing a single-call interface for editing genomes.
#'
#' The \code{donor} parameter controls where replacement sequences come from:
#' \itemize{
#'   \item If \code{donor} is a character vector with one sequence per interval
#'         row, those literal sequences are used directly.
#'   \item If \code{donor} is a single string pointing to an existing directory,
#'         it is treated as a misha database root. Sequences are extracted from
#'         that database at the same coordinates as \code{intervals}.
#' }
#'
#' Perturbations are applied in reverse coordinate order within each chromosome
#' so that earlier coordinates remain valid when later ones are replaced.
#'
#' @param intervals A data.frame with \code{chrom}, \code{start}, \code{end}
#'   columns specifying the regions to replace. Coordinates are 0-based,
#'   half-open (standard misha convention).
#' @param donor Either a character vector of DNA sequences (one per interval
#'   row), or a single string path to a misha database root from which
#'   sequences will be extracted at the same intervals.
#' @param output Path for the output FASTA file.
#' @param genome_fasta Path to the reference FASTA file to edit. If \code{NULL},
#'   the current misha database is exported via \code{gdb.export_fasta}.
#' @param create_trackdb Logical. If \code{TRUE}, creates a misha trackdb from
#'   the output FASTA using \code{gdb.create}.
#' @param trackdb_path Path for the new trackdb. Defaults to
#'   \code{<dirname(output)>/trackdb}.
#' @param line_width Integer. Number of bases per FASTA line. Default: 80.
#' @param overwrite Logical. If \code{TRUE}, overwrite existing output file.
#'   Default: \code{FALSE}.
#'
#' @return Invisibly returns the output FASTA path.
#'
#' @examples
#' \dontrun{
#' gdb.init("/path/to/my/db")
#'
#' # Replace two regions with literal sequences
#' intervals <- data.frame(
#'     chrom = c("chr1", "chr1"),
#'     start = c(100, 500),
#'     end = c(110, 510)
#' )
#' donors <- c("AAAAAAAAAA", "CCCCCCCCCC")
#' ggenome.implant(intervals, donors, output = "/tmp/edited.fa")
#'
#' # Replace regions with sequences from another misha database
#' ggenome.implant(intervals,
#'     donor = "/path/to/donor/db",
#'     output = "/tmp/transplanted.fa"
#' )
#' }
#'
#' @seealso \code{\link{ggenome.transplant}}, \code{\link{gdb.export_fasta}},
#'   \code{\link{gseq.extract}}, \code{\link{gdb.create}}
#' @export
ggenome.implant <- function(intervals, donor, output, genome_fasta = NULL,
                            create_trackdb = TRUE, trackdb_path = NULL,
                            line_width = 80L, overwrite = FALSE) {
    # --- argument validation ---
    if (missing(intervals) || missing(donor) || missing(output)) {
        stop("Usage: ggenome.implant(intervals, donor, output, genome_fasta = NULL, create_trackdb = TRUE, trackdb_path = NULL, line_width = 80, overwrite = FALSE)",
            call. = FALSE
        )
    }

    if (!is.data.frame(intervals)) {
        stop("'intervals' must be a data.frame with chrom, start, end columns", call. = FALSE)
    }
    required_cols <- c("chrom", "start", "end")
    missing_cols <- setdiff(required_cols, names(intervals))
    if (length(missing_cols) > 0L) {
        stop("'intervals' is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    if (nrow(intervals) == 0L) {
        stop("'intervals' must have at least one row", call. = FALSE)
    }

    line_width <- suppressWarnings(as.integer(line_width))
    if (is.na(line_width) || line_width < 1L) {
        stop("'line_width' must be a positive integer", call. = FALSE)
    }

    out_dir <- dirname(output)
    if (!dir.exists(out_dir)) {
        stop(sprintf("Output directory does not exist: %s", out_dir), call. = FALSE)
    }

    if (file.exists(output) && !overwrite) {
        stop(sprintf("Output file already exists: %s. Use overwrite = TRUE to replace it.", output), call. = FALSE)
    }

    # --- resolve donor sequences ---
    donor_is_db <- is.character(donor) && length(donor) == 1L && dir.exists(donor)

    if (donor_is_db) {
        # donor is a misha database root — extract sequences from it
        old_groot <- NULL
        if (exists("GROOT", envir = .misha, inherits = FALSE)) {
            old_groot <- get("GROOT", envir = .misha)
        }
        suppressMessages(gdb.init(donor))
        withr::defer({
            if (!is.null(old_groot) && old_groot != "") {
                suppressMessages(gdb.init(old_groot))
            }
        })
        donor_seqs <- gseq.extract(intervals)
        # Restore root immediately since we need to possibly export reference
        if (!is.null(old_groot) && old_groot != "") {
            suppressMessages(gdb.init(old_groot))
        }
    } else if (is.character(donor)) {
        if (length(donor) != nrow(intervals)) {
            stop(sprintf(
                "Length of 'donor' (%d) must equal number of rows in 'intervals' (%d)",
                length(donor), nrow(intervals)
            ), call. = FALSE)
        }
        donor_seqs <- donor
    } else {
        stop("'donor' must be a character vector of sequences or a path to a misha database", call. = FALSE)
    }

    # --- validate donor sequence lengths match interval widths ---
    interval_widths <- intervals$end - intervals$start
    donor_widths <- nchar(donor_seqs)
    mismatches <- which(donor_widths != interval_widths)
    if (length(mismatches) > 0L) {
        first_mm <- mismatches[1L]
        stop(sprintf(
            "Donor sequence length (%d) does not match interval width (%d) at row %d (%s:%d-%d). Length-changing perturbations are not supported.",
            donor_widths[first_mm], interval_widths[first_mm], first_mm,
            intervals$chrom[first_mm], intervals$start[first_mm], intervals$end[first_mm]
        ), call. = FALSE)
    }

    # --- load reference FASTA ---
    if (is.null(genome_fasta)) {
        .gcheckroot()
        genome_fasta <- tempfile(fileext = ".fa")
        withr::defer(unlink(genome_fasta))
        gdb.export_fasta(genome_fasta)
    } else {
        if (!file.exists(genome_fasta)) {
            stop(sprintf("Reference FASTA file does not exist: %s", genome_fasta), call. = FALSE)
        }
    }

    sequences <- .read_fasta(genome_fasta)

    # --- validate intervals against reference ---
    chroms_in_ref <- names(sequences)
    interval_chroms <- as.character(intervals$chrom)
    missing_chroms <- setdiff(unique(interval_chroms), chroms_in_ref)
    if (length(missing_chroms) > 0L) {
        stop(sprintf(
            "Chromosome(s) not found in reference: %s",
            paste(missing_chroms, collapse = ", ")
        ), call. = FALSE)
    }

    for (i in seq_len(nrow(intervals))) {
        chrom <- interval_chroms[i]
        chrom_len <- nchar(sequences[[chrom]])
        if (intervals$start[i] < 0 || intervals$end[i] > chrom_len) {
            stop(sprintf(
                "Interval out of bounds at row %d: %s:%d-%d (chromosome length: %d)",
                i, chrom, intervals$start[i], intervals$end[i], chrom_len
            ), call. = FALSE)
        }
        if (intervals$start[i] >= intervals$end[i]) {
            stop(sprintf(
                "Invalid interval at row %d: start (%d) must be less than end (%d)",
                i, intervals$start[i], intervals$end[i]
            ), call. = FALSE)
        }
    }

    # --- apply perturbations ---
    # Group by chromosome, sort descending by start within each group
    # so that replacing later positions first preserves earlier coordinates
    order_idx <- order(match(interval_chroms, chroms_in_ref), -intervals$start)

    for (i in order_idx) {
        chrom <- interval_chroms[i]
        # misha coordinates are 0-based; R substr is 1-based
        r_start <- intervals$start[i] + 1L
        r_end <- intervals$end[i]
        seq_str <- sequences[[chrom]]

        prefix <- if (r_start > 1L) substr(seq_str, 1L, r_start - 1L) else ""
        suffix <- if (r_end < nchar(seq_str)) substr(seq_str, r_end + 1L, nchar(seq_str)) else ""
        sequences[[chrom]] <- paste0(prefix, donor_seqs[i], suffix)
    }

    # --- write output ---
    .write_fasta(sequences, output, line_width = line_width)
    .write_fai(output, line_width = line_width)

    # --- create trackdb ---
    if (create_trackdb) {
        if (is.null(trackdb_path)) {
            trackdb_path <- file.path(dirname(output), "trackdb")
        }
        if (dir.exists(trackdb_path) && !overwrite) {
            stop(sprintf("Trackdb directory already exists: %s. Use overwrite = TRUE to replace it.", trackdb_path), call. = FALSE)
        }
        if (dir.exists(trackdb_path)) {
            unlink(trackdb_path, recursive = TRUE)
        }
        suppressMessages(gdb.create(groot = trackdb_path, fasta = output, verbose = FALSE))
    }

    invisible(output)
}


#' Transplant sequences from one genome into another
#'
#' Sugar function that extracts sequences from a source genome at the given
#' intervals and implants them into a target genome. Equivalent to calling
#' \code{ggenome.implant} with \code{donor = source_genome} and
#' \code{genome_fasta = target_genome}.
#'
#' @param intervals A data.frame with \code{chrom}, \code{start}, \code{end}
#'   columns specifying the regions to transplant.
#' @param source_genome Path to the misha database root containing the donor
#'   sequences.
#' @param target_genome Path to the target reference FASTA file. If
#'   \code{NULL}, the current misha database is used.
#' @param output Path for the output FASTA file.
#' @param create_trackdb Logical. If \code{TRUE}, creates a misha trackdb.
#' @param trackdb_path Path for the new trackdb.
#' @param line_width Integer. Number of bases per FASTA line.
#' @param overwrite Logical. If \code{TRUE}, overwrite existing output.
#'
#' @return Invisibly returns the output FASTA path.
#'
#' @examples
#' \dontrun{
#' intervals <- data.frame(
#'     chrom = c("chr1", "chr2"),
#'     start = c(100, 200),
#'     end = c(200, 300)
#' )
#' ggenome.transplant(intervals,
#'     source_genome = "/path/to/donor/db",
#'     target_genome = "/path/to/target.fa",
#'     output = "/tmp/transplanted.fa"
#' )
#' }
#'
#' @seealso \code{\link{ggenome.implant}}, \code{\link{gdb.export_fasta}},
#'   \code{\link{gseq.extract}}
#' @export
ggenome.transplant <- function(intervals, source_genome, target_genome, output,
                               create_trackdb = TRUE, trackdb_path = NULL,
                               line_width = 80L, overwrite = FALSE) {
    if (missing(intervals) || missing(source_genome) || missing(output)) {
        stop("Usage: ggenome.transplant(intervals, source_genome, target_genome, output, create_trackdb = TRUE, trackdb_path = NULL, line_width = 80, overwrite = FALSE)",
            call. = FALSE
        )
    }

    ggenome.implant(
        intervals = intervals,
        donor = source_genome,
        output = output,
        genome_fasta = target_genome,
        create_trackdb = create_trackdb,
        trackdb_path = trackdb_path,
        line_width = line_width,
        overwrite = overwrite
    )
}
