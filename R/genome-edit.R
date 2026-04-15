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


#' Stream-process a FASTA file chromosome by chromosome
#'
#' Reads a FASTA file and calls \code{callback(chrom_name, seq_raw)} for each
#' chromosome, where \code{seq_raw} is the uppercase sequence as a raw vector.
#' Only one chromosome is in memory at a time.
#'
#' @param fasta_path Path to FASTA file
#' @param callback Function taking (chrom_name, seq_raw). Return value ignored.
#' @return Character vector of chromosome names in file order (invisibly).
#' @noRd
.stream_fasta <- function(fasta_path, callback) {
    # Bulk-read all lines (fast), then iterate by chromosome
    lines <- readLines(fasta_path, warn = FALSE)
    header_idx <- which(startsWith(lines, ">"))
    if (length(header_idx) == 0L) {
        stop("No FASTA headers found in: ", fasta_path, call. = FALSE)
    }

    chrom_names <- sub("^>\\s*", "", sub("\\s.*", "", lines[header_idx]))

    for (i in seq_along(header_idx)) {
        seq_start <- header_idx[i] + 1L
        seq_end <- if (i < length(header_idx)) header_idx[i + 1L] - 1L else length(lines)

        if (seq_start > seq_end) {
            seq_raw <- raw(0)
        } else {
            seq_raw <- charToRaw(toupper(paste0(lines[seq_start:seq_end], collapse = "")))
        }
        callback(chrom_names[i], seq_raw)
    }

    invisible(chrom_names)
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
#' gdb.init_examples()
#'
#' # Export the example DB to a reference FASTA
#' ref_fasta <- tempfile(fileext = ".fa")
#' gdb.export_fasta(ref_fasta)
#'
#' # Replace two regions with literal sequences
#' intervals <- data.frame(
#'     chrom = c("chr1", "chr1"),
#'     start = c(100, 200),
#'     end = c(110, 210)
#' )
#' donors <- c("AAAAAAAAAA", "CCCCCCCCCC")
#' out <- tempfile(fileext = ".fa")
#' trackdb <- tempfile()
#' ggenome.implant(intervals, donors,
#'     output = out,
#'     genome_fasta = ref_fasta,
#'     create_trackdb = TRUE,
#'     trackdb_path = trackdb
#' )
#'
#' # Verify the implanted sequences via the new trackdb
#' gdb.init(trackdb)
#' gseq.extract(data.frame(chrom = "chr1", start = 100, end = 110))
#'
#' # Clean up
#' unlink(c(ref_fasta, out, paste0(out, ".fai"), trackdb), recursive = TRUE)
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
        # donor is a misha database root — extract sequences from it.
        # We must restore the original root afterwards, especially when
        # genome_fasta is NULL (which needs the original DB for export).
        old_groot <- NULL
        if (exists("GROOT", envir = .misha, inherits = FALSE)) {
            old_groot <- get("GROOT", envir = .misha)
        }
        if (is.null(old_groot) || old_groot == "") {
            if (is.null(genome_fasta)) {
                stop(
                    "No misha database is initialized and 'genome_fasta' is NULL. ",
                    "Either call gdb.init() first or provide 'genome_fasta'.",
                    call. = FALSE
                )
            }
        }
        suppressMessages(gdb.init(donor))
        donor_seqs <- gseq.extract(intervals)
        # Restore original root immediately
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

    # --- resolve reference FASTA path ---
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

    interval_chroms <- as.character(intervals$chrom)

    # --- C++ fast path: read FASTA, apply perturbations, write output + .fai
    fai_df <- .gcall(
        "C_ggenome_implant",
        genome_fasta, output,
        as.character(intervals$chrom),
        as.integer(intervals$start),
        as.integer(intervals$end),
        toupper(donor_seqs),
        as.integer(line_width)
    )

    # Validate that all interval chromosomes were found in reference
    missing_chroms <- setdiff(unique(interval_chroms), fai_df$name)
    if (length(missing_chroms) > 0L) {
        unlink(output)
        stop(sprintf(
            "Chromosome(s) not found in reference: %s",
            paste(missing_chroms, collapse = ", ")
        ), call. = FALSE)
    }

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
#' gdb.init_examples()
#'
#' # Create a "donor" DB with different sequence (all T's)
#' donor_fasta <- tempfile(fileext = ".fa")
#' cat(">chr1\n", paste(rep("T", 500000), collapse = ""), "\n",
#'     ">chr2\n", paste(rep("T", 300000), collapse = ""), "\n",
#'     file = donor_fasta, sep = ""
#' )
#' donor_db <- tempfile()
#' gdb.create(donor_db, fasta = donor_fasta, verbose = FALSE)
#'
#' # Export the current DB as the target FASTA
#' gdb.init_examples()
#' ref_fasta <- tempfile(fileext = ".fa")
#' gdb.export_fasta(ref_fasta)
#'
#' # Transplant donor sequence into positions 100-200 of chr1
#' intervals <- data.frame(chrom = "chr1", start = 100, end = 200)
#' out <- tempfile(fileext = ".fa")
#' trackdb <- tempfile()
#' ggenome.transplant(intervals,
#'     source_genome = donor_db,
#'     target_genome = ref_fasta,
#'     output = out,
#'     create_trackdb = TRUE,
#'     trackdb_path = trackdb
#' )
#'
#' # Verify: positions 100-200 should now be all T's
#' gdb.init(trackdb)
#' gseq.extract(data.frame(chrom = "chr1", start = 100, end = 200))
#'
#' # Clean up
#' unlink(
#'     c(
#'         donor_fasta, donor_db, ref_fasta, out,
#'         paste0(out, ".fai"), trackdb
#'     ),
#'     recursive = TRUE
#' )
#'
#' @seealso \code{\link{ggenome.implant}}, \code{\link{gdb.export_fasta}},
#'   \code{\link{gseq.extract}}
#' @export
ggenome.transplant <- function(intervals, source_genome, target_genome = NULL, output,
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
