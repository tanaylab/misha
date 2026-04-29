# Importers from common interval file formats (BED, GFF/GTF, VCF) into
# misha 1D intervals data frames.
#
# misha uses 0-based, half-open [start, end) coordinates, like BED.
# GFF/GTF/VCF are 1-based; the importers below subtract 1 from start.

# Internal: read a tabular file, drop common header lines, return a
# stringsAsFactors=FALSE data frame. Used by BED/GFF/VCF importers.
.gread_table_filtered <- function(file, header_pat) {
    bed <- NULL
    if (requireNamespace("data.table", quietly = TRUE)) {
        tryCatch(
            {
                can_grep <- nzchar(Sys.which("grep"))
                if (can_grep) {
                    cmd <- sprintf("grep -vE %s %s", shQuote(header_pat), shQuote(file))
                    bed <- data.table::fread(
                        cmd = cmd,
                        header = FALSE, sep = "\t", quote = "",
                        fill = TRUE, stringsAsFactors = FALSE,
                        data.table = FALSE, showProgress = FALSE
                    )
                } else {
                    bed <- data.table::fread(
                        file,
                        header = FALSE, sep = "\t", quote = "",
                        fill = TRUE, stringsAsFactors = FALSE,
                        data.table = FALSE, showProgress = FALSE
                    )
                    if (!is.null(bed) && nrow(bed) > 0) {
                        v1 <- trimws(bed[[1]])
                        keep <- !grepl(header_pat, v1)
                        bed <- bed[keep, , drop = FALSE]
                    }
                }
                if (is.null(bed) || nrow(bed) == 0) bed <- NULL
            },
            error = function(e) bed <<- NULL
        )
    }
    if (is.null(bed)) {
        bed <- utils::read.table(
            file,
            header = FALSE, sep = "\t", quote = "", comment.char = "",
            fill = TRUE, stringsAsFactors = FALSE, colClasses = "character"
        )
        if (nrow(bed) > 0) {
            v1 <- trimws(bed[[1]])
            keep <- !grepl(header_pat, v1)
            bed <- bed[keep, , drop = FALSE]
        }
    }
    bed
}

# Internal: full-table BED reader (returns N-col data.frame). Used by both
# gtrack.import (BED -> track) and gintervals.import_bed (BED -> intervals).
.gread_bed_table <- function(file) {
    bed <- .gread_table_filtered(file, "^(track|browser|#|$)")
    if (is.null(bed) || nrow(bed) == 0) {
        stop(sprintf("BED file %s appears to be empty or contains no data intervals", file), call. = FALSE)
    }
    if (ncol(bed) < 3) {
        stop(sprintf("BED file %s appears to be malformed (less than 3 columns)", file), call. = FALSE)
    }
    bed
}

# Internal: sort an intervals data frame by (chromid, start), preserving all
# columns. Used by importers because gintervsort strips columns beyond
# chrom/start/end/strand.
.gsort_intervals_df <- function(df) {
    allchroms <- get("ALLGENOME", envir = .misha)[[1]]$chrom
    chromid <- match(as.character(df$chrom), as.character(allchroms))
    df <- df[order(chromid, df$start, df$end), , drop = FALSE]
    rownames(df) <- NULL
    df
}

# Internal: parse a strand string vector to misha numeric strand (1/-1/0).
# Mirrors GInterval::str2strand in C++.
.gparse_strand_vec <- function(s, context = "input") {
    s <- as.character(s)
    out <- rep(NA_integer_, length(s))
    out[s == "+"] <- 1L
    out[s == "-"] <- -1L
    out[s == "." | s == "*" | s == "" | is.na(s)] <- 0L
    bad <- is.na(out)
    if (any(bad)) {
        i <- which(bad)[1]
        stop(sprintf("Invalid strand value \"%s\" at row %d of %s", s[i], i, context), call. = FALSE)
    }
    as.numeric(out)
}


#' Import intervals from a BED file
#'
#' Reads a BED/BED.gz/BED.zip file and returns a misha 1D intervals data
#' frame. Track/browser/comment header lines are skipped automatically.
#' Chromosome names are normalized through the active database's
#' \code{CHROM_ALIAS} mechanism (so \code{chr1} <-> \code{1} works without
#' explicit configuration).
#'
#' BED is already 0-based half-open, so coordinates are taken as-is.
#'
#' @param file path to a BED file (\code{.bed}, \code{.bed.gz}, or
#'   \code{.bed.zip}).
#' @param name if \code{TRUE} and a 4th column exists, include it as
#'   \code{name}.
#' @param score if \code{TRUE} and a 5th (numeric) column exists, include
#'   it as \code{score}.
#' @param strand if \code{TRUE} and a 6th column exists, include it as
#'   \code{strand} (mapped to \code{1}/\code{-1}/\code{0}).
#' @return A 1D intervals data frame, sorted by chrom and start.
#' @seealso \code{\link{gintervals.import_gff}},
#'   \code{\link{gintervals.import_vcf}}.
#' @keywords ~intervals ~import ~BED
#' @export
gintervals.import_bed <- function(file = NULL, name = TRUE, score = TRUE, strand = TRUE) {
    if (is.null(file)) {
        stop("Usage: gintervals.import_bed(file, name = TRUE, score = TRUE, strand = TRUE)", call. = FALSE)
    }
    .gcheckroot()
    if (!file.exists(file)) {
        stop(sprintf("BED file %s does not exist", file), call. = FALSE)
    }

    bed <- .gread_bed_table(file)
    n <- ncol(bed)

    starts <- suppressWarnings(as.numeric(bed[[2]]))
    ends <- suppressWarnings(as.numeric(bed[[3]]))
    if (any(is.na(starts)) || any(is.na(ends))) {
        stop(sprintf("Non-numeric coordinates detected in BED file %s", file), call. = FALSE)
    }

    df <- data.frame(
        chrom = .gchroms(as.character(bed[[1]])),
        start = starts,
        end = ends,
        stringsAsFactors = FALSE
    )

    if (n >= 6 && strand) {
        df$strand <- .gparse_strand_vec(bed[[6]], context = sprintf("BED file %s", file))
    }
    if (n >= 4 && name) {
        df$name <- as.character(bed[[4]])
    }
    if (n >= 5 && score) {
        df$score <- suppressWarnings(as.numeric(as.character(bed[[5]])))
    }

    .gsort_intervals_df(df)
}


#' Import intervals from a GFF/GTF file
#'
#' Reads a GFF3 or GTF file (optionally gzipped) and returns a misha 1D
#' intervals data frame. GFF/GTF are 1-based and inclusive on both ends;
#' coordinates are converted to 0-based half-open by subtracting 1 from
#' \code{start} and leaving \code{end} as-is.
#'
#' Chromosome names are normalized through the active database's
#' \code{CHROM_ALIAS} mechanism.
#'
#' @param file path to a GFF/GTF file.
#' @param feature optional feature-type filter (column 3 of GFF). Pass a
#'   character vector to keep only those types (e.g. \code{"exon"},
#'   \code{c("gene", "transcript")}).
#' @param strand if \code{TRUE}, include the \code{strand} column.
#' @param attrs if \code{TRUE}, include the raw attributes string as
#'   column \code{attrs}. The attribute string is not parsed.
#' @return A 1D intervals data frame with columns chrom, start, end, and
#'   optionally strand, type, source, score, attrs.
#' @seealso \code{\link{gintervals.import_bed}},
#'   \code{\link{gintervals.import_vcf}}.
#' @keywords ~intervals ~import ~GFF ~GTF
#' @export
gintervals.import_gff <- function(file = NULL, feature = NULL, strand = TRUE, attrs = TRUE) {
    if (is.null(file)) {
        stop("Usage: gintervals.import_gff(file, feature = NULL, strand = TRUE, attrs = TRUE)", call. = FALSE)
    }
    .gcheckroot()
    if (!file.exists(file)) {
        stop(sprintf("GFF file %s does not exist", file), call. = FALSE)
    }

    gff <- .gread_table_filtered(file, "^#")
    if (is.null(gff) || nrow(gff) == 0) {
        stop(sprintf("GFF file %s appears to be empty or contains no records", file), call. = FALSE)
    }
    if (ncol(gff) < 8) {
        stop(sprintf("GFF file %s appears to be malformed (expected at least 8 tab-separated columns, got %d)", file, ncol(gff)), call. = FALSE)
    }

    if (!is.null(feature)) {
        keep <- gff[[3]] %in% feature
        gff <- gff[keep, , drop = FALSE]
        if (nrow(gff) == 0) {
            stop(sprintf(
                "No records of feature type(s) %s found in GFF file %s",
                paste(feature, collapse = ", "), file
            ), call. = FALSE)
        }
    }

    starts1 <- suppressWarnings(as.numeric(gff[[4]]))
    ends1 <- suppressWarnings(as.numeric(gff[[5]]))
    if (any(is.na(starts1)) || any(is.na(ends1))) {
        stop(sprintf("Non-numeric coordinates detected in GFF file %s", file), call. = FALSE)
    }

    df <- data.frame(
        chrom = .gchroms(as.character(gff[[1]])),
        start = starts1 - 1,
        end = ends1,
        stringsAsFactors = FALSE
    )

    if (strand) {
        df$strand <- .gparse_strand_vec(gff[[7]], context = sprintf("GFF file %s", file))
    }
    df$source <- as.character(gff[[2]])
    df$type <- as.character(gff[[3]])
    score_num <- suppressWarnings(as.numeric(as.character(gff[[6]])))
    df$score <- score_num
    if (attrs && ncol(gff) >= 9) {
        df$attrs <- as.character(gff[[9]])
    }

    .gsort_intervals_df(df)
}


#' Import intervals from a VCF file
#'
#' Reads a VCF/VCF.gz file and returns a misha 1D intervals data frame
#' with one row per record. VCF is 1-based; \code{start} is set to
#' \code{POS - 1} and \code{end} is set to \code{POS - 1 + nchar(REF)},
#' yielding a 0-based half-open span covering the reference allele.
#'
#' Chromosome names are normalized through the active database's
#' \code{CHROM_ALIAS} mechanism.
#'
#' Multi-allelic records are kept as a single row; the \code{ALT} column
#' contains the original comma-separated string.
#'
#' @param file path to a VCF/VCF.gz file.
#' @param info if \code{TRUE}, include the raw INFO column as
#'   \code{info}. The string is not parsed.
#' @return A 1D intervals data frame with columns chrom, start, end, and
#'   id, ref, alt, qual, filter, optionally info.
#' @seealso \code{\link{gintervals.import_bed}},
#'   \code{\link{gintervals.import_gff}}.
#' @keywords ~intervals ~import ~VCF
#' @export
gintervals.import_vcf <- function(file = NULL, info = TRUE) {
    if (is.null(file)) {
        stop("Usage: gintervals.import_vcf(file, info = TRUE)", call. = FALSE)
    }
    .gcheckroot()
    if (!file.exists(file)) {
        stop(sprintf("VCF file %s does not exist", file), call. = FALSE)
    }

    vcf <- .gread_table_filtered(file, "^#")
    if (is.null(vcf) || nrow(vcf) == 0) {
        stop(sprintf("VCF file %s appears to be empty or contains no records", file), call. = FALSE)
    }
    if (ncol(vcf) < 5) {
        stop(sprintf("VCF file %s appears to be malformed (expected at least 5 tab-separated columns, got %d)", file, ncol(vcf)), call. = FALSE)
    }

    pos <- suppressWarnings(as.numeric(vcf[[2]]))
    if (any(is.na(pos))) {
        stop(sprintf("Non-numeric POS detected in VCF file %s", file), call. = FALSE)
    }
    ref <- as.character(vcf[[4]])
    ref_len <- nchar(ref)
    if (any(ref_len < 1)) {
        stop(sprintf("Empty REF allele detected in VCF file %s", file), call. = FALSE)
    }

    df <- data.frame(
        chrom = .gchroms(as.character(vcf[[1]])),
        start = pos - 1,
        end = pos - 1 + ref_len,
        stringsAsFactors = FALSE
    )
    df$id <- as.character(vcf[[3]])
    df$ref <- ref
    df$alt <- as.character(vcf[[5]])
    if (ncol(vcf) >= 6) df$qual <- suppressWarnings(as.numeric(as.character(vcf[[6]])))
    if (ncol(vcf) >= 7) df$filter <- as.character(vcf[[7]])
    if (info && ncol(vcf) >= 8) df$info <- as.character(vcf[[8]])

    .gsort_intervals_df(df)
}
