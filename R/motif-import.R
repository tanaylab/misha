# Internal helpers for motif import functions

#' De-duplicate motif IDs by appending .1, .2, etc.
#' @noRd
.dedup_motif_ids <- function(ids) {
    if (length(ids) == 0) {
        return(ids)
    }
    dupes <- duplicated(ids)
    if (!any(dupes)) {
        return(ids)
    }
    warning("Duplicate motif IDs found; appending numeric suffixes to disambiguate")
    tab <- table(ids)
    dup_names <- names(tab[tab > 1])
    for (nm in dup_names) {
        idx <- which(ids == nm)
        ids[idx] <- paste0(nm, ".", seq_along(idx))
    }
    ids
}

#' Re-normalize matrix rows that do not sum to 1.0
#' @noRd
.renormalize_rows <- function(mat, tol = 1e-4, context = "") {
    rsums <- rowSums(mat)
    bad <- which(abs(rsums - 1.0) > tol)
    if (length(bad) > 0) {
        warning(
            "Row sums deviate from 1.0 in ", context, "; re-normalizing ",
            length(bad), " row(s)"
        )
        for (i in bad) {
            if (rsums[i] == 0) {
                mat[i, ] <- 0.25
            } else {
                mat[i, ] <- mat[i, ] / rsums[i]
            }
        }
    }
    mat
}


#' Read motifs from a MEME minimal motif format file
#'
#' Parses a MEME minimal motif format file and returns a named list of
#' position probability matrices (PPM). Each matrix has rows corresponding
#' to motif positions and columns \code{A}, \code{C}, \code{G}, \code{T}.
#' The returned matrices are directly usable with \code{\link{gseq.pwm}}.
#'
#' @param file character(1) path to a MEME format file (\code{.meme}, \code{.txt}).
#'
#' @return A named list of numeric matrices. Each matrix has columns
#'   \code{A, C, G, T} and one row per motif position. List names are motif
#'   identifiers. Each matrix carries the following attributes:
#'   \describe{
#'     \item{name}{Motif name / alternate ID (second token on the MOTIF line)}
#'     \item{alength}{Alphabet length (integer, typically 4)}
#'     \item{w}{Motif width (integer, number of positions)}
#'     \item{nsites}{Number of sites used to build the matrix (numeric; \code{NA} if absent)}
#'     \item{E}{E-value (numeric; \code{NA} if absent)}
#'     \item{url}{URL string if present, otherwise \code{NA}}
#'     \item{strand}{Strand specification from the file header (e.g. \code{"+ -"})}
#'     \item{background}{Named numeric vector of background frequencies (\code{c(A=..., C=..., G=..., T=...)}), or \code{NULL} if absent}
#'   }
#'
#' @examples
#' \dontrun{
#' motifs <- gseq.read_meme("JASPAR2024_CORE_vertebrates.meme")
#' names(motifs)
#' m <- motifs[[1]]
#' head(m)
#' attr(m, "name")
#' attr(m, "nsites")
#' }
#'
#' @family motif functions
#' @export
gseq.read_meme <- function(file) {
    if (!is.character(file) || length(file) != 1) {
        stop("'file' must be a single file path string")
    }
    if (!file.exists(file)) {
        stop("File not found: ", file)
    }

    lines <- readLines(file, warn = FALSE)
    # Strip Windows line endings
    lines <- sub("\r$", "", lines)

    # Parse global header fields
    strand <- NA_character_
    background <- NULL

    # Parse strands line
    strand_idx <- grep("^strands:", lines, ignore.case = TRUE)
    if (length(strand_idx) > 0) {
        strand <- trimws(sub("^strands:\\s*", "", lines[strand_idx[1]], ignore.case = TRUE))
    }

    # Parse alphabet line -- validate DNA only
    alpha_idx <- grep("^ALPHABET", lines, ignore.case = TRUE)
    if (length(alpha_idx) > 0) {
        alpha_line <- lines[alpha_idx[1]]
        if (grepl("=", alpha_line)) {
            alpha_val <- trimws(sub("^ALPHABET\\s*=\\s*", "", alpha_line, ignore.case = TRUE))
            if (nchar(alpha_val) > 0 && !grepl("^[ACGTacgt]+$", alpha_val)) {
                stop("Only DNA alphabet (ACGT) is supported")
            }
        }
    }

    # Parse background letter frequencies
    bg_idx <- grep("^Background letter frequencies", lines, ignore.case = TRUE)
    if (length(bg_idx) > 0) {
        # Frequencies may be on the same line or the next line(s)
        bg_text <- ""
        for (j in (bg_idx[1] + 1):min(bg_idx[1] + 2, length(lines))) {
            ln <- trimws(lines[j])
            if (nchar(ln) == 0 || grepl("^MOTIF\\b", ln, ignore.case = TRUE)) break
            bg_text <- paste(bg_text, ln)
        }
        bg_tokens <- strsplit(trimws(bg_text), "\\s+")[[1]]
        if (length(bg_tokens) >= 8) {
            bg_names <- bg_tokens[seq(1, 7, by = 2)]
            bg_vals <- suppressWarnings(as.numeric(bg_tokens[seq(2, 8, by = 2)]))
            if (!any(is.na(bg_vals))) {
                background <- setNames(bg_vals, bg_names)
            }
        }
    }

    # Find MOTIF lines
    motif_starts <- grep("^MOTIF\\b", lines)
    if (length(motif_starts) == 0) {
        stop("No motifs found in ", file)
    }

    # Parse each motif block
    results <- list()
    ids <- character(0)

    for (mi in seq_along(motif_starts)) {
        start_line <- motif_starts[mi]
        # Determine end of this motif block
        end_line <- if (mi < length(motif_starts)) motif_starts[mi + 1] - 1 else length(lines)
        block <- lines[start_line:end_line]

        # Parse MOTIF line: MOTIF <id> [<name>]
        motif_tokens <- strsplit(trimws(block[1]), "\\s+")[[1]]
        motif_id <- motif_tokens[2]
        motif_name <- if (length(motif_tokens) >= 3) paste(motif_tokens[3:length(motif_tokens)], collapse = " ") else NA_character_

        # Find letter-probability matrix header
        lpm_idx <- grep("letter-probability matrix", block, ignore.case = TRUE)
        if (length(lpm_idx) == 0) {
            # Check for log-odds matrix
            lom_idx <- grep("log-odds matrix", block, ignore.case = TRUE)
            if (length(lom_idx) > 0) {
                stop("Log-odds matrices not supported; use letter-probability matrix format")
            }
            warning("Motif '", motif_id, "' has no matrix data; skipping")
            next
        }

        # Parse metadata from the matrix header line
        header_line <- block[lpm_idx[1]]
        alength <- .parse_meme_key(header_line, "alength")
        w <- .parse_meme_key(header_line, "w")
        nsites <- .parse_meme_key(header_line, "nsites")
        E_val <- .parse_meme_key(header_line, "E")

        # Read matrix rows starting after the header line
        mat_start <- lpm_idx[1] + 1
        mat_lines <- character(0)
        for (j in mat_start:length(block)) {
            ln <- trimws(block[j])
            if (nchar(ln) == 0) next
            # Stop at URL line, next MOTIF, or non-numeric line
            if (grepl("^URL\\b", ln, ignore.case = TRUE)) break
            if (grepl("^MOTIF\\b", ln, ignore.case = TRUE)) break
            # Check if this line looks like numbers
            test_vals <- suppressWarnings(as.numeric(strsplit(ln, "\\s+")[[1]]))
            if (any(is.na(test_vals))) break
            mat_lines <- c(mat_lines, ln)
        }

        if (length(mat_lines) == 0) {
            warning("Motif '", motif_id, "' has no matrix data; skipping")
            next
        }

        # Check w consistency
        if (!is.na(w) && length(mat_lines) != as.integer(w)) {
            stop("Expected ", as.integer(w), " rows but found ", length(mat_lines), " for motif '", motif_id, "'")
        }

        # Parse matrix
        mat <- do.call(rbind, lapply(seq_along(mat_lines), function(i) {
            vals <- as.numeric(strsplit(trimws(mat_lines[i]), "\\s+")[[1]])
            if (length(vals) != 4) {
                stop("Expected 4 columns (A,C,G,T) at line ", start_line + lpm_idx[1] + i - 1, " of motif '", motif_id, "'")
            }
            vals
        }))
        colnames(mat) <- c("A", "C", "G", "T")

        # Check for log-odds (negative values)
        if (any(mat < 0, na.rm = TRUE)) {
            stop("Log-odds matrices not supported; use letter-probability matrix format")
        }

        # Check for non-numeric
        if (any(is.na(mat))) {
            stop("Non-numeric value in probability matrix for motif '", motif_id, "'")
        }

        # Re-normalize rows if needed
        mat <- .renormalize_rows(mat, context = paste0("motif '", motif_id, "'"))

        # Parse URL if present
        url_idx <- grep("^URL\\b", block, ignore.case = TRUE)
        url_val <- if (length(url_idx) > 0) {
            trimws(sub("^URL\\s+", "", block[url_idx[1]], ignore.case = TRUE))
        } else {
            NA_character_
        }

        # Set attributes
        attr(mat, "name") <- motif_name
        attr(mat, "alength") <- if (!is.na(alength)) as.integer(alength) else 4L
        attr(mat, "w") <- as.integer(nrow(mat))
        attr(mat, "nsites") <- if (!is.na(nsites)) as.numeric(nsites) else NA_real_
        attr(mat, "E") <- if (!is.na(E_val)) as.numeric(E_val) else NA_real_
        attr(mat, "url") <- url_val
        attr(mat, "strand") <- strand
        attr(mat, "background") <- background

        ids <- c(ids, motif_id)
        results <- c(results, list(mat))
    }

    if (length(results) == 0) {
        stop("No motifs found in ", file)
    }

    ids <- .dedup_motif_ids(ids)
    names(results) <- ids
    results
}


#' Parse a key=value pair from a MEME matrix header line
#' @noRd
.parse_meme_key <- function(line, key) {
    pat <- paste0("\\b", key, "=\\s*([^\\s]+)")
    m <- regmatches(line, regexpr(pat, line, perl = TRUE))
    if (length(m) == 0 || nchar(m) == 0) {
        return(NA_real_)
    }
    val <- sub(paste0("^", key, "=\\s*"), "", m)
    suppressWarnings(as.numeric(val))
}


#' Read motifs from a JASPAR PFM format file
#'
#' Parses a JASPAR Position Frequency Matrix (PFM) file and returns a named
#' list of position probability matrices (PPM). Supports both the standard
#' JASPAR header format (\code{>ID NAME} followed by labeled rows) and the
#' simple 4-row PFM format. Counts are converted to probabilities by dividing
#' each column by its column sum.
#'
#' @param file character(1) path to a JASPAR format file (\code{.jaspar}, \code{.pfm}, \code{.txt}).
#'
#' @return A named list of numeric matrices. Each matrix has columns
#'   \code{A, C, G, T} and one row per motif position. List names are motif
#'   identifiers. Each matrix carries the following attributes:
#'   \describe{
#'     \item{name}{Motif name from the header line}
#'     \item{w}{Motif width (integer)}
#'     \item{nsites}{Total counts per position (numeric; \code{NA} for simple-format files)}
#'     \item{format}{Sub-format detected: \code{"jaspar"} or \code{"simple"}}
#'   }
#'
#' @examples
#' \dontrun{
#' motifs <- gseq.read_jaspar("JASPAR2024_CORE.jaspar")
#' names(motifs)
#' m <- motifs[[1]]
#' head(m)
#' }
#'
#' @family motif functions
#' @export
gseq.read_jaspar <- function(file) {
    if (!is.character(file) || length(file) != 1) {
        stop("'file' must be a single file path string")
    }
    if (!file.exists(file)) {
        stop("File not found: ", file)
    }

    lines <- readLines(file, warn = FALSE)
    lines <- sub("\r$", "", lines)
    # Remove empty lines for easier parsing
    nonempty <- lines[nchar(trimws(lines)) > 0]

    if (length(nonempty) == 0) {
        stop("No motifs found in ", file)
    }

    has_header <- any(grepl("^>", nonempty))

    if (has_header) {
        .parse_jaspar_header(nonempty, file)
    } else {
        .parse_jaspar_simple(nonempty, file)
    }
}


#' Parse JASPAR header format (>ID NAME + labeled rows)
#' @noRd
.parse_jaspar_header <- function(lines, file) {
    header_idx <- grep("^>", lines)
    if (length(header_idx) == 0) {
        stop("No motifs found in ", file)
    }

    results <- list()
    ids <- character(0)

    for (hi in seq_along(header_idx)) {
        start <- header_idx[hi]
        end <- if (hi < length(header_idx)) header_idx[hi + 1] - 1 else length(lines)

        # Parse header: >ID NAME
        hdr <- sub("^>\\s*", "", lines[start])
        hdr_tokens <- strsplit(trimws(hdr), "\\s+")[[1]]
        motif_id <- hdr_tokens[1]
        motif_name <- if (length(hdr_tokens) >= 2) paste(hdr_tokens[-1], collapse = " ") else NA_character_

        # Grab the data rows
        data_lines <- lines[(start + 1):end]
        data_lines <- data_lines[nchar(trimws(data_lines)) > 0]

        if (length(data_lines) != 4) {
            stop("Expected 4 rows (A/C/G/T) for motif '", motif_id, "', got ", length(data_lines))
        }

        # Parse each labeled row
        count_mat <- matrix(NA_real_, nrow = 4, ncol = 0)
        expected_bases <- c("A", "C", "G", "T")
        base_order <- character(4)

        for (ri in 1:4) {
            ln <- trimws(data_lines[ri])
            # Strip brackets and row label
            ln <- gsub("[\\[\\]]", "", ln, perl = TRUE)
            # Extract the base label
            parts <- strsplit(trimws(ln), "\\s+")[[1]]
            # Filter out empty strings from splitting
            parts <- parts[nchar(parts) > 0]
            base_label <- gsub(":$", "", parts[1]) # strip trailing colon
            base_label <- toupper(base_label)
            if (!base_label %in% expected_bases) {
                stop("Unexpected row label '", base_label, "'; expected one of A, C, G, T")
            }
            base_order[ri] <- base_label
            vals <- suppressWarnings(as.numeric(parts[-1]))
            if (any(is.na(vals))) {
                stop("Non-numeric count value in motif '", motif_id, "'")
            }
            if (any(vals < 0)) {
                stop("Negative count values not allowed in JASPAR format")
            }
            if (ri == 1) {
                count_mat <- matrix(NA_real_, nrow = 4, ncol = length(vals))
            } else if (length(vals) != ncol(count_mat)) {
                stop("Rows have different lengths for motif '", motif_id, "'")
            }
            count_mat[ri, ] <- vals
        }
        rownames(count_mat) <- base_order

        # Reorder to A, C, G, T if needed
        count_mat <- count_mat[c("A", "C", "G", "T"), , drop = FALSE]

        # Convert counts to probabilities by column
        col_sums <- colSums(count_mat)
        zero_cols <- which(col_sums == 0)
        if (length(zero_cols) > 0) {
            warning(
                "Position(s) ", paste(zero_cols, collapse = ", "),
                " have zero total counts; using uniform probability"
            )
        }

        prob_mat <- count_mat
        for (j in seq_len(ncol(prob_mat))) {
            if (col_sums[j] == 0) {
                prob_mat[, j] <- 0.25
            } else {
                prob_mat[, j] <- count_mat[, j] / col_sums[j]
            }
        }

        # Transpose: bases-by-positions -> positions-by-bases
        mat <- t(prob_mat)
        colnames(mat) <- c("A", "C", "G", "T")
        rownames(mat) <- NULL

        nsites <- col_sums[1] # total counts per position

        attr(mat, "name") <- motif_name
        attr(mat, "w") <- as.integer(nrow(mat))
        attr(mat, "nsites") <- as.numeric(nsites)
        attr(mat, "format") <- "jaspar"

        ids <- c(ids, motif_id)
        results <- c(results, list(mat))
    }

    ids <- .dedup_motif_ids(ids)
    names(results) <- ids
    results
}


#' Parse JASPAR simple PFM format (4 rows of counts, no header)
#' @noRd
.parse_jaspar_simple <- function(lines, file) {
    if (length(lines) %% 4 != 0) {
        stop("Simple PFM format expects a multiple of 4 non-empty lines, got ", length(lines))
    }

    n_motifs <- length(lines) %/% 4
    results <- list()
    ids <- character(0)

    base_id <- tools::file_path_sans_ext(basename(file))

    for (mi in seq_len(n_motifs)) {
        row_start <- (mi - 1) * 4 + 1
        count_mat <- matrix(NA_real_, nrow = 4, ncol = 0)

        for (ri in 0:3) {
            ln <- trimws(lines[row_start + ri])
            # Strip brackets if present
            ln <- gsub("[\\[\\]]", "", ln, perl = TRUE)
            vals <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]]))
            if (any(is.na(vals))) {
                stop("Non-numeric count value in motif at rows ", row_start, "-", row_start + 3)
            }
            if (any(vals < 0)) {
                stop("Negative count values not allowed in JASPAR format")
            }
            if (ri == 0) {
                count_mat <- matrix(NA_real_, nrow = 4, ncol = length(vals))
            } else if (length(vals) != ncol(count_mat)) {
                stop("Rows have different lengths in simple PFM format")
            }
            count_mat[ri + 1, ] <- vals
        }

        rownames(count_mat) <- c("A", "C", "G", "T")

        # Convert counts to probabilities by column
        col_sums <- colSums(count_mat)
        zero_cols <- which(col_sums == 0)
        if (length(zero_cols) > 0) {
            warning(
                "Position(s) ", paste(zero_cols, collapse = ", "),
                " have zero total counts; using uniform probability"
            )
        }

        prob_mat <- count_mat
        for (j in seq_len(ncol(prob_mat))) {
            if (col_sums[j] == 0) {
                prob_mat[, j] <- 0.25
            } else {
                prob_mat[, j] <- count_mat[, j] / col_sums[j]
            }
        }

        # Transpose: bases-by-positions -> positions-by-bases
        mat <- t(prob_mat)
        colnames(mat) <- c("A", "C", "G", "T")
        rownames(mat) <- NULL

        motif_id <- if (n_motifs == 1) base_id else paste0(base_id, ".", mi)

        attr(mat, "name") <- NA_character_
        attr(mat, "w") <- as.integer(nrow(mat))
        attr(mat, "nsites") <- NA_real_
        attr(mat, "format") <- "simple"

        ids <- c(ids, motif_id)
        results <- c(results, list(mat))
    }

    ids <- .dedup_motif_ids(ids)
    names(results) <- ids
    results
}


#' Read motifs from a HOMER motif format file
#'
#' Parses a HOMER \code{.motif} format file and returns a named list of
#' position probability matrices (PPM). Each matrix has rows corresponding
#' to motif positions and columns \code{A}, \code{C}, \code{G}, \code{T}.
#' The returned matrices are directly usable with \code{\link{gseq.pwm}}.
#'
#' @param file character(1) path to a HOMER motif file (\code{.motif}).
#'
#' @return A named list of numeric matrices. Each matrix has columns
#'   \code{A, C, G, T} and one row per motif position. List names are
#'   derived from the consensus sequence. Each matrix carries the following
#'   attributes:
#'   \describe{
#'     \item{name}{Motif name / description from the header}
#'     \item{consensus}{Consensus sequence from the header}
#'     \item{log_odds_threshold}{Detection threshold (numeric)}
#'     \item{log_p_value}{Log p-value (numeric)}
#'     \item{w}{Motif width (integer)}
#'     \item{source}{\code{"homer"}}
#'   }
#'
#' @examples
#' \dontrun{
#' motifs <- gseq.read_homer("known_motifs.motif")
#' names(motifs)
#' m <- motifs[[1]]
#' head(m)
#' attr(m, "consensus")
#' }
#'
#' @family motif functions
#' @export
gseq.read_homer <- function(file) {
    if (!is.character(file) || length(file) != 1) {
        stop("'file' must be a single file path string")
    }
    if (!file.exists(file)) {
        stop("File not found: ", file)
    }

    lines <- readLines(file, warn = FALSE)
    lines <- sub("\r$", "", lines)

    header_idx <- grep("^>", lines)
    if (length(header_idx) == 0) {
        stop("No motifs found in HOMER file ", file)
    }

    results <- list()
    ids <- character(0)

    for (hi in seq_along(header_idx)) {
        start <- header_idx[hi]
        end <- if (hi < length(header_idx)) header_idx[hi + 1] - 1 else length(lines)

        # Parse header: tab-separated fields
        hdr <- sub("^>", "", lines[start])
        fields <- strsplit(hdr, "\t")[[1]]

        consensus <- if (length(fields) >= 1) trimws(fields[1]) else NA_character_
        motif_name <- if (length(fields) >= 2) trimws(fields[2]) else NA_character_
        log_odds_threshold <- if (length(fields) >= 3) suppressWarnings(as.numeric(fields[3])) else NA_real_
        log_p_value <- if (length(fields) >= 4) suppressWarnings(as.numeric(fields[4])) else NA_real_

        if (length(fields) < 2) {
            warning("Incomplete HOMER header at line ", start)
        }

        # Parse matrix rows
        if (start >= end) {
            warning("Motif '", consensus, "' has no matrix data; skipping")
            next
        }

        mat_lines <- lines[(start + 1):end]
        # Keep only non-empty lines
        mat_lines <- mat_lines[nchar(trimws(mat_lines)) > 0]

        if (length(mat_lines) == 0) {
            warning("Motif '", consensus, "' has no matrix data; skipping")
            next
        }

        mat <- do.call(rbind, lapply(seq_along(mat_lines), function(i) {
            vals <- suppressWarnings(as.numeric(strsplit(trimws(mat_lines[i]), "\\s+")[[1]]))
            if (length(vals) != 4) {
                stop("Expected 4 columns at line ", start + i, " of HOMER file")
            }
            if (any(is.na(vals))) {
                stop("Non-numeric probability value at line ", start + i, " of HOMER file")
            }
            vals
        }))
        colnames(mat) <- c("A", "C", "G", "T")

        # Check for negative values
        if (any(mat < 0, na.rm = TRUE)) {
            stop("Negative probability values not allowed")
        }

        # Re-normalize rows if needed
        mat <- .renormalize_rows(mat, context = paste0("motif '", consensus, "'"))

        attr(mat, "name") <- motif_name
        attr(mat, "consensus") <- consensus
        attr(mat, "log_odds_threshold") <- log_odds_threshold
        attr(mat, "log_p_value") <- log_p_value
        attr(mat, "w") <- as.integer(nrow(mat))
        attr(mat, "source") <- "homer"

        motif_id <- if (!is.na(consensus) && nchar(consensus) > 0) consensus else paste0("motif_", hi)
        ids <- c(ids, motif_id)
        results <- c(results, list(mat))
    }

    if (length(results) == 0) {
        stop("No motifs found in HOMER file ", file)
    }

    ids <- .dedup_motif_ids(ids)
    names(results) <- ids
    results
}
