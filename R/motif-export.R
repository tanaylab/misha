# Motif export: the writers of the formats gseq.read_meme(), gseq.read_jaspar()
# and gseq.read_homer() read. PyMisha's pymisha/motif_export.py writes the same
# bytes for the same input; keep the two in step.

# Whitespace as the writers see it: spelled out, so that R and Python agree.
.MOTIF_WS <- "[ \t\r\n\f\v]+"

#' Trim and collapse whitespace; NA for missing, empty or non-scalar values
#' @noRd
.motif_text <- function(x) {
    if (is.null(x) || length(x) != 1 || is.na(x)) {
        return(NA_character_)
    }
    x <- gsub(.MOTIF_WS, " ", enc2utf8(as.character(x)))
    x <- gsub("^ | $", "", x)
    if (!nzchar(x)) NA_character_ else x
}

#' A finite numeric scalar, or NA
#' @noRd
.motif_number <- function(x) {
    if (is.null(x) || length(x) != 1 || !(is.numeric(x) || is.logical(x))) {
        return(NA_real_)
    }
    x <- as.numeric(x)
    if (!is.finite(x)) NA_real_ else x
}

#' Probabilities: 6 decimals, as MEME writes them
#' @noRd
.motif_fmt_prob <- function(x) sprintf("%.6f", x)

#' Metadata numbers (nsites, E-values): 15 significant digits
#' @noRd
.motif_fmt_number <- function(x) sprintf("%.15g", x)

#' Counts: 6 decimals without trailing zeros, so integer counts print as integers
#' @noRd
.motif_fmt_count <- function(x) {
    x <- sub("0+$", "", sprintf("%.6f", x))
    sub("\\.$", "", x)
}

#' Write lines as UTF-8 with \n line endings on every platform
#' @noRd
.motif_write_lines <- function(lines, file) {
    con <- file(file, open = "wb")
    on.exit(close(con))
    writeBin(charToRaw(paste0(paste(enc2utf8(lines), collapse = "\n"), "\n")), con)
    invisible(file)
}

#' Check the 'file' argument of a writer
#' @noRd
.motif_check_file <- function(file) {
    if (!is.character(file) || length(file) != 1 || is.na(file) || !nzchar(file)) {
        stop("'file' must be a single file path string")
    }
}

#' Recycle a per-motif numeric argument to the number of motifs
#' @noRd
.motif_recycle <- function(x, n, what) {
    if (!is.numeric(x) || !(length(x) %in% c(1L, n)) || any(!is.finite(x))) {
        stop("'", what, "' must be a finite number or one per motif")
    }
    rep_len(as.numeric(x), n)
}

#' One motif as a positions x A,C,G,T probability matrix
#' @noRd
.motif_export_matrix <- function(m, id) {
    acgt <- c("A", "C", "G", "T")
    if (is.data.frame(m)) {
        m <- as.matrix(m)
    }
    if (!is.matrix(m) || !(is.numeric(m) || is.logical(m))) {
        stop("Motif '", id, "' must be a numeric matrix")
    }
    cn <- colnames(m)
    rn <- rownames(m)
    if (!is.null(cn)) cn <- toupper(cn)
    if (!is.null(rn)) rn <- toupper(rn)
    is_acgt <- function(x) length(x) == 4 && !anyDuplicated(x) && all(x %in% acgt)
    if (is_acgt(cn)) {
        m <- m[, match(acgt, cn), drop = FALSE]
    } else if (is_acgt(rn)) {
        m <- t(m[match(acgt, rn), , drop = FALSE])
    } else if (is.null(cn) && ncol(m) == 4) {
        # columns in A, C, G, T order
    } else {
        stop("Motif '", id, "' must have 4 columns A, C, G, T (one row per position)")
    }
    m <- matrix(as.numeric(m), ncol = 4, dimnames = list(NULL, acgt))
    if (nrow(m) == 0) {
        stop("Motif '", id, "' has no positions")
    }
    if (any(!is.finite(m))) {
        stop("Motif '", id, "' has missing or infinite values")
    }
    if (any(m < 0)) {
        stop("Motif '", id, "' has negative values; log-odds matrices cannot be written")
    }
    m[m == 0] <- 0 # no "-0.000000"
    # Summed in a fixed order in double precision (rowSums() sums in long double),
    # so that R and Python agree to the last bit.
    rsums <- m[, 1] + m[, 2] + m[, 3] + m[, 4]
    bad <- which(abs(rsums - 1) > 1e-4)
    if (length(bad) > 0) {
        warning(
            "Row sums deviate from 1.0 in motif '", id, "'; re-normalizing ",
            length(bad), " row(s)"
        )
        for (i in bad) {
            m[i, ] <- if (rsums[i] == 0) 0.25 else m[i, ] / rsums[i]
        }
    }
    m
}

#' The motifs to write: a list of list(id, mat, attrs)
#' @noRd
.motif_export_prepare <- function(motifs) {
    if (is.matrix(motifs) || is.data.frame(motifs)) {
        motifs <- list(motifs)
    }
    if (!is.list(motifs) || length(motifs) == 0) {
        stop("'motifs' must be a matrix or a non-empty list of matrices")
    }
    ids <- names(motifs)
    if (is.null(ids)) {
        ids <- rep("", length(motifs))
    }
    ids <- vapply(ids, .motif_text, character(1), USE.NAMES = FALSE)
    if (any(grepl(" ", ids, fixed = TRUE))) {
        warning("Motif IDs cannot contain whitespace; replacing it with '_'")
        ids <- gsub(" ", "_", ids, fixed = TRUE)
    }
    missing_id <- is.na(ids)
    ids[missing_id] <- paste0("motif_", which(missing_id))
    lapply(seq_along(motifs), function(i) {
        list(
            id = ids[i],
            mat = .motif_export_matrix(motifs[[i]], ids[i]),
            attrs = attributes(motifs[[i]])
        )
    })
}

#' A background vector as c(A=, C=, G=, T=) summing to 1, or NULL if invalid
#' @noRd
.motif_background <- function(bg) {
    acgt <- c("A", "C", "G", "T")
    if (is.list(bg)) {
        bg <- unlist(bg)
    }
    if (!is.numeric(bg) || length(bg) != 4 || any(!is.finite(bg)) || any(bg < 0)) {
        return(NULL)
    }
    nm <- names(bg)
    if (!is.null(nm)) {
        nm <- toupper(nm)
        if (anyDuplicated(nm) || !all(nm %in% acgt)) {
            return(NULL)
        }
        bg <- bg[match(acgt, nm)]
    }
    bg <- as.numeric(bg)
    bg[bg == 0] <- 0
    s <- bg[1] + bg[2] + bg[3] + bg[4]
    if (s == 0) {
        return(NULL)
    }
    if (abs(s - 1) > 1e-4) {
        bg <- bg / s
    }
    stats::setNames(bg, acgt)
}

#' The strands line value ("+ -" or "+"), or NULL if invalid
#' @noRd
.motif_strands <- function(x) {
    x <- .motif_text(x)
    if (is.na(x)) {
        return(NULL)
    }
    x <- gsub(" ", "", x, fixed = TRUE)
    if (x %in% c("+-", "-+")) {
        "+ -"
    } else if (x == "+") {
        "+"
    } else {
        NULL
    }
}

#' The value the motifs agree on for an attribute, after `parse`; NULL if none
#' or if they disagree (with a warning naming `what`)
#' @noRd
.motif_common_attr <- function(ms, name, parse, what) {
    vals <- lapply(ms, function(m) {
        v <- m$attrs[[name]]
        if (is.null(v)) NULL else parse(v)
    })
    vals <- Filter(Negate(is.null), vals)
    if (length(vals) == 0) {
        return(NULL)
    }
    if (!all(vapply(vals, identical, logical(1), vals[[1]]))) {
        warning("Motifs carry different ", what, "; writing the default")
        return(NULL)
    }
    vals[[1]]
}


#' Write motifs to a MEME minimal motif format file
#'
#' Writes position probability matrices to a file in the MEME minimal motif
#' format, which \code{\link{gseq.read_meme}} and the MEME Suite read.
#'
#' Each motif is written as a \code{MOTIF <id> [<name>]} line, a
#' \code{letter-probability matrix:} line with \code{alength}, \code{w} and,
#' when the motif has them, \code{nsites} and \code{E}, the matrix rows and,
#' when the motif has one, a \code{URL} line. The file header declares the DNA
#' alphabet, the strands and the background letter frequencies.
#'
#' Probabilities are written with 6 decimals, as MEME writes them, so a written
#' motif reads back within 5e-7 of the original. \code{nsites} and E-values are
#' written with 15 significant digits.
#'
#' Rows that do not sum to 1 are re-normalized (with a warning), as
#' \code{\link{gseq.read_meme}} does. Whitespace in IDs is replaced by
#' \code{_}; whitespace in names is collapsed to single spaces (other tools may
#' keep only the first word of a name). Motifs with the same ID are written
#' as they are; \code{\link{gseq.read_meme}} adds suffixes to tell them apart.
#'
#' @param motifs a named list of motif matrices, such as returned by
#'   \code{\link{gseq.read_meme}}, \code{\link{gseq.read_jaspar}} or
#'   \code{\link{gseq.read_homer}}, or a single matrix. Each matrix has one row
#'   per position and columns \code{A, C, G, T} (unnamed columns are taken in
#'   that order; a matrix with rows \code{A, C, G, T} is transposed). List names
#'   are the motif IDs; unnamed motifs are called \code{motif_<i>}. The
#'   attributes \code{name}, \code{nsites}, \code{E}, \code{url}, \code{strand}
#'   and \code{background} are written when present.
#' @param file character(1) path of the file to write.
#' @param background background letter frequencies: a numeric vector of 4
#'   values named \code{A, C, G, T} (or unnamed, in that order), normalized to
#'   sum to 1. \code{NULL} (the default) uses the \code{background} attribute
#'   the motifs share (as \code{\link{gseq.read_meme}} returns it), and uniform
#'   frequencies otherwise.
#' @param strands \code{"+ -"} (both strands) or \code{"+"}. \code{NULL} (the
#'   default) uses the \code{strand} attribute the motifs share, and \code{"+ -"}
#'   otherwise.
#'
#' @return \code{file}, invisibly.
#'
#' @examples
#' motifs <- list(
#'     M1 = matrix(c(
#'         0.7, 0.1, 0.1, 0.1,
#'         0.1, 0.1, 0.7, 0.1
#'     ), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
#' )
#' out <- tempfile(fileext = ".meme")
#' gseq.write_meme(motifs, out, background = c(A = 0.3, C = 0.2, G = 0.2, T = 0.3))
#' cat(readLines(out), sep = "\n")
#' gseq.read_meme(out)
#'
#' @family motif functions
#' @export
gseq.write_meme <- function(motifs, file, background = NULL, strands = NULL) {
    .motif_check_file(file)
    ms <- .motif_export_prepare(motifs)

    if (is.null(background)) {
        bg <- .motif_common_attr(ms, "background", .motif_background, "backgrounds")
        if (is.null(bg)) {
            bg <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
        }
    } else {
        bg <- .motif_background(background)
        if (is.null(bg)) {
            stop("'background' must be 4 non-negative frequencies named A, C, G, T (or in that order), not all zero")
        }
    }

    if (is.null(strands)) {
        st <- .motif_common_attr(ms, "strand", .motif_strands, "strands")
        if (is.null(st)) {
            st <- "+ -"
        }
    } else {
        st <- .motif_strands(strands)
        if (is.null(st)) {
            stop("'strands' must be \"+ -\" or \"+\"")
        }
    }

    out <- c(
        "MEME version 4", "",
        "ALPHABET= ACGT", "",
        paste0("strands: ", st), "",
        "Background letter frequencies",
        paste("A", .motif_fmt_prob(bg[1]), "C", .motif_fmt_prob(bg[2]), "G", .motif_fmt_prob(bg[3]), "T", .motif_fmt_prob(bg[4])),
        ""
    )
    for (m in ms) {
        name <- .motif_text(m$attrs[["name"]])
        out <- c(out, if (is.na(name)) paste("MOTIF", m$id) else paste("MOTIF", m$id, name))
        hdr <- paste0("letter-probability matrix: alength= 4 w= ", nrow(m$mat))
        nsites <- .motif_number(m$attrs[["nsites"]])
        if (!is.na(nsites)) {
            hdr <- paste0(hdr, " nsites= ", .motif_fmt_number(nsites))
        }
        e_val <- .motif_number(m$attrs[["E"]])
        if (!is.na(e_val)) {
            hdr <- paste0(hdr, " E= ", .motif_fmt_number(e_val))
        }
        p <- .motif_fmt_prob(m$mat)
        dim(p) <- dim(m$mat)
        out <- c(out, hdr, paste(p[, 1], p[, 2], p[, 3], p[, 4]))
        url <- .motif_text(m$attrs[["url"]])
        if (!is.na(url)) {
            out <- c(out, paste("URL", url))
        }
        out <- c(out, "")
    }
    .motif_write_lines(out, file)
}


#' Write motifs to a JASPAR format file
#'
#' Writes motifs as count matrices in the JASPAR format that
#' \code{\link{gseq.read_jaspar}} reads: a \code{>ID name} header line and
#' one \code{A [ ... ]} row of counts per base. \code{format = "simple"} writes
#' the four rows of counts only, without headers or labels.
#'
#' JASPAR matrices hold counts, so each motif's probabilities are multiplied by
#' its number of sites: the \code{nsites} argument, else the motif's
#' \code{nsites} attribute (as \code{\link{gseq.read_jaspar}} and
#' \code{\link{gseq.read_meme}} return it), else 100. Every position then sums
#' to that number. Counts are written with up to 6 decimals and without
#' trailing zeros, so the integer counts of a JASPAR file are written as
#' integers, and the probabilities read back within 1.5e-6 / \code{nsites}
#' of the original. Some tools, such as the MEME Suite's \code{jaspar2meme},
#' read only integer counts: \code{round_counts = TRUE} rounds the counts, at
#' the cost of precision (probabilities then read back within about
#' 0.5 / \code{nsites}; with an \code{nsites} of a few, whole positions can
#' round to zero, which readers take as uniform).
#'
#' The header names a motif by its \code{name} attribute, or by its ID when it
#' has none (JASPAR headers always carry a name). Rows that do not sum to 1 are
#' re-normalized (with a warning). Whitespace in IDs is replaced by \code{_};
#' whitespace in names is collapsed to single spaces (other tools may keep only
#' the first word of a name). A simple-format file has no IDs:
#' \code{\link{gseq.read_jaspar}} names its motifs after the file.
#'
#' @inheritParams gseq.write_meme
#' @param motifs a named list of motif matrices, such as returned by the
#'   motif readers, or a single matrix (see \code{\link{gseq.write_meme}}). The
#'   attributes \code{name} and \code{nsites} are used when present.
#' @param nsites number of sites to scale probabilities to counts: a positive
#'   number, or one per motif. \code{NULL} (the default) uses each motif's
#'   \code{nsites} attribute, and 100 for motifs without one.
#' @param format \code{"jaspar"} (headers and labeled rows) or \code{"simple"}
#'   (four rows of counts per motif).
#' @param round_counts round counts to integers.
#'
#' @return \code{file}, invisibly.
#'
#' @examples
#' m <- matrix(c(
#'     0.7, 0.1, 0.1, 0.1,
#'     0.1, 0.1, 0.7, 0.1
#' ), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
#' out <- tempfile(fileext = ".jaspar")
#' gseq.write_jaspar(list(M1 = m), out, nsites = 20)
#' cat(readLines(out), sep = "\n")
#' gseq.read_jaspar(out)
#'
#' @family motif functions
#' @export
gseq.write_jaspar <- function(motifs, file, nsites = NULL, format = c("jaspar", "simple"), round_counts = FALSE) {
    .motif_check_file(file)
    format <- match.arg(format)
    if (!is.logical(round_counts) || length(round_counts) != 1 || is.na(round_counts)) {
        stop("'round_counts' must be TRUE or FALSE")
    }
    ms <- .motif_export_prepare(motifs)
    if (!is.null(nsites)) {
        nsites <- .motif_recycle(nsites, length(ms), "nsites")
        if (any(nsites <= 0)) {
            stop("'nsites' must be positive")
        }
    }

    out <- character(0)
    for (i in seq_along(ms)) {
        m <- ms[[i]]
        ns <- if (!is.null(nsites)) nsites[i] else .motif_number(m$attrs[["nsites"]])
        if (is.na(ns) || ns <= 0) {
            ns <- 100
        }
        counts <- t(m$mat * ns)
        if (round_counts) {
            counts <- round(counts) # half to even, as Python's round()
        }
        counts <- .motif_fmt_count(counts)
        counts <- sprintf("%*s", max(nchar(counts)), counts)
        dim(counts) <- c(4, nrow(m$mat))
        rows <- apply(counts, 1, paste, collapse = " ")
        if (i > 1) {
            out <- c(out, "")
        }
        if (format == "jaspar") {
            name <- .motif_text(m$attrs[["name"]])
            if (is.na(name)) {
                name <- m$id
            }
            out <- c(
                out,
                paste0(">", m$id, " ", name),
                paste0(c("A", "C", "G", "T"), "  [ ", rows, " ]")
            )
        } else {
            out <- c(out, rows)
        }
    }
    .motif_write_lines(out, file)
}


#' IUPAC consensus of a probability matrix, by HOMER's rule
#'
#' A position's symbol covers every base within 0.15 of its most probable base.
#' @noRd
.homer_consensus <- function(mat) {
    # Indexed by A=1 + C=2 + G=4 + T=8
    iupac <- c("A", "C", "M", "G", "R", "S", "V", "T", "W", "Y", "H", "K", "D", "B", "N")
    bits <- c(1L, 2L, 4L, 8L)
    paste(vapply(seq_len(nrow(mat)), function(i) {
        row <- mat[i, ]
        iupac[sum(bits[(max(row) - row) < 0.15])]
    }, character(1)), collapse = "")
}

#' The highest log-odds score of a matrix against a uniform background
#' (natural log, as HOMER scores)
#' @noRd
.homer_max_score <- function(mat) {
    s <- 0
    for (i in seq_len(nrow(mat))) {
        s <- s + log(max(mat[i, ]) / 0.25)
    }
    s
}


#' Write motifs to a HOMER motif format file
#'
#' Writes motifs in the HOMER \code{.motif} format that
#' \code{\link{gseq.read_homer}} and HOMER read: a header line
#' \code{>consensus<TAB>name<TAB>log-odds threshold[<TAB>log p-value]} and one
#' tab-separated row of A, C, G, T probabilities per position.
#'
#' The consensus is the motif's \code{consensus} attribute (as
#' \code{\link{gseq.read_homer}} returns it) when it has one of the motif's
#' width, and otherwise is computed by HOMER's rule: each position's IUPAC
#' symbol covers every base whose probability is within 0.15 of the highest.
#' The name is the motif's \code{name} attribute, or its ID when it has none.
#' \code{\link{gseq.read_homer}} names motifs by their consensus, so IDs that
#' are not used as names are not kept.
#'
#' The log-odds threshold (natural log against a uniform background, as HOMER
#' scores) is the \code{threshold} argument, else the motif's
#' \code{log_odds_threshold} attribute, else \code{threshold_frac} times the
#' motif's highest possible score. The log p-value is written when the motif
#' has a \code{log_p_value} attribute.
#'
#' Probabilities, thresholds and log p-values are written with 6 decimals,
#' so a written motif reads back within 5e-7 of the original. Rows that do not
#' sum to 1 are re-normalized (with a warning). Tabs and other whitespace in
#' names are collapsed to single spaces.
#'
#' @inheritParams gseq.write_meme
#' @param motifs a named list of motif matrices, such as returned by the motif
#'   readers, or a single matrix (see \code{\link{gseq.write_meme}}). The
#'   attributes \code{name}, \code{consensus}, \code{log_odds_threshold} and
#'   \code{log_p_value} are used when present.
#' @param threshold log-odds detection threshold: a number, or one per motif.
#'   \code{NULL} (the default) uses each motif's \code{log_odds_threshold}
#'   attribute, and \code{threshold_frac} of its highest score for motifs
#'   without one.
#' @param threshold_frac fraction of a motif's highest log-odds score used as
#'   its threshold when neither \code{threshold} nor the motif gives one.
#'
#' @return \code{file}, invisibly.
#'
#' @examples
#' m <- matrix(c(
#'     0.7, 0.1, 0.1, 0.1,
#'     0.1, 0.1, 0.7, 0.1
#' ), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
#' out <- tempfile(fileext = ".motif")
#' gseq.write_homer(list(M1 = m), out)
#' cat(readLines(out), sep = "\n")
#' gseq.read_homer(out)
#'
#' @family motif functions
#' @export
gseq.write_homer <- function(motifs, file, threshold = NULL, threshold_frac = 0.6) {
    .motif_check_file(file)
    ms <- .motif_export_prepare(motifs)
    if (!is.null(threshold)) {
        threshold <- .motif_recycle(threshold, length(ms), "threshold")
    }
    if (!is.numeric(threshold_frac) || length(threshold_frac) != 1 || !is.finite(threshold_frac)) {
        stop("'threshold_frac' must be a finite number")
    }

    out <- character(0)
    for (i in seq_along(ms)) {
        m <- ms[[i]]
        cons <- .motif_text(m$attrs[["consensus"]])
        if (is.na(cons) || grepl(" ", cons, fixed = TRUE) || nchar(cons) != nrow(m$mat)) {
            cons <- .homer_consensus(m$mat)
        }
        name <- .motif_text(m$attrs[["name"]])
        if (is.na(name)) {
            name <- m$id
        }
        thr <- if (!is.null(threshold)) threshold[i] else .motif_number(m$attrs[["log_odds_threshold"]])
        if (is.na(thr)) {
            thr <- threshold_frac * .homer_max_score(m$mat)
        }
        hdr <- paste0(">", cons, "\t", name, "\t", .motif_fmt_prob(thr))
        log_p <- .motif_number(m$attrs[["log_p_value"]])
        if (!is.na(log_p)) {
            hdr <- paste0(hdr, "\t", .motif_fmt_prob(log_p))
        }
        p <- .motif_fmt_prob(m$mat)
        dim(p) <- dim(m$mat)
        out <- c(out, hdr, paste(p[, 1], p[, 2], p[, 3], p[, 4], sep = "\t"))
    }
    .motif_write_lines(out, file)
}
