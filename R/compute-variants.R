# A numeric matrix with columns lo and hi, one row per sequence vtrack in
# `expr`: a bin [s, e) of that vtrack reads [s + lo, e + hi).
.gvariants_reach <- function(expr) {
    x <- .gcall("C_gvariants_reach", expr, .misha_env())
    if (!length(x)) {
        stop(sprintf("Expression \"%s\" has no sequence virtual track, so no substitution can change it", expr), call. = FALSE)
    }
    matrix(x, ncol = 2, dimnames = list(NULL, c("lo", "hi")))
}

# The value of `expr` on each row, with base `base` (0..3 = A, C, G, T; -1 =
# none) at position `pos` while the row is computed. Rows sorted by chromosome
# and start.
.gvariants_eval <- function(expr, chroms, starts, ends, pos, base) {
    .gcall(
        "C_gextract_variants", expr, as.character(chroms), as.numeric(starts), as.numeric(ends),
        as.numeric(pos), as.integer(base), .misha_env()
    )
}

#' Evaluate a track expression under single-base substitutions
#'
#' For each iterator bin, evaluates \code{expr} once on the reference genome and
#' once for every single-base substitution the bin's sequence virtual tracks
#' read: all three alternative bases at each base in reach (the default), or the
#' variants listed in \code{variants}.
#'
#' A bin \code{[s, e)} of a sequence virtual track reads \code{[s + sshift, e +
#' eshift + W - 1)}, where the shifts come from \code{\link{gvtrack.iterator}}
#' and \code{W - 1} is added when the track extends (\code{extend = TRUE}, motif
#' or model width \code{W}). Two families read further: \code{kmer.count} and
#' \code{kmer.frac} with \code{strand} 0 or -1 start at \code{s + sshift - (k -
#' 1)} for a k-mer of length \code{k} (when \code{extend = TRUE}), and
#' \code{pwm.edit_distance} with \code{max_indels} ends \code{max_indels} bases
#' later, at \code{e + eshift + W - 1 + max_indels}. A substitution is paired with every bin whose reach,
#' the union over the expression's sequence virtual tracks, contains it. Tracks
#' and non-sequence virtual tracks in \code{expr} are read as usual, so they
#' contribute the same value to both columns. The substituted base takes the
#' case of the reference base, so \code{masked.*} tracks see the same mask.
#' Only single-base substitutions are supported.
#'
#' The exact pass runs in one process. For genome-scale runs, split
#' \code{intervals} by chromosome and run the pieces in parallel.
#'
#' Run genome-wide, every bin returns about three rows per base in reach. To
#' find the bins one substitution can move across a threshold, screen first with
#' a bound (for example from \code{potts.max} virtual tracks with
#' \code{max_edits = 1}) and pass the result as \code{intervals}, with
#' \code{score.thresh}.
#'
#' @param expr a track expression with at least one sequence virtual track.
#' @param intervals genomic scope, as in \code{\link{gextract}}.
#' @param variants \code{NULL} for every substitution in reach, or a data frame
#'   with columns \code{chrom}, \code{start} (0-based), \code{ref} and \code{alt},
#'   such as the output of \code{\link{gintervals.import_vcf}}. A variant whose
#'   \code{ref} does not match the genome gets a warning and \code{NA} in
#'   \code{<name>.alt}; \code{<name>.ref} keeps the reference value.
#' @param iterator track expression iterator, as in \code{\link{gextract}}.
#' @param score.thresh if not \code{NULL}, keep only the rows where the reference
#'   value and the substituted value fall on different sides of it.
#' @param colnames name used for the value columns, \code{<name>.ref} and
#'   \code{<name>.alt}; the expression itself by default.
#' @return A data frame with one row per (bin, substitution): the bin's
#'   \code{chrom}, \code{start} and \code{end}; the substitution's \code{pos}
#'   (0-based), \code{ref} and \code{alt}; \code{<name>.ref} and
#'   \code{<name>.alt}; and \code{intervalID} as in \code{\link{gextract}}.
#'   \code{NULL} when there is no row.
#' @seealso \code{\link{gextract}}, \code{\link{gvtrack.create}},
#'   \code{\link{gintervals.import_vcf}}
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' pssm <- matrix(c(
#'     0.7, 0.1, 0.1, 0.1,
#'     0.1, 0.7, 0.1, 0.1,
#'     0.1, 0.1, 0.7, 0.1
#' ), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
#' gvtrack.create("motif", NULL, "pwm.max", params = list(pssm = pssm, bidirect = FALSE))
#' gextract.variants("motif", gintervals(1, 100, 102), iterator = 1)
#'
#' @export
gextract.variants <- function(expr = NULL, intervals = NULL, variants = NULL, iterator = NULL,
                              score.thresh = NULL, colnames = NULL) {
    if (is.null(expr) || is.null(intervals)) {
        stop("Usage: gextract.variants(expr, intervals, variants = NULL, iterator = NULL, score.thresh = NULL, colnames = NULL)", call. = FALSE)
    }
    .gcheckroot()
    expr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
    .gextract_variants_impl(expr, intervals, variants, .iterator, score.thresh, colnames)
}

.gextract_variants_impl <- function(expr, intervals, variants, iterator, score.thresh, colnames,
                                    chunk = NULL) {
    if (length(expr) != 1L) {
        stop("gextract.variants takes a single track expression", call. = FALSE)
    }
    name <- if (is.null(colnames)) expr else colnames
    if (!is.character(name) || length(name) != 1L) {
        stop("colnames must be a single string", call. = FALSE)
    }
    if (!is.null(score.thresh)) {
        score.thresh <- .coerce_score_thresh(score.thresh)
    }
    reach <- .gvariants_reach(expr)
    if (!is.null(variants)) {
        variants <- .gvariants_check(variants)
        # Checked once, against the genome, whether or not any bin reads them.
        genome_ref <- toupper(gseq.extract(data.frame(chrom = variants$chrom, start = variants$pos, end = variants$pos + 1)))
        variants$ok <- genome_ref == variants$ref
        if (!all(variants$ok)) {
            warning(sprintf(
                "%d of %d variants have a ref that does not match the genome; their substituted values are NA",
                sum(!variants$ok), nrow(variants)
            ), call. = FALSE)
        }
    }

    bins <- do.call(gextract, list(expr, intervals = quote(intervals), iterator = iterator, colnames = "value"))
    if (is.null(bins) || !nrow(bins)) {
        return(NULL)
    }
    if (is.null(chunk)) {
        # About 3e6 (bin, substitution) rows per chunk, whatever the reach.
        span <- max(reach[, "hi"]) - min(reach[, "lo"]) + max(bins$end - bins$start)
        span <- max(1, span)
        chunk <- max(1L, as.integer(floor(3e6 / (3 * span))))
    }
    chrom_order <- match(as.character(bins$chrom), as.character(gintervals.all()$chrom))

    out <- list()
    for (from in seq(1L, nrow(bins), by = chunk)) {
        rng <- from:min(from + chunk - 1L, nrow(bins))
        b <- bins[rng, , drop = FALSE]
        p <- .gvariants_bin_positions(b, reach, with_ref = is.null(variants))
        if (is.null(variants)) {
            r <- .gvariants_saturate(p)
            r$ok <- rep(TRUE, nrow(r))
        } else {
            p$chrom <- as.character(b$chrom[p$bin])
            if (!nrow(p)) {
                next
            }
            # Only the variants this chunk can reach go into the join.
            v <- variants[variants$chrom %in% unique(p$chrom) &
                variants$pos >= min(p$pos) & variants$pos <= max(p$pos), , drop = FALSE]
            r <- merge(p, v, by = c("chrom", "pos"))
            r <- r[order(r$bin, r$pos, r$alt), c("bin", "pos", "ref", "alt", "ok"), drop = FALSE]
        }
        if (!nrow(r)) {
            next
        }
        alt_value <- rep(NA_real_, nrow(r))
        e <- which(r$ok)
        if (length(e)) {
            # The engine needs rows sorted by chromosome and start; gextract's
            # bins need not be (overlapping intervals repeat a stretch).
            e <- e[order(chrom_order[rng][r$bin[e]], b$start[r$bin[e]])]
            alt_value[e] <- .gvariants_eval(
                expr, b$chrom[r$bin[e]], b$start[r$bin[e]], b$end[r$bin[e]],
                r$pos[e], match(r$alt[e], c("A", "C", "G", "T")) - 1L
            )
        }
        res <- data.frame(
            chrom = b$chrom[r$bin], start = b$start[r$bin], end = b$end[r$bin],
            pos = r$pos, ref = r$ref, alt = r$alt,
            ref_value = b$value[r$bin], alt_value = alt_value,
            intervalID = b$intervalID[r$bin], stringsAsFactors = FALSE
        )
        if (!is.null(score.thresh)) {
            res <- res[!is.na(res$ref_value) & !is.na(res$alt_value) &
                (res$ref_value >= score.thresh) != (res$alt_value >= score.thresh), , drop = FALSE]
        }
        out[[length(out) + 1L]] <- res
    }

    res <- do.call(rbind, out)
    if (is.null(res) || !nrow(res)) {
        return(NULL)
    }
    rownames(res) <- NULL
    names(res)[names(res) == "ref_value"] <- paste0(name, ".ref")
    names(res)[names(res) == "alt_value"] <- paste0(name, ".alt")
    res
}

# One row per (bin, position) a bin reads, with the genome base (upper case) at
# it (unless `with_ref` is FALSE): `bins` has chrom, start, end; `reach` is
# .gvariants_reach()'s matrix.
.gvariants_bin_positions <- function(bins, reach, with_ref = TRUE) {
    sizes <- gintervals.all()
    chrom_size <- sizes$end[match(as.character(bins$chrom), as.character(sizes$chrom))]
    parts <- lapply(seq_len(nrow(reach)), function(v) {
        from <- pmax(bins$start + reach[v, "lo"], 0)
        to <- pmin(bins$end + reach[v, "hi"], chrom_size)
        n <- pmax(to - from, 0)
        bin <- rep.int(seq_len(nrow(bins)), n)
        data.frame(bin = bin, pos = from[bin] + sequence(n) - 1)
    })
    p <- do.call(rbind, parts)
    p <- p[!duplicated(p$bin * 4294967296 + p$pos), , drop = FALSE]
    p <- p[order(p$bin, p$pos), , drop = FALSE]

    if (!with_ref) {
        return(p)
    }
    span_lo <- pmax(bins$start + min(reach[, "lo"]), 0)
    span_hi <- pmin(bins$end + max(reach[, "hi"]), chrom_size)
    has <- span_hi > span_lo
    seqs <- rep(NA_character_, nrow(bins))
    seqs[has] <- toupper(gseq.extract(data.frame(chrom = bins$chrom[has], start = span_lo[has], end = span_hi[has])))
    off <- p$pos - span_lo[p$bin] + 1
    p$ref <- substr(seqs[p$bin], off, off)
    p
}

# Three rows per A/C/G/T position, one per alternative base.
.gvariants_saturate <- function(p) {
    p <- p[p$ref %in% c("A", "C", "G", "T"), , drop = FALSE]
    alts <- matrix(c("C", "G", "T", "A", "G", "T", "A", "C", "T", "A", "C", "G"), nrow = 3)
    i <- rep(seq_len(nrow(p)), each = 3L)
    data.frame(
        bin = p$bin[i], pos = p$pos[i], ref = p$ref[i],
        alt = as.vector(alts[, match(p$ref, c("A", "C", "G", "T"))]),
        stringsAsFactors = FALSE
    )
}

# chrom (canonical), pos, ref, alt (upper case), one row per variant.
.gvariants_check <- function(variants) {
    if (!is.data.frame(variants) || !all(c("chrom", "start", "ref", "alt") %in% colnames(variants))) {
        stop("variants must be a data frame with columns chrom, start, ref and alt (gintervals.import_vcf() output works)", call. = FALSE)
    }
    ref <- toupper(as.character(variants$ref))
    alt <- toupper(as.character(variants$alt))
    bad <- which(!(ref %in% c("A", "C", "G", "T")) | !(alt %in% c("A", "C", "G", "T")))
    if (length(bad)) {
        stop(sprintf(
            "variants row %d: only single-base substitutions are supported (ref and alt each one of A, C, G, T); got ref \"%s\", alt \"%s\"",
            bad[1], variants$ref[bad[1]], variants$alt[bad[1]]
        ), call. = FALSE)
    }
    same <- which(ref == alt)
    if (length(same)) {
        stop(sprintf("variants row %d: ref and alt are the same base", same[1]), call. = FALSE)
    }
    start <- variants$start
    if (!is.numeric(start) || anyNA(start) || any(start < 0) || any(start != floor(start))) {
        stop("variants: start must be a non-negative whole number (0-based)", call. = FALSE)
    }
    chrom <- as.character(.gchroms(as.character(variants$chrom)))
    sizes <- gintervals.all()
    out_of_range <- which(start >= sizes$end[match(chrom, as.character(sizes$chrom))])
    if (length(out_of_range)) {
        stop(sprintf("variants row %d: start %.0f is past the end of %s", out_of_range[1], start[out_of_range[1]], chrom[out_of_range[1]]), call. = FALSE)
    }
    data.frame(chrom = chrom, pos = as.numeric(start), ref = ref, alt = alt, stringsAsFactors = FALSE)
}
