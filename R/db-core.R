# Core chromosome and genome functions

.gchroms <- function(chroms) {
    if (!is.character(chroms)) {
        chroms <- as.character(chroms)
    }

    allchroms <- get("ALLGENOME", envir = .misha)[[1]]$chrom
    uniq <- unique(chroms)

    # chr-prefix toggling (chr1 <-> 1, chrX <-> X) and mitochondrial aliasing
    # (M <-> MT <-> chrM) are handled here via CHROM_ALIAS, which gsetroot
    # populates from .compute_chrom_aliases at DB-load time. Don't re-implement
    # those toggles below.
    if (exists("CHROM_ALIAS", envir = .misha, inherits = FALSE)) {
        alias_map <- get("CHROM_ALIAS", envir = .misha)
        if (length(alias_map)) {
            alias_names <- names(alias_map)
            # If none of the chromosomes match any alias names, short-circuit
            if (length(alias_names) && !any(uniq %in% alias_names)) {
                if (all(uniq %in% allchroms)) {
                    return(chroms)
                }
            } else {
                mapped <- alias_map[chroms]
                matched <- !is.na(mapped)
                # If there are zero alias hits, and all names are canonical, return early
                if (!any(matched)) {
                    if (all(uniq %in% allchroms)) {
                        return(chroms)
                    }
                }
                if (any(matched)) {
                    chroms[matched] <- mapped[matched]
                }
            }
        } else {
            # Fast path: when no aliases are defined and all names are already canonical,
            # skip the full-length match() to avoid O(N) normalization on large extracts.
            if (all(uniq %in% allchroms)) {
                return(chroms)
            }
        }
    } else {
        # Same fast path when alias map is absent
        if (all(uniq %in% allchroms)) {
            return(chroms)
        }
    }

    indices <- match(chroms, allchroms)

    err.chroms <- chroms[is.na(indices)]
    if (length(err.chroms) > 0) {
        sample_known <- utils::head(as.character(allchroms), 5)
        more <- if (length(allchroms) > 5) sprintf(", ... (%d total)", length(allchroms)) else ""
        stop(sprintf(
            "Chromosome %s does not exist in the database. Known chromosomes: %s%s. To register custom names, set .misha$CHROM_ALIAS.",
            err.chroms[1], paste(sample_known, collapse = ", "), more
        ), call. = FALSE)
    }
    allchroms[indices]
}

.gnormalize_chrom_names <- function(intervals) {
    if (!is.data.frame(intervals)) {
        return(intervals)
    }

    normalize_col <- function(df, col) {
        if (col %in% colnames(df)) {
            df[[col]] <- .gchroms(as.character(df[[col]]))
        }
        df
    }

    intervals <- normalize_col(intervals, "chrom")
    intervals <- normalize_col(intervals, "chrom1")
    intervals <- normalize_col(intervals, "chrom2")
    intervals
}

.is_per_chromosome_db <- function(groot, chromsizes) {
    if (file.exists(file.path(groot, "seq", "genome.idx"))) {
        return(FALSE)
    }

    if (!nrow(chromsizes)) {
        return(FALSE)
    }

    names_chr <- chromsizes$chrom
    no_prefix <- names_chr[!startsWith(names_chr, "chr")]

    if (!length(no_prefix)) {
        return(FALSE)
    }

    if (length(no_prefix) < 0.8 * length(names_chr)) {
        return(FALSE)
    }

    seq_dir <- file.path(groot, "seq")
    sample_names <- head(no_prefix, 5)
    prefixed_exist <- vapply(sample_names, function(name) {
        file.exists(file.path(seq_dir, paste0("chr", name, ".seq")))
    }, logical(1))

    if (!all(prefixed_exist)) {
        return(FALSE)
    }

    unprefixed_exist <- vapply(sample_names, function(name) {
        file.exists(file.path(seq_dir, paste0(name, ".seq")))
    }, logical(1))

    if (any(unprefixed_exist)) {
        return(FALSE)
    }

    TRUE
}

.compute_chrom_aliases <- function(chroms) {
    chroms <- unique(as.character(chroms))
    n <- length(chroms)
    if (!n) {
        return(stats::setNames(character(0), character(0)))
    }

    has_chr_prefix <- startsWith(chroms, "chr")
    unprefixed <- ifelse(has_chr_prefix, substring(chroms, 4), chroms)
    prefixed <- ifelse(has_chr_prefix, chroms, paste0("chr", chroms))
    upper_unpref <- toupper(unprefixed)
    upper_chr <- toupper(chroms)
    is_mito <- upper_unpref %in% c("M", "MT") | upper_chr %in% c("CHRM", "CHRMT")

    # On fragmented assemblies (Phylo447, 2.4M contigs named "CroInd_scaffold_*")
    # the prefixed-alias step generates millions of "chrXXX" entries that no user
    # will ever look up. Skip it when (a) the contig count is large, (b) no
    # canonical name has the "chr" prefix, and (c) no name is mito-pattern -- so
    # Ensembl-style references ("1", "2", "X", "MT") still get their chr-toggle
    # aliases.
    skip_prefixed <- n > 1000L && !any(has_chr_prefix) && !any(is_mito)

    # Pass 1: each canonical chrom maps to itself (already deduplicated above).
    self_names <- chroms
    self_targets <- chroms

    if (skip_prefixed) {
        return(stats::setNames(self_targets, self_names))
    }

    # Pass 2: candidate aliases (unprefixed + prefixed), tagged with source-chrom
    # index so we can reproduce the original loop's "first-seen wins" semantics.
    diff_un <- nzchar(unprefixed) & unprefixed != chroms
    diff_pf <- nzchar(prefixed) & prefixed != chroms
    cand_names <- c(unprefixed[diff_un], prefixed[diff_pf])
    cand_targets <- c(chroms[diff_un], chroms[diff_pf])
    cand_pri <- c(which(diff_un), which(diff_pf))

    # Canonical names always win: drop aliases that collide with a canonical.
    not_canon <- !(cand_names %in% chroms)
    cand_names <- cand_names[not_canon]
    cand_targets <- cand_targets[not_canon]
    cand_pri <- cand_pri[not_canon]

    # First-seen wins: order by source-chrom index, then keep first occurrence.
    o <- order(cand_pri)
    cand_names <- cand_names[o]
    cand_targets <- cand_targets[o]
    keep <- !duplicated(cand_names)
    cand_names <- cand_names[keep]
    cand_targets <- cand_targets[keep]

    out_names <- c(self_names, cand_names)
    out_targets <- c(self_targets, cand_targets)

    # Mitochondrial aliases (M, MT, chrM) all point to the first mito-like canonical.
    # Original loop adds each only if not already present.
    if (any(is_mito)) {
        mt_target <- chroms[which(is_mito)[1]]
        for (mt in c("M", "MT", "chrM")) {
            if (!(mt %in% out_names)) {
                out_names <- c(out_names, mt)
                out_targets <- c(out_targets, mt_target)
            }
        }
    }

    stats::setNames(out_targets, out_names)
}

.store_chrom_aliases <- function(chroms) {
    alias_map <- .compute_chrom_aliases(chroms)
    assign("CHROM_ALIAS", alias_map, envir = .misha)
}

# Helper function for lazy 2D genome generation
.generate_2d_on_demand <- function(intervals, mode = "full") {
    # Preserve the factor levels from the input intervals (includes aliases)
    chrom_levels <- levels(intervals$chrom)

    if (mode == "diagonal") {
        # Only intra-chromosomal pairs (chrom1 == chrom2)
        intervals2d <- cbind(intervals, intervals)
        names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
    } else if (mode == "full") {
        # Full cartesian product of all chromosome pairs
        cartesian <- expand.grid(1:nrow(intervals), 1:nrow(intervals))
        intervals2d <- cbind(intervals[cartesian[, 2], ], intervals[cartesian[, 1], ])
        names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
    } else {
        stop("Unknown 2D generation mode: ", mode, ". Must be 'diagonal' or 'full'")
    }

    # Ensure chrom1 and chrom2 have the same factor levels as intervals$chrom (including aliases)
    intervals2d$chrom1 <- factor(intervals2d$chrom1, levels = chrom_levels)
    intervals2d$chrom2 <- factor(intervals2d$chrom2, levels = chrom_levels)

    rownames(intervals2d) <- 1:nrow(intervals2d)
    intervals2d
}

# Helper function to create deferred 2D placeholder
.create_deferred_2d <- function(n_contigs) {
    intervals2d <- list() # Use empty list instead of NULL so we can set attributes
    attr(intervals2d, "deferred") <- TRUE
    attr(intervals2d, "n_contigs") <- n_contigs
    intervals2d
}

# Helper function to check if 2D is deferred
.is_2d_deferred <- function(intervals2d) {
    is.null(intervals2d) || !is.null(attr(intervals2d, "deferred"))
}
