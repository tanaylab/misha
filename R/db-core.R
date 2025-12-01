# Core chromosome and genome functions

.gchroms <- function(chroms) {
    if (!is.character(chroms)) {
        chroms <- as.character(chroms)
    }

    allchroms <- get("ALLGENOME", envir = .misha)[[1]]$chrom
    uniq <- unique(chroms)

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
        stop(sprintf("Chromosome %s does not exist in the database", err.chroms[1]))
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

    # Pre-compute all string transformations (vectorized)
    has_chr_prefix <- startsWith(chroms, "chr")
    unprefixed <- ifelse(has_chr_prefix, substring(chroms, 4), chroms)
    prefixed <- ifelse(has_chr_prefix, chroms, paste0("chr", chroms))
    upper_unprefixed <- toupper(unprefixed)
    upper_chroms <- toupper(chroms)

    # Build alias mappings efficiently using environment for O(1) existence checks
    seen <- new.env(hash = TRUE, parent = emptyenv())

    # Pre-allocate result vectors (estimate ~4 entries per chromosome)
    max_aliases <- length(chroms) * 4
    alias_names <- character(max_aliases)
    alias_targets <- character(max_aliases)
    idx <- 0

    # Inline alias addition to avoid <<- (CRAN compliance)
    # First pass: add all chromosomes mapping to themselves
    for (i in seq_along(chroms)) {
        chrom <- chroms[i]
        if (!exists(chrom, envir = seen, inherits = FALSE)) {
            seen[[chrom]] <- TRUE
            idx <- idx + 1
            alias_names[idx] <- chrom
            alias_targets[idx] <- chrom
        }
    }

    # Second pass: add aliases using pre-computed vectorized values
    for (i in seq_along(chroms)) {
        chrom <- chroms[i]

        # Add unprefixed alias
        if (unprefixed[i] != chrom) {
            name <- unprefixed[i]
            if (nzchar(name) && !exists(name, envir = seen, inherits = FALSE)) {
                seen[[name]] <- TRUE
                idx <- idx + 1
                alias_names[idx] <- name
                alias_targets[idx] <- chrom
            }
        }

        # Add prefixed alias
        if (prefixed[i] != chrom) {
            name <- prefixed[i]
            if (nzchar(name) && !exists(name, envir = seen, inherits = FALSE)) {
                seen[[name]] <- TRUE
                idx <- idx + 1
                alias_names[idx] <- name
                alias_targets[idx] <- chrom
            }
        }

        # Add mitochondrial aliases
        if (upper_unprefixed[i] %in% c("M", "MT") || upper_chroms[i] %in% c("CHRM", "CHRMT")) {
            for (mt_alias in c("M", "MT", "chrM")) {
                if (!exists(mt_alias, envir = seen, inherits = FALSE)) {
                    seen[[mt_alias]] <- TRUE
                    idx <- idx + 1
                    alias_names[idx] <- mt_alias
                    alias_targets[idx] <- chrom
                }
            }
        }
    }

    # Build final named vector
    alias <- stats::setNames(alias_targets[1:idx], alias_names[1:idx])
    alias
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
