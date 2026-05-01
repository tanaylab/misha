# Core chromosome and genome functions

# Build (or refresh) an environment-based mirror of CHROM_ALIAS for O(1)
# lookup. R's named-vector subscript path uses match(), which rebuilds an
# internal hash on every call — O(M) per call, where M is the alias count.
# For large alias maps (e.g. NCBI assemblies with chrom_aliases.tsv loaded:
# ~12k entries on rabbit), this becomes the dominant cost in .gchroms().
# An environment offers genuine O(1) hashed lookup that survives across calls.
.refresh_chrom_alias_env <- function(full_alias_map = NULL) {
    # Drop any existing mirror first.
    if (exists("CHROM_ALIAS_ENV", envir = .misha, inherits = FALSE)) {
        rm("CHROM_ALIAS_ENV", envir = .misha)
    }
    if (!exists("CHROM_ALIAS", envir = .misha, inherits = FALSE)) {
        return(invisible(NULL))
    }
    # The env mirror exists only to give O(1) lookup for the *extra* aliases
    # (TSV-derived NCBI/GenBank/sequenceName) that don't go into the
    # C++-facing CHROM_ALIAS. If the full map is the same as CHROM_ALIAS
    # (no TSV extras present), skip building the env entirely — .gchroms
    # then uses the original named-vector path with no per-call overhead.
    cpp_map <- get("CHROM_ALIAS", envir = .misha)
    alias_map <- full_alias_map
    if (is.null(alias_map) || length(alias_map) == length(cpp_map)) {
        return(invisible(NULL))
    }
    env <- new.env(
        hash = TRUE, parent = emptyenv(),
        size = max(length(alias_map), 1L)
    )
    if (length(alias_map)) {
        nm <- names(alias_map)
        for (i in seq_along(alias_map)) {
            env[[nm[i]]] <- alias_map[[i]]
        }
    }
    assign("CHROM_ALIAS_ENV", env, envir = .misha)
    invisible(NULL)
}

# Look up a single chrom name in the alias env. Returns the canonical target
# or NA_character_ if not present.
.chrom_alias_lookup <- function(name, env) {
    if (exists(name, envir = env, inherits = FALSE)) {
        get(name, envir = env, inherits = FALSE)
    } else {
        NA_character_
    }
}

.gchroms <- function(chroms) {
    if (!is.character(chroms)) {
        chroms <- as.character(chroms)
    }

    allchroms <- get("ALLGENOME", envir = .misha)[[1]]$chrom
    uniq <- unique(chroms)

    # chr-prefix toggling (chr1 <-> 1, chrX <-> X), mitochondrial aliasing
    # (M <-> MT <-> chrM), and accession aliasing (chr1 <-> NC_*) are all
    # handled via CHROM_ALIAS / CHROM_ALIAS_ENV, populated by
    # .compute_chrom_aliases at DB-load time. The env mirror gives O(1)
    # lookup; falling back to the named vector path keeps backward compat
    # if an external caller manipulates .misha$CHROM_ALIAS directly without
    # refreshing the env.
    has_alias_env <- exists("CHROM_ALIAS_ENV", envir = .misha, inherits = FALSE)
    has_alias_map <- exists("CHROM_ALIAS", envir = .misha, inherits = FALSE)

    if (has_alias_env) {
        env <- get("CHROM_ALIAS_ENV", envir = .misha)
        if (length(env)) {
            # Fast probe: do any unique chroms hit the alias env at all?
            any_hit <- FALSE
            for (u in uniq) {
                if (exists(u, envir = env, inherits = FALSE)) {
                    any_hit <- TRUE
                    break
                }
            }
            if (!any_hit) {
                if (all(uniq %in% allchroms)) {
                    return(chroms)
                }
            } else {
                # Map only the unique values, then expand back to chroms.
                u_targets <- vapply(uniq, .chrom_alias_lookup, character(1), env = env)
                # Build a name->target named character; missing entries stay NA.
                names(u_targets) <- uniq
                mapped <- u_targets[chroms]
                matched <- !is.na(mapped)
                if (any(matched)) {
                    chroms[matched] <- unname(mapped[matched])
                } else if (all(uniq %in% allchroms)) {
                    return(chroms)
                }
            }
        } else if (all(uniq %in% allchroms)) {
            return(chroms)
        }
    } else if (has_alias_map) {
        # Backward-compat path: no env mirror present (e.g. user assigned
        # CHROM_ALIAS manually). Use the named-vector fallback.
        alias_map <- get("CHROM_ALIAS", envir = .misha)
        if (length(alias_map)) {
            alias_names <- names(alias_map)
            if (length(alias_names) && !any(uniq %in% alias_names)) {
                if (all(uniq %in% allchroms)) {
                    return(chroms)
                }
            } else {
                mapped <- alias_map[chroms]
                matched <- !is.na(mapped)
                if (!any(matched) && all(uniq %in% allchroms)) {
                    return(chroms)
                }
                if (any(matched)) chroms[matched] <- mapped[matched]
            }
        } else if (all(uniq %in% allchroms)) {
            return(chroms)
        }
    } else if (all(uniq %in% allchroms)) {
        return(chroms)
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

# Build the full alias map (chr-prefix toggles + MT aliases + optional TSV
# extras like refseqAccession / genbankAccession / sequenceName / chrName).
# The result is consumed in two ways at .store_chrom_aliases time:
#   - The "basic" subset (everything EXCEPT the TSV extras) goes to
#     .misha$CHROM_ALIAS, which the C++ IntervUtils constructor iterates on
#     every .gcall(). Keeping this small bounds per-call overhead.
#   - The full map (including TSV extras) goes to .misha$CHROM_ALIAS_ENV, used
#     by the R-side .gchroms() for query-time alias resolution. The C++ side
#     never sees it.
# Self-mapping entries (chrom -> chrom) are deliberately omitted: canonical
# names already resolve via m_chrom_key / ALLGENOME, so storing them would
# only inflate per-call overhead.
.compute_chrom_aliases <- function(chroms, groot = NULL) {
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
    # Mark where the basic-toggle pass ends so .store_chrom_aliases can split
    # the result into the C++-facing CHROM_ALIAS (basic only) and the R-side
    # CHROM_ALIAS_ENV (full).
    n_basic <- 0L

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
    n_basic <- idx

    # Optional pass: load <groot>/chrom_aliases.tsv (written by gdb.build_genome
    # for NCBI-sourced genomes) and add accession / GenBank / sequenceName /
    # chrName aliases pointing to the canonical name. The TSV's `canonical`
    # column is the on-disk chrom name; aliases that already resolve via the
    # chr-prefix or mitochondrial passes above are skipped.
    if (!is.null(groot) && length(groot) == 1L && nzchar(groot)) {
        tsv_path <- file.path(groot, "chrom_aliases.tsv")
        if (file.exists(tsv_path)) {
            tsv <- tryCatch(
                utils::read.table(tsv_path,
                    sep = "\t", header = TRUE,
                    quote = "", comment.char = "",
                    stringsAsFactors = FALSE,
                    colClasses = "character",
                    na.strings = character(0)
                ),
                error = function(e) NULL
            )
            if (!is.null(tsv) && nrow(tsv) && "canonical" %in% names(tsv)) {
                # Only consider rows whose canonical actually appears among the
                # passed-in chroms; defends against stale TSVs from earlier
                # builds with different naming.
                keep <- tsv$canonical %in% chroms
                tsv <- tsv[keep, , drop = FALSE]
                # Columns whose values become aliases pointing at canonical.
                alias_cols <- intersect(
                    c("refseqAccession", "genbankAccession", "sequenceName", "chrName"),
                    names(tsv)
                )
                for (col in alias_cols) {
                    vals <- tsv[[col]]
                    targets <- tsv$canonical
                    for (j in seq_along(vals)) {
                        name <- vals[j]
                        if (is.na(name) || !nzchar(name)) next
                        if (exists(name, envir = seen, inherits = FALSE)) next
                        seen[[name]] <- TRUE
                        idx <- idx + 1
                        if (idx > length(alias_names)) {
                            # Grow buffers in chunks.
                            alias_names <- c(alias_names, character(length(alias_names)))
                            alias_targets <- c(alias_targets, character(length(alias_targets)))
                        }
                        alias_names[idx] <- name
                        alias_targets[idx] <- targets[j]
                    }
                }
            }
        }
    }

    # Build final named vector. seq_len(0) is integer(0), so an alias-free
    # build yields an empty named character (whereas 1:0 would yield c(1, 0)).
    # Attach `n_basic` as an attribute so .store_chrom_aliases can carve out
    # the C++-facing subset without re-running the algorithm.
    if (idx == 0L) {
        result <- stats::setNames(character(0), character(0))
    } else {
        result <- stats::setNames(alias_targets[seq_len(idx)], alias_names[seq_len(idx)])
    }
    attr(result, "n_basic") <- n_basic
    result
}

.store_chrom_aliases <- function(chroms, groot = NULL) {
    alias_map <- .compute_chrom_aliases(chroms, groot = groot)
    n_basic <- attr(alias_map, "n_basic") %||% length(alias_map)
    if (is.null(n_basic)) n_basic <- length(alias_map)

    # CHROM_ALIAS is what the C++ IntervUtils constructor reads on every
    # .gcall(). Iterating it costs ~2 hash ops per entry, so we only put the
    # basic toggles (chr-prefix + MT) here. TSV-derived extras stay in the
    # R-only env mirror and are invisible to the C++ side.
    if (n_basic > 0L && n_basic <= length(alias_map)) {
        cpp_map <- alias_map[seq_len(n_basic)]
    } else {
        cpp_map <- alias_map
    }
    attr(cpp_map, "n_basic") <- NULL
    assign("CHROM_ALIAS", cpp_map, envir = .misha)
    # The env mirror always carries the full set so .gchroms can resolve
    # accession / GenBank / sequenceName aliases at query time.
    .refresh_chrom_alias_env(full_alias_map = alias_map)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

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
