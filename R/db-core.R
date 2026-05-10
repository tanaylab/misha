# Core chromosome and genome functions

# Build (or refresh) an environment-based mirror of CHROM_ALIAS for O(1)
# lookup. R's named-vector subscript path uses match(), which rebuilds an
# internal hash on every call — O(M) per call, where M is the alias count.
# For large alias maps (e.g. NCBI assemblies with chrom_aliases.tsv loaded:
# ~12k entries on rabbit), this becomes the dominant cost in .gchroms().
# An environment offers genuine O(1) hashed lookup that survives across calls.
.refresh_chrom_alias_env <- function(full_alias_map = NULL,
                                     lazy_threshold = getOption("misha.chrom_alias_lazy_threshold", 1e6)) {
    # Drop any existing mirror first.
    if (exists("CHROM_ALIAS_ENV", envir = .misha, inherits = FALSE)) {
        rm("CHROM_ALIAS_ENV", envir = .misha)
    }
    # Drop any stale pending map; caller will set a fresh one if applicable.
    if (exists("CHROM_ALIAS_PENDING", envir = .misha, inherits = FALSE)) {
        rm("CHROM_ALIAS_PENDING", envir = .misha)
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
    # Defer the env build for very large maps. list2env(as.list(named_char))
    # is ~10us per entry; on the Phylo447 / Crocidura assembly's 7.2M-entry
    # alias map that's ~80s of gsetroot stall. The full map is stashed in
    # CHROM_ALIAS_PENDING and .gchroms() materializes the env on the first
    # query that actually requires alias resolution (i.e. has a non-canonical
    # name). Canonical-only workflows never pay the cost.
    if (length(alias_map) > lazy_threshold) {
        assign("CHROM_ALIAS_PENDING", alias_map, envir = .misha)
        return(invisible(NULL))
    }
    env <- new.env(
        hash = TRUE, parent = emptyenv(),
        size = max(length(alias_map), 1L)
    )
    if (length(alias_map)) {
        list2env(as.list(alias_map), envir = env)
    }
    assign("CHROM_ALIAS_ENV", env, envir = .misha)
    invisible(NULL)
}

# Materialize the alias env from CHROM_ALIAS_PENDING. Called by .gchroms() the
# first time a query needs alias resolution on a DB whose map was deferred at
# gsetroot time.
.materialize_chrom_alias_env <- function() {
    if (!exists("CHROM_ALIAS_PENDING", envir = .misha, inherits = FALSE)) {
        return(invisible(NULL))
    }
    full_map <- get("CHROM_ALIAS_PENDING", envir = .misha)
    env <- new.env(
        hash = TRUE, parent = emptyenv(),
        size = max(length(full_map), 1L)
    )
    if (length(full_map)) {
        list2env(as.list(full_map), envir = env)
    }
    assign("CHROM_ALIAS_ENV", env, envir = .misha)
    rm("CHROM_ALIAS_PENDING", envir = .misha)
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

    # Lazy alias-env materialization: if .refresh_chrom_alias_env deferred the
    # build (huge alias map, e.g. Phylo447), build it now -- but only when the
    # query actually needs alias resolution. Canonical-only queries take a
    # direct match() shortcut and never pay either the env build or the
    # named-vector hash on CHROM_ALIAS.
    if (!has_alias_env &&
        exists("CHROM_ALIAS_PENDING", envir = .misha, inherits = FALSE)) {
        indices <- match(chroms, allchroms)
        if (!anyNA(indices)) {
            # All canonical; allchroms[] yields factor-typed output, matching
            # the function's "alias-resolution" tail return semantics.
            return(allchroms[indices])
        }
        .materialize_chrom_alias_env()
        has_alias_env <- exists("CHROM_ALIAS_ENV", envir = .misha, inherits = FALSE)
    }

    if (has_alias_env) {
        env <- get("CHROM_ALIAS_ENV", envir = .misha)
        # Fast probe: do any unique chroms hit the alias env at all?
        # (We deliberately don't gate on length(env) here -- length() on a
        # 7M-entry env is ~700 ms, while a per-uniq exists() probe on the same
        # env is ~5 us.)
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
            names(u_targets) <- uniq
            mapped <- u_targets[chroms]
            matched <- !is.na(mapped)
            if (any(matched)) {
                chroms[matched] <- unname(mapped[matched])
            } else if (all(uniq %in% allchroms)) {
                return(chroms)
            }
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
# Self-mapping entries (chrom -> chrom) are kept in the alias map so that
# .gchroms() always takes the alias-resolution path for canonical names; that
# path returns chroms typed by ALLGENOME's factor, while the fast "all
# canonical" early-return would otherwise hand back the input character
# vector and silently strip factor-ness from gextract output. The C++
# add_chrom_alias no-ops on entries whose alias already names a chrom, so
# self-entries cost only a no-op lookup on the C++ side.
.compute_chrom_aliases <- function(chroms, groot = NULL) {
    chroms <- unique(as.character(chroms))
    n <- length(chroms)
    if (!n) {
        result <- stats::setNames(character(0), character(0))
        attr(result, "n_basic") <- 0L
        return(result)
    }

    has_chr_prefix <- startsWith(chroms, "chr")
    unprefixed <- ifelse(has_chr_prefix, substring(chroms, 4), chroms)
    prefixed <- ifelse(has_chr_prefix, chroms, paste0("chr", chroms))
    upper_unpref <- toupper(unprefixed)
    upper_chr <- toupper(chroms)
    is_mito <- upper_unpref %in% c("M", "MT") | upper_chr %in% c("CHRM", "CHRMT")

    # Vectorized basic pass: self-aliases, chr-prefix toggles, and MT aliases.
    # Replaces a per-element loop over .compute_chrom_aliases's `seen` env, which
    # took ~90 s on a 2.4M-contig fragmented assembly.
    #
    # On fragmented assemblies (e.g. Phylo447, 2.4M scaffold_*) the prefixed-
    # alias step generates millions of "chrXXX" entries that no user will ever
    # look up. Skip it when (a) contig count is large, (b) no canonical name
    # has the "chr" prefix, and (c) no name is mito-pattern. Ensembl-style
    # references (1/2/X/MT) and any genome with mito stay unaffected.
    skip_prefixed <- n > 1000L && !any(has_chr_prefix) && !any(is_mito)

    if (skip_prefixed) {
        basic_names <- chroms
        basic_targets <- chroms
    } else {
        # Candidate chr-toggle aliases tagged with source-chrom index, so we can
        # reproduce the original loop's "first-seen wins" semantics via
        # order + !duplicated.
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

        basic_names <- c(chroms, cand_names)
        basic_targets <- c(chroms, cand_targets)

        # Mitochondrial: M, MT, chrM all point to the first mito-like canonical.
        # Add only if not already present in basic.
        if (any(is_mito)) {
            mt_target <- chroms[which(is_mito)[1]]
            for (mt in c("M", "MT", "chrM")) {
                if (!(mt %in% basic_names)) {
                    basic_names <- c(basic_names, mt)
                    basic_targets <- c(basic_targets, mt_target)
                }
            }
        }
    }

    n_basic <- length(basic_names)

    # No TSV -> short-circuit before read.table (which would otherwise produce
    # a confusing relative-path read on NULL/empty groot).
    no_tsv <- is.null(groot) || length(groot) != 1L || !nzchar(groot) ||
        !file.exists(file.path(groot, "chrom_aliases.tsv"))
    if (no_tsv) {
        result <- stats::setNames(basic_targets, basic_names)
        attr(result, "n_basic") <- n_basic
        return(result)
    }

    # Optional pass: load <groot>/chrom_aliases.tsv and add aliases pointing at
    # canonical chrom names. Two on-disk shapes are handled:
    #   - Wide (NCBI builds): canonical, refseqAccession, genbankAccession,
    #     sequenceName, chrName, role, length.
    #   - Long (ucsc-hub builds): canonical, alias, source -- one row per
    #     (canonical, alias) pair.
    # The whole pass is vectorized: a per-row add_alias closure with <<- env
    # and vector grows took ~260 s on a 7.2M-row TSV (Phylo447 / Crocidura).
    tsv_path <- file.path(groot, "chrom_aliases.tsv")
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

    extra_names <- character(0)
    extra_targets <- character(0)

    if (!is.null(tsv) && nrow(tsv) && "canonical" %in% names(tsv)) {
        # Defensive: ignore TSV rows whose canonical is no longer among `chroms`
        # (e.g. stale alias file referring to renamed contigs).
        keep_canon <- tsv$canonical %in% chroms
        tsv <- tsv[keep_canon, , drop = FALSE]

        if (nrow(tsv)) {
            if (all(c("alias", "source") %in% names(tsv))) {
                # Long format: each row is one (alias -> canonical) mapping.
                extra_names <- tsv$alias
                extra_targets <- tsv$canonical
            } else {
                # Wide format: stack alias-bearing columns in priority order so
                # earlier columns win first-seen-wins below.
                alias_cols <- intersect(
                    c("refseqAccession", "genbankAccession", "sequenceName", "chrName"),
                    names(tsv)
                )
                for (col in alias_cols) {
                    extra_names <- c(extra_names, tsv[[col]])
                    extra_targets <- c(extra_targets, tsv$canonical)
                }
            }
        }
    }

    if (length(extra_names)) {
        # Drop empty / NA aliases.
        valid <- !is.na(extra_names) & nzchar(extra_names)
        extra_names <- extra_names[valid]
        extra_targets <- extra_targets[valid]

        # Basic entries always win: drop TSV aliases that collide with basic.
        not_basic <- !(extra_names %in% basic_names)
        extra_names <- extra_names[not_basic]
        extra_targets <- extra_targets[not_basic]

        # First-seen wins within the TSV pass (matches the original loop's
        # `seen` env semantics: earlier row / earlier column wins on collisions).
        keep <- !duplicated(extra_names)
        extra_names <- extra_names[keep]
        extra_targets <- extra_targets[keep]
    }

    # n_basic is the length of the "basic" prefix; .store_chrom_aliases uses
    # this to carve out the C++-facing CHROM_ALIAS without re-running.
    result <- stats::setNames(
        c(basic_targets, extra_targets),
        c(basic_names, extra_names)
    )
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
