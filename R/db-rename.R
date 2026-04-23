# Internal helpers for gdb.rename_chroms.

.misha_rename_normalize_mapping <- function(mapping) {
    if (is.data.frame(mapping)) {
        missing_cols <- setdiff(c("old", "new"), colnames(mapping))
        if (length(missing_cols)) {
            stop(sprintf(
                "mapping data.frame is missing columns: %s",
                paste(missing_cols, collapse = ", ")
            ), call. = FALSE)
        }
        old <- as.character(mapping$old)
        new <- as.character(mapping$new)
    } else if (is.character(mapping)) {
        nm <- names(mapping)
        if (is.null(nm) || any(is.na(nm)) || any(!nzchar(nm))) {
            stop("mapping must be a named character vector (names = old chroms)",
                call. = FALSE
            )
        }
        old <- nm
        new <- unname(mapping)
    } else {
        stop("mapping must be a data.frame(old, new) or a named character vector",
            call. = FALSE
        )
    }

    data.frame(old = old, new = new, stringsAsFactors = FALSE)
}

.misha_rename_validate_mapping <- function(mapping, existing) {
    if (nrow(mapping) == 0L) {
        stop("mapping is empty", call. = FALSE)
    }

    bad <- is.na(mapping$old) | !nzchar(mapping$old) |
        is.na(mapping$new) | !nzchar(mapping$new)
    if (any(bad)) {
        stop("mapping contains NA or empty chromosome names", call. = FALSE)
    }

    dup_old <- mapping$old[duplicated(mapping$old)]
    if (length(dup_old)) {
        stop(sprintf(
            "duplicate old chromosome names in mapping: %s",
            paste(unique(dup_old), collapse = ", ")
        ), call. = FALSE)
    }

    dup_new <- mapping$new[duplicated(mapping$new)]
    if (length(dup_new)) {
        stop(sprintf(
            "duplicate new chromosome names in mapping: %s",
            paste(unique(dup_new), collapse = ", ")
        ), call. = FALSE)
    }

    unknown <- setdiff(mapping$old, existing)
    if (length(unknown)) {
        stop(sprintf(
            "old chromosome(s) not present in the database: %s",
            paste(unknown, collapse = ", ")
        ), call. = FALSE)
    }

    collisions <- intersect(mapping$new, existing)
    collisions <- setdiff(collisions, mapping$old)
    if (length(collisions)) {
        stop(sprintf(
            "new chromosome name(s) collides with existing un-mapped chromosomes: %s",
            paste(collisions, collapse = ", ")
        ), call. = FALSE)
    }

    invisible(NULL)
}

.misha_rename_needs_two_phase <- function(mapping) {
    active <- mapping[mapping$old != mapping$new, , drop = FALSE]
    if (nrow(active) == 0L) {
        return(FALSE)
    }
    any(active$new %in% active$old)
}

.misha_rename_remap_factor <- function(f, old, new) {
    if (is.null(f)) {
        return(NULL)
    }
    stopifnot(length(old) == length(new), is.factor(f))
    lv <- levels(f)
    idx <- match(lv, old)
    lv[!is.na(idx)] <- new[idx[!is.na(idx)]]
    levels(f) <- lv
    f
}

.misha_rename_remap_df <- function(df, old, new) {
    if (!is.data.frame(df)) {
        return(df)
    }
    for (col in c("chrom", "chrom1", "chrom2")) {
        if (col %in% colnames(df)) {
            df[[col]] <- .misha_rename_remap_factor(df[[col]], old = old, new = new)
        }
    }
    df
}

# Execute `writer(tmp_path)` to produce a new version of `target`, then
# atomically rename the temp file over the target. On writer error, the
# temp file is removed and the original target is left untouched.
.misha_rename_atomic_rewrite <- function(target, writer) {
    tmp_path <- paste0(target, ".tmp.", Sys.getpid(), ".", as.integer(Sys.time()))
    success <- FALSE
    on.exit(if (!success && file.exists(tmp_path)) unlink(tmp_path), add = TRUE)
    # file(..., "wb") signals EACCES as a warning and returns an invalid
    # connection; promote that to an error so the caller aborts cleanly instead
    # of producing a half-written temp file.
    withCallingHandlers(
        writer(tmp_path),
        warning = function(w) {
            if (grepl("cannot open file|Permission denied", conditionMessage(w))) {
                stop(conditionMessage(w), call. = FALSE)
            }
        }
    )
    if (!file.rename(tmp_path, target)) {
        stop(sprintf("failed to rename %s to %s", tmp_path, target), call. = FALSE)
    }
    success <- TRUE
    invisible(target)
}

# Enumerate every directory whose writability the rename requires.
# Rewrites go through a temp-file + rename pattern, so the CONTAINING directory
# must allow file creation (not just modification of the existing file).
.misha_rename_required_dirs <- function(groot, plan) {
    dirs <- groot  # chrom_sizes.txt rewrite + .rename_interrupted breadcrumb
    if (plan$is_indexed) {
        dirs <- c(dirs, dirname(plan$genome_idx_path))
    } else if (length(plan$seq_renames)) {
        dirs <- c(dirs, unique(dirname(plan$seq_renames_new)))
    }
    dirs <- c(dirs, names(plan$track_dir_renames))
    dirs <- c(dirs, names(plan$interv_dir_renames))
    if (length(plan$meta_rewrites)) {
        dirs <- c(dirs, unique(dirname(plan$meta_rewrites)))
    }
    if (length(plan$single_interv_rewrites)) {
        dirs <- c(dirs, unique(dirname(plan$single_interv_rewrites)))
    }
    unique(dirs)
}

# Probe each directory by creating and removing a zero-byte temp file.
# file.access() is unreliable on NFS/Lustre, so exercise the filesystem.
.misha_rename_check_writable <- function(dirs) {
    bad <- character(0)
    for (d in dirs) {
        probe <- tryCatch(tempfile(tmpdir = d), error = function(e) NULL)
        if (is.null(probe)) {
            bad <- c(bad, d)
            next
        }
        ok <- tryCatch(
            withCallingHandlers(
                {
                    fh <- file(probe, "wb")
                    close(fh)
                    unlink(probe)
                    TRUE
                },
                warning = function(w) {
                    if (grepl("cannot open file|Permission denied", conditionMessage(w))) {
                        stop(conditionMessage(w), call. = FALSE)
                    }
                }
            ),
            error = function(e) FALSE
        )
        if (!isTRUE(ok)) bad <- c(bad, d)
    }
    if (length(bad)) {
        stop(sprintf(
            "gdb.rename_chroms requires write permission on every directory containing a file it rewrites. The following %s not writable:\n  %s",
            if (length(bad) == 1) "directory is" else "directories are",
            paste(bad, collapse = "\n  ")
        ), call. = FALSE)
    }
    invisible(NULL)
}

.misha_rename_rewrite_meta <- function(meta_path, old, new) {
    f <- file(meta_path, "rb")
    meta <- unserialize(f)
    close(f)

    meta$stats <- .misha_rename_remap_df(meta$stats, old = old, new = new)
    meta$zeroline <- .misha_rename_remap_df(meta$zeroline, old = old, new = new)

    .misha_rename_atomic_rewrite(meta_path, function(tmp_path) {
        f <- file(tmp_path, "wb")
        serialize(meta, f)
        close(f)
    })
}

.misha_rename_rewrite_single_interv <- function(path, old, new) {
    f <- file(path, "rb")
    df <- unserialize(f)
    close(f)

    df <- .misha_rename_remap_df(df, old = old, new = new)

    .misha_rename_atomic_rewrite(path, function(tmp_path) {
        f <- file(tmp_path, "wb")
        serialize(df, f)
        close(f)
    })
}

# Build a rename plan describing every filesystem mutation required.
# Returns a list:
#   is_indexed              : TRUE/FALSE (genome.idx + genome.seq present?)
#   seq_renames             : character vector of old per-chrom seq paths (per-chrom DB only)
#   seq_renames_new         : character vector of new per-chrom seq paths (parallel)
#   genome_idx_path         : path or NA
#   track_dir_renames       : named list of data.frame(old, new) per track dir (per-chrom DB only)
#   interv_dir_renames      : named list of data.frame(old, new) per bigset dir (per-chrom DB only)
#   meta_rewrites           : character vector of .meta file paths
#   single_interv_rewrites  : character vector of single-file .interv paths
.misha_rename_build_plan <- function(groot, mapping) {
    seq_dir <- file.path(groot, "seq")
    is_indexed <- file.exists(file.path(seq_dir, "genome.idx")) &&
        file.exists(file.path(seq_dir, "genome.seq"))

    plan <- list(
        is_indexed = is_indexed,
        seq_renames = character(0),
        seq_renames_new = character(0),
        genome_idx_path = if (is_indexed) file.path(seq_dir, "genome.idx") else NA_character_,
        track_dir_renames = list(),
        interv_dir_renames = list(),
        meta_rewrites = character(0),
        single_interv_rewrites = character(0)
    )

    if (!is_indexed) {
        from <- file.path(seq_dir, paste0(mapping$old, ".seq"))
        to <- file.path(seq_dir, paste0(mapping$new, ".seq"))
        missing <- !file.exists(from)
        if (any(missing)) {
            stop(sprintf(
                "Expected per-chromosome sequence file(s) not found: %s. The DB layout may be inconsistent.",
                paste(basename(from[missing]), collapse = ", ")
            ), call. = FALSE)
        }
        plan$seq_renames <- from
        plan$seq_renames_new <- to
    }

    tracks_root <- file.path(groot, "tracks")
    all_track_dirs <- list.files(tracks_root,
        pattern = "\\.track$", recursive = TRUE,
        full.names = TRUE, include.dirs = TRUE
    )
    all_interv <- list.files(tracks_root,
        pattern = "\\.interv$", recursive = TRUE,
        full.names = TRUE, include.dirs = TRUE
    )
    track_info <- file.info(all_track_dirs)$isdir %in% TRUE
    interv_info <- file.info(all_interv)$isdir %in% TRUE
    track_dirs <- all_track_dirs[track_info]
    interv_dirs <- all_interv[interv_info]
    single_interv <- all_interv[!interv_info]

    if (!is_indexed) {
        per_chrom_files_for_dir <- function(d) {
            files <- list.files(d, full.names = FALSE)
            files <- files[!files %in% c(
                ".meta",
                "track.dat", "track.idx",
                "intervals.dat", "intervals.idx",
                "intervals2d.dat", "intervals2d.idx"
            )]
            old_names <- files
            new_names <- files
            idx <- match(old_names, mapping$old)
            new_names[!is.na(idx)] <- mapping$new[idx[!is.na(idx)]]
            is_pair <- grepl("-", files, fixed = TRUE)
            if (any(is_pair)) {
                parts <- strsplit(files[is_pair], "-", fixed = TRUE)
                remapped <- vapply(parts, function(p) {
                    if (length(p) != 2) {
                        return(paste(p, collapse = "-"))
                    }
                    i1 <- match(p[1], mapping$old)
                    i2 <- match(p[2], mapping$old)
                    if (!is.na(i1)) p[1] <- mapping$new[i1]
                    if (!is.na(i2)) p[2] <- mapping$new[i2]
                    paste(p, collapse = "-")
                }, character(1))
                new_names[is_pair] <- remapped
            }
            changed <- old_names != new_names
            if (!any(changed)) {
                return(NULL)
            }
            data.frame(
                old = file.path(d, old_names[changed]),
                new = file.path(d, new_names[changed]),
                stringsAsFactors = FALSE
            )
        }

        for (d in track_dirs) {
            r <- per_chrom_files_for_dir(d)
            if (!is.null(r)) plan$track_dir_renames[[d]] <- r
        }
        for (d in interv_dirs) {
            r <- per_chrom_files_for_dir(d)
            if (!is.null(r)) plan$interv_dir_renames[[d]] <- r
        }
    }

    plan$meta_rewrites <- c(
        file.path(track_dirs, ".meta"),
        file.path(interv_dirs, ".meta")
    )
    plan$meta_rewrites <- plan$meta_rewrites[file.exists(plan$meta_rewrites)]

    plan$single_interv_rewrites <- single_interv

    plan
}

# Execute a rename plan. Assumes validation has already been performed.
# Does NOT handle two-phase (swap) cases; caller must stage swaps first.
.misha_rename_apply_plan <- function(plan, mapping, verbose = FALSE) {
    # Phase 1: sequence data
    if (plan$is_indexed) {
        if (verbose) message("Rewriting genome.idx...")
        .gcall(
            "C_gdb_rewrite_genome_idx",
            plan$genome_idx_path, mapping$old, mapping$new, .misha_env()
        )
    } else if (length(plan$seq_renames)) {
        if (verbose) message(sprintf("Renaming %d seq files...", length(plan$seq_renames)))
        ok <- file.rename(plan$seq_renames, plan$seq_renames_new)
        if (!all(ok)) stop("one or more seq file renames failed", call. = FALSE)
    }

    # Phase 2: per-chromosome track/interv files (only per-chrom DBs)
    if (!plan$is_indexed) {
        for (d in names(plan$track_dir_renames)) {
            r <- plan$track_dir_renames[[d]]
            if (verbose) message(sprintf("Track dir %s: renaming %d files", d, nrow(r)))
            ok <- file.rename(r$old, r$new)
            if (!all(ok)) stop(sprintf("rename failed in %s", d), call. = FALSE)
        }
        for (d in names(plan$interv_dir_renames)) {
            r <- plan$interv_dir_renames[[d]]
            if (verbose) message(sprintf("Interval dir %s: renaming %d files", d, nrow(r)))
            ok <- file.rename(r$old, r$new)
            if (!all(ok)) stop(sprintf("rename failed in %s", d), call. = FALSE)
        }
    }

    # Phase 3: meta rewrites
    for (m in plan$meta_rewrites) {
        if (verbose) message(sprintf("Rewriting %s", m))
        .misha_rename_rewrite_meta(m, old = mapping$old, new = mapping$new)
    }

    # Phase 4: single-file .interv rewrites
    for (p in plan$single_interv_rewrites) {
        if (verbose) message(sprintf("Rewriting %s", p))
        .misha_rename_rewrite_single_interv(p, old = mapping$old, new = mapping$new)
    }

    invisible(NULL)
}

# Apply a plan that may include swaps. If swaps are detected, split the
# mapping into two stages -- (old -> temp) and (temp -> new) -- and re-walk
# the DB once per stage (one extra enumeration is acceptable cost for
# the safety of two-phase rename).
.misha_rename_apply_with_swap <- function(groot, plan, mapping, verbose = FALSE) {
    if (!.misha_rename_needs_two_phase(mapping)) {
        .misha_rename_apply_plan(plan, mapping, verbose = verbose)
        return(invisible(NULL))
    }

    tag <- sprintf(
        "__misha_rename_%s_%d_",
        format(Sys.time(), "%Y%m%d%H%M%S"), Sys.getpid()
    )
    tmp_names <- paste0(tag, mapping$old)

    stage1 <- data.frame(old = mapping$old, new = tmp_names, stringsAsFactors = FALSE)
    plan1 <- .misha_rename_build_plan(groot, stage1)
    .misha_rename_apply_plan(plan1, stage1, verbose = verbose)

    stage2 <- data.frame(old = tmp_names, new = mapping$new, stringsAsFactors = FALSE)
    plan2 <- .misha_rename_build_plan(groot, stage2)
    .misha_rename_apply_plan(plan2, stage2, verbose = verbose)

    invisible(NULL)
}

#' Rename chromosomes in a misha database
#'
#' Performs an in-place rename of chromosomes across an existing misha
#' database. Works on both indexed (\code{seq/genome.idx + genome.seq}) and
#' per-chromosome (\code{seq/*.seq}) database formats.
#'
#' The rename is bijective. Partial mappings are allowed: chromosomes not
#' mentioned in \code{mapping} keep their current names. Swaps and
#' permutations are supported via an internal two-phase rename.
#'
#' The operation rewrites every file that references chromosome names:
#' \code{chrom_sizes.txt}, per-chromosome sequence/track/interval files
#' (per-chromosome DBs only), \code{seq/genome.idx} (indexed DBs only),
#' bigset \code{.meta} files, and single-file \code{.interv} sets. Track
#' and interval \code{.idx}/\code{.dat} files reference chromosomes only
#' via numeric IDs and are left untouched.
#'
#' \strong{Atomicity.} Individual file operations are atomic (temp file +
#' rename). The operation as a whole is \emph{not} atomic across the
#' database. If execution is interrupted, a \code{.rename_interrupted}
#' breadcrumb is left at the database root to help recovery. Back up the
#' database (or at least \code{chrom_sizes.txt}, \code{seq/genome.idx},
#' and all \code{.meta} files) before running this on critical databases.
#'
#' @param groot path to database root. If \code{NULL}, the currently active
#'   database is used.
#' @param mapping a \code{data.frame} with columns \code{old} and \code{new},
#'   or a named character vector (names are old chromosome names, values are
#'   new names). Old names must match the canonical chromosome names as
#'   exposed by \code{ALLGENOME} (for per-chromosome DBs this is the
#'   \code{chr}-prefixed form; \code{chrom_sizes.txt} may store a stripped
#'   form, which is preserved on rewrite).
#' @param force logical; if \code{TRUE}, skip the interactive confirmation.
#' @param dry_run logical; if \code{TRUE}, print the plan and return without
#'   modifying any file.
#' @param verbose logical; if \code{TRUE}, print progress per phase.
#'
#' @return invisible \code{NULL}.
#' @export
gdb.rename_chroms <- function(groot = NULL, mapping = NULL,
                              force = FALSE, dry_run = FALSE, verbose = FALSE) {
    if (is.null(mapping)) {
        stop("Usage: gdb.rename_chroms(groot, mapping, force = FALSE, dry_run = FALSE, verbose = FALSE)",
            call. = FALSE
        )
    }

    was_loaded <- FALSE
    if (is.null(groot)) {
        if (!exists("GROOT", envir = .misha) || is.null(get("GROOT", envir = .misha))) {
            stop("No database is currently active. Provide `groot` or call gdb.init() first.",
                call. = FALSE
            )
        }
        groot <- get("GROOT", envir = .misha)
        was_loaded <- TRUE
    } else {
        if (exists("GROOT", envir = .misha, inherits = FALSE)) {
            cur <- get("GROOT", envir = .misha)
            if (!is.null(cur) && nzchar(cur) &&
                normalizePath(cur, mustWork = FALSE) ==
                    normalizePath(groot, mustWork = FALSE)) {
                was_loaded <- TRUE
            }
        }
    }

    if (!dir.exists(groot)) {
        stop(sprintf("Database directory does not exist: %s", groot), call. = FALSE)
    }

    chrom_sizes_path <- file.path(groot, "chrom_sizes.txt")
    if (!file.exists(chrom_sizes_path)) {
        stop(sprintf("chrom_sizes.txt not found: %s", chrom_sizes_path), call. = FALSE)
    }
    cs <- utils::read.table(chrom_sizes_path,
        sep = "\t", stringsAsFactors = FALSE,
        col.names = c("chrom", "size")
    )

    # Ensure the target DB is loaded so we can use ALLGENOME for validation
    # and name resolution. Capture prior state for restoration.
    prior_groot <- NULL
    if (exists("GROOT", envir = .misha, inherits = FALSE)) {
        pg <- get("GROOT", envir = .misha)
        if (!is.null(pg) && nzchar(pg)) prior_groot <- pg
    }

    if (!was_loaded) {
        suppressMessages(gdb.init(groot))
    }

    allgenome <- get("ALLGENOME", envir = .misha)[[1]]
    existing <- as.character(allgenome$chrom)

    mapping <- .misha_rename_normalize_mapping(mapping)
    .misha_rename_validate_mapping(mapping, existing = existing)

    breadcrumb <- file.path(groot, ".rename_interrupted")
    if (file.exists(breadcrumb) && !force) {
        stop(sprintf(
            "Found %s from a previous interrupted rename. Inspect and remove it, or call with force = TRUE.",
            breadcrumb
        ), call. = FALSE)
    }

    plan <- .misha_rename_build_plan(groot, mapping)

    # Pre-flight: fail fast if we can't write to any directory in the plan.
    # Without this, a partway-through permission error leaves the DB in a
    # half-renamed state.
    .misha_rename_check_writable(.misha_rename_required_dirs(groot, plan))

    n_seq <- length(plan$seq_renames) + (if (plan$is_indexed) 1L else 0L)
    n_track <- sum(vapply(plan$track_dir_renames, nrow, integer(1)))
    n_interv <- sum(vapply(plan$interv_dir_renames, nrow, integer(1)))
    n_meta <- length(plan$meta_rewrites)
    n_single <- length(plan$single_interv_rewrites)

    summary_str <- sprintf(
        "Database: %s\nFormat: %s\nRenames: %d sequence(s), %d track files, %d interval files\nRewrites: %d .meta, %d single-file .interv\nChromosomes affected: %d",
        groot,
        if (plan$is_indexed) "indexed" else "per-chromosome",
        n_seq, n_track, n_interv, n_meta, n_single,
        nrow(mapping)
    )

    if (dry_run) {
        cat(summary_str, "\n", sep = "")
        cat("(dry run -- no changes made)\n")
        return(invisible(NULL))
    }

    if (interactive() && !force) {
        cat(summary_str, "\n", sep = "")
        resp <- readline("Proceed with rename? (yes/no): ")
        if (!(tolower(resp) %in% c("y", "yes"))) {
            message("Cancelled.")
            return(invisible(NULL))
        }
    }

    f <- file(breadcrumb, "wb")
    serialize(list(mapping = mapping, timestamp = Sys.time()), f)
    close(f)

    saved_groot <- if (was_loaded) groot else NULL

    success <- FALSE
    tryCatch({
        .misha_rename_apply_with_swap(groot, plan, mapping, verbose = verbose)

        # Rewrite chrom_sizes.txt while preserving its "form" (stripped vs
        # prefixed). mapping$old / mapping$new are in ALLGENOME canonical form.
        # For each row, detect how cs$chrom relates to the ALLGENOME name and
        # apply the same transformation to the new name.
        cs_new <- cs
        for (i in seq_len(nrow(cs))) {
            cs_name <- cs$chrom[i]
            # Find this row in ALLGENOME by matching strip/prefix variants.
            if (cs_name %in% existing) {
                # Same form (cs uses prefix, ALLGENOME uses prefix).
                allgenome_name <- cs_name
                transform <- identity
            } else if (paste0("chr", cs_name) %in% existing) {
                # cs strips "chr" prefix; ALLGENOME has it.
                allgenome_name <- paste0("chr", cs_name)
                transform <- function(x) sub("^chr", "", x)
            } else {
                next # unknown row -- leave as-is
            }
            mi <- match(allgenome_name, mapping$old)
            if (!is.na(mi)) {
                cs_new$chrom[i] <- transform(mapping$new[mi])
            }
        }
        .misha_rename_atomic_rewrite(chrom_sizes_path, function(tmp_path) {
            utils::write.table(cs_new, tmp_path,
                sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE
            )
        })

        # Invalidate .db.cache so next gdb.init() rescans.
        cache_path <- file.path(groot, ".db.cache")
        if (file.exists(cache_path)) unlink(cache_path)

        success <- TRUE
    }, finally = {
        if (success) {
            unlink(breadcrumb)
        }
        if (was_loaded) {
            # Reload whether we succeeded or failed -- in-memory state is stale.
            suppressMessages(try(gdb.init(saved_groot), silent = TRUE))
        } else if (!is.null(prior_groot) &&
            normalizePath(prior_groot, mustWork = FALSE) !=
                normalizePath(groot, mustWork = FALSE)) {
            # We loaded the DB temporarily; restore whatever was active before.
            suppressMessages(try(gdb.init(prior_groot), silent = TRUE))
        } else {
            # Nothing was loaded before (or prior was this same DB) and we
            # loaded it during this call. Re-init so in-memory state reflects
            # the rename.
            suppressMessages(try(gdb.init(groot), silent = TRUE))
        }
    })

    invisible(NULL)
}
