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
                 call. = FALSE)
        }
        old <- nm
        new <- unname(mapping)
    } else {
        stop("mapping must be a data.frame(old, new) or a named character vector",
             call. = FALSE)
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
    if (nrow(active) == 0L) return(FALSE)
    any(active$new %in% active$old)
}

.misha_rename_remap_factor <- function(f, old, new) {
    if (is.null(f)) return(NULL)
    stopifnot(length(old) == length(new), is.factor(f))
    lv <- levels(f)
    idx <- match(lv, old)
    lv[!is.na(idx)] <- new[idx[!is.na(idx)]]
    levels(f) <- lv
    f
}

.misha_rename_remap_df <- function(df, old, new) {
    if (!is.data.frame(df)) return(df)
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
    writer(tmp_path)
    if (!file.rename(tmp_path, target)) {
        stop(sprintf("failed to rename %s to %s", tmp_path, target), call. = FALSE)
    }
    success <- TRUE
    invisible(target)
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
        to   <- file.path(seq_dir, paste0(mapping$new, ".seq"))
        keep <- file.exists(from)
        plan$seq_renames <- from[keep]
        plan$seq_renames_new <- to[keep]
    }

    tracks_root <- file.path(groot, "tracks")
    all_track_dirs <- list.files(tracks_root, pattern = "\\.track$", recursive = TRUE,
                                 full.names = TRUE, include.dirs = TRUE)
    all_interv    <- list.files(tracks_root, pattern = "\\.interv$", recursive = TRUE,
                                full.names = TRUE, include.dirs = TRUE)
    track_dirs  <- all_track_dirs[file.info(all_track_dirs)$isdir %in% TRUE]
    interv_dirs <- all_interv[file.info(all_interv)$isdir %in% TRUE]
    single_interv <- all_interv[!(file.info(all_interv)$isdir %in% TRUE)]

    if (!is_indexed) {
        per_chrom_files_for_dir <- function(d) {
            files <- list.files(d, full.names = FALSE)
            files <- files[!files %in% c(".meta",
                                         "track.dat", "track.idx",
                                         "intervals.dat", "intervals.idx",
                                         "intervals2d.dat", "intervals2d.idx")]
            old_names <- files
            new_names <- files
            idx <- match(old_names, mapping$old)
            new_names[!is.na(idx)] <- mapping$new[idx[!is.na(idx)]]
            is_pair <- grepl("-", files, fixed = TRUE)
            if (any(is_pair)) {
                parts <- strsplit(files[is_pair], "-", fixed = TRUE)
                remapped <- vapply(parts, function(p) {
                    if (length(p) != 2) return(paste(p, collapse = "-"))
                    i1 <- match(p[1], mapping$old)
                    i2 <- match(p[2], mapping$old)
                    if (!is.na(i1)) p[1] <- mapping$new[i1]
                    if (!is.na(i2)) p[2] <- mapping$new[i2]
                    paste(p, collapse = "-")
                }, character(1))
                new_names[is_pair] <- remapped
            }
            changed <- old_names != new_names
            if (!any(changed)) return(NULL)
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
