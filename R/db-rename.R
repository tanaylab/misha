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
