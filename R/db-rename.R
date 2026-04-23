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
        if (is.null(names(mapping)) || any(!nzchar(names(mapping)))) {
            stop("mapping must be a named character vector (names = old chroms)",
                 call. = FALSE)
        }
        old <- names(mapping)
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
