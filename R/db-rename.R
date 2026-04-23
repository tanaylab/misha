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
