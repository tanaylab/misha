# Internal helpers for export functions
.misha_readr_available <- function() {
    tryCatch(requireNamespace("readr", quietly = TRUE), error = function(...) FALSE)
}

.misha_write_with_readr <- function(...) {
    readr::write_tsv(...)
}
