# Internal helpers for PSSM handling shared across sequence utilities and virtual tracks

.coerce_pssm_matrix <- function(pssm,
                                numeric_msg = "pssm must be a numeric matrix",
                                ncol_msg = "pssm must have columns named A, C, G, T",
                                colnames_msg = "pssm columns must be named A, C, G, T") {
    # Handle data frames by extracting required columns first
    if (is.data.frame(pssm)) {
        cols <- colnames(pssm)
        if (is.null(cols) || !all(c("A", "C", "G", "T") %in% cols)) {
            stop(colnames_msg)
        }
        # Extract only the required columns before converting to matrix
        # This avoids issues with non-numeric columns
        pssm <- as.matrix(pssm[, c("A", "C", "G", "T"), drop = FALSE])
    }

    if (!is.matrix(pssm) || !is.numeric(pssm)) {
        stop(numeric_msg)
    }

    cols <- colnames(pssm)
    if (is.null(cols) || !all(c("A", "C", "G", "T") %in% cols)) {
        stop(colnames_msg)
    }

    pssm[, c("A", "C", "G", "T"), drop = FALSE]
}
