# Internal helpers for PSSM handling shared across sequence utilities and virtual tracks

.coerce_pssm_matrix <- function(pssm,
                                numeric_msg = "pssm must be a numeric matrix",
                                ncol_msg = "pssm must have exactly 4 columns",
                                colnames_msg = "pssm columns must be named A, C, G, T") {
    if (is.data.frame(pssm)) {
        pssm <- as.matrix(pssm)
    }

    if (!is.matrix(pssm) || !is.numeric(pssm)) {
        stop(numeric_msg)
    }

    if (ncol(pssm) != 4L) {
        stop(ncol_msg)
    }

    cols <- colnames(pssm)
    if (is.null(cols) || !setequal(cols, c("A", "C", "G", "T"))) {
        stop(colnames_msg)
    }

    pssm[, c("A", "C", "G", "T"), drop = FALSE]
}
