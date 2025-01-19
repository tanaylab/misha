#' Get reverse complement of DNA sequence
#'
#' Takes a DNA sequence string and returns its reverse complement.
#'
#' @param seq A character vector containing DNA sequences (using A,C,G,T). Ignores other characters and NA values.
#' @return A character vector of the same length as the input, containing the reverse
#'         complement sequences
#' @examples
#' grevcomp("ACTG") # Returns "CAGT"
#' grevcomp(c("ACTG", "GGCC")) # Returns c("CAGT", "GGCC")
#' grevcomp(c("ACTG", NA, "GGCC")) # Returns c("CAGT", NA, "GGCC")
#'
#' @export
grevcomp <- function(seq) {
    if (!is.character(seq)) {
        stop("Sequence must be a character string")
    }
    rev_s <- .Call("C_revcomp", seq)
    if (anyNA(seq)) {
        rev_s[is.na(seq)] <- NA_character_
    }

    if (!is.null(names(seq))) {
        names(rev_s) <- names(seq)
    }
    return(rev_s)
}
