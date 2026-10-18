# A numeric matrix with columns lo and hi, one row per sequence vtrack in
# `expr`: a bin [s, e) of that vtrack reads [s + lo, e + hi).
.gvariants_reach <- function(expr) {
    x <- .gcall("C_gvariants_reach", expr, .misha_env())
    if (!length(x)) {
        stop(sprintf("Expression \"%s\" has no sequence virtual track, so no substitution can change it", expr), call. = FALSE)
    }
    matrix(x, ncol = 2, dimnames = list(NULL, c("lo", "hi")))
}

# The value of `expr` on each row, with base `base` (0..3 = A, C, G, T; -1 =
# none) at position `pos` while the row is computed. Rows sorted by chromosome
# and start.
.gvariants_eval <- function(expr, chroms, starts, ends, pos, base) {
    .gcall(
        "C_gextract_variants", expr, as.character(chroms), as.numeric(starts), as.numeric(ends),
        as.numeric(pos), as.integer(base), .misha_env()
    )
}
