# A numeric matrix with columns lo and hi, one row per sequence vtrack in
# `expr`: a bin [s, e) of that vtrack reads [s + lo, e + hi).
.gvariants_reach <- function(expr) {
    x <- .gcall("C_gvariants_reach", expr, .misha_env())
    if (!length(x)) {
        stop(sprintf("Expression \"%s\" has no sequence virtual track, so no substitution can change it", expr), call. = FALSE)
    }
    matrix(x, ncol = 2, dimnames = list(NULL, c("lo", "hi")))
}
