#' How misha treats NaN values
#'
#' Genomic tracks are sparse: a bin with no underlying data evaluates to `NaN`,
#' not to zero. The query functions do not all treat those bins the same way,
#' and the difference changes the answer rather than just the row count.
#'
#' @section NaN values:
#'
#' A track expression evaluates to `NaN` wherever the iterator produces a bin
#' the track has no data for. What happens next depends on the function:
#'
#' \itemize{
#'   \item \code{\link{gextract}} **keeps** `NaN` rows, so the result has one
#'     row per iterator interval whether or not the track covered it.
#'   \item \code{\link{gsummary}} **counts** them and reports the count as the
#'     "NaN intervals" element, while the statistics themselves are computed
#'     over the non-`NaN` values only.
#'   \item \code{\link{gdist}}, \code{\link{gquantiles}} and
#'     \code{\link{gscreen}} **drop** them: `NaN` bins are not counted into any
#'     distribution bin, do not contribute to a percentile, and never satisfy a
#'     screening condition - including a condition that would be true of every
#'     real value.
#'   \item \code{\link{gsegment}} **spans** them: a `NaN` bin contributes no
#'     evidence to the test that places a boundary, but it still falls inside
#'     whichever segment surrounds it, so the returned segments tile the scope
#'     continuously rather than skipping the gaps.
#' }
#'
#' So on 20 bins of which 7 are `NaN`, `gextract` returns 20 rows,
#' `gsummary` reports 20 total and 7 `NaN`, and `gdist` counts 13; and on a
#' 300 kb scope where 120 of 300 bins are `NaN`, `gsegment` still returns
#' segments covering the full 300 kb.
#'
#' The practical consequence is that `NaN` and zero are different, and
#' collapsing them with `ifelse(is.na(x), 0, x)` turns "no data here" into a
#' measured value of zero. Where that is genuinely what you want, note that it
#' also changes every mean, quantile and distribution computed downstream.
#'
#' @name misha-NaN
#' @keywords internal
NULL
