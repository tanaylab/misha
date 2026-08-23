# Convert intervals + values data.frame to an interval-indexed matrix

Builds a numeric matrix of value columns whose rows are indexed by
intervals. The intervals are carried in `attr(mat, "intervals")` as the
authoritative identity; `rownames(mat)` are display-only (default
`"chrom:start-end"`) and are NEVER parsed back by
[`gintervals.from_mat`](https://tanaylab.github.io/misha/reference/gintervals.from_mat.md).
This avoids the round-trip corruption that occurs when chrom names
contain underscores or other separators.

## Usage

``` r
gintervals.to_mat(df, id_col = NULL, value_cols = NULL, labels = TRUE)
```

## Arguments

- df:

  data.frame with `chrom`, `start`, `end` and zero or more value
  columns. May contain an `intervalID` column (from `gextract`), which
  is kept in the attribute but excluded from the matrix.

- id_col:

  optional column name in `df` whose values become rownames. If `NULL`
  (default), rownames are `"chrom:start-end"`.

- value_cols:

  character vector of column names to use as matrix data. If `NULL`
  (default), auto-detect: all columns except `chrom`, `start`, `end`,
  `intervalID`. Auto-detect errors if any selected column is
  non-numeric; pass `value_cols` explicitly to override.

- labels:

  if `TRUE` (default), set `rownames(mat)` to either `df[[id_col]]` (if
  `id_col` is supplied) or `"chrom:start-end"`. If `FALSE`, leave
  rownames `NULL` - useful in pipelines that don't need the display
  labels and would prefer to skip the construction cost on large inputs.
  When `FALSE`, the `id_col` argument is ignored.

## Value

An `intervs_mat` object: a numeric matrix subclass with the intervals
attached as `attr(., "intervals")`. Supports row/column subsetting (`[`)
and [`rbind()`](https://rdrr.io/r/base/cbind.html) while preserving the
attribute.

## See also

[`gintervals.from_mat`](https://tanaylab.github.io/misha/reference/gintervals.from_mat.md)

## Examples

``` r
df <- data.frame(
    chrom = c("chr1", "chr1", "chr2"),
    start = c(100L, 500L, 200L),
    end   = c(200L, 700L, 400L),
    t1    = c(1.5, 2.5, 3.5),
    t2    = c(10, 20, 30)
)
mat <- gintervals.to_mat(df)
rownames(mat)
#> [1] "chr1:100-200" "chr1:500-700" "chr2:200-400"
# subset preserves intervals:
sub <- mat[c(1, 3), ]
attr(sub, "intervals")
#>   chrom start end
#> 1  chr1   100 200
#> 2  chr2   200 400
# round-trip back to a data.frame:
gintervals.from_mat(sub)
#>   chrom start end  t1 t2
#> 1  chr1   100 200 1.5 10
#> 2  chr2   200 400 3.5 30
```
