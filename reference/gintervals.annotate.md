# Annotates 1D intervals using nearest neighbors

Annotates one-dimensional intervals by finding nearest neighbors in
another set of intervals and adding selected columns from the neighbors
to the original intervals.

## Usage

``` r
gintervals.annotate(
  intervals,
  annotation_intervals,
  annotation_columns = NULL,
  column_names = NULL,
  dist_column = "dist",
  max_dist = Inf,
  na_value = NA,
  maxneighbors = 1,
  tie_method = c("first", "min.start", "min.end"),
  overwrite = FALSE,
  keep_order = TRUE,
  intervals.set.out = NULL,
  ...
)
```

## Arguments

- intervals:

  Intervals to annotate (1D).

- annotation_intervals:

  Source intervals containing annotation data (1D).

- annotation_columns:

  Character vector of column names to copy from `annotation_intervals`.
  If `NULL` (default), all non-basic columns are used, i.e. everything
  beyond the coordinate/strand columns among: chrom, start, end, chrom1,
  start1, end1, chrom2, start2, end2, strand.

- column_names:

  Optional custom names for the annotation columns. If provided, must
  have the same length as `annotation_columns`. Defaults to using the
  original names.

- dist_column:

  Name of the distance column to include. Use `NULL` to omit the
  distance column. Defaults to "dist".

- max_dist:

  Maximum absolute distance. When finite, neighbors with
  `|dist| > max_dist` result in annotation columns being set to
  `na_value` for those rows, while the row itself is retained.

- na_value:

  Value(s) to use for annotations when beyond `max_dist` or when no
  neighbor is found. Can be a single scalar recycled for all columns, or
  a named list/vector supplying per-column values matching
  `column_names`.

- maxneighbors:

  Maximum number of neighbors per interval (duplicates intervals as
  needed). Defaults to 1.

- tie_method:

  Tie-breaking when distances are equal: one of "first" (arbitrary but
  stable), "min.start" (smaller neighbor start first), or "min.end"
  (smaller neighbor end first). Applies when `maxneighbors > 1`.

- overwrite:

  When `FALSE` (default), errors if selected annotation columns would
  overwrite existing columns in `intervals`. When `TRUE`, conflicting
  base columns are replaced by the annotation columns.

- keep_order:

  If `TRUE` (default), preserves the original order of `intervals` rows
  in the output.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

- ...:

  Additional arguments forwarded to `gintervals.neighbors` (e.g.,
  `mindist`, `maxdist`).

## Value

A data frame containing the original intervals plus the requested
annotation columns (and optional distance column). If
`maxneighbors > 1`, rows may be duplicated per input interval to
accommodate multiple neighbors.

## Details

The function wraps and extends `gintervals.neighbors` to provide
convenient column selection/renaming, optional distance inclusion,
distance thresholding with custom NA values, multiple neighbors per
interval, and deterministic tie-breaking. Currently supports 1D
intervals only.

\- When `annotation_columns = NULL`, all non-basic columns present in
`annotation_intervals` are included. - Setting `dist_column = NULL`
omits the distance column. - If no neighbor is found for an interval,
annotation columns are filled with `na_value` and the distance (when
present) is `NA_real_`. - Column name collisions are handled as follows:
when `overwrite=FALSE` a clear error is emitted; when `overwrite=TRUE`,
base columns with the same names are replaced by annotation columns.

## Examples

``` r
# Prepare toy data
intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
ann <- gintervals(1, c(900, 5400), c(950, 5500))
ann$remark <- c("a", "b")
ann$score <- c(10, 20)

# Basic usage with default columns (all non-basic columns)
gintervals.annotate(intervs, ann)
#>   chrom start  end remark score dist
#> 1  chr1  1000 1100      a    10   50
#> 2  chr1  5000 5050      b    20  350

# Select specific columns, with custom names and distance column name
gintervals.annotate(
    intervs, ann,
    annotation_columns = c("remark"),
    column_names = c("ann_remark"),
    dist_column = "ann_dist"
)
#>   chrom start  end ann_remark ann_dist
#> 1  chr1  1000 1100          a       50
#> 2  chr1  5000 5050          b      350

# Distance threshold with scalar NA replacement
gintervals.annotate(
    intervs, ann,
    annotation_columns = c("remark"),
    max_dist = 200,
    na_value = "no_ann"
)
#>   chrom start  end remark dist
#> 1  chr1  1000 1100      a   50
#> 2  chr1  5000 5050 no_ann  350

# Multiple neighbors with deterministic tie-breaking
nbrs <- gintervals.annotate(
    gintervals(1, 1000, 1100),
    {
        x <- gintervals(1, c(800, 1200), c(900, 1300))
        x$label <- c("left", "right")
        x
    },
    annotation_columns = "label",
    maxneighbors = 2,
    tie_method = "min.start"
)
nbrs
#>   chrom start  end label dist
#> 1  chr1  1000 1100  left  100
#> 2  chr1  1000 1100 right  100

# Overwrite existing columns in the base intervals
intervs2 <- intervs
intervs2$remark <- c("orig1", "orig2")
gintervals.annotate(intervs2, ann, annotation_columns = "remark", overwrite = TRUE)
#>   chrom start  end remark dist
#> 1  chr1  1000 1100      a   50
#> 2  chr1  5000 5050      b  350
```
