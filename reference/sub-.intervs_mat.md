# Subset an intervs_mat preserving interval identity

Row subset of an `intervs_mat` subsets `attr(., "intervals")` in
parallel. Column-only subset leaves the intervals attribute unchanged.
Single-row results with `drop = TRUE` return a plain numeric vector
(class dropped).

## Usage

``` r
# S3 method for class 'intervs_mat'
x[i, j, ..., drop = TRUE]
```

## Arguments

- x:

  an `intervs_mat`

- i:

  row selector (logical, integer, or character matching rownames)

- j:

  column selector

- ...:

  unused

- drop:

  if `TRUE` and a single row is selected via `i`, collapse to a named
  vector and drop the class. Defaults to `TRUE` to match base matrix `[`
  for row selection. Column-only subsets (`i` missing) are never dropped
  to a vector, even if `j` selects a single column - this preserves the
  row-interval correspondence.

## Value

An `intervs_mat` (still 2D), or a numeric vector (degenerate).
