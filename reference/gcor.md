# Calculates correlation between track expressions

Calculates correlation between track expressions over iterator bins
inside the supplied genomic scope. Expressions are processed in pairs:
(expr1, expr2), (expr3, expr4), etc. Only bins where both expressions
are not NaN are used.

## Usage

``` r
gcor(
  expr1 = NULL,
  expr2 = NULL,
  ...,
  intervals = NULL,
  iterator = NULL,
  band = NULL,
  method = c("pearson", "spearman", "spearman.exact"),
  details = FALSE,
  names = NULL
)
```

## Arguments

- expr1:

  first track expression

- expr2:

  second track expression

- ...:

  additional track expressions, supplied as pairs (expr3, expr4, ...)

- intervals:

  genomic scope for which the function is applied

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expression.

- band:

  track expression band. If 'NULL' no band is used.

- method:

  correlation method to use. One of 'pearson' (default), 'spearman'
  (approximate, memory-efficient), or 'spearman.exact' (exact, requires
  O(n) memory where n is number of non-NaN pairs).

- details:

  if 'TRUE' returns summary statistics for each pair, otherwise returns
  correlations only. For Pearson, includes n, n.na, mean1, mean2, sd1,
  sd2, cov, cor. For Spearman methods, includes n, n.na, cor.

- names:

  optional names for the pairs. If supplied, length must match the
  number of pairs.

## Value

If 'details' is 'FALSE', a numeric vector of correlations. If 'details'
is 'TRUE', a data frame with summary statistics for each pair.

## See also

[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md),
[`gscreen`](https://tanaylab.github.io/misha/reference/gscreen.md),
[`gsummary`](https://tanaylab.github.io/misha/reference/gsummary.md)

## Examples

``` r

gdb.init_examples()
gcor("dense_track", "sparse_track", intervals = gintervals(1, 0, 10000), iterator = 1000)
#> dense_track~sparse_track 
#>               -0.5135881 

# Spearman correlation (approximate, memory-efficient)
gcor("dense_track", "sparse_track",
    intervals = gintervals(1, 0, 10000),
    iterator = 1000, method = "spearman"
)
#> dense_track~sparse_track 
#>               -0.4103913 

# Exact Spearman correlation
gcor("dense_track", "sparse_track",
    intervals = gintervals(1, 0, 10000),
    iterator = 1000, method = "spearman.exact"
)
#> dense_track~sparse_track 
#>               -0.4103913 

```
