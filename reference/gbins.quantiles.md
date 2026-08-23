# Calculates quantiles of a track expression for bins

Calculates quantiles of a track expression for bins.

## Usage

``` r
gbins.quantiles(
  ...,
  expr = NULL,
  percentiles = 0.5,
  intervals = get("ALLGENOME", envir = .misha),
  include.lowest = FALSE,
  iterator = NULL,
  band = NULL
)
```

## Arguments

- ...:

  pairs of track expressions ('bin_expr') that determines the bins and
  breaks that define the bins. See
  [`gdist`](https://tanaylab.github.io/misha/reference/gdist.md). An
  extra, unpaired trailing argument is taken as 'expr'; it cannot be
  combined with a named 'expr' argument - that combination is ambiguous
  and raises an error.

- expr:

  track expression for which quantiles are calculated

- percentiles:

  an array of percentiles of quantiles in \[0, 1\] range

- intervals:

  genomic scope for which the function is applied.

- include.lowest:

  if 'TRUE', the lowest value of the range determined by breaks is
  included

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expressions.

- band:

  track expression band. If 'NULL' no band is used.

## Value

Multi-dimensional array representing quantiles for each percentile and
bin.

## Details

This function is a binned version of 'gquantiles'. For each iterator
interval the value of 'bin_expr' is calculated and assigned to the
corresponding bin determined by 'breaks'. The quantiles of 'expr' are
calculated then separately for each bin.

The bins can be multi-dimensional depending on the number of
'bin_expr'-'breaks' pairs.

The range of bins is determined by 'breaks' argument. For example:
'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins):
(x1, x2\], (x2, x3\], (x3, x4\].

If 'include.lowest' is 'TRUE' the the lowest value will be included in
the first interval, i.e. in \[x1, x2\].

## See also

[`gquantiles`](https://tanaylab.github.io/misha/reference/gquantiles.md),
[`gintervals.quantiles`](https://tanaylab.github.io/misha/reference/gintervals.quantiles.md),
[`gdist`](https://tanaylab.github.io/misha/reference/gdist.md)

## Examples

``` r

gdb.init_examples()
gbins.quantiles("dense_track", c(0, 0.2, 0.4, 2), "sparse_track",
    percentiles = c(0.2, 0.5),
    intervals = gintervals(1),
    iterator = "dense_track"
)
#>                 0.2   0.5
#> (0,0.2]   0.3600000 0.360
#> (0.2,0.4] 0.4133333 0.475
#> (0.4,2]   0.6480000 0.800
#> attr(,"breaks")
#> attr(,"breaks")[[1]]
#> [1] 0.0 0.2 0.4 2.0
#> 
```
