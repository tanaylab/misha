# Calculates summary statistics of a track expression for bins

Calculates summary statistics of a track expression for bins.

## Usage

``` r
gbins.summary(
  ...,
  expr = NULL,
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

  track expression for which summary statistics is calculated

- intervals:

  genomic scope for which the function is applied

- include.lowest:

  if 'TRUE', the lowest value of the range determined by breaks is
  included

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expressions.

- band:

  track expression band. If 'NULL' no band is used.

## Value

Multi-dimensional array representing summary statistics for each bin.

## Details

This function is a binned version of 'gsummary'. For each iterator
interval the value of 'bin_expr' is calculated and assigned to the
corresponding bin determined by 'breaks'. The summary statistics of
'expr' are calculated then separately for each bin.

The bins can be multi-dimensional depending on the number of
'bin_expr'-'breaks' pairs.

The range of bins is determined by 'breaks' argument. For example:
'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins):
(x1, x2\], (x2, x3\], (x3, x4\].

If 'include.lowest' is 'TRUE' the the lowest value will be included in
the first interval, i.e. in \[x1, x2\].

## See also

[`gsummary`](https://tanaylab.github.io/misha/reference/gsummary.md),
[`gintervals.summary`](https://tanaylab.github.io/misha/reference/gintervals.summary.md),
[`gdist`](https://tanaylab.github.io/misha/reference/gdist.md)

## Examples

``` r

gdb.init_examples()
gbins.summary("dense_track", c(0, 0.2, 0.4, 2), "sparse_track",
    intervals = gintervals(1), iterator = "dense_track"
)
#>           Total intervals NaN intervals       Min  Max       Sum      Mean
#> (0,0.2]              5740          5474 0.3555556 0.65 103.01308 0.3872672
#> (0.2,0.4]             771             0 0.3733333 0.76 373.53933 0.4844868
#> (0.4,2]                35             0 0.5371429 1.24  27.80314 0.7943755
#>              Std dev
#> (0,0.2]   0.04753164
#> (0.2,0.4] 0.07836847
#> (0.4,2]   0.16449928
#> attr(,"breaks")
#> attr(,"breaks")[[1]]
#> [1] 0.0 0.2 0.4 2.0
#> 
```
