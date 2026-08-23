# Calculates quantiles of a track expression for intervals

Calculates quantiles of a track expression for intervals.

## Usage

``` r
gintervals.quantiles(
  expr = NULL,
  percentiles = 0.5,
  intervals = NULL,
  iterator = NULL,
  band = NULL,
  intervals.set.out = NULL
)
```

## Arguments

- expr:

  track expression for which quantiles are calculated

- percentiles:

  an array of percentiles of quantiles in \[0, 1\] range

- intervals:

  set of intervals

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expressions.

- band:

  track expression band. If 'NULL' no band is used.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a set of intervals with additional
columns representing quantiles for each percentile.

## Details

This function calculates quantiles of 'expr' for each interval in
'intervals'.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

## See also

[`gquantiles`](https://tanaylab.github.io/misha/reference/gquantiles.md),
[`gbins.quantiles`](https://tanaylab.github.io/misha/reference/gbins.quantiles.md)

## Examples

``` r

gdb.init_examples()
intervs <- gintervals(c(1, 2), 0, 5000)
gintervals.quantiles("dense_track",
    percentiles = c(0.5, 0.3, 0.9), intervs
)
#>   chrom start  end  0.5   0.3  0.9
#> 1  chr1     0 5000 0.04 0.020 0.16
#> 2  chr2     0 5000 0.08 0.034 0.18
```
