# Calculates summary statistics of track expression

Calculates summary statistics of track expression.

## Usage

``` r
gsummary(expr = NULL, intervals = NULL, iterator = NULL, band = NULL)
```

## Arguments

- expr:

  track expression

- intervals:

  genomic scope for which the function is applied

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expression.

- band:

  track expression band. If 'NULL' no band is used.

## Value

An array that represents summary statistics.

## Details

This function returns summary statistics of a track expression: total
number of bins, total number of bins whose value is NaN, min, max, sum,
mean and standard deviation of the values.

## See also

[`gintervals.summary`](https://tanaylab.github.io/misha/reference/gintervals.summary.md),
[`gbins.summary`](https://tanaylab.github.io/misha/reference/gbins.summary.md)

## Examples

``` r

gdb.init_examples()
gsummary("rects_track")
#> Total intervals   NaN intervals             Min             Max             Sum 
#>         1600.00            0.00            3.00       141495.00     85801800.00 
#>            Mean         Std dev 
#>        53626.12        39867.94 
```
