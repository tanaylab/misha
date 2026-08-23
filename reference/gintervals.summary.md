# Calculates summary statistics of track expression for intervals

Calculates summary statistics of track expression for intervals.

## Usage

``` r
gintervals.summary(
  expr = NULL,
  intervals = NULL,
  iterator = NULL,
  band = NULL,
  intervals.set.out = NULL
)
```

## Arguments

- expr:

  track expression

- intervals:

  set of intervals

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expression.

- band:

  track expression band. If 'NULL' no band is used.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a set of intervals with additional
columns representing summary statistics for each percentile and
interval.

## Details

This function returns summary statistics of a track expression for each
interval 'intervals': total number of bins, total number of bins whose
value is NaN, min, max, sum, mean and standard deviation of the values.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

## See also

[`gsummary`](https://tanaylab.github.io/misha/reference/gsummary.md),
[`gbins.summary`](https://tanaylab.github.io/misha/reference/gbins.summary.md)

## Examples

``` r

gdb.init_examples()
intervs <- gintervals(c(1, 2), 0, 5000)
gintervals.summary("dense_track", intervs)
#>   chrom start  end Total intervals NaN intervals Min  Max      Sum       Mean
#> 1  chr1     0 5000             100             0   0 0.26 5.117778 0.05117778
#> 2  chr2     0 5000             100             0   0 0.34 8.700000 0.08700000
#>      Std dev
#> 1 0.05596612
#> 2 0.07818006
```
