# Expand intervals on both sides

Adds 'expansion' bp to each side of every 1D interval, then clips the
result to the chromosome boundaries (see
[`gintervals.force_range`](https://tanaylab.github.io/misha/reference/gintervals.force_range.md)).

## Usage

``` r
gintervals.expand(intervals = NULL, expansion = NULL)
```

## Arguments

- intervals:

  intervals set

- expansion:

  bp added to each side: one value, or one per interval

## Value

A data frame of the expanded intervals, with the other columns kept.

## See also

[`gintervals.normalize`](https://tanaylab.github.io/misha/reference/gintervals.normalize.md),
[`gintervals.promoters`](https://tanaylab.github.io/misha/reference/gintervals.promoters.md)

## Examples

``` r

gdb.init_examples()
gintervals.expand(gintervals(1, c(1000, 5000), c(2000, 6000)), 100)
#>   chrom start  end
#> 1  chr1   900 2100
#> 2  chr1  4900 6100
```
