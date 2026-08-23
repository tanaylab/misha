# Limits intervals to chromosomal range

Limits intervals to chromosomal range.

## Usage

``` r
gintervals.force_range(intervals = NULL, intervals.set.out = NULL)
```

## Arguments

- intervals:

  intervals

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a data frame representing the
intervals.

## Details

This function enforces the intervals to be within the chromosomal range
\[0, chrom length) by altering the intervals' boundaries. Intervals that
lay entirely outside of the chromosomal range are eliminated. The new
intervals are returned.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

## See also

[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md),
[`gintervals.canonic`](https://tanaylab.github.io/misha/reference/gintervals.canonic.md)

## Examples

``` r

gdb.init_examples()
intervs <- data.frame(
    chrom = "chr1",
    start = c(11000, -100, 10000, 10500),
    end = c(12000, 200, 13000000, 10600)
)
gintervals.force_range(intervs)
#>   chrom start    end
#> 1  chr1 11000  12000
#> 2  chr1     0    200
#> 3  chr1 10000 500000
#> 4  chr1 10500  10600
```
