# Normalize intervals to fixed or variable sizes

This function normalizes intervals by computing their centers and then
expanding them to fixed or variable sizes, while ensuring they don't
cross chromosome boundaries.

## Usage

``` r
gintervals.normalize(intervals = NULL, size = NULL, intervals.set.out = NULL)
```

## Arguments

- intervals:

  intervals set

- size:

  target size(s) for normalized intervals. Can be either:

  - A single positive integer: all intervals normalized to this size

  - A numeric vector: each interval normalized to its corresponding
    size. Vector length must exactly match the number of intervals, OR a
    single interval can be provided with multiple sizes to create
    multiple output intervals (one-to-many expansion).

- intervals.set.out:

  intervals set name where the function result is saved. If NULL, the
  result is returned to the user.

## Value

Normalized intervals set with fixed or variable sizes, or NULL if result
is saved to intervals.set.out

## See also

[`gintervals.force_range`](https://tanaylab.github.io/misha/reference/gintervals.force_range.md)

## Examples

``` r

gdb.init_examples()

# Single size (all intervals normalized to 500bp)
intervs <- gintervals(1, c(1000, 5000), c(2000, 6000))
gintervals.normalize(intervs, 500)
#>   chrom start  end
#> 1  chr1  1250 1750
#> 2  chr1  5250 5750

# Vector of sizes (each interval gets its own size)
intervs <- gintervals(1, c(1000, 3000, 5000), c(2000, 4000, 6000))
gintervals.normalize(intervs, c(500, 1000, 750))
#>   chrom start  end
#> 1  chr1  1250 1750
#> 2  chr1  3000 4000
#> 3  chr1  5125 5875

# One-to-many: single interval with multiple sizes
interv <- gintervals(1, 1000, 2000)
gintervals.normalize(interv, c(500, 1000, 1500))
#>   chrom start  end
#> 1  chr1  1250 1750
#> 2  chr1  1000 2000
#> 3  chr1   750 2250
```
