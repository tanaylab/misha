# Unites two sets of 2D intervals

Unites two sets of two-dimensional intervals.

## Usage

``` r
gintervals.2d.union(intervals1 = NULL, intervals2 = NULL)
```

## Arguments

- intervals1, intervals2:

  two-dimensional intervals (a data frame or the name of an intervals
  set)

## Value

A data frame with the six 2D interval columns, sorted as by
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md).

## Details

Concatenates 'intervals1' and 'intervals2' and returns the combined set
sorted as by
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md).
Overlapping rectangles are not merged: because the union of two
rectangles is not generally a rectangle, the result is simply the
combined set (unlike the 1D
[`gintervals.union`](https://tanaylab.github.io/misha/reference/gintervals.union.md)).

## See also

[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md),
[`gintervals.2d.intersect`](https://tanaylab.github.io/misha/reference/gintervals.2d.intersect.md),
[`gintervals.union`](https://tanaylab.github.io/misha/reference/gintervals.union.md)

## Examples

``` r

gdb.init_examples()
a <- gintervals.2d("chr1", 0, 1000, "chr2", 0, 1000)
b <- gintervals.2d("chr1", 500, 1500, "chr2", 500, 1500)
gintervals.2d.union(a, b)
#>   chrom1 start1 end1 chrom2 start2 end2
#> 1   chr1      0 1000   chr2      0 1000
#> 2   chr1    500 1500   chr2    500 1500
```
