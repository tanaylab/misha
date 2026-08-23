# Intersects two sets of 2D intervals

Intersects two sets of two-dimensional intervals.

## Usage

``` r
gintervals.2d.intersect(intervals1 = NULL, intervals2 = NULL)
```

## Arguments

- intervals1, intervals2:

  two-dimensional intervals (a data frame or the name of an intervals
  set)

## Value

A data frame with the six 2D interval columns, sorted as by
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md),
or 'NULL' if the intersection is empty.

## Details

Returns the pairwise rectangle intersections between 'intervals1' and
'intervals2'. For every pair of rectangles that share the same (chrom1,
chrom2) chromosome pair, the overlapping rectangle is computed by
clipping each axis independently:

- start1 = max(start1 of the two rectangles), end1 = min(end1 ...)

- start2 = max(start2 of the two rectangles), end2 = min(end2 ...)

A rectangle is emitted only when it remains non-empty on both axes
(start1 \< end1 and start2 \< end2). The result is not merged or
canonicalized: because the union of two rectangles is not generally a
rectangle, overlapping outputs are left as-is (unlike the 1D
[`gintervals.intersect`](https://tanaylab.github.io/misha/reference/gintervals.intersect.md)).

## See also

[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md),
[`gintervals.2d.union`](https://tanaylab.github.io/misha/reference/gintervals.2d.union.md),
[`gintervals.intersect`](https://tanaylab.github.io/misha/reference/gintervals.intersect.md)

## Examples

``` r

gdb.init_examples()
a <- gintervals.2d("chr1", 0, 1000, "chr2", 0, 1000)
b <- gintervals.2d("chr1", 500, 1500, "chr2", 500, 1500)
gintervals.2d.intersect(a, b)
#>   chrom1 start1 end1 chrom2 start2 end2
#> 1   chr1    500 1000   chr2    500 1000
```
