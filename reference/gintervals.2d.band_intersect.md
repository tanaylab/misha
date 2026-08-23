# Intersects two-dimensional intervals with a band

Intersects two-dimensional intervals with a band.

## Usage

``` r
gintervals.2d.band_intersect(
  intervals = NULL,
  band = NULL,
  intervals.set.out = NULL
)
```

## Arguments

- intervals:

  two-dimensional intervals

- band:

  track expression band. If 'NULL' no band is used.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a data frame representing the
intervals.

## Details

This function intersects each two-dimensional interval from 'intervals'
with 'band'. If the intersection is not empty, the interval is shrunk to
the minimal rectangle that contains the band and added to the return
value.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

## See also

[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md),
[`gintervals.intersect`](https://tanaylab.github.io/misha/reference/gintervals.intersect.md)

## Examples

``` r

gdb.init_examples()
gintervals.2d.band_intersect(gintervals.2d(1), c(10000, 20000))
#>   chrom1 start1  end1 chrom2 start2   end2
#> 1   chr1  10000 5e+05   chr1      0 490000
```
