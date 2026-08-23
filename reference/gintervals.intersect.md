# Calculates an intersection of two sets of intervals

Calculates an intersection of two sets of intervals.

## Usage

``` r
gintervals.intersect(
  intervals1 = NULL,
  intervals2 = NULL,
  intervals.set.out = NULL
)
```

## Arguments

- intervals1, intervals2:

  set of intervals

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a data frame representing the
intersection of intervals.

## Details

This function returns intervals that represent a genomic space which is
achieved by intersection of 'intervals1' and 'intervals2'.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

## See also

[`gintervals.2d.band_intersect`](https://tanaylab.github.io/misha/reference/gintervals.2d.band_intersect.md),
[`gintervals.diff`](https://tanaylab.github.io/misha/reference/gintervals.diff.md),
[`gintervals.union`](https://tanaylab.github.io/misha/reference/gintervals.union.md),
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md)

## Examples

``` r

gdb.init_examples()

intervs1 <- gscreen("dense_track > 0.15")
intervs2 <- gscreen("dense_track < 0.2")

## 'intervs3' and 'intervs4' are identical
intervs3 <- gintervals.intersect(intervs1, intervs2)
intervs4 <- gscreen("dense_track > 0.15 & dense_track < 0.2")
```
