# Calculates a union of two sets of intervals

Calculates a union of two sets of intervals.

## Usage

``` r
gintervals.union(
  intervals1 = NULL,
  intervals2 = NULL,
  intervals.set.out = NULL
)
```

## Arguments

- intervals1, intervals2:

  set of one-dimensional intervals

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a data frame representing the union of
intervals.

## Details

This function returns intervals that represent a genomic space covered
by either 'intervals1' or 'intervals2'.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

## See also

[`gintervals.intersect`](https://tanaylab.github.io/misha/reference/gintervals.intersect.md),
[`gintervals.diff`](https://tanaylab.github.io/misha/reference/gintervals.diff.md),
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md)

## Examples

``` r

gdb.init_examples()

intervs1 <- gscreen("dense_track > 0.15 & dense_track < 0.18")
intervs2 <- gscreen("dense_track >= 0.18 & dense_track < 0.2")

## 'intervs3' and 'intervs4' are identical
intervs3 <- gintervals.union(intervs1, intervs2)
intervs4 <- gscreen("dense_track > 0.15 & dense_track < 0.2")
```
