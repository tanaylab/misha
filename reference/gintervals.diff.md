# Calculates difference of two intervals sets

Returns difference of two sets of intervals.

## Usage

``` r
gintervals.diff(intervals1 = NULL, intervals2 = NULL, intervals.set.out = NULL)
```

## Arguments

- intervals1, intervals2:

  set of one-dimensional intervals

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a data frame representing the
intervals.

## Details

This function returns a genomic space that is covered by 'intervals1'
but not covered by 'intervals2'.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

## See also

[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.intersect`](https://tanaylab.github.io/misha/reference/gintervals.intersect.md),
[`gintervals.union`](https://tanaylab.github.io/misha/reference/gintervals.union.md)

## Examples

``` r

gdb.init_examples()

intervs1 <- gscreen("dense_track > 0.15")
intervs2 <- gscreen("dense_track < 0.2")

## 'res3' equals to 'res4'
res3 <- gintervals.diff(intervs1, intervs2)
res4 <- gscreen("dense_track >= 0.2")
```
