# Creates a named intervals set

Saves intervals to a named intervals set.

## Usage

``` r
gintervals.save(intervals.set.out = NULL, intervals = NULL)
```

## Arguments

- intervals.set.out:

  name of the new intervals set

- intervals:

  intervals to save

## Value

None.

## Details

This function saves 'intervals' as a named intervals set.

## See also

[`gintervals.rm`](https://tanaylab.github.io/misha/reference/gintervals.rm.md),
[`gintervals.load`](https://tanaylab.github.io/misha/reference/gintervals.load.md),
[`gintervals.exists`](https://tanaylab.github.io/misha/reference/gintervals.exists.md),
[`gintervals.ls`](https://tanaylab.github.io/misha/reference/gintervals.ls.md),
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md)

## Examples

``` r

gdb.init_examples()
intervs <- gintervals(c(1, 2))
gintervals.save("testintervs", intervs)
gintervals.ls()
#> [1] "annotations" "testintervs"
gintervals.rm("testintervs", force = TRUE)
```
