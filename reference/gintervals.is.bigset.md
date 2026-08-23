# Tests for big intervals set

Tests for big intervals set.

## Usage

``` r
gintervals.is.bigset(intervals.set = NULL)
```

## Arguments

- intervals.set:

  name of an intervals set

## Value

'TRUE' if intervals set is big, otherwise 'FALSE'.

## Details

This function tests whether 'intervals.set' is a big intervals set.
Intervals set is big if it is stored in big intervals set format and
given the current limits it cannot be fully loaded into memory.

Memory limit is controlled by 'gmax.data.size' option (see:
'getOption("gmax.data.size")').

## See also

[`gintervals.load`](https://tanaylab.github.io/misha/reference/gintervals.load.md),
[`gintervals.save`](https://tanaylab.github.io/misha/reference/gintervals.save.md),
[`gintervals.exists`](https://tanaylab.github.io/misha/reference/gintervals.exists.md),
[`gintervals.ls`](https://tanaylab.github.io/misha/reference/gintervals.ls.md)

## Examples

``` r

gdb.init_examples()
gintervals.is.bigset("annotations")
#> [1] FALSE
```
