# Tests for a named intervals set existence

Tests for a named intervals set existence.

## Usage

``` r
gintervals.exists(intervals.set = NULL)
```

## Arguments

- intervals.set:

  name of an intervals set

## Value

'TRUE' if a named intervals set exists. Otherwise 'FALSE'.

## Details

This function returns 'TRUE' if a named intervals set exists in Genomic
Database.

## See also

[`gintervals.ls`](https://tanaylab.github.io/misha/reference/gintervals.ls.md),
[`gintervals.load`](https://tanaylab.github.io/misha/reference/gintervals.load.md),
[`gintervals.rm`](https://tanaylab.github.io/misha/reference/gintervals.rm.md),
[`gintervals.save`](https://tanaylab.github.io/misha/reference/gintervals.save.md),
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md)

## Examples

``` r

gdb.init_examples()
gintervals.exists("annotations")
#> [1] TRUE
```
