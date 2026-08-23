# Returns number of intervals per chromosome

Returns number of intervals per chromosome (or chromosome pair).

## Usage

``` r
gintervals.chrom_sizes(intervals = NULL)
```

## Arguments

- intervals:

  intervals set

## Value

Data frame representing number of intervals per chromosome (for 1D
intervals) or chromosome pair (for 2D intervals).

## Details

This function returns number of intervals per chromosome (for 1D
intervals) or chromosome pair (for 2D intervals).

## See also

[`gintervals.load`](https://tanaylab.github.io/misha/reference/gintervals.load.md),
[`gintervals.save`](https://tanaylab.github.io/misha/reference/gintervals.save.md),
[`gintervals.exists`](https://tanaylab.github.io/misha/reference/gintervals.exists.md),
[`gintervals.ls`](https://tanaylab.github.io/misha/reference/gintervals.ls.md),
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md)

## Examples

``` r

gdb.init_examples()
gintervals.chrom_sizes("annotations")
#>   chrom size
#> 1  chr1    2
#> 2  chr2    6
```
