# Keep one interval of each group of overlapping intervals

Groups overlapping intervals as
[`gintervals.mark_overlaps`](https://tanaylab.github.io/misha/reference/gintervals.mark_overlaps.md)
does and keeps one interval of each group.

## Usage

``` r
gintervals.remove_overlaps(
  intervals,
  select = c("first", "last", "random"),
  unify_touching_intervals = TRUE
)
```

## Arguments

- intervals:

  intervals set

- select:

  which interval of each group to keep: "first" or "last" in row order,
  or "random"

- unify_touching_intervals:

  if 'TRUE', touching one-dimensional intervals are unified

## Value

The kept intervals, in their original row order.

## See also

[`gintervals.mark_overlaps`](https://tanaylab.github.io/misha/reference/gintervals.mark_overlaps.md),
[`gintervals.canonic`](https://tanaylab.github.io/misha/reference/gintervals.canonic.md)

## Examples

``` r

gdb.init_examples()
intervs <- data.frame(
    chrom = "chr1",
    start = c(11000, 100, 10000, 10500),
    end = c(12000, 200, 13000, 10600),
    data = c(10, 20, 30, 40)
)
gintervals.remove_overlaps(intervs)
#>   chrom start   end data
#> 1  chr1 11000 12000   10
#> 2  chr1   100   200   20
gintervals.remove_overlaps(intervs, select = "last")
#>   chrom start   end data
#> 2  chr1   100   200   20
#> 4  chr1 10500 10600   40
```
