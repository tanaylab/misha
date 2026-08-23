# Combines several sets of intervals

Combines several sets of intervals into one set.

## Usage

``` r
gintervals.rbind(..., intervals.set.out = NULL)
```

## Arguments

- ...:

  intervals sets to combine

- intervals.set.out:

  intervals set name where the function result is optionally outputted

- intervals:

  intervals set

## Value

If 'intervals.set.out' is 'NULL' a data frame combining intervals sets.

## Details

This function combines several intervals sets into one set. It works in
a similar manner as 'rbind' yet it is faster. Also it supports intervals
sets that are stored in files including the big intervals sets.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. If the format of the output intervals is set to be "big"
(determined implicitly based on the result size and options), the order
of the resulted intervals is altered as they are sorted by chromosome
(or chromosomes pair - for 2D).

## See also

[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md),
[`gintervals.canonic`](https://tanaylab.github.io/misha/reference/gintervals.canonic.md)

## Examples

``` r

gdb.init_examples()

intervs1 <- gextract("sparse_track", gintervals(c(1, 2), 1000, 4000))
intervs2 <- gextract("sparse_track", gintervals(c(2, "X"), 2000, 5000))
gintervals.save("testintervs", intervs2)
gintervals.rbind(intervs1, "testintervs")
#>    chrom start  end sparse_track intervalID
#> 1   chr1  1400 1450    0.5200000          1
#> 2   chr2  2100 2150    0.4800000          2
#> 3   chr2  2300 2350    0.6400000          2
#> 4   chr2  2500 2550    0.4400000          2
#> 5   chr2  2900 2950    0.3600000          2
#> 6   chr2  3250 3300    0.3600000          2
#> 7   chr2  3450 3500    0.4000000          2
#> 8   chr2  3550 3600    0.4000000          2
#> 9   chr2  3700 3750    0.6000000          2
#> 10  chr2  2100 2150    0.4800000          1
#> 11  chr2  2300 2350    0.6400000          1
#> 12  chr2  2500 2550    0.4400000          1
#> 13  chr2  2900 2950    0.3600000          1
#> 14  chr2  3250 3300    0.3600000          1
#> 15  chr2  3450 3500    0.4000000          1
#> 16  chr2  3550 3600    0.4000000          1
#> 17  chr2  3700 3750    0.6000000          1
#> 18  chr2  4300 4350    0.3600000          1
#> 19  chr2  4550 4600    0.5200000          1
#> 20  chr2  4650 4700    0.4800000          1
#> 21  chr2  4900 5000    0.5700000          1
#> 22  chrX  2600 2650    0.3600000          2
#> 23  chrX  3150 3250    0.4200000          2
#> 24  chrX  3500 3650    0.4533333          2
#> 25  chrX  3700 3750    0.4400000          2
#> 26  chrX  3900 3950    0.3600000          2
#> 27  chrX  4100 4150    0.7200000          2
#> 28  chrX  4800 4850    0.3600000          2
#> 29  chrX  4900 5000    0.3800000          2
gintervals.rm("testintervs", force = TRUE)
```
