# Partitions the values of track expression

Converts the values of track expression to intervals that match
corresponding bin.

## Usage

``` r
gpartition(
  expr = NULL,
  breaks = NULL,
  intervals = NULL,
  include.lowest = FALSE,
  iterator = NULL,
  band = NULL,
  intervals.set.out = NULL
)
```

## Arguments

- expr:

  track expression

- breaks:

  breaks that determine the bin

- intervals:

  genomic scope for which the function is applied

- include.lowest:

  if 'TRUE', the lowest value of the range determined by breaks is
  included

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expression.

- band:

  track expression band. If 'NULL' no band is used.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a set of intervals with an additional
column that indicates the corresponding bin index.

## Details

This function converts first the values of track expression into 1-based
bin's index according 'breaks' argument. It returns then the intervals
with the corresponding bin's index.

The range of bins is determined by 'breaks' argument. For example:
'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins):
(x1, x2\], (x2, x3\], (x3, x4\].

If 'include.lowest' is 'TRUE' the the lowest value will be included in
the first interval, i.e. in \[x1, x2\].

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

## See also

[`gscreen`](https://tanaylab.github.io/misha/reference/gscreen.md),
[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md),
[`glookup`](https://tanaylab.github.io/misha/reference/glookup.md),
[`gdist`](https://tanaylab.github.io/misha/reference/gdist.md)

## Examples

``` r

gdb.init_examples()
breaks <- seq(0, 0.2, by = 0.05)
gpartition("dense_track", breaks, gintervals(1, 0, 5000))
#>    chrom start  end bin
#> 1   chr1     0  250   4
#> 2   chr1   300  450   4
#> 3   chr1   450  550   2
#> 4   chr1   550  650   1
#> 5   chr1   850 1000   1
#> 6   chr1  1000 1050   3
#> 7   chr1  1100 1150   1
#> 8   chr1  1250 1300   1
#> 9   chr1  1300 1350   2
#> 10  chr1  1350 1400   1
#> 11  chr1  1450 1500   1
#> 12  chr1  1500 1550   2
#> 13  chr1  1600 1700   2
#> 14  chr1  1700 1750   1
#> 15  chr1  1850 1950   1
#> 16  chr1  2000 2100   3
#> 17  chr1  2150 2200   1
#> 18  chr1  2250 2300   1
#> 19  chr1  2300 2350   2
#> 20  chr1  2400 2450   3
#> 21  chr1  2500 2550   1
#> 22  chr1  2700 2750   1
#> 23  chr1  2750 2850   2
#> 24  chr1  2850 2900   3
#> 25  chr1  2900 2950   1
#> 26  chr1  2950 3000   2
#> 27  chr1  3000 3050   1
#> 28  chr1  3050 3100   2
#> 29  chr1  3100 3150   1
#> 30  chr1  3200 3250   2
#> 31  chr1  3300 3350   1
#> 32  chr1  3400 3500   1
#> 33  chr1  3500 3550   3
#> 34  chr1  3550 3600   1
#> 35  chr1  3650 3700   1
#> 36  chr1  3700 3800   2
#> 37  chr1  3850 3950   1
#> 38  chr1  3950 4000   3
#> 39  chr1  4000 4050   1
#> 40  chr1  4050 4100   2
#> 41  chr1  4100 4200   1
#> 42  chr1  4250 4350   2
#> 43  chr1  4350 4500   1
#> 44  chr1  4500 4550   2
#> 45  chr1  4600 4650   1
#> 46  chr1  4700 4750   3
#> 47  chr1  4750 4800   4
#> 48  chr1  4850 4900   1
#> 49  chr1  4900 4950   2
```
