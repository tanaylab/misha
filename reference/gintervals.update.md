# Updates a named intervals set

Updates a named intervals set.

## Usage

``` r
gintervals.update(
  intervals.set = NULL,
  intervals = "",
  chrom = NULL,
  chrom1 = NULL,
  chrom2 = NULL
)
```

## Arguments

- intervals.set:

  name of an intervals set

- intervals:

  intervals or 'NULL'

- chrom:

  chromosome for 1D intervals set

- chrom1:

  first chromosome for 2D intervals set

- chrom2:

  second chromosome for 2D intervals set

## Value

None.

## Details

This function replaces all intervals of given chromosome (or chromosome
pair) within 'intervals.set' with 'intervals'. Chromosome is specified
by 'chrom' for 1D intervals set or 'chrom1', 'chrom2' for 2D intervals
set.

If 'intervals' is 'NULL' all intervals of given chromosome are removed
from 'intervals.set'.

## See also

[`gintervals.save`](https://tanaylab.github.io/misha/reference/gintervals.save.md),
[`gintervals.load`](https://tanaylab.github.io/misha/reference/gintervals.load.md),
[`gintervals.exists`](https://tanaylab.github.io/misha/reference/gintervals.exists.md),
[`gintervals.ls`](https://tanaylab.github.io/misha/reference/gintervals.ls.md)

## Examples

``` r

gdb.init_examples()
intervs <- gscreen(
    "sparse_track > 0.2",
    gintervals(c(1, 2), 0, 10000)
)
gintervals.save("testintervs", intervs)
gintervals.load("testintervs")
#>    chrom start  end
#> 1   chr1     0   50
#> 2   chr1   100  150
#> 3   chr1   250  300
#> 4   chr1  1400 1450
#> 5   chr1  4750 4800
#> 6   chr1  6350 6400
#> 7   chr1  9900 9950
#> 8   chr2  2100 2150
#> 9   chr2  2300 2350
#> 10  chr2  2500 2550
#> 11  chr2  2900 2950
#> 12  chr2  3250 3300
#> 13  chr2  3450 3500
#> 14  chr2  3550 3600
#> 15  chr2  3700 3750
#> 16  chr2  4300 4350
#> 17  chr2  4550 4600
#> 18  chr2  4650 4700
#> 19  chr2  4900 5100
#> 20  chr2  6250 6300
#> 21  chr2  6450 6500
#> 22  chr2  7100 7200
#> 23  chr2  8150 8200
#> 24  chr2  8700 8800
#> 25  chr2  8950 9050
#> 26  chr2  9400 9450
#> 27  chr2  9550 9650
#> 28  chr2  9750 9850
gintervals.update("testintervs", intervs[intervs$chrom == "chr2", ][1:5, ], chrom = 2)
gintervals.load("testintervs")
#>    chrom start  end
#> 1   chr1     0   50
#> 2   chr1   100  150
#> 3   chr1   250  300
#> 4   chr1  1400 1450
#> 5   chr1  4750 4800
#> 6   chr1  6350 6400
#> 7   chr1  9900 9950
#> 8   chr2  2100 2150
#> 9   chr2  2300 2350
#> 10  chr2  2500 2550
#> 11  chr2  2900 2950
#> 12  chr2  3250 3300
gintervals.update("testintervs", NULL, chrom = 2)
gintervals.load("testintervs")
#>   chrom start  end
#> 1  chr1     0   50
#> 2  chr1   100  150
#> 3  chr1   250  300
#> 4  chr1  1400 1450
#> 5  chr1  4750 4800
#> 6  chr1  6350 6400
#> 7  chr1  9900 9950
gintervals.rm("testintervs", force = TRUE)
```
