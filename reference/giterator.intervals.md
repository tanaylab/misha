# Returns iterator intervals

Returns iterator intervals given track expression, scope, iterator and
band.

## Usage

``` r
giterator.intervals(
  expr = NULL,
  intervals = .misha$ALLGENOME,
  iterator = NULL,
  band = NULL,
  intervals.set.out = NULL,
  interval_relative = FALSE,
  partial_bins = c("clip", "exact", "drop")
)
```

## Arguments

- expr:

  track expression

- intervals:

  genomic scope

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expression.

- band:

  track expression band. If 'NULL' no band is used.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

- interval_relative:

  if TRUE, and iterator is numeric, bins start at each interval's start
  position instead of chromosome position 0. Returns intervalID column.
  Default: FALSE.

- partial_bins:

  how to handle partial bins at interval boundaries when
  interval_relative is TRUE. One of "clip" (default, truncate last bin
  to interval boundary), "exact" or "drop" (only output full-size bins).

## Value

If 'intervals.set.out' is 'NULL' a data frame representing iterator
intervals. When 'interval_relative' is TRUE, includes an 'intervalID'
column.

## Details

This function returns a set of intervals used by the iterator intervals
for the given track expression, genomic scope, iterator and band. Some
functions accept an iterator without accepting a track expression (like
'gtrack.create_pwm_energy'). These functions generate the values for
each iterator interval by themselves. Use set 'expr' to 'NULL' to
simulate the work of these functions.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

When 'interval_relative' is TRUE, bins are aligned to each input
interval's start position rather than chromosome position 0. This mode
requires a numeric iterator (binsize) and returns an additional
'intervalID' column indicating which input interval spawned each bin.

## See also

[`giterator.cartesian_grid`](https://tanaylab.github.io/misha/reference/giterator.cartesian_grid.md)

## Examples

``` r

gdb.init_examples()

## iterator is set implicitly to bin size of 'dense' track
giterator.intervals("dense_track", gintervals(1, 0, 200))
#>   chrom start end
#> 1  chr1     0  50
#> 2  chr1    50 100
#> 3  chr1   100 150
#> 4  chr1   150 200

## iterator = 30
giterator.intervals("dense_track", gintervals(1, 0, 200), 30)
#>   chrom start end
#> 1  chr1     0  30
#> 2  chr1    30  60
#> 3  chr1    60  90
#> 4  chr1    90 120
#> 5  chr1   120 150
#> 6  chr1   150 180
#> 7  chr1   180 200

## iterator is an intervals set named 'annotations'
giterator.intervals("dense_track", .misha$ALLGENOME, "annotations")
#>   chrom start   end
#> 1  chr1    20  2000
#> 2  chr1  2500  2600
#> 3  chr2    20  2000
#> 4  chr2  3000  8000
#> 5  chr2  9000 11000
#> 6  chr2 12000 12001
#> 7  chr2 13000 14000
#> 8  chr2 15000 15500

## iterator is set implicitly to intervals of 'array_track' track
giterator.intervals("array_track", gintervals(1, 0, 200))
#>   chrom start end
#> 1  chr1     0  50
#> 2  chr1   100 150

## iterator is a rectangle 100000 by 50000
giterator.intervals(
    "rects_track",
    gintervals.2d(chroms1 = 1, chroms2 = "chrX"),
    c(100000, 50000)
)
#>    chrom1 start1  end1 chrom2 start2   end2
#> 1    chr1  0e+00 1e+05   chrX      0  50000
#> 2    chr1  1e+05 2e+05   chrX      0  50000
#> 3    chr1  2e+05 3e+05   chrX      0  50000
#> 4    chr1  3e+05 4e+05   chrX      0  50000
#> 5    chr1  4e+05 5e+05   chrX      0  50000
#> 6    chr1  0e+00 1e+05   chrX  50000 100000
#> 7    chr1  1e+05 2e+05   chrX  50000 100000
#> 8    chr1  2e+05 3e+05   chrX  50000 100000
#> 9    chr1  3e+05 4e+05   chrX  50000 100000
#> 10   chr1  4e+05 5e+05   chrX  50000 100000
#> 11   chr1  0e+00 1e+05   chrX 100000 150000
#> 12   chr1  1e+05 2e+05   chrX 100000 150000
#> 13   chr1  2e+05 3e+05   chrX 100000 150000
#> 14   chr1  3e+05 4e+05   chrX 100000 150000
#> 15   chr1  4e+05 5e+05   chrX 100000 150000
#> 16   chr1  0e+00 1e+05   chrX 150000 200000
#> 17   chr1  1e+05 2e+05   chrX 150000 200000
#> 18   chr1  2e+05 3e+05   chrX 150000 200000
#> 19   chr1  3e+05 4e+05   chrX 150000 200000
#> 20   chr1  4e+05 5e+05   chrX 150000 200000

## interval_relative mode: bins aligned to each interval's start
intervs <- gintervals(1, c(100, 500), c(300, 700))
giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE)
#>   chrom start end intervalID
#> 1  chr1   100 150          1
#> 2  chr1   150 200          1
#> 3  chr1   200 250          1
#> 4  chr1   250 300          1
#> 5  chr1   500 550          2
#> 6  chr1   550 600          2
#> 7  chr1   600 650          2
#> 8  chr1   650 700          2
```
