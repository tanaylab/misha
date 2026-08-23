# Creates a cartesian-grid iterator

Creates a cartesian grid two-dimensional iterator that can be used by
any function that accepts an iterator argument.

## Usage

``` r
giterator.cartesian_grid(
  intervals1 = NULL,
  expansion1 = NULL,
  intervals2 = NULL,
  expansion2 = NULL,
  min.band.idx = NULL,
  max.band.idx = NULL
)
```

## Arguments

- intervals1:

  one-dimensional intervals

- expansion1:

  an array of integers that define expansion around intervals1 centers

- intervals2:

  one-dimensional intervals. If 'NULL' then 'intervals2' is considered
  to be equal to 'intervals1'

- expansion2:

  an array of integers that define expansion around intervals2 centers.
  If 'NULL' then 'expansion2' is considered to be equal to 'expansion1'

- min.band.idx, max.band.idx:

  integers that limit iterator intervals to band

## Value

A list containing the definition of cartesian iterator.

## Details

This function creates and returns a cartesian grid two-dimensional
iterator that can be used by any function that accepts an iterator
argument.

Assume 'centers1' and 'centers2' to be the central points of each
interval from 'intervals1' and 'intervals2', and 'C1', 'C2' to be two
points from 'centers1', 'centers2' accordingly. Assume also that the
values in 'expansion1' and 'expansion2' are unique and sorted.

'giterator.cartesian_grid' creates a set of all possible unique and
non-overlapping two-dimensional intervals of form: '(chrom1, start1,
end1, chrom2, start2, end2)'. Each '(chrom1, start1, end1)' is created
by taking a point 'C1' - '(chrom1, coord1)' and converting it to
'start1' and 'end1' such that 'start1 == coord1+E1\[i\]', 'end1 ==
coord1+E1\[i+1\]', where 'E1\[i\]' is one of the sorted 'expansion1'
values. Overlaps between rectangles or expansion beyond the limits of
chromosome are avoided.

'min.band.idx' and 'max.band.idx' parameters control whether a pair of
'C1' and 'C2' is skipped or not. If both of these parameters are not
'NULL' AND if both 'C1' and 'C2' share the same chromosome AND the delta
of indices of 'C1' and 'C2' ('C1 index - C2 index') lays within
'\[min.band.idx, max.band.idx\]' range - only then the pair will be used
to create the intervals. Otherwise 'C1-C2' pair is filtered out. Note:
if 'min.band.idx' and 'max.band.idx' are not 'NULL', i.e. band indices
filtering is applied, then 'intervals2' parameter must be set to 'NULL'.

## See also

[`giterator.intervals`](https://tanaylab.github.io/misha/reference/giterator.intervals.md)

## Examples

``` r

gdb.init_examples()

intervs1 <- gintervals(
    c(1, 1, 2), c(100, 300, 200),
    c(300, 500, 300)
)
intervs2 <- gintervals(
    c(1, 2, 2), c(400, 1000, 3000),
    c(800, 2000, 4000)
)
itr <- giterator.cartesian_grid(
    intervs1, c(-20, 100), intervs2,
    c(-40, -10, 50)
)
giterator.intervals(iterator = itr)
#>    chrom1 start1 end1 chrom2 start2 end2
#> 1    chr1    180  300   chr1    560  590
#> 2    chr1    180  300   chr1    590  650
#> 3    chr1    380  500   chr1    560  590
#> 4    chr1    380  500   chr1    590  650
#> 5    chr1    180  300   chr2   1460 1490
#> 6    chr1    180  300   chr2   1490 1550
#> 7    chr1    180  300   chr2   3460 3490
#> 8    chr1    180  300   chr2   3490 3550
#> 9    chr1    380  500   chr2   1460 1490
#> 10   chr1    380  500   chr2   1490 1550
#> 11   chr1    380  500   chr2   3460 3490
#> 12   chr1    380  500   chr2   3490 3550
#> 13   chr2    230  350   chr1    560  590
#> 14   chr2    230  350   chr1    590  650
#> 15   chr2    230  350   chr2   1460 1490
#> 16   chr2    230  350   chr2   1490 1550
#> 17   chr2    230  350   chr2   3460 3490
#> 18   chr2    230  350   chr2   3490 3550

itr <- giterator.cartesian_grid(intervs1, c(-20, 50, 100))
giterator.intervals(iterator = itr)
#>    chrom1 start1 end1 chrom2 start2 end2
#> 1    chr1    180  250   chr1    180  250
#> 2    chr1    180  250   chr1    250  300
#> 3    chr1    250  300   chr1    180  250
#> 4    chr1    250  300   chr1    250  300
#> 5    chr1    180  250   chr1    380  450
#> 6    chr1    180  250   chr1    450  500
#> 7    chr1    250  300   chr1    380  450
#> 8    chr1    250  300   chr1    450  500
#> 9    chr1    380  450   chr1    180  250
#> 10   chr1    380  450   chr1    250  300
#> 11   chr1    450  500   chr1    180  250
#> 12   chr1    450  500   chr1    250  300
#> 13   chr1    380  450   chr1    380  450
#> 14   chr1    380  450   chr1    450  500
#> 15   chr1    450  500   chr1    380  450
#> 16   chr1    450  500   chr1    450  500
#> 17   chr1    180  250   chr2    230  300
#> 18   chr1    180  250   chr2    300  350
#> 19   chr1    250  300   chr2    230  300
#> 20   chr1    250  300   chr2    300  350
#> 21   chr1    380  450   chr2    230  300
#> 22   chr1    380  450   chr2    300  350
#> 23   chr1    450  500   chr2    230  300
#> 24   chr1    450  500   chr2    300  350
#> 25   chr2    230  300   chr1    180  250
#> 26   chr2    230  300   chr1    250  300
#> 27   chr2    300  350   chr1    180  250
#> 28   chr2    300  350   chr1    250  300
#> 29   chr2    230  300   chr1    380  450
#> 30   chr2    230  300   chr1    450  500
#> 31   chr2    300  350   chr1    380  450
#> 32   chr2    300  350   chr1    450  500
#> 33   chr2    230  300   chr2    230  300
#> 34   chr2    230  300   chr2    300  350
#> 35   chr2    300  350   chr2    230  300
#> 36   chr2    300  350   chr2    300  350

itr <- giterator.cartesian_grid(intervs1, c(-20, 50, 100),
    min.band.idx = -1,
    max.band.idx = 0
)
giterator.intervals(iterator = itr)
#>    chrom1 start1 end1 chrom2 start2 end2
#> 1    chr1    180  250   chr1    180  250
#> 2    chr1    180  250   chr1    250  300
#> 3    chr1    250  300   chr1    180  250
#> 4    chr1    250  300   chr1    250  300
#> 5    chr1    180  250   chr1    380  450
#> 6    chr1    180  250   chr1    450  500
#> 7    chr1    250  300   chr1    380  450
#> 8    chr1    250  300   chr1    450  500
#> 9    chr1    380  450   chr1    380  450
#> 10   chr1    380  450   chr1    450  500
#> 11   chr1    450  500   chr1    380  450
#> 12   chr1    450  500   chr1    450  500
#> 13   chr2    230  300   chr2    230  300
#> 14   chr2    230  300   chr2    300  350
#> 15   chr2    300  350   chr2    230  300
#> 16   chr2    300  350   chr2    300  350
```
