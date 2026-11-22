# Returns track expression values over equal-width intervals as a matrix

Evaluates track expressions in fixed bins over each interval, like
[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md)
with an iterator, and returns the values as a matrix with one row per
interval and one column per bin.

## Usage

``` r
gextract.matrix(..., intervals = NULL, iterator = 1, colnames = NULL)
```

## Arguments

- ...:

  track expressions, as in
  [`gextract`](https://tanaylab.github.io/misha/reference/gextract.md)

- intervals:

  1D intervals of equal width: a data frame or the name of an intervals
  set

- iterator:

  bin size: a positive integer that divides the width of the intervals

- colnames:

  names of the third dimension when there are several expressions. If
  'NULL' the names are the track expressions.

## Value

For one track expression, a numeric matrix with `nrow(intervals)` rows
and `W / iterator` columns. For several, a numeric array of dimensions
`nrow(intervals)`, `W / iterator` and the number of expressions.

## Details

All intervals must have the same width `W`, and `iterator` must divide
it. Row `i` of the result holds the values of row `i` of 'intervals', in
`W / iterator` columns ordered by genomic coordinate, left to right. A
'strand' column is ignored: columns are in genomic order on both
strands. Rows may be unsorted, duplicated, overlapping or on several
chromosomes; each one is extracted on its own.

The bins start at each interval's start, as in
[`giterator.intervals`](https://tanaylab.github.io/misha/reference/giterator.intervals.md)
with `interval_relative = TRUE`: column `j` of row `i` covers
`[start_i + (j - 1) * iterator, start_i + j * iterator)`. When every
start is a multiple of `iterator` these are the bins of gextract's
numeric iterator. As in gextract, an interval that extends past the end
of its chromosome is an error.

Each value is the one gextract returns for that bin as an iterator
interval, NaN included. For a dense track, a bin that overlaps several
of the track's bins gets the plain mean of their non-NaN values, not
weighted by the length of each overlap. The result is allocated once and
filled by calling gextract on consecutive chunks of 'intervals', so
memory use stays close to the size of the result.

## NaN values

A track expression evaluates to `NaN` wherever the iterator produces a
bin the track has no data for. What happens next depends on the
function:

- [`gextract`](https://tanaylab.github.io/misha/reference/gextract.md)
  **keeps** `NaN` rows, so the result has one row per iterator interval
  whether or not the track covered it.

- [`gsummary`](https://tanaylab.github.io/misha/reference/gsummary.md)
  **counts** them and reports the count as the "NaN intervals" element,
  while the statistics themselves are computed over the non-`NaN` values
  only.

- [`gdist`](https://tanaylab.github.io/misha/reference/gdist.md),
  [`gquantiles`](https://tanaylab.github.io/misha/reference/gquantiles.md)
  and [`gscreen`](https://tanaylab.github.io/misha/reference/gscreen.md)
  **drop** them: `NaN` bins are not counted into any distribution bin,
  do not contribute to a percentile, and never satisfy a screening
  condition - including a condition that would be true of every real
  value.

- [`gsegment`](https://tanaylab.github.io/misha/reference/gsegment.md)
  **spans** them: a `NaN` bin contributes no evidence to the test that
  places a boundary, but it still falls inside whichever segment
  surrounds it, so the returned segments tile the scope continuously
  rather than skipping the gaps.

So on 20 bins of which 7 are `NaN`, `gextract` returns 20 rows,
`gsummary` reports 20 total and 7 `NaN`, and `gdist` counts 13; and on a
300 kb scope where 120 of 300 bins are `NaN`, `gsegment` still returns
segments covering the full 300 kb.

The practical consequence is that `NaN` and zero are different, and
collapsing them with `ifelse(is.na(x), 0, x)` turns "no data here" into
a measured value of zero. Where that is genuinely what you want, note
that it also changes every mean, quantile and distribution computed
downstream.

## See also

[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md),
[`giterator.intervals`](https://tanaylab.github.io/misha/reference/giterator.intervals.md)

## Examples

``` r

gdb.init_examples()
intervs <- gintervals(c(1, 2), c(1000, 2000), c(1100, 2100))

## per-base values: a 2 x 100 matrix
m <- gextract.matrix("dense_track", intervals = intervs)
dim(m)
#> [1]   2 100

## a 0/1 mask of the intervals of 'annotations' inside each interval
gvtrack.create("annot_cov", "annotations", func = "coverage")
mask <- gextract.matrix("annot_cov", intervals = intervs)
rowSums(mask)
#> [1] 100   0

## two expressions in 10 bp bins: a 2 x 10 x 2 array
a <- gextract.matrix("dense_track", "sparse_track",
    intervals = intervs, iterator = 10,
    colnames = c("dense", "sparse")
)
a[, , "dense"]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,] 0.10 0.10 0.10 0.10 0.10 0.00 0.00 0.00 0.00  0.00
#> [2,] 0.08 0.08 0.08 0.08 0.08 0.14 0.14 0.14 0.14  0.14

## 32 bp bins counted from each start: a 2 x 3 matrix
gextract.matrix("dense_track",
    intervals = gintervals(c(1, 2), c(1005, 2017), c(1101, 2113)),
    iterator = 32
)
#>      [,1] [,2] [,3]
#> [1,] 0.10 0.05 0.01
#> [2,] 0.08 0.11 0.19
```
