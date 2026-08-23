# Finds intervals that match track expression

Finds all intervals where track expression is 'TRUE'.

## Usage

``` r
gscreen(
  expr = NULL,
  intervals = NULL,
  iterator = NULL,
  band = NULL,
  intervals.set.out = NULL
)
```

## Arguments

- expr:

  logical track expression

- intervals:

  genomic scope for which the function is applied

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expression.

- band:

  track expression band. If 'NULL' no band is used.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a set of intervals that match track
expression.

## Details

This function finds all intervals where track expression's value is
'TRUE'.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

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
  and `gscreen` **drop** them: `NaN` bins are not counted into any
  distribution bin, do not contribute to a percentile, and never satisfy
  a screening condition - including a condition that would be true of
  every real value.

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

[`gsegment`](https://tanaylab.github.io/misha/reference/gsegment.md),
[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md)

## Examples

``` r

gdb.init_examples()
gscreen("dense_track > 0.2 & sparse_track < 0.4",
    iterator = "dense_track"
)
#>    chrom  start    end
#> 1   chr1  34850  34950
#> 2   chr1  47500  47550
#> 3   chr1  65100  65200
#> 4   chr1  65850  65950
#> 5   chr1  66050  66100
#> 6   chr1  92650  92700
#> 7   chr1  95350  95400
#> 8   chr1 100650 100700
#> 9   chr1 139050 139100
#> 10  chr1 152600 152650
#> 11  chr1 245500 245550
#> 12  chr1 341100 341150
#> 13  chr1 371450 371500
#> 14  chr1 373000 373100
#> 15  chr1 383550 383600
#> 16  chr1 459400 459450
#> 17  chr2  19200  19250
#> 18  chr2  22100  22150
#> 19  chr2  34250  34300
#> 20  chr2  34500  34550
#> 21  chr2  37400  37450
#> 22  chr2  50700  50750
#> 23  chr2  64050  64100
#> 24  chr2  69200  69250
#> 25  chr2  83550  83600
#> 26  chr2 100800 100850
#> 27  chr2 101450 101500
#> 28  chr2 139600 139650
#> 29  chr2 168600 168650
#> 30  chr2 233050 233100
#> 31  chr2 245200 245250
#> 32  chr2 282300 282350
#> 33  chrX  85150  85200
```
