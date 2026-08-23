# Divides track expression into segments

Divides the values of track expression into segments by using Wilcoxon
test.

## Usage

``` r
gsegment(
  expr = NULL,
  minsegment = NULL,
  maxpval = 0.05,
  onetailed = TRUE,
  intervals = NULL,
  iterator = NULL,
  intervals.set.out = NULL
)
```

## Arguments

- expr:

  track expression

- minsegment:

  minimal segment size

- maxpval:

  maximal P-value that separates two adjacent segments

- onetailed:

  if 'TRUE', Wilcoxon test is performed one tailed, otherwise two tailed

- intervals:

  genomic scope for which the function is applied

- iterator:

  track expression iterator of "fixed bin" type. If 'NULL' iterator is
  determined implicitly based on track expression.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a set of intervals where each interval
represents a segment.

## Details

This function divides the values of track expression into segments,
where each segment size is at least of 'minsegment' size and the P-value
of comparing the segment with the first 'minsegment' values from the
next segment is at most 'maxpval'. Comparison is done using Wilcoxon
(also known as Mann-Whitney) test.

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
  and [`gscreen`](https://tanaylab.github.io/misha/reference/gscreen.md)
  **drop** them: `NaN` bins are not counted into any distribution bin,
  do not contribute to a percentile, and never satisfy a screening
  condition - including a condition that would be true of every real
  value.

- `gsegment` **spans** them: a `NaN` bin contributes no evidence to the
  test that places a boundary, but it still falls inside whichever
  segment surrounds it, so the returned segments tile the scope
  continuously rather than skipping the gaps.

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

[`gscreen`](https://tanaylab.github.io/misha/reference/gscreen.md),
[`gwilcox`](https://tanaylab.github.io/misha/reference/gwilcox.md)

## Examples

``` r

gdb.init_examples()
gsegment("dense_track", 5000, 0.0001)
#>    chrom  start    end
#> 1   chr1      0   9200
#> 2   chr1   9200  31550
#> 3   chr1  31550  41350
#> 4   chr1  41350  65600
#> 5   chr1  65600  75800
#> 6   chr1  75800  99150
#> 7   chr1  99150 121800
#> 8   chr1 121800 129250
#> 9   chr1 129250 150650
#> 10  chr1 150650 157850
#> 11  chr1 157850 162900
#> 12  chr1 162900 224200
#> 13  chr1 224200 246200
#> 14  chr1 246200 309950
#> 15  chr1 309950 326350
#> 16  chr1 326350 346700
#> 17  chr1 346700 352200
#> 18  chr1 352200 362000
#> 19  chr1 362000 369400
#> 20  chr1 369400 379300
#> 21  chr1 379300 384800
#> 22  chr1 384800 426600
#> 23  chr1 426600 433950
#> 24  chr1 433950 500000
#> 25  chr2      0  28400
#> 26  chr2  28400  33600
#> 27  chr2  33600  38950
#> 28  chr2  38950  45800
#> 29  chr2  45800  50800
#> 30  chr2  50800  75600
#> 31  chr2  75600  85150
#> 32  chr2  85150  90600
#> 33  chr2  90600  99600
#> 34  chr2  99600 116150
#> 35  chr2 116150 122550
#> 36  chr2 122550 143400
#> 37  chr2 143400 154050
#> 38  chr2 154050 184350
#> 39  chr2 184350 194450
#> 40  chr2 194450 199600
#> 41  chr2 199600 207250
#> 42  chr2 207250 224400
#> 43  chr2 224400 229400
#> 44  chr2 229400 246000
#> 45  chr2 246000 276550
#> 46  chr2 276550 300000
#> 47  chrX      0  11300
#> 48  chrX  11300  16350
#> 49  chrX  16350  90000
#> 50  chrX  90000  95000
#> 51  chrX  95000 106200
#> 52  chrX 106200 116800
#> 53  chrX 116800 127850
#> 54  chrX 127850 135050
#> 55  chrX 135050 153450
#> 56  chrX 153450 200000
```
