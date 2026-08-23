# Calculates quantiles of a track expression

Calculates the quantiles of a track expression for the given
percentiles.

## Usage

``` r
gquantiles(
  expr = NULL,
  percentiles = 0.5,
  intervals = get("ALLGENOME", envir = .misha),
  iterator = NULL,
  band = NULL
)
```

## Arguments

- expr:

  track expression

- percentiles:

  an array of percentiles of quantiles in \[0, 1\] range

- intervals:

  genomic scope for which the function is applied

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expression.

- band:

  track expression band. If 'NULL' no band is used.

## Value

An array that represent quantiles.

## Details

This function calculates the quantiles for the given percentiles.

If data size exceeds the limit (see: 'getOption(gmax.data.size)'), the
data is randomly sampled to fit the limit. A warning message is
generated. Call [`set.seed()`](https://rdrr.io/r/base/Random.html)
before this function to make the sample reproducible.

Note: this function is capable to run in multitasking mode. Sampling may
vary according to the extent of multitasking. Since multitasking depends
on the number of available CPU cores, running the function on two
different machines might give different results. Please switch off
multitasking if you want to achieve identical results on any machine.
For more information regarding multitasking please refer "User Manual".

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
  `gquantiles` and
  [`gscreen`](https://tanaylab.github.io/misha/reference/gscreen.md)
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

[`gbins.quantiles`](https://tanaylab.github.io/misha/reference/gbins.quantiles.md),
[`gintervals.quantiles`](https://tanaylab.github.io/misha/reference/gintervals.quantiles.md),
[`gdist`](https://tanaylab.github.io/misha/reference/gdist.md)

## Examples

``` r

gdb.init_examples()
gquantiles("dense_track", c(0.1, 0.6, 0.8), gintervals(c(1, 2)))
#>  0.1  0.6  0.8 
#> 0.02 0.10 0.14 
```
