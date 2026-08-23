# Track-parallel gextract via mclapply

Splits `tracks` into chunks across at most `getOption("gmax.processes")`
worker processes. Each worker runs gextract on its track subset with
`gmultitasking=FALSE` so misha's tile-parallel multitask doesn't nest.
Results are merged column-wise: interval/intervalID columns come from
the first worker, value columns are cbind'd from each.

## Usage

``` r
.gextract_track_parallel(
  intervals,
  tracks,
  colnames,
  iterator,
  band,
  file,
  intervals.set.out,
  envir
)
```
