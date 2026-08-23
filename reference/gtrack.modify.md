# Modifies track contents

Modifies 'Dense' track contents.

## Usage

``` r
gtrack.modify(track = NULL, expr = NULL, intervals = NULL)
```

## Arguments

- track:

  track name

- expr:

  track expression

- intervals:

  genomic scope for which track is modified

## Value

None.

## Details

This function modifies the contents of a 'Dense' track by the values of
'expr'. 'intervals' argument controls which portion of the track is
modified. The iterator policy is set internally to the bin size of the
track.

## See also

[`gtrack.create`](https://tanaylab.github.io/misha/reference/gtrack.create.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md)

## Examples

``` r

gdb.init_examples()
intervs <- gintervals(1, 300, 800)
gextract("dense_track", intervs)
#>    chrom start end dense_track intervalID
#> 1   chr1   300 350        0.16          1
#> 2   chr1   350 400        0.16          1
#> 3   chr1   400 450        0.16          1
#> 4   chr1   450 500        0.06          1
#> 5   chr1   500 550        0.06          1
#> 6   chr1   550 600        0.02          1
#> 7   chr1   600 650        0.04          1
#> 8   chr1   650 700        0.00          1
#> 9   chr1   700 750        0.00          1
#> 10  chr1   750 800        0.00          1
gtrack.modify("dense_track", "dense_track * 2", intervs)
gextract("dense_track", intervs)
#>    chrom start end dense_track intervalID
#> 1   chr1   300 350        0.32          1
#> 2   chr1   350 400        0.32          1
#> 3   chr1   400 450        0.32          1
#> 4   chr1   450 500        0.12          1
#> 5   chr1   500 550        0.12          1
#> 6   chr1   550 600        0.04          1
#> 7   chr1   600 650        0.08          1
#> 8   chr1   650 700        0.00          1
#> 9   chr1   700 750        0.00          1
#> 10  chr1   750 800        0.00          1
gtrack.modify("dense_track", "dense_track / 2", intervs)
```
