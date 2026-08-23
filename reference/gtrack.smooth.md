# Creates a new track from smoothed values of track expression

Creates a new track from smoothed values of track expression.

## Usage

``` r
gtrack.smooth(
  track = NULL,
  description = NULL,
  expr = NULL,
  winsize = NULL,
  weight_thr = 0,
  smooth_nans = FALSE,
  alg = "LINEAR_RAMP",
  iterator = NULL
)
```

## Arguments

- track:

  track name

- description:

  a character string description

- expr:

  track expression

- winsize:

  size of smoothing window

- weight_thr:

  smoothing weight threshold

- smooth_nans:

  if 'FALSE' track value is always set to 'NaN' if central window value
  is 'NaN', otherwise it is calculated from the rest of non 'NaN' values

- alg:

  smoothing algorithm - "MEAN" or "LINEAR_RAMP"

- iterator:

  track expression iterator of 'Fixed bin' type

## Value

None.

## Details

This function creates a new 'Dense' track named 'track'. The values of
the track are results of smoothing the values of 'expr'.

Each track value at coordinate 'C' is determined by smoothing non 'NaN'
values of 'expr' over the window around 'C'. The window size is
controlled by 'winsize' and is given in coordinate units (not in number
of bins), defining the total regions to be considered when smoothing (on
both sides of the central point). Two different algorithms can be used
for smoothing:

"MEAN" - an arithmetic average.

"LINEAR_RAMP" - a weighted arithmetic average, where the weights
linearly decrease as the distance from the center of the window
increases.

'weight_thr' determines the function behavior when some of the values in
the window are missing or 'NaN' (missing values may occur at the edges
of each chromosome when the window covers an area beyond chromosome
boundaries). 'weight_thr' sets the weight sum threshold below which
smoothing algorithm returns 'NaN' rather than a smoothing value based on
non 'NaN' values in the window.

'smooth_nans' controls what would be the smoothed value if the central
value in the window is 'NaN'. If 'smooth_nans' is 'FALSE' then the
smoothed value is set to 'NaN' regardless of 'weight_thr' parameter.
Otherwise it is calculated normally.

'description' is added as a track attribute.

Iterator policy must be of "fixed bin" type.

## See also

[`gtrack.create`](https://tanaylab.github.io/misha/reference/gtrack.create.md),
[`gtrack.2d.create`](https://tanaylab.github.io/misha/reference/gtrack.2d.create.md),
[`gtrack.create_sparse`](https://tanaylab.github.io/misha/reference/gtrack.create_sparse.md),
[`gtrack.modify`](https://tanaylab.github.io/misha/reference/gtrack.modify.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md)

## Examples

``` r

gdb.init_examples()
gtrack.smooth("smoothed_track", "Test track", "dense_track", 500)
gextract("dense_track", "smoothed_track", gintervals(1, 0, 1000))
#>    chrom start  end dense_track smoothed_track intervalID
#> 1   chr1     0   50   0.1777778     0.17079365          1
#> 2   chr1    50  100   0.1600000     0.17034188          1
#> 3   chr1   100  150   0.1800000     0.17037037          1
#> 4   chr1   150  200   0.1600000     0.16949494          1
#> 5   chr1   200  250   0.1600000     0.16615874          1
#> 6   chr1   250  300   0.2000000     0.16049382          1
#> 7   chr1   300  350   0.1600000     0.14888889          1
#> 8   chr1   350  400   0.1600000     0.13444445          1
#> 9   chr1   400  450   0.1600000     0.11555555          1
#> 10  chr1   450  500   0.0600000     0.09277777          1
#> 11  chr1   500  550   0.0600000     0.07111111          1
#> 12  chr1   550  600   0.0200000     0.05055555          1
#> 13  chr1   600  650   0.0400000     0.03500000          1
#> 14  chr1   650  700   0.0000000     0.02277778          1
#> 15  chr1   700  750   0.0000000     0.01611111          1
#> 16  chr1   750  800   0.0000000     0.01666667          1
#> 17  chr1   800  850   0.0000000     0.01888889          1
#> 18  chr1   850  900   0.0200000     0.02333333          1
#> 19  chr1   900  950   0.0400000     0.02722222          1
#> 20  chr1   950 1000   0.0400000     0.03000000          1
gtrack.rm("smoothed_track", force = TRUE)
```
