# Find PWM motif hits in genomic intervals

Returns, for each interval, the motif windows that reach `score.thresh`:
exactly the windows a `pwm.count` virtual track with the same parameters
counts over that interval, as intervals with their strand and score. It
replaces running
[`gscreen`](https://tanaylab.github.io/misha/reference/gscreen.md) over
a 1 bp iterator.

## Usage

``` r
gintervals.pwm_hits(
  intervals = NULL,
  pssm = NULL,
  score.thresh = NULL,
  bidirect = TRUE,
  strand = 1L,
  prior = 0.01,
  extend = TRUE
)
```

## Arguments

- intervals:

  genomic intervals to scan

- pssm:

  numeric matrix or data frame with columns A, C, G, T; each row holds
  the probabilities or counts of the four bases at one motif position.

- score.thresh:

  single number; windows scoring at or above it are hits.

- bidirect:

  logical; if TRUE (default) both strands are scanned.

- strand:

  `1` or `-1`: the strand scanned when `bidirect = FALSE`.

- prior:

  numeric; pseudocount added to the frequencies (default 0.01).

- extend:

  logical; if TRUE (default) a window anchored near an interval's end
  may read bases past it.

## Value

An intervals data frame with columns `chrom`, `start`, `end` (the
window), `strand`, `score` and `intervalID` (the row of `intervals` the
window is anchored in), sorted by `intervalID` and `start`; `NULL` when
nothing passes. Overlapping input intervals report a window once per
interval it is anchored in.

## Details

A window is anchored at its first (leftmost) base, which must lie in the
interval; with `extend = TRUE` it may run up to w - 1 bases past the
interval's end (w being the motif length), as the PWM virtual tracks do.
Scores are natural-log likelihoods including the prior, computed as the
`pwm` virtual tracks compute them (in single precision), and a window
holding `N` or any other non-ACGT base is not scored. Windows that
cannot reach the threshold are skipped without being scored, which makes
a strict threshold fast; the hits are the same as scoring every window.

Strand: under `bidirect = TRUE` a window is one hit, whose `score` is
the log-sum-exp of its forward and reverse-complement scores (the value
`pwm.count` thresholds), and whose `strand` is the strand that scores
higher on its own (`1` on a tie), as `pwm.max.pos` names it. Under
`bidirect = FALSE` only the strand named by `strand` is read and
`strand` is that strand. A `strand` column of `intervals` is ignored, as
by the virtual tracks. `start` and `end` are genomic coordinates of the
window whatever the strand.

## See also

[`gseq.pwm_hits`](https://tanaylab.github.io/misha/reference/gseq.pwm_hits.md),
[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)
for the `pwm.count` virtual track.

## Examples

``` r

gdb.init_examples()
pssm <- matrix(
    c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7
    ),
    ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T"))
)
intervs <- gintervals(1, 0, 10000)
hits <- gintervals.pwm_hits(intervs, pssm, score.thresh = -3)
head(hits)
#>   chrom start end strand     score intervalID
#> 1  chr1   176 180      1 -2.698482          1
#> 2  chr1   504 508      1 -2.698482          1
#> 3  chr1   561 565      1 -2.698482          1
#> 4  chr1   607 611      1 -2.698482          1
#> 5  chr1   664 668      1 -2.698482          1
#> 6  chr1   674 678      1 -2.698482          1

# the windows pwm.count counts
gvtrack.create("motif_count", NULL, "pwm.count", pssm = pssm, score.thresh = -3)
gextract("motif_count", intervs, iterator = intervs)
#>   chrom start   end motif_count intervalID
#> 1  chr1     0 10000         269          1
```
