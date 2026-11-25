# Find PWM motif hits in DNA sequences

Returns every window of each sequence whose PWM score reaches
`score.thresh`: the windows `gseq.pwm(mode = "count")` counts, with
their positions, strand and score.

## Usage

``` r
gseq.pwm_hits(
  seqs,
  pssm,
  score.thresh,
  bidirect = TRUE,
  strand = 1L,
  start_pos = NULL,
  end_pos = NULL,
  extend = FALSE,
  prior = 0.01
)
```

## Arguments

- seqs:

  character vector of DNA sequences (A/C/G/T/N; case-insensitive)

- pssm:

  numeric matrix or data frame with columns named A, C, G, T (additional
  columns are allowed and will be ignored). Each row holds the
  probabilities or counts of the four bases at one motif position;
  log-odds matrices, `NA` values and rows summing to 0 are rejected.

- score.thresh:

  single number; windows scoring at or above it are hits.

- bidirect:

  logical; if TRUE, scans both strands (default: TRUE)

- strand:

  integer; `1` for the forward strand, `-1` for the reverse (default:
  `1`). Used only when `bidirect = FALSE`, and required there - `0` is
  refused, since exactly one strand is read and `0` cannot name it.
  Under `bidirect = TRUE` both strands are read and this is ignored.

- start_pos:

  integer or NULL; 1-based inclusive start of ROI (default: 1)

- end_pos:

  integer or NULL; 1-based inclusive end of ROI (default: sequence
  length)

- extend:

  logical or integer; how far past `end_pos` a window may run, as in
  [`gseq.pwm`](https://tanaylab.github.io/misha/reference/gseq.pwm.md).

- prior:

  numeric; pseudocount added to frequencies (default: 0.01). Set to 0
  for no pseudocounts.

## Value

A data frame with one row per hit and columns `seq_idx` (1-based index
into `seqs`), `start` and `end` (1-based, inclusive window positions on
the sequence), `strand` (`1` or `-1`) and `score`, ordered by `seq_idx`
and `start`. It has no rows when nothing passes.

## Details

Windows are scored exactly as
[`gseq.pwm`](https://tanaylab.github.io/misha/reference/gseq.pwm.md)
scores them in `"count"` mode with its default `neutral_chars` and no
spatial weights: a natural-log likelihood (including the prior), with a
window holding `N` or any other non-ACGT character left out. Gap
characters are not skipped - they are non-ACGT characters like any
other - so on a sequence without gaps the number of hits per sequence
equals `gseq.pwm(mode = "count")`'s count. Windows that cannot reach the
threshold are skipped without being scored, which makes a strict
threshold fast; the hits are the same as scoring every window.

Strand: under `bidirect = TRUE` a window is one hit, whose `score` is
the log-sum-exp of its forward and reverse-complement scores (the value
thresholded, as in `gseq.pwm` and `pwm.count`), and whose `strand` is
the strand that scores higher on its own (`1` on a tie), as
`gseq.pwm(mode = "pos")` names it. Under `bidirect = FALSE` only the
strand named by `strand` is read and `strand` is that strand. `start`
and `end` are the window's positions on the input sequence whatever the
strand.

## See also

[`gseq.pwm`](https://tanaylab.github.io/misha/reference/gseq.pwm.md),
[`gintervals.pwm_hits`](https://tanaylab.github.io/misha/reference/gintervals.pwm_hits.md)
for hits over genomic intervals.

## Examples

``` r
pssm <- matrix(
    c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7
    ),
    ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T"))
)
seqs <- c("TTACGTAAACGTTT", "GGGGGGGG")
gseq.pwm_hits(seqs, pssm, score.thresh = -3)
#>   seq_idx start end strand     score
#> 1       1     3   6      1 -0.833697
#> 2       1     9  12      1 -0.833697

# the same windows gseq.pwm(mode = "count") counts
gseq.pwm(seqs, pssm, mode = "count", score.thresh = -3)
#> [1] 2 0
```
