# Score DNA sequences with a PWM over a region of interest

Scores full DNA sequences using a Position Weight Matrix (PWM) over a
specified region of interest (ROI). The ROI is defined by `start_pos`
and `end_pos` (1-based, inclusive), with optional extension controlled
by `extend`. All reported positions are on the full input sequence.

## Usage

``` r
gseq.pwm(
  seqs,
  pssm,
  mode = c("lse", "max", "pos", "count"),
  bidirect = TRUE,
  strand = 1L,
  score.thresh = NULL,
  start_pos = NULL,
  end_pos = NULL,
  extend = FALSE,
  spat.factor = NULL,
  spat.bin = 1L,
  spat.min = NULL,
  spat.max = NULL,
  return_strand = FALSE,
  skip_gaps = TRUE,
  gap_chars = c("-", "."),
  neutral_chars = c("N", "n", "*"),
  neutral_chars_policy = c("average", "log_quarter", "na"),
  prior = 0.01
)
```

## Arguments

- seqs:

  character vector of DNA sequences (A/C/G/T/N; case-insensitive)

- pssm:

  numeric matrix or data frame with columns named A, C, G, T (additional
  columns are allowed and will be ignored)

- mode:

  character; one of "lse", "max", "pos", or "count"

- bidirect:

  logical; if TRUE, scans both strands (default: TRUE)

- strand:

  integer; `1` for the forward strand, `-1` for the reverse (default:
  `1`). Used only when `bidirect = FALSE`, and required there - `0` is
  refused, since exactly one strand is read and `0` cannot name it.
  Under `bidirect = TRUE` both strands are read and this is ignored.

- score.thresh:

  single number; windows scoring at or above this value are counted.
  Required when `mode="count"` and ignored otherwise. PWM scores are
  log-likelihoods, so the usable range depends on the PSSM, the prior
  and any spatial weights - calibrate with `mode="max"` before choosing
  one.

- start_pos:

  integer or NULL; 1-based inclusive start of ROI (default: 1)

- end_pos:

  integer or NULL; 1-based inclusive end of ROI (default: sequence
  length)

- extend:

  logical or integer; extension of allowed window starts (default:
  FALSE)

- spat.factor:

  numeric vector; spatial weighting factors (optional)

- spat.bin:

  integer; bin size for spatial weighting

- spat.min:

  numeric; start of scanning window

- spat.max:

  numeric; end of scanning window

- return_strand:

  logical; if TRUE and `mode="pos"`, returns data.frame with `pos` and
  `strand` columns

- skip_gaps:

  logical; if TRUE, treat gap characters as holes and skip them while
  scanning. Windows are w consecutive non-gap bases (default: TRUE)

- gap_chars:

  character vector; which characters count as gaps (default: c("-",
  "."))

- neutral_chars:

  character vector; bases treated as unknown and scored with the average
  log probability per position (default: c("N", "n", "\*"))

- neutral_chars_policy:

  character string; how to treat neutral characters. One of `"average"`
  (default; use the column's mean log-probability), `"log_quarter"`
  (always use `log(1/4)`), or `"na"` (return NA when a neutral character
  is encountered in the scanning window).

- prior:

  numeric; pseudocount added to frequencies (default: 0.01). Set to 0
  for no pseudocounts.

## Value

Numeric vector (for "lse"/"max"/"count" modes), integer vector (for
"pos" mode), or data.frame with `pos` and `strand` columns (for "pos"
mode with `return_strand=TRUE`). Returns NA when no valid windows exist.

## Details

This function scores DNA sequences directly without requiring a genomics
database. For detailed documentation on PWM scoring modes, parameters,
and spatial weighting, see
[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)
(functions "pwm", "pwm.max", "pwm.max.pos", "pwm.count").

The ROI (region of interest) is defined by `start_pos` and `end_pos`.
The `extend` parameter controls whether motif matches can extend beyond
the ROI boundaries.

When `skip_gaps=TRUE`, characters specified in `gap_chars` are treated
as gaps. Windows are defined as w consecutive non-gap bases. All
positions (`pos`) are reported as 1-based indices on the original full
sequence (including gaps). `start_pos` and `end_pos` are interpreted as
physical coordinates on the full sequence.

Neutral characters (`neutral_chars`, default `c("N", "n", "*")`) are
treated as unknown bases in both orientations. Each neutral contributes
the mean log-probability of the corresponding PSSM column, yielding
identical penalties on forward and reverse strands without hard-coded
background scores.

Under `bidirect = TRUE` the forward and reverse-complement matches are
two readings of the same window, and every mode but `"pos"` combines
them at each window by log-sum-exp, after any spatial weights are
applied: `"lse"` then log-sum-exps those unions across the windows,
`"max"` takes the largest of them, and `"count"` thresholds each one
once, so a window is at most a single hit however well its reverse
complement matches it. That is the same union the `pwm`, `pwm.max` and
`pwm.count` virtual tracks take, and over the same windows `gseq.pwm()`
agrees with them to the precision the two share - the virtual tracks
accumulate in single precision and `gseq.pwm()` in double, so expect
agreement to about `1e-5` relative, not to the last bit. `mode = "pos"`
instead takes the better of the two strands, because it has to name
one - what the `pwm.max.pos` virtual track does - so `"max"` and `"pos"`
can select different windows, and the value from `"max"` is not
necessarily the score at the position from `"pos"`. The asymmetry is the
`pwm` family's own (`pwm.max` combines the strands, `pwm.max.pos` does
not) and is mirrored here rather than fixed on one side. With
`bidirect = FALSE` only the strand named by `strand` is read and there
is nothing to combine.

## See also

[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)
for detailed PWM parameter documentation

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a PSSM (position-specific scoring matrix) with frequency values
pssm <- matrix(
    c(
        0.7, 0.1, 0.1, 0.1, # Position 1: mostly A
        0.1, 0.7, 0.1, 0.1, # Position 2: mostly C
        0.1, 0.1, 0.7, 0.1, # Position 3: mostly G
        0.1, 0.1, 0.1, 0.7 # Position 4: mostly T
    ),
    ncol = 4, byrow = TRUE
)
colnames(pssm) <- c("A", "C", "G", "T")

# Example sequences
seqs <- c("ACGTACGTACGT", "GGGGACGTCCCC", "TTTTTTTTTTT")

# Score sequences using log-sum-exp (default mode)
gseq.pwm(seqs, pssm, mode = "lse")

# Get maximum score
gseq.pwm(seqs, pssm, mode = "max")

# Find position of best match
gseq.pwm(seqs, pssm, mode = "pos")

# Find position with strand information
gseq.pwm(seqs, pssm, mode = "pos", bidirect = TRUE, return_strand = TRUE)

# Count matches above threshold. score.thresh is mandatory in "count" mode
# and has no default: PWM scores are log-likelihoods, so calibrate it
# against the score range of your own PSSM.
range(gseq.pwm(seqs, pssm, mode = "max"))
gseq.pwm(seqs, pssm, mode = "count", score.thresh = -3)

# Score only a region of interest
gseq.pwm(seqs, pssm, mode = "max", start_pos = 3, end_pos = 10)

# Allow matches to extend beyond ROI boundaries
gseq.pwm(seqs, pssm, mode = "count", score.thresh = -3, start_pos = 5, end_pos = 8, extend = TRUE)

# Spatial weighting example: higher weight in the center
spatial_weights <- c(0.5, 1.0, 2.0, 1.0, 0.5)
gseq.pwm(seqs, pssm,
    mode = "lse",
    spat.factor = spatial_weights,
    spat.bin = 2
)
} # }
```
