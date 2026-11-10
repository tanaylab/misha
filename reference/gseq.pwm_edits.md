# Show optimal edits to reach a PWM score threshold

For each input sequence (or genomic interval), finds the optimal motif
window and the specific base changes needed to reach `score.thresh`.
Returns a long-format data frame with one row per edit.

## Usage

``` r
gseq.pwm_edits(
  seqs,
  pssm,
  score.thresh,
  max_edits = NULL,
  max_indels = NULL,
  bidirect = TRUE,
  prior = 0.01,
  score.min = NULL,
  score.max = NULL,
  extend = TRUE,
  strand = 1L,
  direction = "above"
)
```

## Arguments

- seqs:

  character vector of DNA sequences, OR a data frame of genomic
  intervals (with columns chrom, start, end). When intervals are
  provided, sequences are extracted automatically via `gseq.extract`,
  and the windows scanned are those the `pwm.edit_distance` virtual
  track scans for the same interval. A `strand` column is ignored, as by
  the virtual track: the strand scanned is set by `bidirect` and
  `strand`.

- pssm:

  numeric matrix or data frame with columns A, C, G, T. Each row is a
  motif position and holds the probabilities or counts of the four
  bases; log-odds matrices, `NA` values and rows summing to 0 are
  rejected.

- score.thresh:

  single number; target PWM log-likelihood score to reach.

- max_edits:

  integer or NULL; maximum number of edits to search. NULL means no cap.
  Default NULL.

- max_indels:

  integer or NULL; maximum number of insertions and deletions allowed
  (default NULL, substitutions only). When \> 0, a banded
  Needleman-Wunsch DP is used. Edits are reported with `edit_type`
  "sub", "ins", or "del". An inserted base is the best its motif column
  allows (the column maximum for `"above"`, the minimum for `"below"`)
  and is scored as such.

- bidirect:

  logical; scan both strands? Default TRUE.

- prior:

  numeric; pseudocount for PSSM frequencies. Default 0.01.

- score.min:

  numeric or NULL; skip windows with PWM score below this. Default NULL
  (no filter). The score is the window's ungapped score (the PSSM-length
  window from its start), also with `max_indels > 0`; a window without
  one (a non-ACGT base among its first PSSM-length bases) is skipped
  when a filter is set.

- score.max:

  numeric or NULL; skip windows with PWM score above this. Default NULL
  (no filter). Same score as `score.min`.

- extend:

  logical or integer; extend sequence for boundary motifs. Default TRUE.
  For interval input it means what it means for the `pwm.edit_distance`
  virtual track: `TRUE` or a positive integer lets a window start at
  every position of the interval and read bases past its end; `FALSE` or
  0 keeps every window (deletions included) inside the interval. For a
  sequence, its end acts as a chromosome end: with `max_indels > 0`,
  `TRUE` lets a window start wherever at least `max(1, w - max_indels)`
  bases remain (insertions supply the rest), as the virtual track does
  at a chromosome end, and an integer lets a window run at most that
  many bases past the end.

- strand:

  integer; which strand to scan when `bidirect=FALSE`. 1=forward,
  -1=reverse. Default 1.

- direction:

  character; direction of the edit distance query. `"above"` (default)
  finds minimum edits to bring score above `score.thresh`; `"below"`
  finds minimum edits to bring score below `score.thresh`. With
  `"below"` and `bidirect = TRUE`, as in the `pwm.edit_distance` virtual
  track, both strand readings of a window must end up below the
  threshold, so the strand needing more edits decides `n_edits` (a
  strand that cannot be brought below makes the window unreachable); the
  reported edits and `strand` are those of that strand's reading, and
  the edits are what disrupts that reading.

## Value

A data frame (long format) with one row per edit, containing columns:

- seq_idx:

  Index into input sequences/intervals (1-based)

- strand:

  +1 (forward) or -1 (reverse strand)

- window_start:

  1-based position of optimal window within sequence

- score_before:

  PWM score before edits

- score_after:

  PWM score after all edits

- n_edits:

  Total number of edits needed

- edit_num:

  Which edit this row represents (1, 2, ...)

- motif_col:

  1-based position within the motif where the edit occurs

- ref:

  Current base at this position

- alt:

  Suggested replacement base

- gain:

  Score improvement from this individual edit

- window_seq:

  Motif-length sequence at the optimal window (as seen by PSSM,
  reverse-complemented if on reverse strand)

- mutated_seq:

  Same sequence with all edits applied

When intervals are provided, additional columns `chrom`, `start`, `end`
are included.

Sequences already above the threshold produce a single row with
`n_edits = 0`. Unreachable sequences are omitted from the output.

## Details

This is an investigation tool: use it on a small set of positions (e.g.,
from `gscreen`) to see what mutations would activate latent binding
sites.

## See also

[`gseq.pwm`](https://tanaylab.github.io/misha/reference/gseq.pwm.md),
[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)

## Examples

``` r

gdb.init_examples()

# Simple PSSM
pssm <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0),
    nrow = 2,
    dimnames = list(NULL, c("A", "C", "G", "T"))
)

# What edits are needed?
gseq.pwm_edits("CCGTACGT", pssm, score.thresh = -0.5, prior = 0)
#>   seq_idx strand window_start score_before score_after n_edits edit_num
#> 1       1     -1            1         -Inf           0       1        1
#>   motif_col ref alt gain edit_type window_seq mutated_seq
#> 1         1   G   A  Inf       sub         GG          AG
```
