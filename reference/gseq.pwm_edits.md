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
  provided, sequences are extracted automatically via `gseq.extract`.

- pssm:

  numeric matrix or data frame with columns A, C, G, T. Each row is a
  motif position.

- score.thresh:

  single number; target PWM log-likelihood score to reach.

- max_edits:

  integer or NULL; maximum number of edits to search. NULL means no cap.
  Default NULL.

- max_indels:

  integer or NULL; maximum number of insertions and deletions allowed
  (default NULL, substitutions only). When \> 0, a banded
  Needleman-Wunsch DP is used. Edits are reported with `edit_type`
  "sub", "ins", or "del".

- bidirect:

  logical; scan both strands? Default TRUE.

- prior:

  numeric; pseudocount for PSSM frequencies. Default 0.01.

- score.min:

  numeric or NULL; skip windows with PWM score below this. Default NULL
  (no filter).

- score.max:

  numeric or NULL; skip windows with PWM score above this. Default NULL
  (no filter).

- extend:

  logical or integer; extend sequence for boundary motifs. Default TRUE.

- strand:

  integer; which strand to scan when `bidirect=FALSE`. 1=forward,
  -1=reverse. Default 1.

- direction:

  character; direction of the edit distance query. `"above"` (default)
  finds minimum edits to bring score above `score.thresh`; `"below"`
  finds minimum edits to bring score below `score.thresh`.

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
