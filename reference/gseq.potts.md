# Score sequences under a Potts (pairwise energy) model

Scores each sequence under a pairwise energy model - the \`potts\`
family's sequence-level entry point, the way
[`gseq.pwm`](https://tanaylab.github.io/misha/reference/gseq.pwm.md) is
the \`pwm\` family's. The score of one width-`W` window is
`intercept + sum_i e[i, x_i] + sum_k J_k[x_i, x_j]`, and a sequence
longer than `W` is reduced over every window by `mode`.

## Usage

``` r
gseq.potts(
  seqs,
  model,
  mode = c("lse", "max", "pos", "count"),
  bidirect = TRUE,
  strand = 0L,
  score.thresh = NULL,
  extend = FALSE
)
```

## Arguments

- seqs:

  character vector of sequences (case-insensitive).

- model:

  the energy model: a list with `e` (a `W x 4` numeric matrix, columns
  `A`, `C`, `G`, `T`), `J` (a list of `4x4` matrices, or an `npair x 16`
  matrix), `pairs` (an `npair x 2` matrix of 1-based position pairs,
  lower position first) and `intercept`. Extra elements named `width`,
  `pair_strength`, `attr` or `link` are accepted and ignored, so an
  object with exactly these fields can be passed verbatim.

- mode:

  `"lse"` for the log-sum-exp over windows, `"max"` for the best window,
  `"pos"` for the 1-based position of the best window (signed by strand
  when `bidirect = TRUE`), `"count"` for the number of windows scoring
  at least `score.thresh`.

- bidirect:

  if `TRUE` (default) both strands are read. The two strands are
  combined at each window by log-sum-exp for `"lse"`, `"max"` and
  `"count"`, and by the maximum for `"pos"`, which has to name a strand.
  [`gseq.pwm`](https://tanaylab.github.io/misha/reference/gseq.pwm.md)
  splits the two the same way.

- strand:

  used only when `bidirect = FALSE`, and required there: `1` for the
  forward strand, `-1` for the reverse. The default `0` is accepted only
  under `bidirect = TRUE`, where both strands are read and it is
  ignored.

- score.thresh:

  required for `mode = "count"` and ignored otherwise. A Potts score is
  an energy whose usable range depends entirely on the model, so there
  is no default; read one off `mode = "max"` over your own sequences.

- extend:

  accepted for signature parity with
  [`gseq.pwm`](https://tanaylab.github.io/misha/reference/gseq.pwm.md)
  and ignored - a bare sequence has no genome to extend into.

## Value

a numeric vector, one value per sequence.

## Details

The parameterisation is not unique: adding a constant to a position's
four `e` entries, or a row or column shift to a `J` block, changes no
score, so `e` is comparable across positions only in the zero-sum gauge,
where `rowSums(e)` is 0 and every `J` block has zero row and column
sums.

A Potts carries energies, not probabilities, so there is no `prior` and
no fallback for an ambiguous base: a window containing any non-ACGT base
is not scored and is left out of the reduction. A sequence that has
windows but none of them scorable returns `NA`, and `0` for
`mode = "count"`: it counted, and found none. A sequence shorter than
the model, or `NA`, returns `NA` for every mode including `"count"` -
there was no window to count.

Under `bidirect = TRUE`, `mode = "max"` combines the two strands at each
window by log-sum-exp, matching the `potts` and `potts.max` virtual
tracks - and
[`gseq.pwm`](https://tanaylab.github.io/misha/reference/gseq.pwm.md)
with `mode = "max"`, which combines them the same way.

## See also

[`gseq.pwm`](https://tanaylab.github.io/misha/reference/gseq.pwm.md),
[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)
for the `potts` virtual track functions.

## Examples

``` r
# A W=2 model with a real coupling. On the per-position energies alone the
# model prefers "AC" (1 + 1), but the pair term penalises A-then-C and
# rewards G-then-T, so "GT" wins instead. Zero J out and the score collapses
# to intercept + sum_i e[i, x_i], which is a PWM and not a Potts.
e <- matrix(c(
    1, 0, 0, 0,
    0, 1, 0, 0
), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
J <- list(matrix(c(
    0, -3, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 3,
    0, 0, 0, 0
), 4, 4, byrow = TRUE))
pairs <- matrix(c(1L, 2L), ncol = 2)
model <- list(e = e, J = J, pairs = pairs, intercept = 0)

# -1 and 3: the coupling picks the winner, not the per-position energies
gseq.potts(c("AC", "GT"), model, mode = "max", bidirect = FALSE, strand = 1)
#> [1] -1  3

# the reverse strand scores the reverse complement, so the two swap
gseq.potts(c("AC", "GT"), model, mode = "max", bidirect = FALSE, strand = -1)
#> [1]  3 -1

# bidirect unions the strands by LOG-SUM-EXP, not by maximum. "CA" scores 0
# on both strands, so it comes back as log(2) = 0.693 rather than 0 - and a
# window with an N is left out of the reduction entirely.
gseq.potts(c("AC", "GT", "CA", "AN"), model, mode = "max", bidirect = TRUE)
#> [1] 3.0181499 3.0181499 0.6931472        NA

gseq.potts("ACGTACAC", model,
    mode = "count", score.thresh = 1.5, bidirect = FALSE, strand = 1
)
#> [1] 1
```
