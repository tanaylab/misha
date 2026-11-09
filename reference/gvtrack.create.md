# Creates a new virtual track

Creates a new virtual track.

## Usage

``` r
gvtrack.create(
  vtrack = NULL,
  src = NULL,
  func = NULL,
  params = NULL,
  dim = NULL,
  sshift = NULL,
  eshift = NULL,
  filter = NULL,
  ...
)
```

## Arguments

- vtrack:

  virtual track name

- src:

  source (track/intervals). NULL for PWM functions. For value-based
  tracks, provide a data frame with columns `chrom`, `start`, `end`, and
  one numeric value column. The data frame functions as an in-memory
  sparse track and supports all track-based summarizer functions.
  Intervals must not overlap.

- func:

  function name (see above)

- params:

  function parameters (see above)

- dim:

  use 'NULL' or '0' for 1D iterators. '1' converts 2D iterator to
  (chrom1, start1, end1) , '2' converts 2D iterator to (chrom2, start2,
  end2)

- sshift:

  shift of 'start' coordinate

- eshift:

  shift of 'end' coordinate

- filter:

  genomic mask to apply. Can be:

  - A data.frame with columns 'chrom', 'start', 'end' (intervals to
    mask)

  - A character string naming an intervals set

  - A character string naming a track (must be intervals-type track)

  - A list of any combination of the above (all will be unified)

  - NULL to clear the filter

- ...:

  additional named parameters for functions that accept them this way
  instead of `params` - currently the PWM, k-mer and edit-distance
  families, and `neighbor.count`. Use one style or the other: supplying
  `params` together with named `...` arguments in the same call is
  ambiguous (the named arguments would be ignored) and raises an error.
  Parameter names are validated against the accepted names for `func`,
  whether they arrive through `params` or through `...`; an unrecognised
  name, or any name passed via `...` for a `func` that takes no
  additional named parameters, raises an error naming the offending
  argument and the accepted ones. `masked.count` / `masked.frac` are the
  exception: they warn and ignore the extra arguments.

## Value

None.

## Details

This function creates a new virtual track named 'vtrack' with the given
source, function and parameters. 'src' can be either a track, intervals
(1D or 2D), or a data frame with intervals and a numeric value column
(value-based track). The tables below summarize the supported
combinations.

**Value-based tracks** Value-based tracks are data frames containing
genomic intervals with associated numeric values. They function as
in-memory sparse tracks without requiring track creation in the
database. To create a value-based track, provide a data frame with
columns `chrom`, `start`, `end`, and one numeric value column (any name
is acceptable). Value-based tracks support all track-based summarizer
functions (e.g., `avg`, `min`, `max`, `sum`, `weighted.sum`, `lse`,
`stddev`, `quantile`, `nearest`, `exists`, `size`, `first`, `last`,
`sample`, and position functions), but do not support overlapping
intervals. They behave like sparse tracks in aggregation: values are
aggregated using count-based averaging (each interval contributes
equally regardless of length), not coverage-based averaging;
`weighted.sum` weighs each interval by its overlap with the iterator
interval.

**Track-based summarizers**

|  |  |  |  |
|----|----|----|----|
| Source | func | params | Description |
| Track | avg | NULL | Average track value in the iterator interval. |
| Track (1D) | exists | vals (optional) | Returns 1 if any value exists (or specific vals if provided), 0 otherwise. |
| Track (1D) | first | NULL | First value in the iterator interval. |
| Track (1D) | last | NULL | Last value in the iterator interval. |
| Track | max | NULL | Maximum track value in the iterator interval. |
| Track | min | NULL | Minimum track value in the iterator interval. |
| Dense / Sparse / Array track | nearest | NULL | Average value inside the iterator; for sparse tracks with no samples in the interval, falls back to the closest sample outside the interval (by genomic distance). |
| Track (1D) | sample | NULL | Uniformly sampled source value from the iterator interval. |
| Track (1D) | size | NULL | Number of non-NaN values in the iterator interval. |
| Dense / Sparse / Array track | stddev | NULL | Unbiased standard deviation of values in the iterator interval. |
| Dense / Sparse / Array track | sum | NULL | Sum of the values of the bins (dense) or intervals (sparse, array) that the iterator interval touches. Each counts in full, also when the interval covers only part of it: a 312 bp window (iterator 32, `sshift = -140`, `eshift = 140`) over a 20 bp track touches 16 or 17 bins depending on where it starts, so on a track that is 1 everywhere `sum` alternates between 16 and 17. Use `weighted.sum` to weigh each bin by its overlap. |
| Dense / Sparse / Array track | weighted.sum | NULL | Sum of value times overlap: each bin or interval adds its value times the number of its bp inside the iterator interval (after `sshift`/`eshift`, clipped to the chromosome), the definition [`gtrack.create_dense`](https://tanaylab.github.io/misha/reference/gtrack.create_dense.md) uses for `func = "weighted.sum"`. The 312 bp window above gives 312 on a track that is 1 everywhere; on an interval aligned with the bins it equals `sum` times the bin size. NaN values are skipped; `NaN` when no value overlaps the interval, as for `sum`. |
| Dense / Sparse / Array track | lse | NULL | Log-sum-exp of track values in the iterator interval: log(sum(exp(x))). Numerically stable. |
| Dense / Sparse / Array track | quantile | Percentile in \[0, 1\] | Quantile of values in the iterator interval. |
| Dense track | global.percentile | NULL | Percentile of the interval average relative to the full-track distribution. |
| Dense track | global.percentile.max | NULL | Percentile of the interval maximum relative to the full-track distribution. |
| Dense track | global.percentile.min | NULL | Percentile of the interval minimum relative to the full-track distribution. |

**Track position summarizers**

|  |  |  |  |
|----|----|----|----|
| Source | func | params | Description |
| Track (1D) | first.pos.abs | NULL | Absolute genomic coordinate of the first value. |
| Track (1D) | first.pos.relative | NULL | Zero-based position (relative to interval start) of the first value. |
| Track (1D) | last.pos.abs | NULL | Absolute genomic coordinate of the last value. |
| Track (1D) | last.pos.relative | NULL | Zero-based position (relative to interval start) of the last value. |
| Track (1D) | max.pos.abs | NULL | Absolute genomic coordinate of the maximum value inside the iterator interval. |
| Track (1D) | max.pos.relative | NULL | Zero-based position (relative to interval start) of the maximum value. |
| Track (1D) | min.pos.abs | NULL | Absolute genomic coordinate of the minimum value inside the iterator interval. |
| Track (1D) | min.pos.relative | NULL | Zero-based position (relative to interval start) of the minimum value. |
| Track (1D) | sample.pos.abs | NULL | Absolute genomic coordinate of a uniformly sampled value. |
| Track (1D) | sample.pos.relative | NULL | Zero-based position (relative to interval start) of a uniformly sampled value. |

For `max.pos.relative`, `min.pos.relative`, `first.pos.relative`,
`last.pos.relative`, `sample.pos.relative`, iterator modifiers
(including `sshift` / `eshift` and 1D projections generated via
`gvtrack.iterator`) are applied before the position is reported. In
other words, the returned coordinate is always 0-based and measured from
the start of the iterator interval after all modifier adjustments.

**Interval-based summarizers**

|  |  |  |  |
|----|----|----|----|
| Source | func | params | Description |
| 1D intervals | distance | Minimal distance from center (default 0) | Signed distance using normalized formula when inside intervals, distance to edge when outside; see notes below for exact formula. |
| 1D intervals | distance.center | NULL | Distance from iterator center to the closest interval center, `NA` if outside all intervals. |
| 1D intervals | distance.edge | NULL | Edge-to-edge distance from iterator interval to closest source interval (like `gintervals.neighbors`); see notes below for strand handling. |
| 1D intervals | coverage | NULL | Fraction of iterator length covered by source intervals (after unifying overlaps). |
| 1D intervals | neighbor.count | Max distance (\>= 0) | Number of source intervals whose edge-to-edge distance from the iterator interval is within params (no unification). |

**2D track summarizers**

|  |  |  |  |
|----|----|----|----|
| Source | func | params | Description |
| 2D track | area | NULL | Area covered by intersections of track rectangles with the iterator interval. |
| 2D track | weighted.sum | NULL | Weighted sum of values where each weight equals the intersection area. |

**Motif (PWM) summarizers**

|  |  |  |  |
|----|----|----|----|
| Source | func | Key params | Description |
| NULL (sequence) | pwm | pssm, bidirect, prior, extend, spat\_\* | Log-sum-exp score of motif likelihoods across all anchors inside the iterator interval. |
| NULL (sequence) | pwm.max | pssm, bidirect, prior, extend, spat\_\* | Maximum log-likelihood score among all anchors (per-position union across strands). |
| NULL (sequence) | pwm.max.pos | pssm, bidirect, prior, extend, spat\_\* | 1-based position of the best-scoring anchor (signed by strand when `bidirect = TRUE`); coordinates are always relative to the iterator interval after any [`gvtrack.iterator()`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.md) shifts/extensions. |
| NULL (sequence) | pwm.count | pssm, score.thresh (required), bidirect, prior, extend, strand, spat\_\* | Count of anchors scoring \>= `score.thresh` (per-position union). |

**Potts (pairwise energy) summarizers**

|  |  |  |  |
|----|----|----|----|
| Source | func | Key params | Description |
| NULL (sequence) | potts | e, J, pairs, intercept, bidirect, extend, strand | Log-sum-exp of the pairwise-energy score over all anchors inside the iterator interval. |
| NULL (sequence) | potts.max | e, J, pairs, intercept, bidirect, extend, strand, max_edits | Highest-scoring anchor; with `max_edits = 1`, the highest score reachable with at most one substitution. |
| NULL (sequence) | potts.max.pos | e, J, pairs, intercept, bidirect, extend, strand | 1-based position of the best anchor, signed by strand when `bidirect = TRUE`. |
| NULL (sequence) | potts.count | e, J, pairs, intercept, score.thresh (required), bidirect, extend, strand | Number of anchors scoring at least `score.thresh`. |

**Edit distance summarizers**

|  |  |  |  |
|----|----|----|----|
| Source | func | Key params | Description |
| NULL (sequence) | pwm.edit_distance | pssm, score.thresh, max_edits, max_indels, score.min, score.max, direction, bidirect, prior, extend, strand | Minimum number of edits (substitutions + indels) needed to reach `score.thresh` across all windows. |
| NULL (sequence) | pwm.edit_distance.pos | pssm, score.thresh, max_edits, max_indels, score.min, score.max, direction, bidirect, prior, extend, strand | 1-based position of the window achieving minimum edit distance (signed by strand when bidirect). |
| NULL (sequence) | pwm.max.edit_distance | pssm, score.thresh, max_edits, max_indels, score.min, score.max, direction, bidirect, prior, extend, strand | Edit distance at the best-scoring PWM window (same location as pwm.max/pwm.max.pos). |

**Mutation count summarizers**

|  |  |  |  |
|----|----|----|----|
| Source | func | Key params | Description |
| NULL (sequence) | pwm.n_mutations | pssm, score.thresh, score.min, score.max, direction, bidirect, prior, extend, strand | For each window, the number of single-base substitutions that each on its own bring the window to `score.thresh` (0 when the window already satisfies it, none counted when no single substitution does); the interval value is the largest of these over its windows (and strands), and NA when no window has a value. An interval that already holds a match therefore returns more than 0 when another of its windows is one substitution away. |

**LSE edit distance summarizers**

|  |  |  |  |
|----|----|----|----|
| Source | func | Key params | Description |
| NULL (sequence) | pwm.edit_distance.lse | pssm, score.thresh, max_edits, score.min, score.max, direction, bidirect, prior, extend, strand | Minimum number of substitution edits needed to bring the LSE score (log-sum-exp of per-window PWM scores) to `score.thresh`. Exhaustive for k\<=2; for k\>=3 a heuristic (beam) search, whose count can exceed the true minimum and whose NA can mean that it found no solution rather than that none exists. |
| NULL (sequence) | pwm.edit_distance.lse.pos | pssm, score.thresh, max_edits, score.min, score.max, direction, bidirect, prior, extend, strand | 1-based position of an edit in the optimal edit set for raising the LSE score. Exactly "the most impactful single edit" only when one edit suffices; when the optimum needs two or more edits this is one member of that set, and which member is reported depends on scan orientation. |

**K-mer summarizers**

|  |  |  |  |
|----|----|----|----|
| Source | func | Key params | Description |
| NULL (sequence) | kmer.count | kmer, extend, strand | Number of k-mer occurrences whose anchor lies inside the iterator interval. |
| NULL (sequence) | kmer.frac | kmer, extend, strand | Fraction of possible anchors within the interval that match the k-mer. |

**Masked sequence summarizers**

|  |  |  |  |
|----|----|----|----|
| Source | func | Key params | Description |
| NULL (sequence) | masked.count | NULL | Number of masked (lowercase) base pairs in the iterator interval. |
| NULL (sequence) | masked.frac | NULL | Fraction of base pairs in the iterator interval that are masked (lowercase). |

The sections below provide additional notes for motif, potts, interval,
k-mer, and masked sequence functions.

**Motif (PWM) notes**

- `pssm`: Position-specific scoring matrix (matrix or data frame) with
  columns `A`, `C`, `G`, `T`; extra columns are ignored. Each row holds
  the probabilities or counts of the four bases at one motif position
  (rows are normalized to sum to 1); log-odds or log-probability
  matrices are rejected, as are `NA` values and rows that sum to 0.

- `bidirect`: When TRUE (default), both strands are scanned and combined
  per genomic start (per-position union). The `strand` argument is
  ignored. When FALSE, only the strand specified by `strand` is scanned.

- `prior`: Pseudocount added to frequencies (default 0.01). Set to 0 to
  disable.

- `extend`: Extends the fetched sequence so boundary-anchored motifs
  retain full context (default TRUE). The END coordinate is padded by
  motif_length - 1 for all strand modes. Anchors must still start inside
  the iterator, so with `extend = FALSE` an iterator interval narrower
  than the model holds no anchor at all and every value comes back `NA`;
  under the default `extend = TRUE` the END padding supplies the missing
  context and such an interval is scored.

- A window containing `N` or any other character that is not A, C, G or
  T (in either case) is not scored, on either strand: it is left out, as
  an anchor outside the interval would be. An interval with no scorable
  window gives `-Inf` for `pwm` and `pwm.max`, `NA` for `pwm.max.pos`
  and `0` for `pwm.count`.

- `strand`: Used only when `bidirect = FALSE`; 1 scans the forward
  strand, -1 scans the reverse strand. The `*.pos` funcs report the same
  thing for both strands: the 1-based position of the first base of the
  match in forward-strand orientation.

- `score.thresh`: Threshold for `pwm.count`, and mandatory for it -
  there is no default, and it must be a single value. Anchors with
  log-likelihood \>= `score.thresh` are counted; only one count per
  genomic start. PWM scores are log-likelihoods, so the usable range
  depends on the PSSM, the prior and any spatial weights; calibrate with
  a `pwm` or `pwm.max` virtual track over the same intervals before
  choosing one. `pwm`, `pwm.max` and `pwm.max.pos` accept `score.thresh`
  but ignore it.

- Spatial weighting (`spat_factor`, `spat_bin`, `spat_min`, `spat_max`):
  optional position-dependent weights applied in log-space. Provide a
  positive numeric vector `spat_factor`; `spat_bin` (integer \> 0)
  defines bin width; `spat_min`/`spat_max` (1-based positions counted
  from the iterator interval start) restrict the scan to anchors at or
  after `spat_min` whose window ends at or before `spat_max`, and the
  spatial bins are counted from `spat_min`.

- `pwm.max.pos`: Positions are reported 1-based relative to the iterator
  interval start (after any
  [`gvtrack.iterator()`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.md)
  shifts); `spat_min`/`spat_max` restrict which anchors are scanned but
  do not move this origin. Values are signed when `bidirect = TRUE`
  (positive for forward, negative for reverse). Ties are broken by scan
  order, which is strand-dependent: with `strand = 1` (and under
  `bidirect = TRUE`) the most 5' tied anchor wins and the forward strand
  wins a tie at the same coordinate, but `bidirect = FALSE, strand = -1`
  scans the reverse-complemented sequence and so keeps the most 3' tied
  anchor in forward coordinates.

**Potts (pairwise energy) notes**

- The score of one window is
  `intercept + sum_i e[i, x_i] + sum_k J_k[x_i, x_j]`: `e` is a `W x 4`
  matrix of per-position energies (columns `A`, `C`, `G`, `T`), and each
  row of `pairs` names a coupled pair of positions (1-based, lower
  position first) whose `4 x 4` coupling table is the matching entry of
  `J`, indexed with the lower position on the rows. The parameterisation
  is not unique: adding a constant to a position's four `e` entries, or
  a row or column shift to a `J` block, changes no score, so `e` is
  comparable across positions only in the zero-sum gauge, where
  `rowSums(e)` is 0 and every `J` block has zero row and column sums.

- **No `prior` and no `spat_*`.** A Potts carries energies, not
  probabilities, so there is no pseudocount and no spatial weighting. A
  window containing any non-ACGT base is not scored and is left out of
  the reduction.

- **"Counted nothing" and "nothing to count" are different answers.** An
  interval that has anchors but none of them scorable (for example, all
  `N`) gives `NaN` for `potts`, `potts.max` and `potts.max.pos`, and `0`
  for `potts.count`: it counted, and found none. An interval with no
  anchor at all gives `NaN` for all four, `potts.count` included - there
  was nothing to count, and a `0` would be indistinguishable from the
  line above. That covers an interval narrower than the model with
  `extend = FALSE`; a chromosome end with the default `extend = TRUE`,
  where the end padding is clipped at the contig boundary and the
  fetched sequence comes back shorter than the model (`pwm.count`
  reports `0` rather than `NA` in this one case); and a
  [`gvtrack.filter`](https://tanaylab.github.io/misha/reference/gvtrack.filter.md)
  that leaves no unmasked part, or only parts narrower than the model.

- `strand`: Used only when `bidirect = FALSE`; it is clamped to 1 when
  `bidirect = TRUE`, as in the `pwm` family.

- **The strand union, and the wart.** `potts`, `potts.max` and
  `potts.count` combine the two strands at each anchor by log-sum-exp;
  `potts.max.pos` takes the maximum, because it has to name a strand. So
  `potts.max` and `potts.max.pos` can select **different anchors**. This
  matches the `pwm` family exactly (`pwm.max` combines by log-sum-exp,
  `pwm.max.pos` does not) and is inherited on purpose rather than fixed
  on one side.

- **Reproducing a model fitted elsewhere.** Because the union is
  log-sum-exp, `bidirect = TRUE` does not reproduce a model whose own
  definition takes the **maximum** over the two strands: the two differ
  by up to `log 2`, and by exactly `log 2` wherever the strands tie. For
  those, declare two virtual tracks over the same model with
  `bidirect = FALSE` and `strand = 1` / `strand = -1`, and combine them
  yourself with `pmax`. The same applies when several models have to
  share ONE strand decision - a core plus flanking windows, say. Give
  each window a strand-locked pair, place it with
  [`gvtrack.iterator`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.md),
  and settle the strand once with an `ifelse` on the pair that decides
  it; letting each window take its own maximum computes a different and
  strictly larger quantity.

- `e`/`J`/`pairs`/`intercept` can be supplied as a whole fitted model: a
  fitted Potts model has exactly the fields this family needs plus
  `width`, `pair_strength`, `attr` and `link`; those four are accepted
  and ignored, so the model can be passed verbatim as `params`, and
  `width` is cross-checked against `nrow(e)` when present.

- `score.thresh`: Mandatory for `potts.count` - there is no default,
  since a Potts score is an energy whose usable range depends on the
  model - and ignored by the rest of the family.

- `max_edits` (`potts.max` only): 0 (default) or 1. With 1, each anchor
  scores the best of its own window and every window one substitution
  away, and the vtrack returns the best anchor. It is exact, costs about
  4 to 30 times a plain `potts.max` (more for wide models and with
  `bidirect = TRUE`), and is meant as a genome-wide bound before
  [`gextract.variants`](https://tanaylab.github.io/misha/reference/gextract.variants.md).
  A window with a non-ACGT base is still not scored.

**Edit distance notes**

The edit distance functions answer "how many single-base changes are
needed to reach a target PWM score?" Each genomic window (same length as
the PSSM) gets an edit distance: the minimum number of substitutions
(and optionally indels) required to bring its PWM log-likelihood score
to the target threshold. The virtual track returns the minimum edit
distance across all windows in the iterator interval.

**Direction.** The `direction` parameter controls which side of the
threshold to target:

- `"above"` (default): minimum edits to make the score **reach or
  exceed** `score.thresh`. Use this to ask "how close is this sequence
  to becoming a motif match?" A window already scoring above the
  threshold needs 0 edits.

- `"below"`: minimum edits to make the score **fall below**
  `score.thresh`. Use this to ask "how fragile is this motif match — how
  many mutations would disrupt it?" A window already scoring below the
  threshold needs 0 edits.

**Strand handling differs by direction.** When `bidirect = TRUE`, both
strands are scanned at each genomic position. The way their edit
distances are combined depends on the direction:

- `direction = "above"`: takes the **minimum** across both strands. A
  motif match on *either* strand is sufficient (e.g., for TF binding),
  so the easier strand wins.

- `direction = "below"`: takes the **maximum** across both strands. A
  single genomic substitution changes both strands simultaneously (the
  forward base and its reverse complement). To truly disrupt a binding
  site you must bring *both* strands below the threshold, so the harder
  strand determines the answer. For example, if the forward strand needs
  3 edits to go below the threshold and the reverse needs 1, the result
  is 3 — after just 1 edit the forward strand would still have a match.
  A strand that cannot be brought below within `max_edits` makes the
  window NA; only a strand that `score.min` / `score.max` leave out (not
  a current match) is ignored. `pwm.n_mutations` likewise counts the
  substitutions that bring both strands below the threshold.

Then, across all windows in the iterator interval, the minimum is
reported (the position where it is easiest to create or disrupt a
match).

The LSE variants (`pwm.edit_distance.lse`, `pwm.edit_distance.lse.pos`)
score the whole interval rather than one window. With `bidirect = TRUE`
their target is the strand-combined log-sum-exp, the value the `pwm`
virtual track returns, and each substitution changes both strands at
once, so the two strands are searched together in either direction.
`pwm.edit_distance.lse.pos` is then signed by the strand whose windows
the reported edit changes more.

**Parameters:**

- `score.thresh`: The target PWM log-likelihood score that the edit
  distance algorithm tries to reach; mandatory, and it must be a single
  value. For `direction = "above"`, this is the score to reach or
  exceed; for `direction = "below"`, this is the score to fall below.

- `score.min`, `score.max`: Optional numeric pre-filters on the
  **current** (unedited) PWM score. Windows scoring outside the
  \[`score.min`, `score.max`\] range are skipped (edit distance returns
  NA). Default NULL (no filter). The score is the window's ungapped
  score (the PSSM-length window from its start), also with
  `max_indels > 0`; a window whose first PSSM-length bases are not all
  A, C, G or T has no such score and is skipped when a filter is set.
  These are independent of `score.thresh` — they control *which windows
  to evaluate*, not the target score.

  **Typical usage with `direction = "below"`**: use `score.min` to
  restrict to windows that currently have a motif match. For example,
  with `score.thresh = -18` and `score.min = -14`, only windows scoring
  \>= -14 are considered, and the edit distance tells you how many edits
  to push each one below -18. Without `score.min`, windows already
  scoring below -18 will return 0 edits.

  For LSE variants, the filters apply to the aggregate LSE score, not
  individual windows.

- `max_edits`: Optional positive integer. Cap the search depth: if
  reaching the threshold requires more than `max_edits` edits, NA is
  returned for that window. Default NULL (exact computation, no cap).

- `max_indels`: Optional non-negative integer (default 0). When \> 0,
  allows insertions and deletions in addition to substitutions, using a
  banded Needleman-Wunsch DP. Each insertion or deletion counts as one
  edit. An inserted base scores the best its motif column allows (the
  column maximum for `"above"`, the minimum for `"below"`), and a
  deleted base scores nothing. A window may then hold fewer bases than
  the PSSM (down to PSSM length minus `max_indels`, the rest inserted):
  such windows are scanned where the A, C, G, T run they start in ends,
  before an `N` or other non-ACGT base, and at a chromosome end (with
  `extend = TRUE`). Typical values are 1-2. Only supported by the
  non-LSE variants (`pwm.edit_distance`, `pwm.edit_distance.pos`,
  `pwm.max.edit_distance`).

- `direction`: `"above"` (default) or `"below"`. See Direction above.

**Spatial weighting** enables position-dependent weighting for modeling
positional biases. Bins are 0-indexed from the start of the scanned
range: the iterator interval start, or position `spat_min` of it when
`spat_min` is set. When using
[`gvtrack.iterator()`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.md)
shifts (e.g., `sshift = -50`, `eshift = 50`), the interval is the
shifted one, so bins index from the expanded scan window start, not the
original interval. Both strands use the same bin at each genomic
position. Positions beyond the last bin reuse the final bin's weight. If
the window size is not divisible by `spat_bin`, the last bin is shorter
(e.g., scanning 500 bp with 40 bp bins yields bins 0-11 of 40 bp plus
bin 12 of 20 bp). Use `spat_min` and `spat_max` to restrict scanning to
a range divisible by `spat_bin` if needed.

PWM parameters can be supplied either as a single list (`params`) or via
named arguments (see examples), but not both in the same call.

**Interval distance notes**

`distance`: Given the center 'C' of the current iterator interval,
returns 'DC \* X/2' where 'DC' is the normalized distance to the center
of the interval that contains 'C', and 'X' is the value of the parameter
(default: 0). If no interval contains 'C', the result is 'D + X/2' where
'D' is the distance between 'C' and the edge of the closest interval.

`distance.center`: Given the center 'C' of the current iterator
interval, returns `NaN` if 'C' is outside of all intervals, otherwise
returns the distance between 'C' and the center of the closest interval.

`distance.edge`: Computes edge-to-edge distance from the iterator
interval to the closest source interval, using the same calculation as
`gintervals.neighbors`. Returns 0 for overlapping intervals. Distance
sign depends on the strand column of source intervals; returns unsigned
(absolute) distance if no strand column exists. Returns `NA` if no
source intervals exist on the current chromosome.

For `distance` and `distance.center`, distance can be positive or
negative depending on the position of the coordinate relative to the
interval and the strand (-1 or 1) of the interval. Distance is always
positive if `strand = 0` or if the strand column is missing. The result
is `NA` if no intervals exist for the current chromosome.

**Difference between distance functions:** The `distance` function
measures from the *center* of the iterator interval (a single coordinate
point) to the closest *edge* of source intervals when outside, or
returns a normalized distance within the interval when inside. The
`distance.center` function measures from the center of the iterator
interval to the *center* of source intervals. The `distance.edge`
function measures *edge-to-edge* distance between intervals, exactly
like `gintervals.neighbors`. Use `distance.edge` when you need the same
distance computation as `gintervals.neighbors` within a virtual track
context.

**K-mer notes**

- `kmer`: DNA sequence (case-insensitive) to count.

- `extend`: If TRUE (default), counts kmers whose anchor lies in the
  interval even if the kmer extends beyond it; when FALSE, only kmers
  fully contained in the interval are considered.

- `strand`: 1 counts matches on the forward strand, -1 counts matches on
  the reverse-complement strand, 0 counts both strands and sums them
  (default). Unlike the PWM functions - where the both-strand mode is
  set by `bidirect` and `strand` must be 1 or -1 - the k-mer functions
  accept `strand = 0` directly. For palindromic kmers (equal to their
  own reverse complement) use 1 or -1 to avoid double counting;
  `strand = 0` issues a warning in that case.

- For `kmer.frac` the denominator is the number of possible anchor
  positions in the interval: every anchor with `extend = TRUE` (minus
  the last k-1 at a chromosome end), the first L-k+1 of an interval of
  length L with `extend = FALSE`; with `strand = 0` positions on both
  strands are counted, so the fraction stays in `[0, 1]` and remains
  comparable to the single-strand fraction.

K-mer parameters can be supplied as a list or via named arguments (see
examples), but not both in the same call.

Modify iterator behavior with 'gvtrack.iterator' or
'gvtrack.iterator.2d'.

## See also

[`gvtrack.info`](https://tanaylab.github.io/misha/reference/gvtrack.info.md),
[`gvtrack.iterator`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.md),
[`gvtrack.iterator.2d`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.2d.md),
[`gvtrack.array.slice`](https://tanaylab.github.io/misha/reference/gvtrack.array.slice.md),
[`gvtrack.ls`](https://tanaylab.github.io/misha/reference/gvtrack.ls.md),
[`gvtrack.rm`](https://tanaylab.github.io/misha/reference/gvtrack.rm.md)

[`gvtrack.iterator`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.md),
[`gvtrack.iterator.2d`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.2d.md),
[`gvtrack.filter`](https://tanaylab.github.io/misha/reference/gvtrack.filter.md)

## Examples

``` r

gdb.init_examples()

gvtrack.create("vtrack1", "dense_track", "max")
gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
gextract("dense_track", "vtrack1", "vtrack2",
    gintervals(1, 0, 10000),
    iterator = 1000
)
#>    chrom start   end dense_track vtrack1 vtrack2 intervalID
#> 1   chr1     0  1000  0.08988889    0.20    0.06          1
#> 2   chr1  1000  2000  0.04100000    0.26    0.02          1
#> 3   chr1  2000  3000  0.04100000    0.10    0.04          1
#> 4   chr1  3000  4000  0.04100000    0.14    0.04          1
#> 5   chr1  4000  5000  0.04300000    0.18    0.04          1
#> 6   chr1  5000  6000  0.04000000    0.10    0.04          1
#> 7   chr1  6000  7000  0.05800000    0.22    0.04          1
#> 8   chr1  7000  8000  0.02400000    0.06    0.02          1
#> 9   chr1  8000  9000  0.03300000    0.10    0.03          1
#> 10  chr1  9000 10000  0.06900000    0.18    0.07          1

gvtrack.create("vtrack3", "dense_track", "global.percentile")
gvtrack.create("vtrack4", "annotations", "distance")
gdist(
    "vtrack3", seq(0, 1, l = 10), "vtrack4",
    seq(-500, 500, 200)
)
#>                     (-500,-300] (-300,-100] (-100,100] (100,300] (300,500]
#> (0,0.111111]                  2           5         77         7         8
#> (0.111111,0.222222]           0           0          0         0         0
#> (0.222222,0.333333]           5           7         36         9         7
#> (0.333333,0.444444]           5           2         32         3         1
#> (0.444444,0.555556]           2           2         24         2         2
#> (0.555556,0.666667]           4           1         29         6         2
#> (0.666667,0.777778]           0           2         29         2         4
#> (0.777778,0.888889]           0           1         35         0         1
#> (0.888889,1]                  2           3         18         2         3
#> attr(,"breaks")
#> attr(,"breaks")[[1]]
#>  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
#>  [8] 0.7777778 0.8888889 1.0000000
#> 
#> attr(,"breaks")[[2]]
#> [1] -500 -300 -100  100  300  500
#> 

gvtrack.create("cov", "annotations", "coverage")
gextract("cov", gintervals(1, 0, 1000), iterator = 100)
#>    chrom start  end cov intervalID
#> 1   chr1     0  100 0.8          1
#> 2   chr1   100  200 1.0          1
#> 3   chr1   200  300 1.0          1
#> 4   chr1   300  400 1.0          1
#> 5   chr1   400  500 1.0          1
#> 6   chr1   500  600 1.0          1
#> 7   chr1   600  700 1.0          1
#> 8   chr1   700  800 1.0          1
#> 9   chr1   800  900 1.0          1
#> 10  chr1   900 1000 1.0          1

pssm <- matrix(
    c(
        0.7, 0.1, 0.1, 0.1, # Example PSSM
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.7, 0.1
    ),
    ncol = 4, byrow = TRUE
)
colnames(pssm) <- c("A", "C", "G", "T")
gvtrack.create(
    "motif_score", NULL, "pwm",
    list(pssm = pssm, bidirect = TRUE, prior = 0.01)
)
gvtrack.create("max_motif_score", NULL, "pwm.max",
    pssm = pssm, bidirect = TRUE, prior = 0.01
)
gvtrack.create("max_motif_pos", NULL, "pwm.max.pos",
    pssm = pssm
)
gextract(
    c(
        "dense_track", "motif_score", "max_motif_score",
        "max_motif_pos"
    ),
    gintervals(1, 0, 10000),
    iterator = 500
)
#>    chrom start   end dense_track motif_score max_motif_score max_motif_pos
#> 1   chr1     0   500   0.1577778  -1.7766492       -5.326688          -146
#> 2   chr1   500  1000   0.0220000  -0.9113593       -4.131331           178
#> 3   chr1  1000  1500   0.0540000  -1.3073429       -4.131331           -32
#> 4   chr1  1500  2000   0.0280000  -1.5331372       -4.131331           203
#> 5   chr1  2000  2500   0.0420000  -1.4151751       -4.151339           127
#> 6   chr1  2500  3000   0.0400000  -0.9677328       -2.289690           135
#> 7   chr1  3000  3500   0.0320000  -1.3787566       -4.151339            22
#> 8   chr1  3500  4000   0.0500000  -1.2041818       -4.151339            -4
#> 9   chr1  4000  4500   0.0360000  -1.0488112       -2.289690           -39
#> 10  chr1  4500  5000   0.0500000  -0.9755265       -4.131331           -55
#> 11  chr1  5000  5500   0.0460000  -1.4299539       -4.131331           -11
#> 12  chr1  5500  6000   0.0340000  -1.2019985       -4.131331            62
#> 13  chr1  6000  6500   0.0800000  -1.6382387       -4.151339           152
#> 14  chr1  6500  7000   0.0360000  -1.1377999       -4.131331            20
#> 15  chr1  7000  7500   0.0220000  -1.0761303       -4.131331             8
#> 16  chr1  7500  8000   0.0260000  -1.0666244       -4.131331          -174
#> 17  chr1  8000  8500   0.0320000  -0.9999073       -2.289690            75
#> 18  chr1  8500  9000   0.0340000  -0.9541373       -2.289690           189
#> 19  chr1  9000  9500   0.0580000  -0.8253099       -2.289690          -450
#> 20  chr1  9500 10000   0.0800000  -1.5154095       -4.151339          -248
#>    intervalID
#> 1           1
#> 2           1
#> 3           1
#> 4           1
#> 5           1
#> 6           1
#> 7           1
#> 8           1
#> 9           1
#> 10          1
#> 11          1
#> 12          1
#> 13          1
#> 14          1
#> 15          1
#> 16          1
#> 17          1
#> 18          1
#> 19          1
#> 20          1

# Potts (pairwise energy) examples - a small hand-written W=2 model whose
# coupling does the work: the per-position energies prefer "AC", the pair
# term prefers "GT". A zero J would make this a PWM.
potts_e <- matrix(c(
    1, 0, 0, 0,
    0, 1, 0, 0
), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
potts_J <- list(matrix(c(
    0, -3, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 3,
    0, 0, 0, 0
), 4, 4, byrow = TRUE))
potts_pairs <- matrix(c(1L, 2L), ncol = 2)
potts_model <- list(e = potts_e, J = potts_J, pairs = potts_pairs, intercept = 0)

gvtrack.create("potts_max", NULL, "potts.max", params = potts_model)
gextract("potts_max", gintervals(1, 0, 1000), iterator = 200)
#>   chrom start  end potts_max intervalID
#> 1  chr1     0  200   3.01815          1
#> 2  chr1   200  400   3.01815          1
#> 3  chr1   400  600   3.01815          1
#> 4  chr1   600  800   3.01815          1
#> 5  chr1   800 1000   3.01815          1

gvtrack.create("potts_count", NULL, "potts.count",
    params = c(potts_model, list(score.thresh = 1.5))
)
gextract("potts_count", gintervals(1, 0, 1000), iterator = 200)
#>   chrom start  end potts_count intervalID
#> 1  chr1     0  200          34          1
#> 2  chr1   200  400          32          1
#> 3  chr1   400  600          41          1
#> 4  chr1   600  800          60          1
#> 5  chr1   800 1000          56          1

# Reproducing a model that MAXIMISES over the two strands rather than
# log-sum-exping them: declare each strand on its own and combine with pmax.
gvtrack.create("potts_p", NULL, "potts.max",
    params = c(potts_model, list(bidirect = FALSE, strand = 1))
)
gvtrack.create("potts_m", NULL, "potts.max",
    params = c(potts_model, list(bidirect = FALSE, strand = -1))
)
gextract("pmax(potts_p, potts_m)", gintervals(1, 0, 1000), iterator = 200)
#>   chrom start  end pmax(potts_p, potts_m) intervalID
#> 1  chr1     0  200                      3          1
#> 2  chr1   200  400                      3          1
#> 3  chr1   400  600                      3          1
#> 4  chr1   600  800                      3          1
#> 5  chr1   800 1000                      3          1

# Kmer counting examples
gvtrack.create("cg_count", NULL, "kmer.count", kmer = "CG", strand = 1)
gvtrack.create("cg_frac", NULL, "kmer.frac", kmer = "CG", strand = 1)
gextract(c("cg_count", "cg_frac"), gintervals(1, 0, 10000), iterator = 1000)
#>    chrom start   end cg_count cg_frac intervalID
#> 1   chr1     0  1000       73   0.073          1
#> 2   chr1  1000  2000       27   0.027          1
#> 3   chr1  2000  3000       23   0.023          1
#> 4   chr1  3000  4000       15   0.015          1
#> 5   chr1  4000  5000       27   0.027          1
#> 6   chr1  5000  6000       25   0.025          1
#> 7   chr1  6000  7000       18   0.018          1
#> 8   chr1  7000  8000       25   0.025          1
#> 9   chr1  8000  9000       20   0.020          1
#> 10  chr1  9000 10000       24   0.024          1

gvtrack.create("at_pos", NULL, "kmer.count", kmer = "AT", strand = 1)
gvtrack.create("at_neg", NULL, "kmer.count", kmer = "AT", strand = -1)
gvtrack.create("at_both", NULL, "kmer.count", kmer = "AT", strand = 0)
#> Warning: kmer sequence 'AT' is palindromic, please set strand to 1 or -1 to avoid double counting
gextract(c("at_pos", "at_neg", "at_both"), gintervals(1, 0, 10000), iterator = 1000)
#>    chrom start   end at_pos at_neg at_both intervalID
#> 1   chr1     0  1000      4      4       8          1
#> 2   chr1  1000  2000     47     47      94          1
#> 3   chr1  2000  3000     28     28      56          1
#> 4   chr1  3000  4000     33     33      66          1
#> 5   chr1  4000  5000     24     24      48          1
#> 6   chr1  5000  6000     19     19      38          1
#> 7   chr1  6000  7000     34     34      68          1
#> 8   chr1  7000  8000     31     31      62          1
#> 9   chr1  8000  9000     40     40      80          1
#> 10  chr1  9000 10000     32     32      64          1

# GC content
gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
gextract("g_frac + c_frac", gintervals(1, 0, 10000),
    iterator = 1000,
    colnames = "gc_content"
)
#>    chrom start   end gc_content intervalID
#> 1   chr1     0  1000      0.629          1
#> 2   chr1  1000  2000      0.522          1
#> 3   chr1  2000  3000      0.603          1
#> 4   chr1  3000  4000      0.570          1
#> 5   chr1  4000  5000      0.591          1
#> 6   chr1  5000  6000      0.622          1
#> 7   chr1  6000  7000      0.557          1
#> 8   chr1  7000  8000      0.626          1
#> 9   chr1  8000  9000      0.587          1
#> 10  chr1  9000 10000      0.565          1

# Masked base pair counting
gvtrack.create("masked_count", NULL, "masked.count")
gvtrack.create("masked_frac", NULL, "masked.frac")
gextract(c("masked_count", "masked_frac"), gintervals(1, 0, 10000), iterator = 1000)
#>    chrom start   end masked_count masked_frac intervalID
#> 1   chr1     0  1000         1000       1.000          1
#> 2   chr1  1000  2000          604       0.604          1
#> 3   chr1  2000  3000            0       0.000          1
#> 4   chr1  3000  4000            0       0.000          1
#> 5   chr1  4000  5000            0       0.000          1
#> 6   chr1  5000  6000           91       0.091          1
#> 7   chr1  6000  7000           37       0.037          1
#> 8   chr1  7000  8000            0       0.000          1
#> 9   chr1  8000  9000          143       0.143          1
#> 10  chr1  9000 10000          190       0.190          1

# Combined with GC content (unmasked regions only)
gvtrack.create("gc", NULL, "kmer.frac", kmer = "G")
gextract("gc * (1 - masked_frac)",
    gintervals(1, 0, 10000),
    iterator = 1000,
    colnames = "gc_unmasked"
)
#>    chrom start   end gc_unmasked intervalID
#> 1   chr1     0  1000   0.0000000          1
#> 2   chr1  1000  2000   0.1033560          1
#> 3   chr1  2000  3000   0.3015000          1
#> 4   chr1  3000  4000   0.2850000          1
#> 5   chr1  4000  5000   0.2955000          1
#> 6   chr1  5000  6000   0.2826990          1
#> 7   chr1  6000  7000   0.2681955          1
#> 8   chr1  7000  8000   0.3130000          1
#> 9   chr1  8000  9000   0.2515295          1
#> 10  chr1  9000 10000   0.2288250          1

# Value-based track examples
# Create a data frame with intervals and numeric values
intervals_with_values <- data.frame(
    chrom = "chr1",
    start = c(100, 300, 500),
    end = c(200, 400, 600),
    score = c(10, 20, 30)
)
# Use as value-based sparse track (functions like sparse track)
gvtrack.create("value_track", intervals_with_values, "avg")
gvtrack.create("value_track_max", intervals_with_values, "max")
gextract(c("value_track", "value_track_max"),
    gintervals(1, 0, 10000),
    iterator = 1000
)
#>    chrom start   end value_track value_track_max intervalID
#> 1   chr1     0  1000          20              30          1
#> 2   chr1  1000  2000         NaN             NaN          1
#> 3   chr1  2000  3000         NaN             NaN          1
#> 4   chr1  3000  4000         NaN             NaN          1
#> 5   chr1  4000  5000         NaN             NaN          1
#> 6   chr1  5000  6000         NaN             NaN          1
#> 7   chr1  6000  7000         NaN             NaN          1
#> 8   chr1  7000  8000         NaN             NaN          1
#> 9   chr1  8000  9000         NaN             NaN          1
#> 10  chr1  9000 10000         NaN             NaN          1

# Spatial PWM examples
# Create a PWM with higher weight in the center of intervals
pssm <- matrix(
    c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.1, 0.7
    ),
    ncol = 4, byrow = TRUE
)
colnames(pssm) <- c("A", "C", "G", "T")

# Spatial factors: low weight at edges, high in center
# For 200bp intervals with 40bp bins: bins 0, 40, 80, 120, 160
spatial_weights <- c(0.5, 1.0, 2.0, 1.0, 0.5)

gvtrack.create(
    "spatial_pwm", NULL, "pwm",
    list(
        pssm = pssm,
        bidirect = TRUE,
        spat_factor = spatial_weights,
        spat_bin = 40L
    )
)

# Compare with non-spatial PWM
gvtrack.create(
    "regular_pwm", NULL, "pwm",
    list(pssm = pssm, bidirect = TRUE)
)

gextract(c("spatial_pwm", "regular_pwm"),
    gintervals(1, 0, 10000),
    iterator = 200
)
#>    chrom start   end  spatial_pwm regular_pwm intervalID
#> 1   chr1     0   200 -0.180401430 -0.14769733          1
#> 2   chr1   200   400 -0.313298523 -0.27115873          1
#> 3   chr1   400   600 -0.115642570 -0.13521937          1
#> 4   chr1   600   800  0.140714884  0.14071478          1
#> 5   chr1   800  1000  0.483687818  0.59744269          1
#> 6   chr1  1000  1200 -0.338178247 -0.32684132          1
#> 7   chr1  1200  1400 -0.103174701 -0.01625533          1
#> 8   chr1  1400  1600  0.025532642 -0.01764396          1
#> 9   chr1  1600  1800  0.125026315  0.08757108          1
#> 10  chr1  1800  2000 -0.131757185 -0.01831027          1
#> 11  chr1  2000  2200  0.408821940  0.39407605          1
#> 12  chr1  2200  2400 -0.135600701 -0.14926806          1
#> 13  chr1  2400  2600  0.066557795  0.08941976          1
#> 14  chr1  2600  2800 -0.076511867 -0.10243648          1
#> 15  chr1  2800  3000 -0.156360179 -0.18365006          1
#> 16  chr1  3000  3200 -0.178761199 -0.09572130          1
#> 17  chr1  3200  3400 -0.045557301 -0.02312547          1
#> 18  chr1  3400  3600 -0.104294769 -0.18529195          1
#> 19  chr1  3600  3800 -0.078706115 -0.10243648          1
#> 20  chr1  3800  4000 -0.040662788  0.05026907          1
#> 21  chr1  4000  4200 -0.219249174 -0.20749611          1
#> 22  chr1  4200  4400  0.037397459  0.21399346          1
#> 23  chr1  4400  4600  0.211919934  0.19952914          1
#> 24  chr1  4600  4800 -0.177553684 -0.11072104          1
#> 25  chr1  4800  5000  0.171405658  0.04122563          1
#> 26  chr1  5000  5200 -0.102811173 -0.21086201          1
#> 27  chr1  5200  5400 -0.168331370 -0.21505040          1
#> 28  chr1  5400  5600 -0.224751011 -0.28187045          1
#> 29  chr1  5600  5800 -0.261454135 -0.22261991          1
#> 30  chr1  5800  6000  0.098980054  0.02286761          1
#> 31  chr1  6000  6200 -0.152791902 -0.16430125          1
#> 32  chr1  6200  6400  0.009540267 -0.05986627          1
#> 33  chr1  6400  6600 -0.003594513 -0.01348384          1
#> 34  chr1  6600  6800  0.154161200  0.10841974          1
#> 35  chr1  6800  7000  0.536666572  0.53904897          1
#> 36  chr1  7000  7200 -0.154778585 -0.17233904          1
#> 37  chr1  7200  7400  0.234409571  0.21766372          1
#> 38  chr1  7400  7600  0.004531178  0.02222776          1
#> 39  chr1  7600  7800 -0.358309329 -0.28187045          1
#> 40  chr1  7800  8000 -0.037823148 -0.06488629          1
#> 41  chr1  8000  8200 -0.123648830 -0.06131684          1
#> 42  chr1  8200  8400 -0.065596923 -0.05483555          1
#> 43  chr1  8400  8600  0.300766498  0.43346265          1
#> 44  chr1  8600  8800  0.233623773  0.16052382          1
#> 45  chr1  8800  9000  0.043815937  0.05156659          1
#> 46  chr1  9000  9200  0.534710526  0.34400457          1
#> 47  chr1  9200  9400 -0.122908115 -0.13060568          1
#> 48  chr1  9400  9600  0.076439612  0.26928109          1
#> 49  chr1  9600  9800 -0.048389841  0.02025047          1
#> 50  chr1  9800 10000 -0.175184920 -0.17636883          1

# Using spatial parameters with iterator shifts
gvtrack.create(
    "spatial_extended", NULL, "pwm.max",
    pssm = pssm,
    spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5),
    spat_bin = 40L
)
# Scan window will be 280bp (100bp + 2*90bp)
gvtrack.iterator("spatial_extended", sshift = -90, eshift = 90)
gextract("spatial_extended", gintervals(1, 0, 10000), iterator = 100)
#>     chrom start   end spatial_extended intervalID
#> 1    chr1     0   100       -2.0053344          1
#> 2    chr1   100   200       -2.0053344          1
#> 3    chr1   200   300       -2.6984816          1
#> 4    chr1   300   400       -3.6469755          1
#> 5    chr1   400   500       -2.0053344          1
#> 6    chr1   500   600       -1.7821908          1
#> 7    chr1   600   700       -1.7821908          1
#> 8    chr1   700   800       -1.7821908          1
#> 9    chr1   800   900       -1.5268441          1
#> 10   chr1   900  1000        0.0825938          1
#> 11   chr1  1000  1100       -0.8336970          1
#> 12   chr1  1100  1200       -2.0053344          1
#> 13   chr1  1200  1300       -1.7821908          1
#> 14   chr1  1300  1400       -2.0053344          1
#> 15   chr1  1400  1500       -2.0053344          1
#> 16   chr1  1500  1600       -1.7821908          1
#> 17   chr1  1600  1700       -2.0053344          1
#> 18   chr1  1700  1800       -1.7821908          1
#> 19   chr1  1800  1900       -1.7821908          1
#> 20   chr1  1900  2000       -2.0053344          1
#> 21   chr1  2000  2100       -0.8336970          1
#> 22   chr1  2100  2200       -0.1405498          1
#> 23   chr1  2200  2300       -1.5268441          1
#> 24   chr1  2300  2400       -1.7821908          1
#> 25   chr1  2400  2500       -1.7821908          1
#> 26   chr1  2500  2600       -1.7821908          1
#> 27   chr1  2600  2700       -1.7821908          1
#> 28   chr1  2700  2800       -2.0053344          1
#> 29   chr1  2800  2900       -2.0053344          1
#> 30   chr1  2900  3000       -2.0053344          1
#> 31   chr1  3000  3100       -2.0053344          1
#> 32   chr1  3100  3200       -1.7821908          1
#> 33   chr1  3200  3300       -1.7821908          1
#> 34   chr1  3300  3400       -2.0053344          1
#> 35   chr1  3400  3500       -1.7821908          1
#> 36   chr1  3500  3600       -1.7821908          1
#> 37   chr1  3600  3700       -2.0053344          1
#> 38   chr1  3700  3800       -1.7821908          1
#> 39   chr1  3800  3900       -1.7821908          1
#> 40   chr1  3900  4000       -1.7821908          1
#> 41   chr1  4000  4100       -1.7821908          1
#> 42   chr1  4100  4200       -0.8336970          1
#> 43   chr1  4200  4300       -0.1405498          1
#> 44   chr1  4300  4400       -1.5268441          1
#> 45   chr1  4400  4500       -1.7821908          1
#> 46   chr1  4500  4600       -1.7821908          1
#> 47   chr1  4600  4700       -1.7821908          1
#> 48   chr1  4700  4800       -1.7821908          1
#> 49   chr1  4800  4900       -1.7821908          1
#> 50   chr1  4900  5000       -1.7821908          1
#> 51   chr1  5000  5100       -2.0053344          1
#> 52   chr1  5100  5200       -1.7821908          1
#> 53   chr1  5200  5300       -2.0053344          1
#> 54   chr1  5300  5400       -2.0053344          1
#> 55   chr1  5400  5500       -2.0053344          1
#> 56   chr1  5500  5600       -2.0053344          1
#> 57   chr1  5600  5700       -1.7821908          1
#> 58   chr1  5700  5800       -2.0053344          1
#> 59   chr1  5800  5900       -2.0053344          1
#> 60   chr1  5900  6000       -1.7821908          1
#> 61   chr1  6000  6100       -2.0053344          1
#> 62   chr1  6100  6200       -1.7821908          1
#> 63   chr1  6200  6300       -1.7821908          1
#> 64   chr1  6300  6400       -1.7821908          1
#> 65   chr1  6400  6500       -1.7821908          1
#> 66   chr1  6500  6600       -1.7821908          1
#> 67   chr1  6600  6700       -1.7821908          1
#> 68   chr1  6700  6800       -1.7821908          1
#> 69   chr1  6800  6900       -1.7821908          1
#> 70   chr1  6900  7000       -1.7821908          1
#> 71   chr1  7000  7100       -1.7821908          1
#> 72   chr1  7100  7200       -2.0053344          1
#> 73   chr1  7200  7300       -1.7821908          1
#> 74   chr1  7300  7400       -1.7821908          1
#> 75   chr1  7400  7500       -1.7821908          1
#> 76   chr1  7500  7600       -1.7821908          1
#> 77   chr1  7600  7700       -2.6984816          1
#> 78   chr1  7700  7800       -2.0053344          1
#> 79   chr1  7800  7900       -1.7821908          1
#> 80   chr1  7900  8000       -1.7821908          1
#> 81   chr1  8000  8100       -1.7821908          1
#> 82   chr1  8100  8200       -2.0053344          1
#> 83   chr1  8200  8300       -2.0053344          1
#> 84   chr1  8300  8400       -1.7821908          1
#> 85   chr1  8400  8500       -1.7821908          1
#> 86   chr1  8500  8600       -0.1405498          1
#> 87   chr1  8600  8700       -0.1405498          1
#> 88   chr1  8700  8800       -2.0053344          1
#> 89   chr1  8800  8900       -1.7821908          1
#> 90   chr1  8900  9000       -1.7821908          1
#> 91   chr1  9000  9100       -0.1405498          1
#> 92   chr1  9100  9200       -0.1405498          1
#> 93   chr1  9200  9300       -2.0053344          1
#> 94   chr1  9300  9400       -0.1405498          1
#> 95   chr1  9400  9500       -0.1405498          1
#> 96   chr1  9500  9600       -1.7821908          1
#> 97   chr1  9600  9700       -1.7821908          1
#> 98   chr1  9700  9800       -1.7821908          1
#> 99   chr1  9800  9900       -2.0053344          1
#> 100  chr1  9900 10000       -2.0053344          1

# Using spat_min/spat_max to restrict scanning to a window
# For 500bp intervals, scan only positions 30-470 (440bp window)
gvtrack.create(
    "window_pwm", NULL, "pwm",
    pssm = pssm,
    bidirect = TRUE,
    spat_min = 30, # 1-based position
    spat_max = 470 # 1-based position
)
gextract("window_pwm", gintervals(1, 0, 10000), iterator = 500)
#>    chrom start   end window_pwm intervalID
#> 1   chr1     0   500  0.5735636          1
#> 2   chr1   500  1000  1.0993481          1
#> 3   chr1  1000  1500  0.6392645          1
#> 4   chr1  1500  2000  0.8320074          1
#> 5   chr1  2000  2500  0.8932184          1
#> 6   chr1  2500  3000  0.6802177          1
#> 7   chr1  3000  3500  0.6951929          1
#> 8   chr1  3500  4000  0.7342752          1
#> 9   chr1  4000  4500  0.8058116          1
#> 10  chr1  4500  5000  0.8210369          1
#> 11  chr1  5000  5500  0.5609006          1
#> 12  chr1  5500  6000  0.6062215          1
#> 13  chr1  6000  6500  0.6843127          1
#> 14  chr1  6500  7000  1.0550103          1
#> 15  chr1  7000  7500  0.7978342          1
#> 16  chr1  7500  8000  0.7205425          1
#> 17  chr1  8000  8500  0.7662502          1
#> 18  chr1  8500  9000  1.0681919          1
#> 19  chr1  9000  9500  1.0339508          1
#> 20  chr1  9500 10000  0.7142706          1

# Combining spatial weighting with window restriction
# Scan positions 50-450 with spatial weights favoring the center
gvtrack.create(
    "window_spatial_pwm", NULL, "pwm",
    pssm = pssm,
    bidirect = TRUE,
    spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5, 1.0, 0.5, 0.5),
    spat_bin = 40L,
    spat_min = 50,
    spat_max = 450
)
gextract("window_spatial_pwm", gintervals(1, 0, 10000), iterator = 500)
#>    chrom start   end window_spatial_pwm intervalID
#> 1   chr1     0   500          0.6662636          1
#> 2   chr1   500  1000          1.0076169          1
#> 3   chr1  1000  1500          0.6896258          1
#> 4   chr1  1500  2000          0.9325638          1
#> 5   chr1  2000  2500          0.9136859          1
#> 6   chr1  2500  3000          0.7101438          1
#> 7   chr1  3000  3500          0.8078886          1
#> 8   chr1  3500  4000          0.7714550          1
#> 9   chr1  4000  4500          0.9492338          1
#> 10  chr1  4500  5000          0.7719348          1
#> 11  chr1  5000  5500          0.5994145          1
#> 12  chr1  5500  6000          0.6014583          1
#> 13  chr1  6000  6500          0.7587047          1
#> 14  chr1  6500  7000          1.0583307          1
#> 15  chr1  7000  7500          0.9049057          1
#> 16  chr1  7500  8000          0.6299827          1
#> 17  chr1  8000  8500          0.7886875          1
#> 18  chr1  8500  9000          1.1660967          1
#> 19  chr1  9000  9500          0.9576747          1
#> 20  chr1  9500 10000          0.7346352          1

# --- Edit distance examples ---
# (Using the same 4-position PSSM defined above)

# First, check the PWM score distribution to pick sensible thresholds
gvtrack.create("pwm_score", NULL, "pwm.max", pssm = pssm)
head(gextract("pwm_score", gintervals(1, 500, 520), iterator = 1))
#>   chrom start end pwm_score intervalID
#> 1  chr1   500 501 -6.428051          1
#> 2  chr1   501 502 -6.428051          1
#> 3  chr1   502 503 -8.292835          1
#> 4  chr1   503 504 -8.292835          1
#> 5  chr1   504 505 -2.698482          1
#> 6  chr1   505 506 -6.428051          1
# Scores range from ~-8.3 (poor) to ~-0.8 (strong match)

# --- direction = "above": how many edits to CREATE a motif match? ---
# "How many substitutions until this window scores >= -5?"
gvtrack.create("edist_above", NULL, "pwm.edit_distance",
    pssm = pssm, score.thresh = -5
)
# At 1bp: each position gets its own edit distance (0 if already matching)
head(gextract(c("pwm_score", "edist_above"),
    gintervals(1, 500, 520),
    iterator = 1
))
#>   chrom start end pwm_score edist_above intervalID
#> 1  chr1   500 501 -6.428051           2          1
#> 2  chr1   501 502 -6.428051           2          1
#> 3  chr1   502 503 -8.292835           3          1
#> 4  chr1   503 504 -8.292835           3          1
#> 5  chr1   504 505 -2.698482           0          1
#> 6  chr1   505 506 -6.428051           2          1

# --- direction = "below": how many edits to DISRUPT a motif match? ---
#
# score.thresh is the target: "push the score below this value"
# score.min is a pre-filter: "only consider windows currently scoring above this"
# These are independent — score.min controls WHICH windows to evaluate,
# score.thresh controls WHAT score to target.
#
# Example: among windows scoring >= -5 (strong matches),
# how many edits to push below -7 (disrupt the match)?
gvtrack.create("edist_below", NULL, "pwm.edit_distance",
    pssm = pssm, score.thresh = -7, score.min = -5,
    direction = "below"
)
head(gextract(c("pwm_score", "edist_above", "edist_below"),
    gintervals(1, 500, 520),
    iterator = 1
))
#>   chrom start end pwm_score edist_above edist_below intervalID
#> 1  chr1   500 501 -6.428051           2         NaN          1
#> 2  chr1   501 502 -6.428051           2         NaN          1
#> 3  chr1   502 503 -8.292835           3         NaN          1
#> 4  chr1   503 504 -8.292835           3         NaN          1
#> 5  chr1   504 505 -2.698482           0           2          1
#> 6  chr1   505 506 -6.428051           2         NaN          1
# score >= -5 (e.g. -2.7): gets edit distance (2 edits to go below -7)
# score < -5 (e.g. -6.4): NA (filtered out by score.min)

# Without score.min, ALL windows are evaluated — windows already
# below -7 return 0 edits:
gvtrack.create("edist_below_all", NULL, "pwm.edit_distance",
    pssm = pssm, score.thresh = -7, direction = "below"
)
head(gextract(c("pwm_score", "edist_below", "edist_below_all"),
    gintervals(1, 500, 520),
    iterator = 1
))
#>   chrom start end pwm_score edist_below edist_below_all intervalID
#> 1  chr1   500 501 -6.428051         NaN               0          1
#> 2  chr1   501 502 -6.428051         NaN               0          1
#> 3  chr1   502 503 -8.292835         NaN               0          1
#> 4  chr1   503 504 -8.292835         NaN               0          1
#> 5  chr1   504 505 -2.698482           2               2          1
#> 6  chr1   505 506 -6.428051         NaN               0          1
# score >= -5: same result as edist_below (2 edits)
# score between -7 and -5 (e.g. -6.4): edist_below=NA, edist_below_all=0
# score < -7: both return 0 (already below threshold)

# --- max_edits: cap the search depth ---
# Return NA if more than 2 substitutions are needed
gvtrack.create("edist_max2", NULL, "pwm.edit_distance",
    pssm = pssm, score.thresh = -5, max_edits = 2
)
# Positions needing 3+ edits now return NA instead of 3
head(gextract(c("pwm_score", "edist_above", "edist_max2"),
    gintervals(1, 500, 520),
    iterator = 1
))
#>   chrom start end pwm_score edist_above edist_max2 intervalID
#> 1  chr1   500 501 -6.428051           2          2          1
#> 2  chr1   501 502 -6.428051           2          2          1
#> 3  chr1   502 503 -8.292835           3        NaN          1
#> 4  chr1   503 504 -8.292835           3        NaN          1
#> 5  chr1   504 505 -2.698482           0          0          1
#> 6  chr1   505 506 -6.428051           2          2          1

# --- Aggregation over intervals ---
# With a larger iterator, the minimum edit distance across the
# interval is returned.
gextract(c("pwm_score", "edist_above", "edist_below_all"),
    gintervals(1, 0, 2000),
    iterator = 200
)
#>    chrom start  end pwm_score edist_above edist_below_all intervalID
#> 1   chr1     0  200 -2.698482           0               0          1
#> 2   chr1   200  400 -4.563266           1               0          1
#> 3   chr1   400  600 -2.698482           0               0          1
#> 4   chr1   600  800 -2.698482           0               0          1
#> 5   chr1   800 1000 -0.833697           0               0          1
#> 6   chr1  1000 1200 -2.698482           0               0          1
#> 7   chr1  1200 1400 -2.698482           0               0          1
#> 8   chr1  1400 1600 -2.698482           0               0          1
#> 9   chr1  1600 1800 -2.698482           0               0          1
#> 10  chr1  1800 2000 -2.698482           0               0          1
```
