# Score the genome under a trained gsynth model

Writes a misha fixed-bin dense track whose value at each output bin is
the summed natural-log conditional probability of the reference sequence
under the trained Markov model: \$\$T(a) = \sum\_{i=a}^{a+r-1} \log
P_M(\mathrm{seq}\[i\] \mid \mathrm{seq}\[i-k..i-1\], b(i))\$\$ where
\\b(i)\\ is the model's stratum bin (constant within each
`model$iterator`-bp window). The first \\k\\ bp of every chromosome are
NA (no upstream context available); under default policies, any output
bin containing an NA per-bp contribution is NA.

## Usage

``` r
gsynth.score(
  model,
  track,
  description = NULL,
  intervals = NULL,
  mask = NULL,
  resolution = NULL,
  sparse_policy = c("NA", "uniform"),
  n_policy = c("NA", "uniform"),
  overwrite = FALSE
)
```

## Arguments

- model:

  A `gsynth.model` object (from
  [`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md)).

- track:

  Name of the misha track to create.

- description:

  Optional track description.

- intervals:

  Intervals to score. Defaults to
  [`gintervals.all()`](https://tanaylab.github.io/misha/reference/gintervals.all.md).
  Best results when interval starts are aligned to multiples of
  `model$iterator`; otherwise the first stratum window is shorter than
  `model$iterator` and its bin label may differ from training.

- mask:

  Optional intervals to NA-out in the output (e.g. repeats). Intersects
  with `intervals` per bp; every output bin containing a masked bp
  becomes `NaN`.

- resolution:

  Output bin size in bp. Defaults to `model$iterator`. `1` produces a
  per-bp track; any positive integer is allowed.

- sparse_policy:

  How to score positions whose stratum bin is sparse in the model:
  `"NA"` (default) propagates NA; `"uniform"` contributes `log(1/4)` per
  bp.

- n_policy:

  How to score positions whose k-mer *context* contains an N: `"NA"`
  (default) or `"uniform"` (`log(1/4)` per bp). The predicted base
  itself is always NA when N — the model has no \\\log P\\ for non-ACGT
  bases.

- overwrite:

  If `TRUE`, replace an existing track of the same name.

## Value

Invisibly `NULL`. Side effect: creates the named misha track. Output
bins are written as `NaN` where the sum is undefined (out of intervals,
or any per-bp NA).

## Details

Two models scored against the same reference give a windowed log Bayes
factor as a track expression:


        gsynth.score(genome_model, "genome_score")
        gsynth.score(cre_model,    "cre_score")
        gextract("cre_score - genome_score", iterator=1000, ...)

## See also

[`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md),
[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md)
