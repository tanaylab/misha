# Forbid a k-mer pattern in a trained gsynth model

Returns a new `gsynth.model` whose samples are guaranteed not to contain
`pattern` as a substring (subject to the seeding caveat below).
Analytically equivalent to rejection sampling the output, implemented by
zeroing every transition that would produce the pattern and
renormalizing per state-row.

## Usage

``` r
gsynth.forbid_kmer(model, pattern, check = TRUE)
```

## Arguments

- model:

  A `gsynth.model` from
  [`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md).

- pattern:

  Character scalar, uppercase DNA (`ACGT` only), with
  `nchar(pattern) <= model$k + 1`. Patterns longer than one transition
  cannot be forbidden locally and error.

- check:

  Logical. If `TRUE` (default), print a short summary of how many
  transitions and how many bins were affected.

## Value

A new `gsynth.model` with modified `model_data$counts` and
`model_data$cdf`. The original model is not mutated.

## Details

Useful for building CpG-null, motif-null, or repeat-class-null synthetic
backgrounds from a standard
[`gsynth.train()`](https://tanaylab.github.io/misha/reference/gsynth.train.md)
model without retraining.

**Seeding caveat.**
[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md)
initializes the first `k` bases of each sampling interval by uniform
random draw, so those seed bases may themselves contain `pattern`. If
the seed lands on a state k-mer that already contains `pattern` as a
substring, every possible next base would extend that occurrence and
thus be forbidden; such "trapped" states fall back to uniform sampling
(not the forbid'd CDF) until the pattern slides out of the state window.
The guarantee applies to the Markov-sampled bases downstream of the
trap-escape window, not to the first few bases of the interval. Expected
residual per interval is small but nonzero; for strict pattern-free
output, pass `mask_copy` to
[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md)
to seed from a known pattern-free reference, or scrub residuals after
sampling.

## See also

[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md),
[`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# CpG-null synthetic background: train on the genome, then forbid CG.
model <- gsynth.train(
    list(expr = "gc_vt", breaks = seq(0, 1, 0.05)),
    intervals = gintervals.all(),
    iterator = 200
)
model_no_cg <- gsynth.forbid_kmer(model, "CG")
seqs <- gsynth.sample(model_no_cg,
    output_format = "vector",
    intervals = some_regions, seed = 42
)

# Motif-null background: forbid a 4-mer TF consensus substring.
model_no_ebox <- gsynth.forbid_kmer(model, "CACG")
} # }
```
