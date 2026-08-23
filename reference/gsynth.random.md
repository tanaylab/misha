# Generate random genome sequences

Generates random DNA sequences based on nucleotide probabilities without
using a trained Markov model. Each nucleotide is sampled independently
according to the specified probabilities.

## Usage

``` r
gsynth.random(
  intervals = NULL,
  output_path = NULL,
  output_format = c("misha", "fasta", "vector"),
  nuc_probs = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
  mask_copy = NULL,
  preserve_n = TRUE,
  seed = NULL,
  n_samples = 1,
  iterator = 1
)
```

## Arguments

- intervals:

  Genomic intervals to sample. If NULL, uses all chromosomes.

- output_path:

  Path to the output file (ignored when output_format = "vector")

- output_format:

  Output format:

  - "misha": .seq binary format (default)

  - "fasta": FASTA text format

  - "vector": Return sequences as a character vector (does not write to
    file)

- nuc_probs:

  Nucleotide probabilities. Can be specified as:

  - A named vector: `c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)`

  - An unnamed vector in A, C, G, T order: `c(0.3, 0.2, 0.2, 0.3)`

  Probabilities are automatically normalized to sum to 1. Default is
  uniform (0.25 each).

- mask_copy:

  Optional intervals to copy from the original genome instead of random
  sampling. Use this to preserve specific regions exactly as they appear
  in the reference.

- preserve_n:

  Logical; default `TRUE`. When `TRUE`, positions whose original
  reference is `N` (or `n`) are written to the output verbatim rather
  than filled with a random ACGT base. Same semantics as in
  [`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md);
  `mask_copy` intervals take precedence.

- seed:

  Random seed for reproducibility. If NULL, uses current random state.

- n_samples:

  Number of samples to generate per interval. Default is 1.

- iterator:

  Iterator for position resolution. Default is 1 (base-pair resolution).
  Larger values may speed up processing but are typically not needed for
  random sampling.

## Value

When output_format is "misha" or "fasta", returns invisible NULL and
writes the random sequences to output_path. When output_format is
"vector", returns a character vector of sequences (length = n_intervals
\* n_samples).

## Details

Unlike
[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md)
which uses a trained Markov model to generate sequences that preserve
k-mer statistics, `gsynth.random` generates purely random sequences
where each nucleotide is sampled independently. This is useful for
generating baseline random sequences or sequences with specific GC
content.

**Nucleotide ordering:** When using an unnamed vector for `nuc_probs`,
the order is A, C, G, T. Named vectors can be in any order.

## See also

[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md),
[`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md)

## Examples

``` r
gdb.init_examples()

# Generate random sequences with uniform nucleotide probabilities
seqs <- gsynth.random(
    intervals = gintervals(1, 0, 1000),
    output_format = "vector",
    seed = 42
)
#> Setting up random sampling positions...
#> Generating random sequences (1 samples per interval)...
#> Generated 1 random sequence(s)

# Generate GC-rich sequences (60% GC)
gc_rich <- gsynth.random(
    intervals = gintervals(1, 0, 1000),
    output_format = "vector",
    nuc_probs = c(A = 0.2, C = 0.3, G = 0.3, T = 0.2),
    seed = 42
)
#> Setting up random sampling positions...
#> Generating random sequences (1 samples per interval)...
#> Generated 1 random sequence(s)

# Generate AT-rich sequences
at_rich <- gsynth.random(
    intervals = gintervals(1, 0, 1000),
    output_format = "vector",
    nuc_probs = c(A = 0.35, C = 0.15, G = 0.15, T = 0.35),
    seed = 42
)
#> Setting up random sampling positions...
#> Generating random sequences (1 samples per interval)...
#> Generated 1 random sequence(s)
```
