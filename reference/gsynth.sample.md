# Sample a synthetic genome from a trained Markov model

Generates a synthetic genome by sampling from a trained stratified
Markov model. The generated genome preserves the k-mer statistics of the
original genome within each stratification bin.

## Usage

``` r
gsynth.sample(
  model,
  output_path = NULL,
  output_format = c("misha", "fasta", "vector"),
  mask_copy = NULL,
  preserve_n = TRUE,
  seed = NULL,
  intervals = NULL,
  n_samples = 1,
  bin_merge = NULL,
  cell_merge = NULL
)
```

## Arguments

- model:

  A gsynth.model object from
  [`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md)

- output_path:

  Path to the output file (ignored when output_format = "vector")

- output_format:

  Output format:

  - "misha": .seq binary format (default)

  - "fasta": FASTA text format

  - "vector": Return sequences as a character vector (does not write to
    file)

- mask_copy:

  Optional intervals to copy from the original genome instead of
  sampling. Use this to preserve specific regions (e.g., repeats,
  regulatory elements) exactly as they appear in the reference. Regions
  not in mask_copy will be sampled using the Markov model. Note:
  mask_copy intervals should be non-overlapping and sorted by start
  position within each chromosome. Overlapping intervals may result in
  only the first overlapping region being copied, with subsequent
  overlaps skipped due to cursor advancement during sequential
  processing.

- preserve_n:

  Logical; default `TRUE`. When `TRUE`, positions whose original
  reference is `N` (or `n`) are written to the output verbatim instead
  of being filled in with a sampled ACGT base. Case is preserved (so
  soft-masked `n` round-trips). `mask_copy` regions take precedence:
  inside a `mask_copy` interval the original byte is copied regardless.
  Set to `FALSE` to restore the pre-5.6.22 behavior of treating every
  position as something to sample.

- seed:

  Random seed for reproducibility. If NULL, uses current random state.

- intervals:

  Genomic intervals to sample. If NULL, uses all chromosomes.

- n_samples:

  Number of samples to generate per interval. Default is 1. When
  n_samples \> 1 and output_format = "fasta", headers include
  "\_sampleN". When output_format = "vector", returns n_samples \*
  n_intervals sequences.

- bin_merge:

  Optional list of bin merge specifications to apply during sampling,
  one per dimension (length must equal model\$n_dims). Each element
  should be:

  - A list of merge specifications (same format as in
    [`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md):
    each spec is `list(from = ..., to = ...)`)

  - Or NULL to use the bin mapping from training for that dimension

  This allows merging sparse bins at sampling time without re-training.
  Example for a 2D model:


      bin_merge = list(
        # Dimension 1: merge bins with values >= 0.8 to bin [0.7, 0.8)
        list(list(from = c(0.8, Inf), to = c(0.7, 0.8))),
        # Dimension 2: use training-time bin_map (no override)
        NULL
      )
             

- cell_merge:

  Optional per-*joint-cell* redirect table, applied after `bin_merge`
  and after any sparse-bin uniform fallback. Each entry redirects one
  specific training cell to another specific training cell whose
  Cartesian position cannot be expressed as a per-axis merge. Format: a
  list where each element is
  `list(from = c(v1, v2, ...), to = c(v1, v2, ...))`, with one
  representative value per dimension; or a data frame with columns
  `from_1, from_2, ..., to_1, to_2, ...`. Internally resolved via
  [`gsynth.cell_merge`](https://tanaylab.github.io/misha/reference/gsynth.cell_merge.md);
  at sampling time each source cell's CDF is replaced with the target
  cell's CDF (a pointer-level copy inside `cdf_list` — no matrix
  duplication and no change to the C++ hot path).

## Value

When output_format is "misha" or "fasta", returns invisible NULL and
writes the synthetic genome to output_path. When output_format is
"vector", returns a character vector of sequences (length = n_intervals
\* n_samples).

## Details

**FASTA index (.fai):** When `output_format = "fasta"`, the function
also writes a samtools-compatible `.fai` file alongside the FASTA (byte
offsets are tracked during the write loop, so this is effectively free).
The index is suitable for any downstream tool that expects a
samtools-indexed reference.

**N bases during sampling:** By default (`preserve_n = TRUE`) positions
whose original reference is N are copied verbatim into the output, so
gaps and centromeres remain N rather than being fabricated as ACGT.
Where the sampler still has to fill an N-adjacent k-mer context (the
first k bp inside an interval, or any k-mer window containing a
preserved N), it falls back to uniform random base selection until a
valid context is re-established. Similarly, if a bin has no learned
statistics (sparse bin with NA CDF), uniform sampling is used for that
position. Pass `preserve_n = FALSE` to recover the pre-5.6.22 behavior
of sampling every position regardless of the reference.

**Sparse bins:** If the model has sparse bins (from `min_obs` during
training), a warning is issued when sampling regions that fall into
these bins. Consider using `bin_merge` to redirect sparse bins to
well-populated ones.

## See also

[`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md),
[`gsynth.save`](https://tanaylab.github.io/misha/reference/gsynth.save.md)

## Examples

``` r
gdb.init_examples()

# Create virtual tracks
gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
gvtrack.create("cg_frac", NULL, "kmer.frac", kmer = "CG")
#> Warning: kmer sequence 'CG' is palindromic, please set strand to 1 or -1 to avoid double counting
gvtrack.create("masked_frac", NULL, "masked.frac")

# Define repeat mask (regions to preserve from original)
repeats <- gscreen("masked_frac > 0.5",
    intervals = gintervals.all(),
    iterator = 100
)

# Train model (excluding repeats from training)
model <- gsynth.train(
    list(expr = "g_frac + c_frac", breaks = seq(0, 1, 0.025)),
    list(expr = "cg_frac", breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.2)),
    mask = repeats,
    iterator = 200,
    min_obs = 1000
)
#> Extracting track values...
#> Training Markov model...
#> Warning: 108 bin(s) had zero observations; their prior fell back to uniform 1/4.
#> Warning: 124 out of 200 bins have fewer than 1000 observations and are marked as NA.
#> Use bin_merge to merge sparse bins before sampling, or they will use uniform sampling.
#> Trained Markov-5 model: 835,310 6-mers across 200 bins (2 dimensions)
#>   Note: 124 bins have < 1000 observations (marked as NA)

# Sample with mask_copy to preserve repeats from original genome
temp_dir <- tempdir()
synthetic_genome_file <- file.path(temp_dir, "synthetic_genome.fa")
gsynth.sample(model, synthetic_genome_file,
    output_format = "fasta",
    mask_copy = repeats,
    seed = 60427,
    bin_merge = list(
        list(list(from = 0.7, to = c(0.675, 0.7))),
        list(list(from = 0.04, to = c(0.03, 0.04)))
    )
)
#> Extracting track values for sampling...
#> Warning: 16 sparse bins (with < 1000 observations) will be used during sampling.
#> These bins will use uniform base distribution. Consider using bin_merge to merge them.
#> Sampling synthetic genome...
#> Synthetic genome written to: /tmp/RtmpMHZs6a/synthetic_genome.fa
```
