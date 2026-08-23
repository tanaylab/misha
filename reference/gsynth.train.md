# Train a stratified Markov model from genome sequences

Computes a Markov model of order `k` (default 5) optionally stratified
by bins of one or more track expressions (e.g., GC content and CG
dinucleotide frequency). This model can be used to generate synthetic
genomes that preserve the k-mer statistics of the original genome within
each stratification bin. When called with no dimension specifications,
trains a single unstratified model.

## Usage

``` r
gsynth.train(
  ...,
  mask = NULL,
  intervals = NULL,
  iterator = NULL,
  pseudocount = 1,
  min_obs = 0,
  k = 5L,
  prior = "marginal"
)
```

## Arguments

- ...:

  Zero or more dimension specifications. Each specification is a list
  containing:

  expr

  :   Track expression for this dimension (required)

  breaks

  :   Numeric vector of bin boundaries for this dimension (required)

  bin_merge

  :   Optional list of merge specifications for merging sparse bins.
      Each specification is a named list with 'from' and 'to' elements.

  If no dimensions are provided, trains an unstratified model with a
  single bin.

- mask:

  Optional intervals to exclude from training. Regions in the mask will
  not contribute to k-mer counts. Can be computed using
  [`gscreen()`](https://tanaylab.github.io/misha/reference/gscreen.md).

- intervals:

  Genomic intervals to process. If NULL, uses all chromosomes.

- iterator:

  Iterator for track evaluation, determines the resolution at which
  track values are computed.

- pseudocount:

  Total Dirichlet concentration alpha used in the smoothed posterior
  `P(a|c,b) = (N + alpha * pi_a(b)) / (sum_a N + alpha)`. Default is 1.

- min_obs:

  Minimum number of observations ((k+1)-mers) required per bin. Bins
  with fewer observations will be marked as NA (not learned) and a
  warning will be issued. Default is 0 (no minimum). During sampling, NA
  bins will fall back to uniform sampling unless merged via `bin_merge`.

- k:

  Integer Markov order (1–10). Default is 5, which models 6-mer (context
  of length 5 plus the emitted base) transition probabilities. Higher
  values capture longer-range sequence dependencies but require
  exponentially more memory (\\4^k\\ context states).

- prior:

  Per-base Dirichlet prior \\\pi_a(b)\\. One of:

  - `"marginal"` (default): per-bin empirical base composition learned
    from the trainer's own counts (post bin-merge). Bins with zero
    observations fall back to uniform with a warning.

  - `"global"`: a single base composition pooled over all bins,
    broadcast to every bin.

  - `NULL` or `"uniform"`: uniform prior (1/4 per base) – the pre-5.6.21
    fallback.

  - Length-4 numeric (optionally named A, C, G, T): user-supplied global
    \\\pi\\, broadcast.

  - `n_bins x 4` numeric matrix: user-supplied per-bin \\\pi\\.

  Together with `pseudocount`, this defines the Dirichlet posterior. To
  reproduce the pre-5.6.21 Laplace-add-one behavior, pass
  `prior = NULL, pseudocount = 4`.

## Value

A `gsynth.model` object containing:

- k:

  Markov order used for training

- num_kmers:

  Number of context states (\\4^k\\)

- n_dims:

  Number of stratification dimensions

- dim_specs:

  List of dimension specifications (expr, breaks, num_bins, bin_map)

- dim_sizes:

  Vector of bin counts per dimension

- total_bins:

  Total number of bins (product of dim_sizes)

- total_kmers:

  Total number of valid (k+1)-mers counted

- per_bin_kmers:

  Number of (k+1)-mers counted per bin

- total_masked:

  Number of positions skipped due to mask

- total_n:

  Number of positions skipped due to N bases

- model_data:

  Internal model data (counts and CDFs)

## Details

**Strand symmetry:** The training process counts both the forward strand
(k+1)-mer and its reverse complement for each position, ensuring
strand-symmetric transition probabilities. This means the reported
total_kmers is approximately double the number of genomic positions
processed.

**N bases:** Positions where the (k+1)-mer contains any N (unknown)
bases are skipped during training and counted in `total_n`. The model
only learns from valid A/C/G/T sequences.

## See also

[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md),
[`gsynth.save`](https://tanaylab.github.io/misha/reference/gsynth.save.md),
[`gsynth.load`](https://tanaylab.github.io/misha/reference/gsynth.load.md),
[`gsynth.bin_map`](https://tanaylab.github.io/misha/reference/gsynth.bin_map.md)

## Examples

``` r
gdb.init_examples()

# Create virtual tracks for stratification
gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
gvtrack.create("cg_frac", NULL, "kmer.frac", kmer = "CG")
#> Warning: kmer sequence 'CG' is palindromic, please set strand to 1 or -1 to avoid double counting
gvtrack.create("masked_frac", NULL, "masked.frac")

# Define repeat mask
repeats <- gscreen("masked_frac > 0.5",
    intervals = gintervals.all(),
    iterator = 100
)

# Train unstratified model (no stratification)
model_0d <- gsynth.train(
    mask = repeats,
    intervals = gintervals.all(),
    iterator = 200
)
#> Setting up iterator positions...
#> Training Markov model...
#> Trained unstratified Markov-5 model: 835,310 6-mers (no stratification)

# Train model with 2D stratification (GC content and CG dinucleotide)
model <- gsynth.train(
    list(
        expr = "g_frac + c_frac",
        breaks = seq(0, 1, 0.025),
        bin_merge = list(list(from = 0.7, to = c(0.675, 0.7)))
    ),
    list(
        expr = "cg_frac",
        breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.2),
        bin_merge = list(list(from = 0.04, to = c(0.03, 0.04)))
    ),
    mask = repeats,
    intervals = gintervals.all(),
    iterator = 200
)
#> Extracting track values...
#> Training Markov model...
#> Warning: 127 bin(s) had zero observations; their prior fell back to uniform 1/4.
#> Trained Markov-5 model: 835,310 6-mers across 200 bins (2 dimensions)
```
