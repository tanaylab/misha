# Evaluate a track expression under single-base substitutions

For each iterator bin, evaluates `expr` once on the reference genome and
once for every single-base substitution the bin's sequence virtual
tracks read: all three alternative bases at each base in reach (the
default), or the variants listed in `variants`.

## Usage

``` r
gextract.variants(
  expr = NULL,
  intervals = NULL,
  variants = NULL,
  iterator = NULL,
  score.thresh = NULL,
  colnames = NULL
)
```

## Arguments

- expr:

  a track expression with at least one sequence virtual track.

- intervals:

  genomic scope, as in
  [`gextract`](https://tanaylab.github.io/misha/reference/gextract.md).

- variants:

  `NULL` for every substitution in reach, or a data frame with columns
  `chrom`, `start` (0-based), `ref` and `alt`, such as the output of
  [`gintervals.import_vcf`](https://tanaylab.github.io/misha/reference/gintervals.import_vcf.md).
  A variant whose `ref` does not match the genome gets a warning and
  `NA` in `<name>.alt`; `<name>.ref` keeps the reference value.

- iterator:

  track expression iterator, as in
  [`gextract`](https://tanaylab.github.io/misha/reference/gextract.md).

- score.thresh:

  if not `NULL`, keep only the rows where the reference value and the
  substituted value fall on different sides of it.

- colnames:

  name used for the value columns, `<name>.ref` and `<name>.alt`; the
  expression itself by default.

## Value

A data frame with one row per (bin, substitution): the bin's `chrom`,
`start` and `end`; the substitution's `pos` (0-based), `ref` and `alt`;
`<name>.ref` and `<name>.alt`; and `intervalID` as in
[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md).
`NULL` when there is no row.

## Details

A bin `[s, e)` of a sequence virtual track reads
`[s + sshift, e + eshift + W - 1)`, where the shifts come from
[`gvtrack.iterator`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.md)
and `W - 1` is added when the track extends (`extend = TRUE`, motif or
model width `W`). Two families read further: `kmer.count` and
`kmer.frac` with `strand` 0 or -1 start at `s + sshift - (k - 1)` for a
k-mer of length `k` (when `extend = TRUE`), and `pwm.edit_distance` with
`max_indels` ends `max_indels` bases later, at
`e + eshift + W - 1 + max_indels`. A substitution is paired with every
bin whose reach, the union over the expression's sequence virtual
tracks, contains it. Tracks and non-sequence virtual tracks in `expr`
are read as usual, so they contribute the same value to both columns.
The substituted base takes the case of the reference base, so `masked.*`
tracks see the same mask. Only single-base substitutions are supported.

The exact pass runs in one process. For genome-scale runs, split
`intervals` by chromosome and run the pieces in parallel.

Run genome-wide, every bin returns about three rows per base in reach.
To find the bins one substitution can move across a threshold, screen
first with a bound (for example from `potts.max` virtual tracks with
`max_edits = 1`) and pass the result as `intervals`, with
`score.thresh`.

## See also

[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md),
[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md),
[`gintervals.import_vcf`](https://tanaylab.github.io/misha/reference/gintervals.import_vcf.md)

## Examples

``` r

gdb.init_examples()
pssm <- matrix(c(
    0.7, 0.1, 0.1, 0.1,
    0.1, 0.7, 0.1, 0.1,
    0.1, 0.1, 0.7, 0.1
), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
gvtrack.create("motif", NULL, "pwm.max", params = list(pssm = pssm, bidirect = FALSE))
gextract.variants("motif", gintervals(1, 100, 102), iterator = 1)
#>    chrom start end pos ref alt motif.ref motif.alt intervalID
#> 1   chr1   100 101 100   C   A -4.874702 -3.009918          1
#> 2   chr1   100 101 100   C   G -4.874702 -4.874702          1
#> 3   chr1   100 101 100   C   T -4.874702 -4.874702          1
#> 4   chr1   100 101 101   C   A -4.874702 -6.739487          1
#> 5   chr1   100 101 101   C   G -4.874702 -6.739487          1
#> 6   chr1   100 101 101   C   T -4.874702 -6.739487          1
#> 7   chr1   100 101 102   T   A -4.874702 -4.874702          1
#> 8   chr1   100 101 102   T   C -4.874702 -4.874702          1
#> 9   chr1   100 101 102   T   G -4.874702 -3.009918          1
#> 10  chr1   101 102 101   C   A -6.739487 -4.874702          1
#> 11  chr1   101 102 101   C   G -6.739487 -6.739487          1
#> 12  chr1   101 102 101   C   T -6.739487 -6.739487          1
#> 13  chr1   101 102 102   T   A -6.739487 -6.739487          1
#> 14  chr1   101 102 102   T   C -6.739487 -4.874702          1
#> 15  chr1   101 102 102   T   G -6.739487 -6.739487          1
#> 16  chr1   101 102 103   A   C -6.739487 -6.739487          1
#> 17  chr1   101 102 103   A   G -6.739487 -4.874702          1
#> 18  chr1   101 102 103   A   T -6.739487 -6.739487          1
```
