# Import intervals from a BED file

Reads a BED file, plain or gzipped, and returns a misha 1D intervals
data frame. Track/browser/comment header lines are skipped
automatically. Chromosome names are normalized through the active
database's `CHROM_ALIAS` mechanism (so `chr1` \<-\> `1` works without
explicit configuration).

## Usage

``` r
gintervals.import_bed(file = NULL, name = TRUE, score = TRUE, strand = TRUE)
```

## Arguments

- file:

  path to a BED file (`.bed` or `.bed.gz`). Zip archives are not
  supported - unzip them first.

- name:

  if `TRUE` and a 4th column exists, include it as `name`.

- score:

  if `TRUE` and a 5th (numeric) column exists, include it as `score`.

- strand:

  if `TRUE` and a 6th column exists, include it as `strand` (mapped to
  `1`/`-1`/`0`).

## Value

A 1D intervals data frame, sorted by chrom and start.

## Details

BED is already 0-based half-open, so coordinates are taken as-is.

## See also

[`gintervals.import_gff`](https://tanaylab.github.io/misha/reference/gintervals.import_gff.md),
[`gintervals.import_vcf`](https://tanaylab.github.io/misha/reference/gintervals.import_vcf.md).
