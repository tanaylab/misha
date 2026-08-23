# Import intervals from a VCF file

Reads a VCF/VCF.gz file and returns a misha 1D intervals data frame with
one row per record. VCF is 1-based; `start` is set to `POS - 1` and
`end` is set to `POS - 1 + nchar(REF)`, yielding a 0-based half-open
span covering the reference allele.

## Usage

``` r
gintervals.import_vcf(file = NULL, info = TRUE)
```

## Arguments

- file:

  path to a VCF/VCF.gz file.

- info:

  if `TRUE`, include the raw INFO column as `info`. The string is not
  parsed.

## Value

A 1D intervals data frame with columns chrom, start, end, and id, ref,
alt, qual, filter, optionally info.

## Details

Chromosome names are normalized through the active database's
`CHROM_ALIAS` mechanism.

Multi-allelic records are kept as a single row; the `ALT` column
contains the original comma-separated string.

## See also

[`gintervals.import_bed`](https://tanaylab.github.io/misha/reference/gintervals.import_bed.md),
[`gintervals.import_gff`](https://tanaylab.github.io/misha/reference/gintervals.import_gff.md).
