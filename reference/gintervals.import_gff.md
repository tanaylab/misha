# Import intervals from a GFF/GTF file

Reads a GFF3 or GTF file (optionally gzipped) and returns a misha 1D
intervals data frame. GFF/GTF are 1-based and inclusive on both ends;
coordinates are converted to 0-based half-open by subtracting 1 from
`start` and leaving `end` as-is.

## Usage

``` r
gintervals.import_gff(file = NULL, feature = NULL, strand = TRUE, attrs = TRUE)
```

## Arguments

- file:

  path to a GFF/GTF file.

- feature:

  optional feature-type filter (column 3 of GFF). Pass a character
  vector to keep only those types (e.g. `"exon"`,
  `c("gene", "transcript")`).

- strand:

  if `TRUE`, include the `strand` column.

- attrs:

  if `TRUE`, include the raw attributes string as column `attrs`. The
  attribute string is not parsed.

## Value

A 1D intervals data frame with columns chrom, start, end, and optionally
strand, type, source, score, attrs.

## Details

Chromosome names are normalized through the active database's
`CHROM_ALIAS` mechanism.

## See also

[`gintervals.import_bed`](https://tanaylab.github.io/misha/reference/gintervals.import_bed.md),
[`gintervals.import_vcf`](https://tanaylab.github.io/misha/reference/gintervals.import_vcf.md).
