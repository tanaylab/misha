# Read motifs from a JASPAR PFM format file

Parses a JASPAR Position Frequency Matrix (PFM) file and returns a named
list of position probability matrices (PPM). Supports both the standard
JASPAR header format (`>ID NAME` followed by labeled rows) and the
simple 4-row PFM format. Counts are converted to probabilities by
dividing each column by its column sum.

## Usage

``` r
gseq.read_jaspar(file)
```

## Arguments

- file:

  character(1) path to a JASPAR format file (`.jaspar`, `.pfm`, `.txt`).

## Value

A named list of numeric matrices. Each matrix has columns `A, C, G, T`
and one row per motif position. List names are motif identifiers. Each
matrix carries the following attributes:

- name:

  Motif name from the header line

- w:

  Motif width (integer)

- nsites:

  Total counts per position (numeric; `NA` for simple-format files)

- format:

  Sub-format detected: `"jaspar"` or `"simple"`

## See also

Other motif functions:
[`gseq.read_homer()`](https://tanaylab.github.io/misha/reference/gseq.read_homer.md),
[`gseq.read_meme()`](https://tanaylab.github.io/misha/reference/gseq.read_meme.md)

## Examples

``` r
if (FALSE) { # \dontrun{
motifs <- gseq.read_jaspar("JASPAR2024_CORE.jaspar")
names(motifs)
m <- motifs[[1]]
head(m)
} # }
```
