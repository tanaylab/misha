# Read motifs from a MEME minimal motif format file

Parses a MEME minimal motif format file and returns a named list of
position probability matrices (PPM). Each matrix has rows corresponding
to motif positions and columns `A`, `C`, `G`, `T`. The returned matrices
are directly usable with
[`gseq.pwm`](https://tanaylab.github.io/misha/reference/gseq.pwm.md).

## Usage

``` r
gseq.read_meme(file)
```

## Arguments

- file:

  character(1) path to a MEME format file (`.meme`, `.txt`).

## Value

A named list of numeric matrices. Each matrix has columns `A, C, G, T`
and one row per motif position. List names are motif identifiers. Each
matrix carries the following attributes:

- name:

  Motif name / alternate ID (second token on the MOTIF line)

- alength:

  Alphabet length (integer, typically 4)

- w:

  Motif width (integer, number of positions)

- nsites:

  Number of sites used to build the matrix (numeric; `NA` if absent)

- E:

  E-value (numeric; `NA` if absent)

- url:

  URL string if present, otherwise `NA`

- strand:

  Strand specification from the file header (e.g. `"+ -"`)

- background:

  Named numeric vector of background frequencies
  (`c(A=..., C=..., G=..., T=...)`), or `NULL` if absent

## See also

Other motif functions:
[`gseq.read_homer()`](https://tanaylab.github.io/misha/reference/gseq.read_homer.md),
[`gseq.read_jaspar()`](https://tanaylab.github.io/misha/reference/gseq.read_jaspar.md)

## Examples

``` r
if (FALSE) { # \dontrun{
motifs <- gseq.read_meme("JASPAR2024_CORE_vertebrates.meme")
names(motifs)
m <- motifs[[1]]
head(m)
attr(m, "name")
attr(m, "nsites")
} # }
```
