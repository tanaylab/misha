# Read motifs from a HOMER motif format file

Parses a HOMER `.motif` format file and returns a named list of position
probability matrices (PPM). Each matrix has rows corresponding to motif
positions and columns `A`, `C`, `G`, `T`. The returned matrices are
directly usable with
[`gseq.pwm`](https://tanaylab.github.io/misha/reference/gseq.pwm.md).

## Usage

``` r
gseq.read_homer(file)
```

## Arguments

- file:

  character(1) path to a HOMER motif file (`.motif`).

## Value

A named list of numeric matrices. Each matrix has columns `A, C, G, T`
and one row per motif position. List names are derived from the
consensus sequence. Each matrix carries the following attributes:

- name:

  Motif name / description from the header

- consensus:

  Consensus sequence from the header

- log_odds_threshold:

  Detection threshold (numeric)

- log_p_value:

  Log p-value (numeric)

- w:

  Motif width (integer)

- source:

  `"homer"`

## See also

Other motif functions:
[`gseq.read_jaspar()`](https://tanaylab.github.io/misha/reference/gseq.read_jaspar.md),
[`gseq.read_meme()`](https://tanaylab.github.io/misha/reference/gseq.read_meme.md)

## Examples

``` r
if (FALSE) { # \dontrun{
motifs <- gseq.read_homer("known_motifs.motif")
names(motifs)
m <- motifs[[1]]
head(m)
attr(m, "consensus")
} # }
```
