# Write motifs to a JASPAR format file

Writes motifs as count matrices in the JASPAR format that
[`gseq.read_jaspar`](https://tanaylab.github.io/misha/reference/gseq.read_jaspar.md)
reads: a `>ID name` header line and one `A [ ... ]` row of counts per
base. `format = "simple"` writes the four rows of counts only, without
headers or labels.

## Usage

``` r
gseq.write_jaspar(
  motifs,
  file,
  nsites = NULL,
  format = c("jaspar", "simple"),
  round_counts = TRUE
)
```

## Arguments

- motifs:

  a named list of motif matrices, such as returned by the motif readers,
  or a single matrix (see
  [`gseq.write_meme`](https://tanaylab.github.io/misha/reference/gseq.write_meme.md)).
  The attributes `name` and `nsites` are used when present.

- file:

  character(1) path of the file to write.

- nsites:

  number of sites to scale probabilities to counts: a positive number,
  or one per motif. `NULL` (the default) uses each motif's `nsites`
  attribute, and 100 for motifs without one.

- format:

  `"jaspar"` (headers and labeled rows) or `"simple"` (four rows of
  counts per motif).

- round_counts:

  round counts to whole numbers (the default), or write them with up to
  6 decimals (`FALSE`).

## Value

`file`, invisibly.

## Details

JASPAR matrices hold counts, so each motif's probabilities are
multiplied by its number of sites: the `nsites` argument, else the
motif's `nsites` attribute (as
[`gseq.read_jaspar`](https://tanaylab.github.io/misha/reference/gseq.read_jaspar.md)
and
[`gseq.read_meme`](https://tanaylab.github.io/misha/reference/gseq.read_meme.md)
return it), else 100. Every position then sums to that number. Counts
are rounded to whole numbers (halves to even), as JASPAR files hold them
and as tools such as the MEME Suite's `jaspar2meme` and universalmotif
require, so probabilities read back within about 0.5 / `nsites`; the
integer counts of a JASPAR file are written back unchanged. With an
`nsites` of a few, whole positions can round to zero, which readers take
as uniform. `round_counts = FALSE` writes counts with up to 6 decimals
and without trailing zeros instead, so that probabilities read back
within 1.5e-6 / `nsites` of the original.

The header names a motif by its `name` attribute, or by its ID when it
has none (JASPAR headers always carry a name). Rows that do not sum to 1
are re-normalized (with a warning). Whitespace in IDs is replaced by
`_`; whitespace in names is collapsed to single spaces (other tools may
keep only the first word of a name). A simple-format file has no IDs:
[`gseq.read_jaspar`](https://tanaylab.github.io/misha/reference/gseq.read_jaspar.md)
names its motifs after the file.

## See also

Other motif functions:
[`gseq.read_homer()`](https://tanaylab.github.io/misha/reference/gseq.read_homer.md),
[`gseq.read_jaspar()`](https://tanaylab.github.io/misha/reference/gseq.read_jaspar.md),
[`gseq.read_meme()`](https://tanaylab.github.io/misha/reference/gseq.read_meme.md),
[`gseq.write_homer()`](https://tanaylab.github.io/misha/reference/gseq.write_homer.md),
[`gseq.write_meme()`](https://tanaylab.github.io/misha/reference/gseq.write_meme.md)

## Examples

``` r
m <- matrix(c(
    0.7, 0.1, 0.1, 0.1,
    0.1, 0.1, 0.7, 0.1
), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
out <- tempfile(fileext = ".jaspar")
gseq.write_jaspar(list(M1 = m), out, nsites = 20)
cat(readLines(out), sep = "\n")
#> >M1 M1
#> A  [ 14  2 ]
#> C  [  2  2 ]
#> G  [  2 14 ]
#> T  [  2  2 ]
gseq.read_jaspar(out)
#> $M1
#>        A   C   G   T
#> [1,] 0.7 0.1 0.1 0.1
#> [2,] 0.1 0.1 0.7 0.1
#> attr(,"name")
#> [1] "M1"
#> attr(,"w")
#> [1] 2
#> attr(,"nsites")
#> [1] 20
#> attr(,"format")
#> [1] "jaspar"
#> 
```
