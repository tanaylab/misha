# Write motifs to a MEME minimal motif format file

Writes position probability matrices to a file in the MEME minimal motif
format, which
[`gseq.read_meme`](https://tanaylab.github.io/misha/reference/gseq.read_meme.md)
and the MEME Suite read.

## Usage

``` r
gseq.write_meme(motifs, file, background = NULL, strands = NULL)
```

## Arguments

- motifs:

  a named list of motif matrices, such as returned by
  [`gseq.read_meme`](https://tanaylab.github.io/misha/reference/gseq.read_meme.md),
  [`gseq.read_jaspar`](https://tanaylab.github.io/misha/reference/gseq.read_jaspar.md)
  or
  [`gseq.read_homer`](https://tanaylab.github.io/misha/reference/gseq.read_homer.md),
  or a single matrix. Each matrix has one row per position and columns
  `A, C, G, T` (unnamed columns are taken in that order; a matrix with
  rows `A, C, G, T` is transposed). List names are the motif IDs;
  unnamed motifs are called `motif_<i>`. The attributes `name`,
  `nsites`, `E`, `url`, `strand` and `background` are written when
  present.

- file:

  character(1) path of the file to write.

- background:

  background letter frequencies: a numeric vector of 4 values named
  `A, C, G, T` (or unnamed, in that order), normalized to sum to 1.
  `NULL` (the default) uses the `background` attribute the motifs share
  (as
  [`gseq.read_meme`](https://tanaylab.github.io/misha/reference/gseq.read_meme.md)
  returns it), and uniform frequencies otherwise.

- strands:

  `"+ -"` (both strands) or `"+"`. `NULL` (the default) uses the
  `strand` attribute the motifs share, and `"+ -"` otherwise.

## Value

`file`, invisibly.

## Details

Each motif is written as a `MOTIF <id> [<name>]` line, a
`letter-probability matrix:` line with `alength`, `w` and, when the
motif has them, `nsites` and `E`, the matrix rows and, when the motif
has one, a `URL` line. The file header declares the DNA alphabet, the
strands and the background letter frequencies.

Probabilities are written with 6 decimals, as MEME writes them, so a
written motif reads back within 5e-7 of the original. `nsites` and
E-values are written with 15 significant digits.

Rows that do not sum to 1 are re-normalized (with a warning), as
[`gseq.read_meme`](https://tanaylab.github.io/misha/reference/gseq.read_meme.md)
does. Whitespace in IDs is replaced by `_`; whitespace in names is
collapsed to single spaces (other tools may keep only the first word of
a name). Motifs with the same ID are written as they are;
[`gseq.read_meme`](https://tanaylab.github.io/misha/reference/gseq.read_meme.md)
adds suffixes to tell them apart.

## See also

Other motif functions:
[`gseq.read_homer()`](https://tanaylab.github.io/misha/reference/gseq.read_homer.md),
[`gseq.read_jaspar()`](https://tanaylab.github.io/misha/reference/gseq.read_jaspar.md),
[`gseq.read_meme()`](https://tanaylab.github.io/misha/reference/gseq.read_meme.md),
[`gseq.write_homer()`](https://tanaylab.github.io/misha/reference/gseq.write_homer.md),
[`gseq.write_jaspar()`](https://tanaylab.github.io/misha/reference/gseq.write_jaspar.md)

## Examples

``` r
motifs <- list(
    M1 = matrix(c(
        0.7, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1
    ), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
)
out <- tempfile(fileext = ".meme")
gseq.write_meme(motifs, out, background = c(A = 0.3, C = 0.2, G = 0.2, T = 0.3))
cat(readLines(out), sep = "\n")
#> MEME version 4
#> 
#> ALPHABET= ACGT
#> 
#> strands: + -
#> 
#> Background letter frequencies
#> A 0.300000 C 0.200000 G 0.200000 T 0.300000
#> 
#> MOTIF M1
#> letter-probability matrix: alength= 4 w= 2
#> 0.700000 0.100000 0.100000 0.100000
#> 0.100000 0.100000 0.700000 0.100000
#> 
gseq.read_meme(out)
#> $M1
#>        A   C   G   T
#> [1,] 0.7 0.1 0.1 0.1
#> [2,] 0.1 0.1 0.7 0.1
#> attr(,"name")
#> [1] NA
#> attr(,"alength")
#> [1] 4
#> attr(,"w")
#> [1] 2
#> attr(,"nsites")
#> [1] NA
#> attr(,"E")
#> [1] NA
#> attr(,"url")
#> [1] NA
#> attr(,"strand")
#> [1] "+ -"
#> attr(,"background")
#>   A   C   G   T 
#> 0.3 0.2 0.2 0.3 
#> 
```
