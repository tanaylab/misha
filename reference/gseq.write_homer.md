# Write motifs to a HOMER motif format file

Writes motifs in the HOMER `.motif` format that
[`gseq.read_homer`](https://tanaylab.github.io/misha/reference/gseq.read_homer.md)
and HOMER read: a header line
`>consensus<TAB>name<TAB>log-odds threshold[<TAB>log p-value]` and one
tab-separated row of A, C, G, T probabilities per position.

## Usage

``` r
gseq.write_homer(motifs, file, threshold = NULL, threshold_frac = 0.6)
```

## Arguments

- motifs:

  a named list of motif matrices, such as returned by the motif readers,
  or a single matrix (see
  [`gseq.write_meme`](https://tanaylab.github.io/misha/reference/gseq.write_meme.md)).
  The attributes `name`, `consensus`, `log_odds_threshold` and
  `log_p_value` are used when present.

- file:

  character(1) path of the file to write.

- threshold:

  log-odds detection threshold: a number, or one per motif. `NULL` (the
  default) uses each motif's `log_odds_threshold` attribute, and
  `threshold_frac` of its highest score for motifs without one.

- threshold_frac:

  fraction of a motif's highest log-odds score used as its threshold
  when neither `threshold` nor the motif gives one.

## Value

`file`, invisibly.

## Details

The consensus is the motif's `consensus` attribute (as
[`gseq.read_homer`](https://tanaylab.github.io/misha/reference/gseq.read_homer.md)
returns it) when it has one of the motif's width, and otherwise is
computed by HOMER's rule: each position's IUPAC symbol covers every base
whose probability is within 0.15 of the highest. The name is the motif's
`name` attribute, or its ID when it has none.
[`gseq.read_homer`](https://tanaylab.github.io/misha/reference/gseq.read_homer.md)
names motifs by their consensus, so IDs that are not used as names are
not kept.

The log-odds threshold (natural log against a uniform background, as
HOMER scores) is the `threshold` argument, else the motif's
`log_odds_threshold` attribute, else `threshold_frac` times the motif's
highest possible score. The log p-value is written when the motif has a
`log_p_value` attribute.

Probabilities, thresholds and log p-values are written with 6 decimals,
so a written motif reads back within 5e-7 of the original. Rows that do
not sum to 1 are re-normalized (with a warning). Tabs and other
whitespace in names are collapsed to single spaces.

## See also

Other motif functions:
[`gseq.read_homer()`](https://tanaylab.github.io/misha/reference/gseq.read_homer.md),
[`gseq.read_jaspar()`](https://tanaylab.github.io/misha/reference/gseq.read_jaspar.md),
[`gseq.read_meme()`](https://tanaylab.github.io/misha/reference/gseq.read_meme.md),
[`gseq.write_jaspar()`](https://tanaylab.github.io/misha/reference/gseq.write_jaspar.md),
[`gseq.write_meme()`](https://tanaylab.github.io/misha/reference/gseq.write_meme.md)

## Examples

``` r
m <- matrix(c(
    0.7, 0.1, 0.1, 0.1,
    0.1, 0.1, 0.7, 0.1
), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T")))
out <- tempfile(fileext = ".motif")
gseq.write_homer(list(M1 = m), out)
cat(readLines(out), sep = "\n")
#> >AG   M1  1.235543
#> 0.700000 0.100000    0.100000    0.100000
#> 0.100000 0.100000    0.700000    0.100000
gseq.read_homer(out)
#> $AG
#>        A   C   G   T
#> [1,] 0.7 0.1 0.1 0.1
#> [2,] 0.1 0.1 0.7 0.1
#> attr(,"name")
#> [1] "M1"
#> attr(,"consensus")
#> [1] "AG"
#> attr(,"log_odds_threshold")
#> [1] 1.235543
#> attr(,"log_p_value")
#> [1] NA
#> attr(,"w")
#> [1] 2
#> attr(,"source")
#> [1] "homer"
#> 
```
