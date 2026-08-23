# Complement DNA sequence

Takes a DNA sequence string and returns its complement (without
reversing).

## Usage

``` r
gseq.comp(seq)
```

## Arguments

- seq:

  A character vector containing DNA sequences (using A,C,G,T). Preserves
  case and handles NA values.

## Value

A character vector of the same length as the input, containing the
complemented sequences

## See also

[`gseq.revcomp`](https://tanaylab.github.io/misha/reference/gseq.revcomp.md),
[`gseq.rev`](https://tanaylab.github.io/misha/reference/gseq.rev.md)

## Examples

``` r
gseq.comp("ACTG") # Returns "TGAC"
#> [1] "TGAC"
gseq.comp(c("ACTG", "GGCC")) # Returns c("TGAC", "CCGG")
#> [1] "TGAC" "CCGG"
gseq.comp(c("ACTG", NA, "GGCC")) # Returns c("TGAC", NA, "CCGG")
#> [1] "TGAC" NA     "CCGG"
```
