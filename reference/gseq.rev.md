# Reverse DNA sequence

Takes a DNA sequence string and returns its reverse (without
complementing).

## Usage

``` r
gseq.rev(seq)
```

## Arguments

- seq:

  A character vector containing DNA sequences. Preserves case and
  handles NA values.

## Value

A character vector of the same length as the input, containing the
reversed sequences

## See also

[`gseq.revcomp`](https://tanaylab.github.io/misha/reference/gseq.revcomp.md),
[`gseq.comp`](https://tanaylab.github.io/misha/reference/gseq.comp.md)

## Examples

``` r
gseq.rev("ACTG") # Returns "GTCA"
#> [1] "GTCA"
gseq.rev(c("ACTG", "GGCC")) # Returns c("GTCA", "CCGG")
#> [1] "GTCA" "CCGG"
gseq.rev(c("ACTG", NA, "GGCC")) # Returns c("GTCA", NA, "CCGG")
#> [1] "GTCA" NA     "CCGG"
```
