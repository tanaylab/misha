# Get reverse complement of DNA sequence

Alias for
[`grevcomp`](https://tanaylab.github.io/misha/reference/grevcomp.md).
Takes a DNA sequence string and returns its reverse complement.

## Usage

``` r
gseq.revcomp(seq)
```

## Arguments

- seq:

  A character vector containing DNA sequences (using A,C,G,T). Ignores
  other characters and NA values.

## Value

A character vector of the same length as the input, containing the
reverse complement sequences

## See also

[`grevcomp`](https://tanaylab.github.io/misha/reference/grevcomp.md),
[`gseq.rev`](https://tanaylab.github.io/misha/reference/gseq.rev.md),
[`gseq.comp`](https://tanaylab.github.io/misha/reference/gseq.comp.md)

## Examples

``` r
gseq.revcomp("ACTG") # Returns "CAGT"
#> [1] "CAGT"
gseq.revcomp(c("ACTG", "GGCC")) # Returns c("CAGT", "GGCC")
#> [1] "CAGT" "GGCC"
```
