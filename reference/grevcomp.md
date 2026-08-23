# Get reverse complement of DNA sequence

Takes a DNA sequence string and returns its reverse complement.

## Usage

``` r
grevcomp(seq)
```

## Arguments

- seq:

  A character vector containing DNA sequences (using A,C,G,T). Ignores
  other characters and NA values.

## Value

A character vector of the same length as the input, containing the
reverse complement sequences

## Examples

``` r
grevcomp("ACTG") # Returns "CAGT"
#> [1] "CAGT"
grevcomp(c("ACTG", "GGCC")) # Returns c("CAGT", "GGCC")
#> [1] "CAGT" "GGCC"
grevcomp(c("ACTG", NA, "GGCC")) # Returns c("CAGT", NA, "GGCC")
#> [1] "CAGT" NA     "GGCC"
```
