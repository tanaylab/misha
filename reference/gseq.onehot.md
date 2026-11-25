# One-hot encoding of DNA sequences

Returns the one-hot encoding of genomic intervals of equal width, read
directly from the genome sequence files, or of character strings of
equal length.

## Usage

``` r
gseq.onehot(x = NULL, pad = FALSE, dtype = c("raw", "integer", "double"))
```

## Arguments

- x:

  intervals of equal width (a data frame or the name of an intervals
  set), or a character vector of sequences of equal length

- pad:

  for intervals: if TRUE, an interval may run past either end of its
  chromosome and the positions outside it are all 0, as the "N" of
  `gseq.extract(pad = TRUE)`. If FALSE (the default), such an interval
  is an error.

- dtype:

  storage type of the result: "raw" (one byte per value, the default),
  "integer" or "double"

## Value

An array of dimensions `c(N, 4, L)`.

## Details

The result is an array of dimensions `c(N, 4, L)`: `x[n, c, p]` is 1
when position `p` of sequence `n` is base `c`, and 0 otherwise. The
channels are A, C, G and T, in that order (the names of the second
dimension). Lowercase (soft-masked) bases are encoded like uppercase
ones. Any other character (N, IUPAC ambiguity codes, gaps) and, with
`pad = TRUE`, every position outside the chromosome has all four
channels 0.

The axes mean the same as in pymisha's `gseq_onehot`, whose numpy array
has shape `(N, 4, L)` as well: element `[n, c, p]` is the same value in
both. (The memory order differs: R arrays are column-major, numpy's are
row-major. Converting to a torch tensor in either language keeps the
indexing, i.e. `(batch, channels, length)`.)

For intervals with a 'strand' column, a row with strand -1 (or "-") is
reverse complemented, as in
[`gseq.extract`](https://tanaylab.github.io/misha/reference/gseq.extract.md):
its encoding is that of the reverse complement of the forward sequence.
The encoding is computed in C++ straight from the sequence store; no
character strings are made.

The array takes `N * 4 * L` bytes for `dtype = "raw"`, 4 times that for
"integer" and 8 times for "double".

## See also

[`gseq.extract`](https://tanaylab.github.io/misha/reference/gseq.extract.md),
[`gextract.windows`](https://tanaylab.github.io/misha/reference/gextract.windows.md)

## Examples

``` r

gdb.init_examples()
x <- gseq.onehot(gintervals(1, c(10000, 20000), c(10010, 20010)))
dim(x)
#> [1]  2  4 10
x[1, , ]
#>   [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> A   00   00   01   00   00   00   00   00   00    01
#> C   01   01   00   00   00   00   00   01   01    00
#> G   00   00   00   01   01   01   01   00   00    00
#> T   00   00   00   00   00   00   00   00   00    00

gseq.onehot(c("ACGTN", "acgtn"), dtype = "integer")[2, , ]
#>   [,1] [,2] [,3] [,4] [,5]
#> A    1    0    0    0    0
#> C    0    1    0    0    0
#> G    0    0    1    0    0
#> T    0    0    0    1    0
```
