# Convert an interval-indexed matrix back to an intervals + values data.frame

Inverse of
[`gintervals.to_mat`](https://tanaylab.github.io/misha/reference/gintervals.to_mat.md).
Recovers intervals from `attr(mat, "intervals")` if `mat` is an
`intervs_mat`, or from the explicit `intervals` argument if `mat` is a
plain matrix. Never parses rownames.

## Usage

``` r
gintervals.from_mat(mat, intervals = NULL)
```

## Arguments

- mat:

  an `intervs_mat` (produced by `gintervals.to_mat`) or a plain matrix
  (then `intervals` must be supplied).

- intervals:

  data.frame with `chrom`, `start`, `end` columns, required when `mat`
  is a plain matrix. Must satisfy `nrow(intervals) == nrow(mat)`;
  alignment is strictly positional. Must NOT be supplied when `mat` is
  already an `intervs_mat`.

## Value

A data.frame with `chrom`, `start`, `end` (plus `intervalID` if it was
present in the original input), followed by the value columns.

## See also

[`gintervals.to_mat`](https://tanaylab.github.io/misha/reference/gintervals.to_mat.md)

## Examples

``` r
df <- data.frame(
    chrom = c("chr1", "chr2"),
    start = c(100L, 200L),
    end   = c(200L, 400L),
    t1    = c(1.0, 2.0)
)
mat <- gintervals.to_mat(df)
identical(gintervals.from_mat(mat)$t1, df$t1)
#> [1] TRUE

# plain matrix path:
plain <- unclass(mat)
attr(plain, "intervals") <- NULL
gintervals.from_mat(plain, intervals = df[, c("chrom", "start", "end")])
#>   chrom start end t1
#> 1  chr1   100 200  1
#> 2  chr2   200 400  2
```
