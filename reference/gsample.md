# Returns samples from the values of track expression

Returns a sample of the specified size from the values of track
expression.

## Usage

``` r
gsample(expr = NULL, n = NULL, intervals = NULL, iterator = NULL, band = NULL)
```

## Arguments

- expr:

  track expression

- n:

  a number of items to choose

- intervals:

  genomic scope for which the function is applied

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expression.

- band:

  track expression band. If 'NULL' no band is used.

## Value

An array that represent quantiles.

## Details

This function returns a sample of the specified size from the values of
track expression. If 'n' is less than the total number of values, the
data is randomly sampled. Call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before this function
to make the sample reproducible.

If 'n' is higher than the total number of values, all values are
returned (yet reshuffled).

## See also

[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md)

## Examples

``` r

gdb.init_examples()
gsample("sparse_track", 10)
#>  [1] 0.44 0.40 0.36 0.44 0.36 0.36 0.36 0.44 0.40 0.46
```
