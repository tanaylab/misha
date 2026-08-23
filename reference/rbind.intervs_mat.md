# Row-bind intervs_mat objects concatenating their intervals

Concatenates the rows of the matrices AND the rows of their
`"intervals"` attributes. If any input is not an `intervs_mat`, falls
back to plain matrix `rbind` (the result is a base matrix with no
intervals attribute).

## Usage

``` r
# S3 method for class 'intervs_mat'
rbind(..., deparse.level = 1)
```

## Arguments

- ...:

  `intervs_mat` objects (and/or other matrix-like inputs).

- deparse.level:

  passed to base `rbind`.

## Value

An `intervs_mat` if all inputs were `intervs_mat`; otherwise a plain
matrix.

## See also

\[gintervals.to_mat()\]
