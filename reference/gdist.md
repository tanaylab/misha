# Calculates distribution of track expressions

Calculates distribution of track expressions' values over the given set
of bins.

## Usage

``` r
gdist(
  ...,
  intervals = NULL,
  include.lowest = FALSE,
  iterator = NULL,
  band = NULL,
  dataframe = FALSE,
  names = NULL
)
```

## Arguments

- ...:

  pairs of 'expr', 'breaks' where 'expr' is a track expression and the
  breaks determine the bin. An extra, unpaired trailing argument is
  taken as the genomic scope (equivalent to 'intervals'); it cannot be
  combined with a named 'intervals' argument - that combination is
  ambiguous and raises an error.

- intervals:

  genomic scope for which the function is applied

- include.lowest:

  if 'TRUE', the lowest value of the range determined by breaks is
  included

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expressions.

- band:

  track expression band. If 'NULL' no band is used.

- dataframe:

  return a data frame instead of an N-dimensional vector.

- names:

  names for track expressions in the returned dataframe (only relevant
  when `dataframe == TRUE`)

## Value

N-dimensional vector where N is the number of 'expr'-'breaks' pairs. If
`dataframe == TRUE` - a data frame with a column for each track
expression and an additional column 'n' with counts.

## Details

This function calculates the distribution of values of the numeric track
expressions over the given set of bins.

The range of bins is determined by 'breaks' argument. For example:
'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins):
(x1, x2\], (x2, x3\], (x3, x4\].

If 'include.lowest' is 'TRUE' the the lowest value will be included in
the first interval, i.e. in \[x1, x2\]

'gdist' can work with any number of dimensions. If more than one
'expr'-'breaks' pair is passed, the result is a multidimensional vector,
and an individual value can be accessed by \[i1,i2,...,iN\] notation,
where 'i1' is the first track and 'iN' is the last track expression.

## See also

[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md)

## Examples

``` r

gdb.init_examples()

## calculate the distribution of dense_track for bins:
## (0, 0.2], (0.2, 0.5] and (0.5, 1]
gdist("dense_track", c(0, 0.2, 0.5, 1))
#>   (0,0.2] (0.2,0.5]   (0.5,1] 
#>     12463      1518        10 
#> attr(,"breaks")
#> attr(,"breaks")[[1]]
#> [1] 0.0 0.2 0.5 1.0
#> 

## calculate two-dimensional distribution:
## dense_track vs. sparse_track
gdist("dense_track", seq(0, 1, by = 0.1), "sparse_track",
    seq(0, 2, by = 0.2),
    iterator = 100
)
#>           (0,0.2] (0.2,0.4] (0.4,0.6] (0.6,0.8] (0.8,1] (1,1.2] (1.2,1.4]
#> (0,0.1]         0        58       114         8       3       0         0
#> (0.1,0.2]       0       377       828        38       0       0         0
#> (0.2,0.3]       0         5       294        68      21       0         0
#> (0.3,0.4]       0         2        13        20       4       2         1
#> (0.4,0.5]       0         0         0         0       1       0         0
#> (0.5,0.6]       0         0         0         1       0       0         0
#> (0.6,0.7]       0         0         0         0       0       0         0
#> (0.7,0.8]       0         0         0         0       0       0         0
#> (0.8,0.9]       0         0         0         0       0       0         0
#> (0.9,1]         0         0         0         0       0       0         0
#>           (1.4,1.6] (1.6,1.8] (1.8,2]
#> (0,0.1]           0         0       0
#> (0.1,0.2]         0         0       0
#> (0.2,0.3]         0         0       0
#> (0.3,0.4]         0         0       0
#> (0.4,0.5]         0         0       0
#> (0.5,0.6]         0         0       0
#> (0.6,0.7]         0         0       0
#> (0.7,0.8]         0         0       0
#> (0.8,0.9]         0         0       0
#> (0.9,1]           0         0       0
#> attr(,"breaks")
#> attr(,"breaks")[[1]]
#>  [1] 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
#> 
#> attr(,"breaks")[[2]]
#>  [1] 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0
#> 
```
