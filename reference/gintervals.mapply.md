# Applies a function to values of track expressions

Applies a function to values of track expressions for each interval.

## Usage

``` r
gintervals.mapply(
  FUN = NULL,
  ...,
  intervals = NULL,
  enable.gapply.intervals = FALSE,
  iterator = NULL,
  band = NULL,
  intervals.set.out = NULL,
  colnames = "value"
)
```

## Arguments

- FUN:

  function to apply, found via 'match.fun'

- ...:

  track expressions whose values are used as arguments for 'FUN'

- intervals:

  intervals for which track expressions are calculated

- enable.gapply.intervals:

  if 'TRUE', then a variable 'GAPPLY.INTERVALS' is available

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expressions.

- band:

  track expression band. If 'NULL' no band is used.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

- colnames:

  name of the column that contains the return values of 'FUN'. Default
  is "value".

## Value

If 'intervals.set.out' is 'NULL' a data frame representing intervals
with an additional column that contains the return values of 'FUN'. The
name of this additional column is specified by the 'colnames' parameter.

## Details

This function evaluates track expressions for each interval from
'intervals'. The resulted vectors are passed then as arguments to 'FUN'.

If the intervals are one-dimensional and have an additional column named
'strand' whose value is '-1', the values of the track expression are
placed to the vector in reverse order.

The current interval index (1-based) is stored in 'GAPPLY.INTERVID'
variable that is available during the execution of 'gintervals.mapply'.
There is no guarantee about the order in which the intervals are
processed. Do not rely on any specific order and use
'GITERATOR.INTERVID' variable to detect the current interval id.

If 'enable.gapply.intervals' is 'TRUE', an additional variable
'GAPPLY.INTERVALS' is defined during the execution of
'gintervals.mapply'. This variable stores the current iterator intervals
prior to track expression evaluation. Please note that setting
'enable.gapply.intervals' to 'TRUE' might severely affect the run-time
of the function.

Note: all the changes made in R environment by 'FUN' will be void if
multitasking mode is switched on. One should also refrain from
performing any other operations in 'FUN' that might be not "thread-safe"
such as updating files, etc. Please switch off multitasking
('options(gmultitasking = FALSE)') if you wish to perform such
operations.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

## See also

[`mapply`](https://rdrr.io/r/base/mapply.html)

## Examples

``` r

gdb.init_examples()
gintervals.mapply(
    max, "dense_track",
    gintervals(c(1, 2), 0, 10000)
)
#>   chrom start   end value
#> 1  chr1     0 10000  0.26
#> 2  chr2     0 10000  0.44
gintervals.mapply(
    function(x, y) {
        max(x + y)
    }, "dense_track",
    "sparse_track", gintervals(c(1, 2), 0, 10000),
    iterator = "sparse_track"
)
#>   chrom start   end value
#> 1  chr1     0 10000  0.78
#> 2  chr2     0 10000  0.96
# Using custom column name
gintervals.mapply(
    max, "dense_track",
    gintervals(c(1, 2), 0, 10000),
    colnames = "max_value"
)
#>   chrom start   end max_value
#> 1  chr1     0 10000      0.26
#> 2  chr2     0 10000      0.44
```
