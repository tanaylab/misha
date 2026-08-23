# Returns values from 'Array' track

Returns values from 'Array' track.

## Usage

``` r
gtrack.array.extract(
  track = NULL,
  slice = NULL,
  intervals = NULL,
  file = NULL,
  intervals.set.out = NULL
)
```

## Arguments

- track:

  track name

- slice:

  a vector of column names or column indices or 'NULL'

- intervals:

  genomic scope for which the function is applied

- file:

  file name where the function result is to be saved. If 'NULL' result
  is returned to the user.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'file' and 'intervals.set.out' are 'NULL' a set of intervals with
additional columns for 'Array' track column values and 'intervalID'.

## Details

This function returns the column values of an 'Array' track in the
genomic scope specified by 'intervals'. 'slice' parameter determines
which columns should appear in the result. The columns can be indicated
by their names or their indices. If 'slice' is 'NULL' the values of all
track columns are returned.

The order inside the result might not be the same as the order of
intervals. An additional column 'intervalID' is added to the return
value. Use this column to refer to the index of the original interval
from the supplied 'intervals'.

If 'file' parameter is not 'NULL' the result is saved to a tab-delimited
text file (without 'intervalID' column) rather than returned to the
user. This can be especially useful when the result is too big to fit
into the physical memory. The resulted file can be used as an input for
'gtrack.array.import' function.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Similarly to 'file' parameter 'intervals.set.out' can be useful to
overcome the limits of the physical memory.

## See also

[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md),
[`gtrack.array.get_colnames`](https://tanaylab.github.io/misha/reference/gtrack.array.get_colnames.md),
[`gtrack.array.import`](https://tanaylab.github.io/misha/reference/gtrack.array.import.md)

## Examples

``` r

gdb.init_examples()
gtrack.array.extract(
    "array_track", c("col3", "col5"),
    gintervals(1, 0, 2000)
)
#>   chrom start  end col3 col5 intervalID
#> 1  chr1     0   50  NaN  NaN          1
#> 2  chr1   100  150  NaN  105          1
#> 3  chr1   250  300  203  NaN          1
#> 4  chr1  1400 1450  303  305          1
```
