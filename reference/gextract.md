# Returns evaluated track expression

Returns the result of track expressions evaluation for each of the
iterator intervals.

## Usage

``` r
gextract(
  ...,
  intervals = NULL,
  colnames = NULL,
  iterator = NULL,
  band = NULL,
  file = NULL,
  intervals.set.out = NULL,
  intervals_join = c("id", "intervals", "none")
)
```

## Arguments

- ...:

  track expression

- intervals:

  genomic scope for which the function is applied

- colnames:

  sets the columns names in the returned value. If 'NULL' names are set
  to track expression.

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expressions.

- band:

  track expression band. If 'NULL' no band is used.

- file:

  file name where the function result is optionally outputted in
  tab-delimited format

- intervals.set.out:

  intervals set name where the function result is optionally outputted

- intervals_join:

  how the output relates to the input intervals data frame. `"id"`
  (default) appends an `intervalID` integer column carrying the 1-based
  row index of the originating input interval. `"intervals"` drops
  `intervalID` and instead attaches every column of the input intervals
  data frame (coords + metadata) to each output row, suffixing names
  that collide with existing output columns with `"1"`. `"none"` drops
  `intervalID` and attaches nothing. The `"intervals"` mode is only
  supported when the result is returned in memory; combining it with
  `file` or `intervals.set.out` raises an error.

## Value

If 'file' and 'intervals.set.out' are 'NULL' a set of intervals with an
additional column for each of the track expressions and 'columnID'
column.

## Details

This function returns the result of track expressions evaluation for
each of the iterator intervals. The returned value is a set of intervals
with an additional column for each of the track expressions. This value
can be used as an input for any other function that accepts intervals.
If the intervals inside 'intervals' argument overlap gextract returns
the overlapped coordinate more than once.

The order inside the result might not be the same as the order of
intervals. An additional column 'intervalID' is added to the return
value. Use this column to refer to the index of the original interval
from the supplied 'intervals'.

If 'file' parameter is not 'NULL' the result is outputted to a
tab-delimited text file (without 'intervalID' column) rather than
returned to the user. This can be especially useful when the result is
too big to fit into the physical memory. The resulted file can be used
as an input for 'gtrack.import' or 'gtrack.array.import' functions.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Similarly to 'file' parameter 'intervals.set.out' can be useful to
overcome the limits of the physical memory.

'colnames' parameter controls the names of the columns that contain the
evaluated expressions. By default the column names match the track
expressions.

## See also

[`gtrack.array.extract`](https://tanaylab.github.io/misha/reference/gtrack.array.extract.md),
[`gsample`](https://tanaylab.github.io/misha/reference/gsample.md),
[`gtrack.import`](https://tanaylab.github.io/misha/reference/gtrack.import.md),
[`gtrack.array.import`](https://tanaylab.github.io/misha/reference/gtrack.array.import.md),
[`glookup`](https://tanaylab.github.io/misha/reference/glookup.md),
[`gpartition`](https://tanaylab.github.io/misha/reference/gpartition.md),
[`gdist`](https://tanaylab.github.io/misha/reference/gdist.md)

## Examples

``` r

gdb.init_examples()

## get values of 'dense_track' for [0, 400), chrom 1
gextract("dense_track", gintervals(1, 0, 400))
#>   chrom start end dense_track intervalID
#> 1  chr1     0  50   0.1777778          1
#> 2  chr1    50 100   0.1600000          1
#> 3  chr1   100 150   0.1800000          1
#> 4  chr1   150 200   0.1600000          1
#> 5  chr1   200 250   0.1600000          1
#> 6  chr1   250 300   0.2000000          1
#> 7  chr1   300 350   0.1600000          1
#> 8  chr1   350 400   0.1600000          1

## get values of 'rects_track' (a 2D track) for a 2D interval
gextract(
    "rects_track",
    gintervals.2d("chr1", 0, 4000, "chr2", 2000, 5000)
)
#> NULL
```
