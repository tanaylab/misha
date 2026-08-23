# Defines rules for a single value calculation of a virtual 'Array' track

Defines how a single value within an interval is achieved for a virtual
track based on 'Array' track.

## Usage

``` r
gvtrack.array.slice(vtrack = NULL, slice = NULL, func = "avg", params = NULL)
```

## Arguments

- vtrack:

  virtual track name

- slice:

  a vector of column names or column indices or 'NULL'

- func, params:

  see below

## Value

None.

## Details

A track (regular or virtual) used in a track expression is expected to
return one value for each track interval. 'Array' tracks store multiple
values per interval (one for each 'column') and hence if used in a track
expression one must define the way of how a single value should be
deduced from several ones.

By default if an 'Array' track is used in a track expressions, its
interval value would be the average of all column values that are not
NaN. 'gvtrack.array.slice' allows to select specific columns and to
specify the function applied to their values.

'slice' parameter allows to choose the columns. Columns can be indicated
by their names or their indices. If 'slice' is 'NULL' the non-NaN values
of all track columns are used.

'func' parameter determines the function applied to the columns' values.
Use the following table for a reference of all valid functions and
parameters combinations:

*func = "avg", params = NULL*  
Average of columns' values.

*func = "max", params = NULL*  
Maximum of columns' values.

*func = "min", params = NULL*  
Minimum of columns' values.

*func = "stdev", params = NULL*  
Unbiased standard deviation of columns' values.

*func = "sum", params = NULL*  
Sum of columns' values.

*func = "quantile", params = \[Percentile in the range of \[0, 1\]\]*  
Quantile of columns' values.

## See also

[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md),
[`gtrack.array.get_colnames`](https://tanaylab.github.io/misha/reference/gtrack.array.get_colnames.md),
[`gtrack.array.extract`](https://tanaylab.github.io/misha/reference/gtrack.array.extract.md)

## Examples

``` r

gdb.init_examples()
gvtrack.create("vtrack1", "array_track")
gvtrack.array.slice("vtrack1", c("col2", "col4"), "max")
gextract("vtrack1", gintervals(1, 0, 1000))
#>   chrom start end vtrack1 intervalID
#> 1  chr1     0  50       4          1
#> 2  chr1   100 150     102          1
#> 3  chr1   250 300     NaN          1
```
