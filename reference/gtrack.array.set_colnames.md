# Sets column names of array track

Sets column names of array track.

## Usage

``` r
gtrack.array.set_colnames(track = NULL, names = NULL)
```

## Arguments

- track:

  track name

- names:

  vector of column names

## Value

None.

## Details

This sets the column names of an array track.

## See also

[`gtrack.array.get_colnames`](https://tanaylab.github.io/misha/reference/gtrack.array.get_colnames.md),
[`gtrack.array.extract`](https://tanaylab.github.io/misha/reference/gtrack.array.extract.md),
[`gvtrack.array.slice`](https://tanaylab.github.io/misha/reference/gvtrack.array.slice.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md)

## Examples

``` r

old.names <- gtrack.array.get_colnames("array_track")
new.names <- paste("modified", old.names, sep = "_")
gtrack.array.set_colnames("array_track", new.names)
gtrack.array.get_colnames("array_track")
#>  [1] "modified_col0" "modified_col1" "modified_col2" "modified_col3"
#>  [5] "modified_col4" "modified_col5" "modified_col6" "modified_col7"
#>  [9] "modified_col8" "modified_col9"
gtrack.array.set_colnames("array_track", old.names)
gtrack.array.get_colnames("array_track")
#>  [1] "col0" "col1" "col2" "col3" "col4" "col5" "col6" "col7" "col8" "col9"
```
