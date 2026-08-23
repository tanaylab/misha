# Returns column names of array track

Returns column names of array track.

## Usage

``` r
gtrack.array.get_colnames(track = NULL)
```

## Arguments

- track:

  track name

## Value

A character vector with column names.

## Details

This function returns the column names of an array track.

## See also

[`gtrack.array.set_colnames`](https://tanaylab.github.io/misha/reference/gtrack.array.set_colnames.md),
[`gtrack.array.extract`](https://tanaylab.github.io/misha/reference/gtrack.array.extract.md),
[`gvtrack.array.slice`](https://tanaylab.github.io/misha/reference/gvtrack.array.slice.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md)

## Examples

``` r

gtrack.array.get_colnames("array_track")
#>  [1] "col0" "col1" "col2" "col3" "col4" "col5" "col6" "col7" "col8" "col9"
```
