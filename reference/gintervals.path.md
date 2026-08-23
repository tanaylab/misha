# Returns the path on disk of an interval set

Returns the path on disk of an interval set.

## Usage

``` r
gintervals.path(intervals.set = NULL)
```

## Arguments

- intervals.set:

  name of an interval set or a vector of interval set names

## Value

A character vector containing the full paths to the interval sets on
disk.

## Details

This function returns the actual file system path where an interval set
is stored. The function works with a single interval set name or a
vector of names.

## See also

[`gintervals.exists`](https://tanaylab.github.io/misha/reference/gintervals.exists.md),
[`gintervals.ls`](https://tanaylab.github.io/misha/reference/gintervals.ls.md),
[`gtrack.path`](https://tanaylab.github.io/misha/reference/gtrack.path.md)

## Examples

``` r

gdb.init_examples()
gintervals.path("annotations")
#> [1] "/tmp/Rtmp4gKoGb/trackdb/test/tracks/annotations.interv"
gintervals.path(c("annotations", "coding"))
#> [1] "/tmp/Rtmp4gKoGb/trackdb/test/tracks/annotations.interv"
#> [2] "/tmp/Rtmp4gKoGb/trackdb/test/tracks/coding.interv"     
```
