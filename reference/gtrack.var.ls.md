# Returns a list of track variables for a track

Returns a list of track variables for a track.

## Usage

``` r
gtrack.var.ls(
  track = NULL,
  pattern = "",
  ignore.case = FALSE,
  perl = FALSE,
  fixed = FALSE,
  useBytes = FALSE
)
```

## Arguments

- track:

  track name

- pattern, ignore.case, perl, fixed, useBytes:

  see 'grep'

## Value

An array that contains the names of track variables.

## Details

This function returns a list of track variables of a track that match
the pattern (see 'grep'). If called without any arguments all track
variables of a track are returned.

## See also

[`grep`](https://rdrr.io/r/base/grep.html),
[`gtrack.var.get`](https://tanaylab.github.io/misha/reference/gtrack.var.get.md),
[`gtrack.var.set`](https://tanaylab.github.io/misha/reference/gtrack.var.set.md),
[`gtrack.var.rm`](https://tanaylab.github.io/misha/reference/gtrack.var.rm.md)

## Examples

``` r

gdb.init_examples()
gtrack.var.ls("sparse_track")
#> character(0)
gtrack.var.set("sparse_track", "test_var1", 1:10)
gtrack.var.set("sparse_track", "test_var2", "v")
gtrack.var.ls("sparse_track")
#> [1] "test_var1" "test_var2"
gtrack.var.ls("sparse_track", pattern = "2")
#> [1] "test_var2"
gtrack.var.rm("sparse_track", "test_var1")
gtrack.var.rm("sparse_track", "test_var2")
```
