# Returns information about a track

Returns information about a track.

## Usage

``` r
gtrack.info(track = NULL, validate = FALSE)
```

## Arguments

- track:

  track name

- validate:

  if TRUE, validates the track index file integrity (for indexed
  tracks). Default: FALSE

## Value

A list that contains track properties

## Details

Returns information about the track (type, dimensions, size in bytes,
etc.). The fields in the returned value vary depending on the type of
the track.

## See also

[`gtrack.exists`](https://tanaylab.github.io/misha/reference/gtrack.exists.md),
[`gtrack.ls`](https://tanaylab.github.io/misha/reference/gtrack.ls.md)

## Examples

``` r

gdb.init_examples()
gtrack.info("dense_track")
#> $type
#> [1] "dense"
#> 
#> $dimensions
#> [1] 1
#> 
#> $size.in.bytes
#> [1] 80012
#> 
#> $format
#> [1] "per-chromosome"
#> 
#> $bin.size
#> [1] 50
#> 
gtrack.info("rects_track")
#> $type
#> [1] "rectangles"
#> 
#> $dimensions
#> [1] 2
#> 
#> $size.in.bytes
#> [1] 104844
#> 
#> $format
#> [1] "per-chromosome"
#> 
```
