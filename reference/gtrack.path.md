# Returns the path on disk of a track

Returns the path on disk of a track.

## Usage

``` r
gtrack.path(track = NULL)
```

## Arguments

- track:

  track name or a vector of track names

## Value

A character vector containing the full paths to the tracks on disk.

## Details

This function returns the actual file system path where a track is
stored. The function works with a single track name or a vector of track
names.

## See also

[`gtrack.exists`](https://tanaylab.github.io/misha/reference/gtrack.exists.md),
[`gtrack.ls`](https://tanaylab.github.io/misha/reference/gtrack.ls.md),
[`gintervals.path`](https://tanaylab.github.io/misha/reference/gintervals.path.md)

## Examples

``` r

gdb.init_examples()
gtrack.path("dense_track")
#> [1] "/tmp/RtmpCsKjkH/trackdb/test/tracks/dense_track.track"
gtrack.path(c("dense_track", "sparse_track"))
#> [1] "/tmp/RtmpCsKjkH/trackdb/test/tracks/dense_track.track" 
#> [2] "/tmp/RtmpCsKjkH/trackdb/test/tracks/sparse_track.track"
```
