# Returns the database paths that contain track(s)

Returns all database paths that contain a version of a track.

## Usage

``` r
gtrack.dbs(track = NULL, dataframe = FALSE)
```

## Arguments

- track:

  track name or a vector of track names

- dataframe:

  return a data frame with columns `track` and `db` instead of a named
  character vector.

## Value

A named character vector of database paths for each track. If
`dataframe` is TRUE, returns a data frame with columns `track` and `db`,
with multiple rows per track when it appears in multiple databases.

## Details

When datasets are loaded, a track may exist in multiple locations
(working database and/or datasets). This function computes on-demand and
returns all such paths, which is useful for debugging when using
`force=TRUE` with
[`gdataset.load()`](https://tanaylab.github.io/misha/reference/gdataset.load.md).

## See also

[`gtrack.dataset`](https://tanaylab.github.io/misha/reference/gtrack.dataset.md),
[`gtrack.exists`](https://tanaylab.github.io/misha/reference/gtrack.exists.md),
[`gtrack.ls`](https://tanaylab.github.io/misha/reference/gtrack.ls.md),
[`gdataset.ls`](https://tanaylab.github.io/misha/reference/gdataset.ls.md)

## Examples

``` r

gdb.init_examples()
gtrack.dbs("dense_track")
#>                    dense_track 
#> "/tmp/RtmptwYWUn/trackdb/test" 
gtrack.dbs(gtrack.ls(), dataframe = TRUE)
#>                 track                           db
#> 1         array_track /tmp/RtmptwYWUn/trackdb/test
#> 2         dense_track /tmp/RtmptwYWUn/trackdb/test
#> 3         rects_track /tmp/RtmptwYWUn/trackdb/test
#> 4        sparse_track /tmp/RtmptwYWUn/trackdb/test
#> 5 subdir.dense_track2 /tmp/RtmptwYWUn/trackdb/test
```
