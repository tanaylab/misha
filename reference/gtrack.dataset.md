# Returns the database/dataset path for a track

Returns the path of the database or dataset containing a track.

## Usage

``` r
gtrack.dataset(track = NULL)
```

## Arguments

- track:

  track name or a vector of track names

## Value

Character vector of database/dataset paths. Returns NA for non-existent
tracks.

## Details

When datasets are loaded, tracks can come from either the working
database or from loaded datasets. This function returns the source path
for each track.

## See also

[`gtrack.dbs`](https://tanaylab.github.io/misha/reference/gtrack.dbs.md),
[`gtrack.exists`](https://tanaylab.github.io/misha/reference/gtrack.exists.md),
[`gtrack.ls`](https://tanaylab.github.io/misha/reference/gtrack.ls.md),
[`gdataset.ls`](https://tanaylab.github.io/misha/reference/gdataset.ls.md)

## Examples

``` r

gdb.init_examples()
gtrack.dataset("dense_track")
#> [1] "/tmp/RtmpigTFOh/trackdb/test"
```
