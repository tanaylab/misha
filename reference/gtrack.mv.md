# Renames or moves a track

Renames a track or moves it to a different namespace within the same
database.

## Usage

``` r
gtrack.mv(src = NULL, dest = NULL)
```

## Arguments

- src:

  source track name

- dest:

  destination track name

## Value

None.

## Details

This function renames a track or moves it to a different namespace
(directory) within the same database. The track cannot be moved to a
different database. Use
[`gtrack.copy`](https://tanaylab.github.io/misha/reference/gtrack.copy.md)
followed by
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md)
if you need to move a track between databases.

## See also

[`gtrack.copy`](https://tanaylab.github.io/misha/reference/gtrack.copy.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.exists`](https://tanaylab.github.io/misha/reference/gtrack.exists.md),
[`gtrack.ls`](https://tanaylab.github.io/misha/reference/gtrack.ls.md)

## Examples

``` r

gdb.init_examples()
gtrack.create_sparse("test_track", "Test", gintervals(1, 0, 100), 1)
gtrack.mv("test_track", "renamed_track")
gtrack.exists("renamed_track")
#> [1] TRUE
gtrack.rm("renamed_track", force = TRUE)
```
