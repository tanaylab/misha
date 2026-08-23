# Deletes a track

Deletes a track.

## Usage

``` r
gtrack.rm(track = NULL, force = FALSE, db = NULL)
```

## Arguments

- track:

  track name

- force:

  if 'TRUE', suppresses user confirmation of a named track removal

- db:

  optional database path to delete the track from when multiple
  databases are connected

## Value

None.

## Details

This function deletes a track from the Genomic Database. By default
'gtrack.rm' requires the user to interactively confirm the deletion. Set
'force' to 'TRUE' to suppress the user prompt.

## See also

[`gtrack.exists`](https://tanaylab.github.io/misha/reference/gtrack.exists.md),
[`gtrack.ls`](https://tanaylab.github.io/misha/reference/gtrack.ls.md),
[`gtrack.create`](https://tanaylab.github.io/misha/reference/gtrack.create.md),
[`gtrack.2d.create`](https://tanaylab.github.io/misha/reference/gtrack.2d.create.md),
[`gtrack.create_sparse`](https://tanaylab.github.io/misha/reference/gtrack.create_sparse.md),
[`gtrack.smooth`](https://tanaylab.github.io/misha/reference/gtrack.smooth.md)

## Examples

``` r

gdb.init_examples()
gtrack.create("new_track", "Test track", "2 * dense_track")
gtrack.exists("new_track")
#> [1] TRUE
gtrack.rm("new_track", force = TRUE)
gtrack.exists("new_track")
#> [1] FALSE
```
