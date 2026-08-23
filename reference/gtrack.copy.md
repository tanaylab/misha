# Copies one or more tracks

Creates a copy of an existing track, optionally to a different database.
Transparently handles format mismatches (per-chromosome vs indexed) and
chromosome-order differences between source and destination databases.

## Usage

``` r
gtrack.copy(src = NULL, dest = NULL, db = NULL, overwrite = FALSE)
```

## Arguments

- src:

  source track name(s). Either a single name or a character vector of
  names.

- dest:

  destination name. If `src` is a single name, this is the destination
  track name (defaults to `src`). If `src` is a vector, `dest` is
  treated as a namespace prefix (e.g. `"ns"` produces `"ns.track1"`,
  `"ns.track2"`, ...). NULL keeps each track's name.

- db:

  destination database root. Must be the current `GROOT` or a member of
  `GDATASETS`. NULL means the current working directory db.

- overwrite:

  if TRUE, replace an existing destination track.

## Value

invisibly, the character vector of created track names.

## Details

Chromosomes that exist in the source database but not in the destination
are dropped with a warning.

For 2D tracks (rectangles, points), cross-database copy requires
identical chromosome order in source and destination. Format conversion
(per-chromosome to indexed) is supported only when the destination is
the active database.

## See also

[`gtrack.mv`](https://tanaylab.github.io/misha/reference/gtrack.mv.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md)

## Examples

``` r

gdb.init_examples()
gtrack.copy("dense_track", "dense_track_copy")
gtrack.exists("dense_track_copy")
#> [1] TRUE
gtrack.rm("dense_track_copy", force = TRUE)
```
