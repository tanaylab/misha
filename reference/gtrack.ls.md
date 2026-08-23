# Returns a list of track names

Returns a list of track names in Genomic Database.

## Usage

``` r
gtrack.ls(
  ...,
  pattern = NULL,
  db = NULL,
  ignore.case = FALSE,
  perl = FALSE,
  fixed = FALSE,
  useBytes = FALSE
)
```

## Arguments

- ...:

  additional arguments of either form 'pattern' (matched against track
  names) or 'attribute = pattern' (matched against a track attribute
  value)

- pattern:

  optional pattern to match against track names, equivalent to passing
  it positionally (e.g. `gtrack.ls("den*")` and
  `gtrack.ls(pattern = "den*")` are the same). Note that this means a
  track attribute literally named "pattern" can no longer be filtered by
  name via `...`; use `gtrack.ls()` followed by
  [`gtrack.attr.get`](https://tanaylab.github.io/misha/reference/gtrack.attr.get.md)
  instead.

- db:

  optional database path to filter tracks. If specified, only tracks
  from that database are returned.

- ignore.case, perl, fixed, useBytes:

  see 'grep'

## Value

An array that contains the names of tracks that match the supplied
patterns.

## Details

This function returns a list of tracks whose name or track attribute
value match a pattern (see 'grep'). If called without any arguments all
tracks are returned.

If pattern is specified without a track attribute (i.e. in the form of
'pattern', either positionally or via the 'pattern' argument) then
filtering is applied to the track names. If pattern is supplied with a
track attribute (i.e. in the form of 'name = pattern') then track
attribute is matched against the pattern.

Multiple patterns are applied one after another. The resulted list of
tracks should match all the patterns.

When multiple databases are connected, the 'db' parameter can be used to
filter tracks to only those from a specific database.

## See also

[`grep`](https://rdrr.io/r/base/grep.html),
[`gtrack.exists`](https://tanaylab.github.io/misha/reference/gtrack.exists.md),
[`gtrack.create`](https://tanaylab.github.io/misha/reference/gtrack.create.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.dataset`](https://tanaylab.github.io/misha/reference/gtrack.dataset.md)

## Examples

``` r

gdb.init_examples()

# get all track names
gtrack.ls()
#> [1] "array_track"         "dense_track"         "rects_track"        
#> [4] "sparse_track"        "subdir.dense_track2"

# get track names that match the pattern "den*"
gtrack.ls("den*")
#> [1] "dense_track"         "subdir.dense_track2"

# equivalent, using the 'pattern' argument
gtrack.ls(pattern = "den*")
#> [1] "dense_track"         "subdir.dense_track2"

# get track names whose "created.by" attribute match the pattern
# "create_sparse"
gtrack.ls(created.by = "create_sparse")
#> [1] "sparse_track"

# get track names whose names match the pattern "den*" and whose
# "created.by" attribute match the pattern "track"
gtrack.ls("den*", created.by = "track")
#> [1] "subdir.dense_track2"
```
