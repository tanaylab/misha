# Creates a track from a track expression

Creates a track from a track expression.

## Usage

``` r
gtrack.create(
  track = NULL,
  description = NULL,
  expr = NULL,
  iterator = NULL,
  band = NULL
)
```

## Arguments

- track:

  track name

- description:

  a character string description

- expr:

  track expression

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expression.

- band:

  track expression band. If 'NULL' no band is used.

## Value

None.

## Details

This function creates a new track named track. The values of the track
are determined by evaluation of 'expr' - a numeric track expression. The
type of the new track is determined by the type of the iterator. 'Fixed
bin', 'Sparse' or 'Rectangles' track can be created accordingly.
'description' is added as a track attribute.

When multiple databases are connected via
[`gsetroot`](https://tanaylab.github.io/misha/reference/gdb.init.md),
the track is created in the current working directory (.misha\$GWD),
which defaults to the last connected database. Use
[`gdir.cd`](https://tanaylab.github.io/misha/reference/gdir.cd.md) with
an absolute path to change where new tracks are created.

## See also

[`gtrack.2d.create`](https://tanaylab.github.io/misha/reference/gtrack.2d.create.md),
[`gtrack.create_sparse`](https://tanaylab.github.io/misha/reference/gtrack.create_sparse.md),
[`gtrack.smooth`](https://tanaylab.github.io/misha/reference/gtrack.smooth.md),
[`gtrack.modify`](https://tanaylab.github.io/misha/reference/gtrack.modify.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md)

## Examples

``` r

gdb.init_examples()

## Creates a new track that is a sum of values from 'dense' and
## 2 * non-nan values of 'sparse' track. The new track type is
## Dense with a bin size that equals to '70'.
gtrack.create("mixed_track", "Test track",
    "dense_track +
              replace(sparse_track, is.nan(sparse_track), 0) * 2",
    iterator = 70
)
gtrack.info("mixed_track")
#> $type
#> [1] "dense"
#> 
#> $dimensions
#> [1] 1
#> 
#> $size.in.bytes
#> [1] 57160
#> 
#> $format
#> [1] "per-chromosome"
#> 
#> $bin.size
#> [1] 70
#> 
gtrack.rm("mixed_track", force = TRUE)
```
