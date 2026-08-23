# Creates a 'Sparse' track from intervals and values

Creates a 'Sparse' track from intervals and values.

## Usage

``` r
gtrack.create_sparse(
  track = NULL,
  description = NULL,
  intervals = NULL,
  values = NULL
)
```

## Arguments

- track:

  track name

- description:

  a character string description

- intervals:

  a set of one-dimensional intervals

- values:

  an array of numeric values - one for each interval, in the same order
  as the rows of `intervals`. If `NULL`, the `value` column of
  `intervals` is used.

## Value

None.

## Details

This function creates a new 'Sparse' track with values at given
intervals. 'description' is added as a track attribute.

When multiple databases are connected via
[`gsetroot`](https://tanaylab.github.io/misha/reference/gdb.init.md),
the track is created in the current working directory (.misha\$GWD),
which defaults to the last connected database. Use
[`gdir.cd`](https://tanaylab.github.io/misha/reference/gdir.cd.md) with
an absolute path to change where new tracks are created.

`values` is matched to `intervals` row by row, in the order the
intervals are passed; `intervals` need not be sorted. Note however that
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md)
returns its result sorted in the canonical chromosome order, so building
`intervals` with
[`gintervals()`](https://tanaylab.github.io/misha/reference/gintervals.md)
while keeping `values` in the original order will bind values to the
wrong intervals. Keep the values in a `value` column of the intervals
data frame and omit `values` to make such a mismatch impossible.

## See also

[`gtrack.create`](https://tanaylab.github.io/misha/reference/gtrack.create.md),
[`gtrack.2d.create`](https://tanaylab.github.io/misha/reference/gtrack.2d.create.md),
[`gtrack.smooth`](https://tanaylab.github.io/misha/reference/gtrack.smooth.md),
[`gtrack.modify`](https://tanaylab.github.io/misha/reference/gtrack.modify.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md)

## Examples

``` r

gdb.init_examples()
intervs <- gintervals.load("annotations")
gtrack.create_sparse(
    "test_sparse", "Test track", intervs,
    1:dim(intervs)[1]
)
gextract("test_sparse", .misha$ALLGENOME)
#>   chrom start   end test_sparse intervalID
#> 1  chr1    20  2000           1          1
#> 2  chr1  2500  2600           2          1
#> 3  chr2    20  2000           3          2
#> 4  chr2  3000  8000           4          2
#> 5  chr2  9000 11000           5          2
#> 6  chr2 12000 12001           6          2
#> 7  chr2 13000 14000           7          2
#> 8  chr2 15000 15500           8          2
gtrack.rm("test_sparse", force = TRUE)
```
