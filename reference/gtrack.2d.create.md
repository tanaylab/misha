# Creates a 'Rectangles' track from intervals and values

Creates a 'Rectangles' track from intervals and values.

## Usage

``` r
gtrack.2d.create(
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

  a set of two-dimensional intervals

- values:

  an array of numeric values - one for each interval

## Value

None.

## Details

This function creates a new 'Rectangles' (two-dimensional) track with
values at given intervals. 'description' is added as a track attribute.

When multiple databases are connected via
[`gsetroot`](https://tanaylab.github.io/misha/reference/gdb.init.md),
the track is created in the current working directory (.misha\$GWD),
which defaults to the last connected database. Use
[`gdir.cd`](https://tanaylab.github.io/misha/reference/gdir.cd.md) with
an absolute path to change where new tracks are created.

## See also

[`gtrack.create`](https://tanaylab.github.io/misha/reference/gtrack.create.md),
[`gtrack.create_sparse`](https://tanaylab.github.io/misha/reference/gtrack.create_sparse.md),
[`gtrack.smooth`](https://tanaylab.github.io/misha/reference/gtrack.smooth.md),
[`gtrack.modify`](https://tanaylab.github.io/misha/reference/gtrack.modify.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md),
[`gtrack.attr.get`](https://tanaylab.github.io/misha/reference/gtrack.attr.get.md)

## Examples

``` r

gdb.init_examples()
intervs1 <- gintervals.2d(
    1, (1:4) * 200, (1:4) * 200 + 100,
    1, (1:4) * 300, (1:4) * 300 + 200
)
intervs2 <- gintervals.2d(
    "X", (7:10) * 100, (7:10) * 100 + 50,
    2, (1:4) * 200, (1:4) * 200 + 130
)
intervs <- rbind(intervs1, intervs2)
gtrack.2d.create(
    "test_rects", "Test 2d track", intervs,
    runif(dim(intervs)[1], 1, 100)
)
gextract("test_rects", .misha$ALLGENOME)
#>   chrom1 start1 end1 chrom2 start2 end2 test_rects intervalID
#> 1   chr1    200  300   chr1    300  500   13.16927          1
#> 2   chr1    400  500   chr1    600  800   93.13802          1
#> 3   chr1    600  700   chr1    900 1100   94.41634          1
#> 4   chr1    800  900   chr1   1200 1400   76.89293          1
#> 5   chrX    700  750   chr2    200  330   82.59595          8
#> 6   chrX    800  850   chr2    400  530   39.15188          8
#> 7   chrX    900  950   chr2    600  730   70.77190          8
#> 8   chrX   1000 1050   chr2    800  930   20.99844          8
gtrack.rm("test_rects", force = TRUE)
```
