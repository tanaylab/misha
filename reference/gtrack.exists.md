# Tests for a track existence

Tests for a track existence.

## Usage

``` r
gtrack.exists(track = NULL)
```

## Arguments

- track:

  track name

## Value

'TRUE' if a track exists. Otherwise 'FALSE'.

## Details

This function returns 'TRUE' if a track exists in Genomic Database.

## See also

[`gtrack.ls`](https://tanaylab.github.io/misha/reference/gtrack.ls.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md),
[`gtrack.create`](https://tanaylab.github.io/misha/reference/gtrack.create.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md)

## Examples

``` r

gdb.init_examples()
gtrack.exists("dense_track")
#> [1] TRUE
```
