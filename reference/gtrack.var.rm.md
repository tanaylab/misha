# Deletes a track variable

Deletes a track variable.

## Usage

``` r
gtrack.var.rm(track = NULL, var = NULL)
```

## Arguments

- track:

  track name

- var:

  track variable name

## Value

None.

## Details

This function deletes a track variable.

## See also

[`gtrack.var.get`](https://tanaylab.github.io/misha/reference/gtrack.var.get.md),
[`gtrack.var.set`](https://tanaylab.github.io/misha/reference/gtrack.var.set.md),
[`gtrack.var.ls`](https://tanaylab.github.io/misha/reference/gtrack.var.ls.md)

## Examples

``` r

gdb.init_examples()
gtrack.var.set("sparse_track", "test_var1", 1:10)
gtrack.var.set("sparse_track", "test_var2", "v")
gtrack.var.ls("sparse_track")
#> [1] "test_var1" "test_var2"
gtrack.var.rm("sparse_track", "test_var1")
gtrack.var.rm("sparse_track", "test_var2")
gtrack.var.ls("sparse_track")
#> character(0)
```
