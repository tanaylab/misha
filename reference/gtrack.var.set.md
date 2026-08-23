# Assigns value to a track variable

Assigns value to a track variable.

## Usage

``` r
gtrack.var.set(track = NULL, var = NULL, value = NULL)
```

## Arguments

- track:

  track name

- var:

  track variable name

- value:

  value

## Value

None.

## Details

This function creates a track variable and assigns 'value' to it. If the
track variable already exists its value is overwritten.

## See also

[`gtrack.var.get`](https://tanaylab.github.io/misha/reference/gtrack.var.get.md),
[`gtrack.var.ls`](https://tanaylab.github.io/misha/reference/gtrack.var.ls.md),
[`gtrack.var.rm`](https://tanaylab.github.io/misha/reference/gtrack.var.rm.md)

## Examples

``` r

gdb.init_examples()
gtrack.var.set("sparse_track", "test_var", 1:10)
gtrack.var.get("sparse_track", "test_var")
#>  [1]  1  2  3  4  5  6  7  8  9 10
gtrack.var.rm("sparse_track", "test_var")
```
