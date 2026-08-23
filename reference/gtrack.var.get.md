# Returns value of a track variable

Returns value of a track variable.

## Usage

``` r
gtrack.var.get(track = NULL, var = NULL)
```

## Arguments

- track:

  track name

- var:

  track variable name

## Value

Track variable value.

## Details

This function returns the value of a track variable. If the variable
does not exist an error is reported.

## See also

[`gtrack.var.set`](https://tanaylab.github.io/misha/reference/gtrack.var.set.md),
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
