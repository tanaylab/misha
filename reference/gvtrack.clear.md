# Deletes all virtual tracks

Deletes all virtual tracks of the current working directory.

## Usage

``` r
gvtrack.clear()
```

## Value

None.

## Details

This function removes every virtual track defined for the current
working directory in one call. After it returns
[`gvtrack.ls`](https://tanaylab.github.io/misha/reference/gvtrack.ls.md)
reports no virtual tracks. It is a convenient way to reset virtual-track
state between analyses. Use
[`gvtrack.rm`](https://tanaylab.github.io/misha/reference/gvtrack.rm.md)
to remove a single virtual track.

## See also

[`gvtrack.rm`](https://tanaylab.github.io/misha/reference/gvtrack.rm.md),
[`gvtrack.ls`](https://tanaylab.github.io/misha/reference/gvtrack.ls.md),
[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)

## Examples

``` r

gdb.init_examples()
gvtrack.create("vtrack1", "dense_track", "max")
gvtrack.create("vtrack2", "dense_track", "avg")
gvtrack.ls()
#> [1] "g_frac"      "c_frac"      "cg_frac"     "masked_frac" "vtrack1"    
#> [6] "vtrack2"    
gvtrack.clear()
gvtrack.ls()
#> NULL
```
