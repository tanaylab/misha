# Returns the definition of a virtual track

Returns the definition of a virtual track.

## Usage

``` r
gvtrack.info(vtrack = NULL)
```

## Arguments

- vtrack:

  virtual track name

## Value

Internal representation of a virtual track.

## Details

This function returns the internal representation of a virtual track.

## See also

[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)

## Examples

``` r

gdb.init_examples()
gvtrack.create("vtrack1", "dense_track", "max")
gvtrack.info("vtrack1")
#> $src
#> [1] "dense_track"
#> 
#> $func
#> [1] "max"
#> 
```
