# Returns value of a track attribute

Returns value of a track attribute.

## Usage

``` r
gtrack.attr.get(track = NULL, attr = NULL)
```

## Arguments

- track:

  track name

- attr:

  attribute name

## Value

Track attribute value.

## Details

This function returns the value of a track attribute. If the attribute
does not exist an empty sting is returned.

## See also

[`gtrack.attr.import`](https://tanaylab.github.io/misha/reference/gtrack.attr.import.md),
[`gtrack.attr.set`](https://tanaylab.github.io/misha/reference/gtrack.attr.set.md)

## Examples

``` r

gdb.init_examples()
gtrack.attr.set("sparse_track", "test_attr", "value")
gtrack.attr.get("sparse_track", "test_attr")
#> [1] "value"
gtrack.attr.set("sparse_track", "test_attr", "")
```
