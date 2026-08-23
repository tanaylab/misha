# Assigns value to a track attribute

Assigns value to a track attribute.

## Usage

``` r
gtrack.attr.set(track = NULL, attr = NULL, value = NULL)
```

## Arguments

- track:

  track name

- attr:

  attribute name

- value:

  value

## Value

None.

## Details

This function creates a track attribute and assigns 'value' to it. If
the attribute already exists its value is overwritten.

If 'value' is an empty string the attribute is removed.

Error is reported on an attempt to modify a value of a read-only
attribute.

## See also

[`gtrack.attr.get`](https://tanaylab.github.io/misha/reference/gtrack.attr.get.md),
[`gtrack.attr.import`](https://tanaylab.github.io/misha/reference/gtrack.attr.import.md),
[`gtrack.var.set`](https://tanaylab.github.io/misha/reference/gtrack.var.set.md),
[`gdb.get_readonly_attrs`](https://tanaylab.github.io/misha/reference/gdb.get_readonly_attrs.md)

## Examples

``` r

gdb.init_examples()
gtrack.attr.set("sparse_track", "test_attr", "value")
gtrack.attr.get("sparse_track", "test_attr")
#> [1] "value"
gtrack.attr.set("sparse_track", "test_attr", "")
```
