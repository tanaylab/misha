# Returns value of an interval set attribute

Returns value of an interval set attribute.

## Usage

``` r
gintervals.attr.get(intervals.set = NULL, attr = NULL)
```

## Arguments

- intervals.set:

  interval set name

- attr:

  attribute name

## Value

Interval set attribute value (character string).

## Details

This function returns the value of an interval set attribute. If the
attribute does not exist an empty string is returned.

## See also

[`gintervals.attr.set`](https://tanaylab.github.io/misha/reference/gintervals.attr.set.md),
[`gintervals.attr.export`](https://tanaylab.github.io/misha/reference/gintervals.attr.export.md),
[`gintervals.attr.import`](https://tanaylab.github.io/misha/reference/gintervals.attr.import.md)

## Examples

``` r

gdb.init_examples()
gintervals.attr.set("annotations", "test_attr", "value")
gintervals.attr.get("annotations", "test_attr")
#> [1] "value"
gintervals.attr.set("annotations", "test_attr", "")
```
