# Assigns value to an interval set attribute

Assigns value to an interval set attribute.

## Usage

``` r
gintervals.attr.set(intervals.set = NULL, attr = NULL, value = NULL)
```

## Arguments

- intervals.set:

  interval set name

- attr:

  attribute name

- value:

  value (character string)

## Value

None.

## Details

This function creates an interval set attribute and assigns 'value' to
it. If the attribute already exists its value is overwritten.

If 'value' is an empty string the attribute is removed.

## See also

[`gintervals.attr.get`](https://tanaylab.github.io/misha/reference/gintervals.attr.get.md),
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
