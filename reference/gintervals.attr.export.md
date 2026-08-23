# Returns interval set attributes values

Returns interval set attributes values.

## Usage

``` r
gintervals.attr.export(intervals.set = NULL, attrs = NULL)
```

## Arguments

- intervals.set:

  a vector of interval set names or 'NULL'

- attrs:

  a vector of attribute names or 'NULL'

## Value

A data frame containing interval set attribute values.

## Details

This function returns a data frame that contains interval set attribute
values. Column names of the data frame consist of the attribute names,
row names contain the interval set names.

The list of required interval sets is specified by 'intervals.set'
argument. If 'intervals.set' is 'NULL' the attribute values of all
existing interval sets are returned.

Likewise the list of required attributes is controlled by 'attrs'
argument. If 'attrs' is 'NULL' all attribute values of the specified
interval sets are returned. The columns are sorted then by "popularity"
of an attribute, i.e. the number of interval sets containing this
attribute. This sorting is not applied if 'attrs' is not 'NULL'.

Empty character string in a table cell marks a non-existing attribute.

## See also

[`gintervals.attr.import`](https://tanaylab.github.io/misha/reference/gintervals.attr.import.md),
[`gintervals.attr.get`](https://tanaylab.github.io/misha/reference/gintervals.attr.get.md),
[`gintervals.attr.set`](https://tanaylab.github.io/misha/reference/gintervals.attr.set.md)

## Examples

``` r

gdb.init_examples()
gintervals.attr.set("annotations", "test_attr", "value1")
gintervals.attr.export()
#>             test_attr
#> annotations    value1
gintervals.attr.export(intervals.set = "annotations")
#>             test_attr
#> annotations    value1
gintervals.attr.export(attrs = "test_attr")
#>             test_attr
#> annotations    value1
gintervals.attr.set("annotations", "test_attr", "")
```
