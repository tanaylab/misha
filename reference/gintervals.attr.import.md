# Imports interval set attributes values

Imports interval set attributes values.

## Usage

``` r
gintervals.attr.import(table = NULL, remove.others = FALSE)
```

## Arguments

- table:

  a data frame containing attribute values

- remove.others:

  specifies what to do with the attributes that are not in the table

## Value

None.

## Details

This function imports attribute values contained in a data frame
'table'. The format of a table is similar to the one returned by
'gintervals.attr.export'. The values of the table must be character
strings. Column names of the table should specify the attribute names,
while row names should contain the interval set names.

The specified attributes of the specified interval sets are modified. If
an attribute value is an empty string this attribute is removed from the
interval set.

If 'remove.others' is 'TRUE' all attributes that do not appear in the
table are removed, otherwise they are preserved unchanged.

## See also

[`gintervals.attr.export`](https://tanaylab.github.io/misha/reference/gintervals.attr.export.md),
[`gintervals.attr.set`](https://tanaylab.github.io/misha/reference/gintervals.attr.set.md),
[`gintervals.attr.get`](https://tanaylab.github.io/misha/reference/gintervals.attr.get.md)

## Examples

``` r

gdb.init_examples()
t <- data.frame(
    myattr = "val1", row.names = "annotations",
    stringsAsFactors = FALSE
)
gintervals.attr.import(t)
gintervals.attr.export(attrs = "myattr")
#>             myattr
#> annotations   val1

# roll-back the changes
t$myattr <- ""
gintervals.attr.import(t)
```
