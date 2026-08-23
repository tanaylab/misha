# Imports track attributes values

Imports track attributes values.

## Usage

``` r
gtrack.attr.import(table = NULL, remove.others = FALSE)
```

## Arguments

- table:

  a data frame containing attribute values

- remove.others:

  specifies what to do with the attributes that are not in the table

## Value

None.

## Details

This function makes imports attribute values contained in a data frame
'table'. The format of a table is similar to the one returned by
'gtrack.attr.export'. The values of the table must be character strings.
Column names of the table should specify the attribute names, while row
names should contain the track names.

The specified attributes of the specified tracks are modified. If an
attribute value is an empty string this attribute is removed from the
track.

If 'remove.others' is 'TRUE' all non-readonly attributes that do not
appear in the table are removed, otherwise they are preserved unchanged.

Error is reported on an attempt to modify a value of a read-only
attribute.

## See also

`gtrack.attr.import`,
[`gtrack.attr.set`](https://tanaylab.github.io/misha/reference/gtrack.attr.set.md),
[`gtrack.attr.get`](https://tanaylab.github.io/misha/reference/gtrack.attr.get.md),
[`gdb.get_readonly_attrs`](https://tanaylab.github.io/misha/reference/gdb.get_readonly_attrs.md)

## Examples

``` r

gdb.init_examples()
t <- gtrack.attr.export()
t$newattr <- as.character(1:dim(t)[1])
gtrack.attr.import(t)
gtrack.attr.export(attrs = "newattr")
#>                     newattr
#> array_track               1
#> dense_track               2
#> rects_track               3
#> sparse_track              4
#> subdir.dense_track2       5

# roll-back the changes
t$newattr <- ""
gtrack.attr.import(t)
```
