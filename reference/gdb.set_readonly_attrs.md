# Sets read-only track attributes

Sets read-only track attributes.

## Usage

``` r
gdb.set_readonly_attrs(attrs)
```

## Arguments

- attrs:

  a vector of read-only attributes names or 'NULL'

## Value

None.

## Details

This function sets the list of read-only track attributes. The specified
attributes may or may not already exist in the tracks.

If 'attrs' is 'NULL' the list of read-only attributes is emptied.

## See also

[`gdb.get_readonly_attrs`](https://tanaylab.github.io/misha/reference/gdb.get_readonly_attrs.md),
[`gtrack.attr.get`](https://tanaylab.github.io/misha/reference/gtrack.attr.get.md),
[`gtrack.attr.set`](https://tanaylab.github.io/misha/reference/gtrack.attr.set.md)
