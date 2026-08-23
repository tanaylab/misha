# Returns a list of read-only track attributes

Returns a list of read-only track attributes.

## Usage

``` r
gdb.get_readonly_attrs()
```

## Value

A list of read-only track attributes.

## Details

This function returns a list of read-only track attributes. These
attributes are not allowed to be modified or deleted.

If no attributes are marked as read-only a 'NULL' is returned.

## See also

[`gdb.set_readonly_attrs`](https://tanaylab.github.io/misha/reference/gdb.set_readonly_attrs.md),
[`gtrack.attr.get`](https://tanaylab.github.io/misha/reference/gtrack.attr.get.md),
[`gtrack.attr.set`](https://tanaylab.github.io/misha/reference/gtrack.attr.set.md)
