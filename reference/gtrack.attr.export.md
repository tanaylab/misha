# Returns track attributes values

Returns track attributes values.

## Usage

``` r
gtrack.attr.export(tracks = NULL, attrs = NULL)
```

## Arguments

- tracks:

  a vector of track names or 'NULL'

- attrs:

  a vector of attribute names or 'NULL'

## Value

A data frame containing track attributes values.

## Details

This function returns a data frame that contains track attributes
values. Column names of the data frame consist of the attribute names,
row names contain the track names.

The list of required tracks is specified by 'tracks' argument. If
'tracks' is 'NULL' the attribute values of all existing tracks are
returned.

Likewise the list of required attributes is controlled by 'attrs'
argument. If 'attrs' is 'NULL' all attribute values of the specified
tracks are returned. The columns are also sorted then by "popularity" of
an attribute, i.e. the number of tracks containing this attribute. This
sorting is not applied if 'attrs' is not 'NULL'.

Empty character string in a table cell marks a non-existing attribute.

## See also

[`gtrack.attr.import`](https://tanaylab.github.io/misha/reference/gtrack.attr.import.md),
[`gtrack.attr.get`](https://tanaylab.github.io/misha/reference/gtrack.attr.get.md),
[`gtrack.attr.set`](https://tanaylab.github.io/misha/reference/gtrack.attr.set.md)

## Examples

``` r

gdb.init_examples()
gtrack.attr.export()
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           created.by
#> array_track                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   gtrack.array.import("array_track", description, src = c("/home/aviezerl/temp/array1"))
#> dense_track                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    immaculate conception
#> rects_track                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              gtrack.2d.import(rects_track, description, c("/home/aviezerl/temp/rects1"))
#> sparse_track        gtrack.create_sparse(sparse_track, description, structure(list(chrom = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, , b[, 4])
#> subdir.dense_track2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  gtrack.create(subdir.dense_track2, description, dense_track * 2, iterator=NULL)
#>                                 created.date  description
#> array_track         Thu Sep  7 11:11:36 2023  array track
#> dense_track         Tue Dec  9 18:15:21 2014             
#> rects_track         Thu Sep  7 11:07:56 2023  rects track
#> sparse_track        Thu Sep  7 11:05:10 2023 sparse track
#> subdir.dense_track2 Wed Jun 11 10:38:22 2014             
gtrack.attr.export(tracks = c("sparse_track", "dense_track"))
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    created.by
#> sparse_track gtrack.create_sparse(sparse_track, description, structure(list(chrom = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, , b[, 4])
#> dense_track                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             immaculate conception
#>                          created.date  description
#> sparse_track Thu Sep  7 11:05:10 2023 sparse track
#> dense_track  Tue Dec  9 18:15:21 2014             
gtrack.attr.export(attrs = "created.by")
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           created.by
#> array_track                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   gtrack.array.import("array_track", description, src = c("/home/aviezerl/temp/array1"))
#> dense_track                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    immaculate conception
#> rects_track                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              gtrack.2d.import(rects_track, description, c("/home/aviezerl/temp/rects1"))
#> sparse_track        gtrack.create_sparse(sparse_track, description, structure(list(chrom = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, , b[, 4])
#> subdir.dense_track2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  gtrack.create(subdir.dense_track2, description, dense_track * 2, iterator=NULL)
```
