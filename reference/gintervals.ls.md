# Returns a list of named intervals sets

Returns a list of named intervals sets in Genomic Database.

## Usage

``` r
gintervals.ls(
  pattern = "",
  db = NULL,
  ignore.case = FALSE,
  perl = FALSE,
  fixed = FALSE,
  useBytes = FALSE
)
```

## Arguments

- pattern, ignore.case, perl, fixed, useBytes:

  see 'grep'

- db:

  optional database path to filter intervals. If specified, only
  interval sets from that database are returned.

## Value

An array that contains the names of intervals sets.

## Details

This function returns a list of named intervals sets that match the
pattern (see 'grep'). If called without any arguments all named
intervals sets are returned.

When multiple databases are connected, the 'db' parameter can be used to
filter intervals to only those from a specific database.

## See also

[`grep`](https://rdrr.io/r/base/grep.html),
[`gintervals.exists`](https://tanaylab.github.io/misha/reference/gintervals.exists.md),
[`gintervals.load`](https://tanaylab.github.io/misha/reference/gintervals.load.md),
[`gintervals.save`](https://tanaylab.github.io/misha/reference/gintervals.save.md),
[`gintervals.rm`](https://tanaylab.github.io/misha/reference/gintervals.rm.md),
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md),
[`gintervals.dataset`](https://tanaylab.github.io/misha/reference/gintervals.dataset.md)

## Examples

``` r

gdb.init_examples()
gintervals.ls()
#> [1] "annotations"
gintervals.ls(pattern = "annot*")
#> [1] "annotations"
```
