# Returns the database/dataset path for interval sets

Returns the path of the database or dataset containing an interval set.

## Usage

``` r
gintervals.dataset(intervals = NULL)
```

## Arguments

- intervals:

  interval set name or a vector of interval set names

## Value

A character vector containing the database paths for each interval set.
Returns NA for interval sets that don't exist in any connected database.

## Details

When datasets are loaded, interval sets can come from either the working
database or from loaded datasets. This function returns the source path
for each interval set.

## See also

[`gintervals.dbs`](https://tanaylab.github.io/misha/reference/gintervals.dbs.md),
[`gintervals.exists`](https://tanaylab.github.io/misha/reference/gintervals.exists.md),
[`gintervals.ls`](https://tanaylab.github.io/misha/reference/gintervals.ls.md),
[`gdataset.ls`](https://tanaylab.github.io/misha/reference/gdataset.ls.md)

## Examples

``` r

gdb.init_examples()
gintervals.dataset("annotations1")
#> [1] NA
```
