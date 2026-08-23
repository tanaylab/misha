# Returns all database paths containing an interval set

Returns all database paths that contain a version of an interval set.

## Usage

``` r
gintervals.dbs(intervals = NULL, dataframe = FALSE)
```

## Arguments

- intervals:

  interval set name

- dataframe:

  return a data frame with columns `intervals` and `db`

## Value

A named character vector of database paths. If `dataframe` is TRUE,
returns a data frame with columns `intervals` and `db`.

## Details

When datasets are loaded, an interval set may exist in multiple
locations. This function computes on-demand and returns all such paths.

## See also

[`gintervals.dataset`](https://tanaylab.github.io/misha/reference/gintervals.dataset.md),
[`gintervals.ls`](https://tanaylab.github.io/misha/reference/gintervals.ls.md),
[`gdataset.ls`](https://tanaylab.github.io/misha/reference/gdataset.ls.md)

## Examples

``` r

gdb.init_examples()
gintervals.dbs("annotations1")
#> annotations1 
#>           NA 
```
