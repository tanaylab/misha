# Deletes a named intervals set

Deletes a named intervals set.

## Usage

``` r
gintervals.rm(intervals.set = NULL, force = FALSE, db = NULL)
```

## Arguments

- intervals.set:

  name of an intervals set

- force:

  if 'TRUE', suppresses user confirmation of a named intervals set
  removal

- db:

  optional database path. When multiple databases are connected, this
  specifies which database to delete the intervals set from. If NULL
  (the default), the intervals set is deleted from the working database
  (GROOT).

## Value

None.

## Details

This function deletes a named intervals set from the Genomic Database.
By default 'gintervals.rm' requires the user to interactively confirm
the deletion. Set 'force' to 'TRUE' to suppress the user prompt.

## See also

[`gintervals.save`](https://tanaylab.github.io/misha/reference/gintervals.save.md),
[`gintervals.exists`](https://tanaylab.github.io/misha/reference/gintervals.exists.md),
[`gintervals.ls`](https://tanaylab.github.io/misha/reference/gintervals.ls.md),
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md)

## Examples

``` r

gdb.init_examples()
intervs <- gintervals(c(1, 2))
gintervals.save("testintervs", intervs)
gintervals.ls()
#> [1] "annotations" "testintervs"
gintervals.rm("testintervs", force = TRUE)
gintervals.ls()
#> [1] "annotations"
```
