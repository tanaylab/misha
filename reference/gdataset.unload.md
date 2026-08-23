# Unload a dataset from the namespace

Removes all tracks and intervals from a previously loaded dataset. If a
track was shadowing another, the shadowed track becomes visible again.

## Usage

``` r
gdataset.unload(path, validate = FALSE)
```

## Arguments

- path:

  Path to a previously loaded dataset

- validate:

  If TRUE, error if path is not currently loaded; otherwise silently
  no-op

## Value

Invisible NULL

## See also

[`gdataset.load`](https://tanaylab.github.io/misha/reference/gdataset.load.md),
[`gdataset.ls`](https://tanaylab.github.io/misha/reference/gdataset.ls.md)

## Examples

``` r

dataset_path <- gdataset.example_path()
gdataset.load(dataset_path)
gdataset.unload(dataset_path, validate = TRUE)
```
