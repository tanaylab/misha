# Load a dataset into the namespace

Loads tracks and intervals from a dataset directory, making them
available for analysis alongside the working database.

## Usage

``` r
gdataset.load(path, force = FALSE, verbose = FALSE)
```

## Arguments

- path:

  Path to a dataset or misha database directory

- force:

  If TRUE, ignore name collisions (working db wins; for
  dataset-to-dataset, later-loaded wins)

- verbose:

  If TRUE, print loaded track/interval names and summary counts

## Value

Invisibly returns a list with:

- tracks:

  Number of visible tracks loaded

- intervals:

  Number of visible intervals loaded

- shadowed_tracks:

  Number of tracks shadowed by collisions

- shadowed_intervals:

  Number of intervals shadowed by collisions

## See also

[`gdataset.unload`](https://tanaylab.github.io/misha/reference/gdataset.unload.md),
[`gdataset.save`](https://tanaylab.github.io/misha/reference/gdataset.save.md),
[`gdataset.ls`](https://tanaylab.github.io/misha/reference/gdataset.ls.md)

## Examples

``` r

dataset_path <- gdataset.example_path()
gdataset.load(dataset_path)
gdataset.unload(dataset_path)
```
