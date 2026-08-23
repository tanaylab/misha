# Save a dataset

Creates a new dataset directory containing selected tracks and/or
intervals from the working database.

## Usage

``` r
gdataset.save(
  path,
  description,
  tracks = NULL,
  intervals = NULL,
  symlinks = FALSE,
  copy_seq = FALSE
)
```

## Arguments

- path:

  Destination directory (must not exist)

- description:

  Required description for metadata

- tracks:

  Character vector of track names to include

- intervals:

  Character vector of interval set names to include

- symlinks:

  If TRUE, create symlinks to tracks/intervals instead of copying

- copy_seq:

  If TRUE, copy seq/ directory instead of symlinking

## Value

Invisible path

## See also

[`gdataset.load`](https://tanaylab.github.io/misha/reference/gdataset.load.md),
[`gdataset.info`](https://tanaylab.github.io/misha/reference/gdataset.info.md)

## Examples

``` r

gdb.init_examples()
example_intervs <- gintervals(1, 0, 10000)
gintervals.save("example_dataset_intervals", example_intervs)
gtrack.create(
    "example_dataset_track",
    "Example dataset track",
    "dense_track",
    iterator = "example_dataset_intervals"
)
dataset_path <- tempfile("misha_dataset_")
gdataset.save(
    path = dataset_path,
    description = "Example dataset",
    tracks = "example_dataset_track",
    intervals = "example_dataset_intervals"
)
gtrack.rm("example_dataset_track", force = TRUE)
gintervals.rm("example_dataset_intervals", force = TRUE)
```
