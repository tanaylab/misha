# Get dataset information

Returns metadata and contents of a dataset.

## Usage

``` r
gdataset.info(path)
```

## Arguments

- path:

  Path to any dataset (loaded or not)

## Value

List with dataset information

## See also

[`gdataset.ls`](https://tanaylab.github.io/misha/reference/gdataset.ls.md),
[`gdataset.load`](https://tanaylab.github.io/misha/reference/gdataset.load.md)

## Examples

``` r

dataset_path <- gdataset.example_path()
gdataset.info(dataset_path)
#> $description
#> [1] "Example dataset created on the fly"
#> 
#> $author
#> [1] "runner"
#> 
#> $created
#> [1] "2026-08-25T11:47:11Z"
#> 
#> $original_db
#> [1] "/tmp/RtmpQPXCab/trackdb/test"
#> 
#> $misha_version
#> [1] "5.11.24"
#> 
#> $track_count
#> [1] 1
#> 
#> $interval_count
#> [1] 1
#> 
#> $genome
#> [1] "a2f87ab337da0fdc0a930472e56f4587"
#> 
#> $is_loaded
#> [1] FALSE
#> 
```
