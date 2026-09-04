# List working database and loaded datasets

Returns a list of the working database and all loaded datasets.

## Usage

``` r
gdataset.ls(dataframe = FALSE)
```

## Arguments

- dataframe:

  If FALSE, return character vector; if TRUE, return data frame

## Value

Character vector of paths or data frame with detailed information

## See also

[`gdataset.load`](https://tanaylab.github.io/misha/reference/gdataset.load.md),
[`gdataset.info`](https://tanaylab.github.io/misha/reference/gdataset.info.md)

## Examples

``` r

dataset_path <- gdataset.example_path()
gdataset.load(dataset_path)
gdataset.ls()
#> [1] "/tmp/RtmptwYWUn/trackdb/test"              
#> [2] "/tmp/RtmptwYWUn/misha_dataset_1bb74d2e1c4b"
gdataset.unload(dataset_path)
```
