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
#> [1] "/tmp/Rtmpr1ocL6/trackdb/test"              
#> [2] "/tmp/Rtmpr1ocL6/misha_dataset_1ad566b3257d"
gdataset.unload(dataset_path)
```
