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
#> [1] "/tmp/RtmplYc5bX/trackdb/test"              
#> [2] "/tmp/RtmplYc5bX/misha_dataset_1bb175eb3c0e"
gdataset.unload(dataset_path)
```
