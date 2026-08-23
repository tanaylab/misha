# Create an example dataset on the fly

Creates a small dataset in a temporary directory using the built-in
example database. This function has side effects: it calls
[`gdb.init_examples`](https://tanaylab.github.io/misha/reference/gdb.init_examples.md)
which resets the working database, and it creates then deletes temporary
tracks ('example_dataset_track') and intervals
('example_dataset_intervals') in that database.

## Usage

``` r
gdataset.example_path()
```

## Value

Path to the created dataset directory (in a temporary location)

## Details

This function performs the following steps:

1.  Calls
    [`gdb.init_examples()`](https://tanaylab.github.io/misha/reference/gdb.init_examples.md)
    to set the working database

2.  Removes any existing 'example_dataset_track' and
    'example_dataset_intervals'

3.  Creates temporary track and intervals in the example database

4.  Saves them to a new dataset in a temporary directory

5.  Removes the temporary track and intervals from the example database

This is primarily intended for use in examples and tests. Users should
be aware that calling this function will change their current working
database.

## See also

[`gdataset.save`](https://tanaylab.github.io/misha/reference/gdataset.save.md),
[`gdataset.load`](https://tanaylab.github.io/misha/reference/gdataset.load.md),
[`gdb.init_examples`](https://tanaylab.github.io/misha/reference/gdb.init_examples.md)

## Examples

``` r

dataset_path <- gdataset.example_path()
gdataset.load(dataset_path)
gdataset.unload(dataset_path)
```
