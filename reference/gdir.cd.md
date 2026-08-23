# Changes current working directory in Genomic Database

Changes current working directory in Genomic Database.

## Usage

``` r
gdir.cd(dir = NULL)
```

## Arguments

- dir:

  directory path

## Value

None.

## Details

This function changes the current working directory in Genomic Database
(not to be confused with shell's current working directory). The list of
database objects - tracks, intervals, track variables - is rescanned
recursively under 'dir'. Object names are updated with the respect to
the new current working directory. Example: a track named 'subdir.dense'
will be referred as 'dense' once current working directory is set to
'subdir'. All virtual tracks are removed.

## See also

[`gdb.init`](https://tanaylab.github.io/misha/reference/gdb.init.md),
[`gdir.cwd`](https://tanaylab.github.io/misha/reference/gdir.cwd.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md),
[`gdir.rm`](https://tanaylab.github.io/misha/reference/gdir.rm.md)

## Examples

``` r

gdb.init_examples()
gdir.cd("subdir")
gtrack.ls()
#> [1] "dense_track2"
gdir.cd("..")
gtrack.ls()
#> [1] "array_track"         "dense_track"         "rects_track"        
#> [4] "sparse_track"        "subdir.dense_track2"
```
