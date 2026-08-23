# Creates a new directory in Genomic Database

Creates a new directory in Genomic Database.

## Usage

``` r
gdir.create(dir = NULL, showWarnings = TRUE, mode = "0777")
```

## Arguments

- dir:

  directory path

- showWarnings:

  see 'dir.create'

- mode:

  see 'dir.create'

## Value

None.

## Details

This function creates a new directory in Genomic Database. Creates only
the last element in the specified path.

## Note

A new directory cannot be created within an existing track directory.

## See also

[`dir.create`](https://rdrr.io/r/base/files2.html),
[`gdb.init`](https://tanaylab.github.io/misha/reference/gdb.init.md),
[`gdir.cwd`](https://tanaylab.github.io/misha/reference/gdir.cwd.md),
[`gdir.rm`](https://tanaylab.github.io/misha/reference/gdir.rm.md)
