# Deletes a directory from Genomic Database

Deletes a directory from Genomic Database.

## Usage

``` r
gdir.rm(dir = NULL, recursive = FALSE, force = FALSE)
```

## Arguments

- dir:

  directory path

- recursive:

  if 'TRUE', the directory is deleted recursively

- force:

  if 'TRUE', suppresses user confirmation of tracks/intervals removal

## Value

None.

## Details

This function deletes a directory from Genomic Database. If 'recursive'
is 'TRUE', the directory is deleted with all the files/directories it
contains. If the directory contains tracks or intervals, the user is
prompted to confirm the deletion. Set 'force' to 'TRUE' to suppress the
prompt.

## See also

[`gdb.init`](https://tanaylab.github.io/misha/reference/gdb.init.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md),
[`gdir.cd`](https://tanaylab.github.io/misha/reference/gdir.cd.md),
[`gdir.cwd`](https://tanaylab.github.io/misha/reference/gdir.cwd.md)
