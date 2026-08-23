# Create a linked database with symlinks to a parent database

Creates a new database directory structure with symbolic links to the
parent database's seq/ directory and chrom_sizes.txt file.

## Usage

``` r
gdb.create_linked(path, parent)
```

## Arguments

- path:

  Path for the new linked database

- parent:

  Path to the parent database (with seq and chrom_sizes.txt)

## Value

Invisible TRUE on success

## Details

This is useful for creating a writable database that shares sequence
data with a read-only main database. The new database can be set as the
working database via
[`gsetroot()`](https://tanaylab.github.io/misha/reference/gdb.init.md),
and then the parent database can be loaded as a dataset via
[`gdataset.load()`](https://tanaylab.github.io/misha/reference/gdataset.load.md).

## See also

[`gsetroot`](https://tanaylab.github.io/misha/reference/gdb.init.md),
[`gdb.create`](https://tanaylab.github.io/misha/reference/gdb.create.md),
[`gdataset.load`](https://tanaylab.github.io/misha/reference/gdataset.load.md),
[`gdataset.ls`](https://tanaylab.github.io/misha/reference/gdataset.ls.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create linked database sharing sequence data with main database
gdb.create_linked("~/my_tracks", parent = "/shared/genomics/hg38")

# Set linked database as working database and load parent as dataset
gsetroot("~/my_tracks")
gdataset.load("/shared/genomics/hg38")
} # }
```
