# Downloads files from FTP server

Downloads multiple files from FTP server

## Usage

``` r
gwget(url = NULL, path = NULL)
```

## Arguments

- url:

  URL of FTP server

- path:

  directory path where the downloaded files are stored

## Value

An array of file names that have been downloaded.

## Details

This function downloads files from FTP server given by 'url'. The
address in 'url' can contain wildcards to download more than one file at
once. Files are downloaded to a directory given by 'path' argument. If
'path' is 'NULL', file are downloaded into 'GROOT/downloads'.

## See also

[`gtrack.import_set`](https://tanaylab.github.io/misha/reference/gtrack.import_set.md)

## Examples

``` r
gdb.init_examples()
# \donttest{
outdir <- tempdir()
gwget("ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/md5sum.txt", path = outdir)
#> Downloading ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/md5sum.txt
# }
```
