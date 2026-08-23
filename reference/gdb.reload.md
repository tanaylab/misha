# Reloads database from the disk

Reloads database from disk: list of tracks, intervals, etc.

## Usage

``` r
gdb.reload(rescan = TRUE)
```

## Arguments

- rescan:

  indicates whether the file structure should be rescanned

## Value

No return value, called for side effects.

## Details

Reloads Genomic Database from disk: list of tracks, intervals, etc. Use
this function if you manually add tracks or if for any reason the
database becomes corrupted. If 'rescan' is 'TRUE', the list of tracks
and intervals is achieved by rescanning directory structure under the
current current working directory. Otherwise 'gdb.reload' attempts to
use the cached list that resides in 'GROOT/.db.cache' file.

## See also

[`gdb.init`](https://tanaylab.github.io/misha/reference/gdb.init.md),
[`gdb.create`](https://tanaylab.github.io/misha/reference/gdb.create.md),
[`gdir.cd`](https://tanaylab.github.io/misha/reference/gdir.cd.md),
