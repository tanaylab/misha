# Initializes connection with Genomic Database

Initializes connection with Genomic Database: loads the list of tracks,
intervals, etc.

## Usage

``` r
gdb.init(groot = NULL, dir = NULL, rescan = FALSE)

gsetroot(groot = NULL, dir = NULL, rescan = FALSE)
```

## Arguments

- groot:

  the root directory of the Genomic Database

- dir:

  the current working directory inside the Genomic Database

- rescan:

  indicates whether the file structure should be rescanned

## Value

None.

## Details

'gdb.init' initializes the connection with the Genomic Database. It is
typically called first prior to any other function. When the package is
attached it internally calls to 'gdb.init.examples' which opens the
connection with the database located at 'PKGDIR/trackdb/test' directory,
where 'PKGDIR' is the directory where the package is installed.

The current working directory inside the Genomic Database is set to
'dir'. If 'dir' is 'NULL', the current working directory is set to
'GROOT/tracks'.

If 'rescan' is 'TRUE', the list of tracks and intervals is achieved by
rescanning directory structure under the current current working
directory. Otherwise 'gdb.init' attempts to use the cached list that
resides in 'groot/.db.cache' file.

Upon completion the connection is established with the database. If
auto-completion mode is switched on (see 'gset_input_method') the list
of tracks and intervals sets is loaded and added as variables to the
global environment allowing auto-completion of object names with \<TAB\>
key. Also a few variables are defined at an environment called `.misha`,
and can be accessed using `.misha$variable`, e.g. `.misha$ALLGENOME`.
These variables should not be modified by user.

|  |  |
|----|----|
| GROOT | Root directory of Genomic Database |
| GWD | Current working directory inside Genomic Database |
| GTRACKS | List of all available tracks |
| GINTERVS | List of all available intervals |
| GVTRACKS | List of all available virtual tracks |
| ALLGENOME | List of all chromosomes and their sizes |
| GITERATOR.INTERVALS | A set of iterator intervals for which the track expression is evaluated |

When option 'gmulticontig.indexed_format' is set to TRUE, the function
loads a database with "indexed" track format.

In a per-chromosome database whose 'chrom_sizes.txt' names lack the
"chr" prefix of its '.seq' files, the chromosomes get that prefix, and
their order (in ALLGENOME and by chromosome id) is that of the sorted
names: the order R's [`order()`](https://rdrr.io/r/base/order.html)
gives in an en_US.UTF-8 session of an R built with ICU, whatever the
session's locale. A name with a character outside ASCII sorts after the
ASCII ones, byte by byte, which may differ from that order.

## See also

[`gdb.reload`](https://tanaylab.github.io/misha/reference/gdb.reload.md),
[`gdb.create`](https://tanaylab.github.io/misha/reference/gdb.create.md),
[`gdir.cd`](https://tanaylab.github.io/misha/reference/gdir.cd.md),
[`gtrack.ls`](https://tanaylab.github.io/misha/reference/gtrack.ls.md),
[`gintervals.ls`](https://tanaylab.github.io/misha/reference/gintervals.ls.md),
[`gvtrack.ls`](https://tanaylab.github.io/misha/reference/gvtrack.ls.md)
