# Creates an array track from array tracks or files

Creates an array track from array tracks or files.

## Usage

``` r
gtrack.array.import(track = NULL, description = NULL, ...)
```

## Arguments

- track:

  name of the newly created track

- description:

  a character string description

- ...:

  array track or name of a tab-delimited file

## Value

None.

## Details

This function creates a new 'Array' track from one or more "sources".
Each source can be either another 'Array' track or a tab-delimited file
that contains one-dimensional intervals and column values that should be
added to the newly created track. One can learn about the exact format
of the file by running 'gtrack.array.extract' or 'gextract' functions
with a 'file' parameter and inspecting the output file.

There might be more than one source used to create the new track. In
that case the new track will contain the columns from all the sources.
The equally named columns are merged. Intervals that appear in one
source but not in the other are added and the values for the missing
columns are set to NaN. Intervals with all NaN values are not added.
Partial overlaps between two intervals from different sources are
forbidden.

'description' is added as a track attribute.

## See also

[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md),
[`gtrack.array.extract`](https://tanaylab.github.io/misha/reference/gtrack.array.extract.md),
[`gtrack.array.set_colnames`](https://tanaylab.github.io/misha/reference/gtrack.array.set_colnames.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md)

## Examples

``` r

f1 <- tempfile()
gextract("sparse_track", gintervals(1, 5000, 20000), file = f1)
f2 <- tempfile()
gtrack.array.extract("array_track", c("col2", "col3", "col4"),
    gintervals(1, 0, 20000),
    file = f2
)
f3 <- tempfile()
gtrack.array.extract("array_track", c("col1", "col3"),
    gintervals(1, 0, 20000),
    file = f3
)

gtrack.array.import("test_track1", "Test array track 1", f1, f2)
gtrack.array.extract("test_track1", NULL, .misha$ALLGENOME)
#>    chrom start   end sparse_track col2 col3 col4 intervalID
#> 1   chr1     0    50          NaN    2  NaN    4          1
#> 2   chr1   100   150          NaN  102  NaN  NaN          1
#> 3   chr1   250   300          NaN  NaN  203  NaN          1
#> 4   chr1  1400  1450          NaN  302  303  304          1
#> 5   chr1  4750  4800          NaN  402  NaN  404          1
#> 6   chr1  6350  6400         0.44  502  NaN  NaN          1
#> 7   chr1  9900  9950         0.36  NaN  603  604          1
#> 8   chr1 10400 10450         0.44  702  NaN  NaN          1
#> 9   chr1 10550 10600         0.36  NaN  803  804          1
#> 10  chr1 12150 12200         0.40  NaN  903  NaN          1
#> 11  chr1 12650 12700         0.40 1002 1003 1004          1
#> 12  chr1 13550 13600         0.36 1102 1103  NaN          1
#> 13  chr1 14300 14350         0.36 1202 1203 1204          1
#> 14  chr1 15100 15150         0.40 1302 1303 1304          1
#> 15  chr1 15200 15250         0.40 1402 1403 1404          1
#> 16  chr1 16500 16550         0.36  NaN 1503  NaN          1
#> 17  chr1 16750 16800         0.40 1602  NaN 1604          1
#> 18  chr1 16950 17000         0.48 1702  NaN 1704          1
#> 19  chr1 17050 17100         0.44 1802 1803 1804          1
#> 20  chr1 17250 17300         0.56 1902 1903 1904          1
#> 21  chr1 17750 17800         0.48  NaN 2003 2004          1
#> 22  chr1 18000 18050         0.36 2102 2103 2104          1
#> 23  chr1 18250 18300         0.48  NaN 2203 2204          1

gtrack.array.import(
    "test_track2", "Test array track 2",
    "test_track1", f3
)
gtrack.array.extract("test_track2", NULL, .misha$ALLGENOME)
#>    chrom start   end sparse_track col2 col3 col4 col1 intervalID
#> 1   chr1     0    50          NaN    2  NaN    4  NaN          1
#> 2   chr1   100   150          NaN  102  NaN  NaN  NaN          1
#> 3   chr1   250   300          NaN  NaN  203  NaN  201          1
#> 4   chr1  1400  1450          NaN  302  303  304  301          1
#> 5   chr1  4750  4800          NaN  402  NaN  404  401          1
#> 6   chr1  6350  6400         0.44  502  NaN  NaN  501          1
#> 7   chr1  9900  9950         0.36  NaN  603  604  NaN          1
#> 8   chr1 10400 10450         0.44  702  NaN  NaN  701          1
#> 9   chr1 10550 10600         0.36  NaN  803  804  NaN          1
#> 10  chr1 12150 12200         0.40  NaN  903  NaN  901          1
#> 11  chr1 12650 12700         0.40 1002 1003 1004 1001          1
#> 12  chr1 13550 13600         0.36 1102 1103  NaN 1101          1
#> 13  chr1 14300 14350         0.36 1202 1203 1204 1201          1
#> 14  chr1 15100 15150         0.40 1302 1303 1304 1301          1
#> 15  chr1 15200 15250         0.40 1402 1403 1404 1401          1
#> 16  chr1 16500 16550         0.36  NaN 1503  NaN  NaN          1
#> 17  chr1 16750 16800         0.40 1602  NaN 1604  NaN          1
#> 18  chr1 16950 17000         0.48 1702  NaN 1704  NaN          1
#> 19  chr1 17050 17100         0.44 1802 1803 1804  NaN          1
#> 20  chr1 17250 17300         0.56 1902 1903 1904 1901          1
#> 21  chr1 17750 17800         0.48  NaN 2003 2004  NaN          1
#> 22  chr1 18000 18050         0.36 2102 2103 2104 2101          1
#> 23  chr1 18250 18300         0.48  NaN 2203 2204  NaN          1

gtrack.rm("test_track1", TRUE)
gtrack.rm("test_track2", TRUE)
unlink(c(f1, f2, f3))
```
