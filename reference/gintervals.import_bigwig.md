# Read the records of a bigWig file that overlap intervals

Reads the records of a bigWig file over a set of intervals without
importing the file into a track. Only the data blocks that the file's
index places over the intervals are read and decompressed, so reading a
few thousand regions of a large bigWig is fast.

## Usage

``` r
gintervals.import_bigwig(file, intervals = NULL)
```

## Arguments

- file:

  path to a bigWig file.

- intervals:

  intervals to read the records of. `NULL` reads all the records of the
  file.

## Value

A data frame with columns chrom, start, end and value.

## Details

Each record that overlaps `intervals` is returned once, whole (not
clipped to the intervals), also when it overlaps several of them, and
the records are sorted by chrom and start. The full-resolution records
are read (bedGraph, variableStep and fixedStep sections alike), never
the zoom levels; zero-length records are skipped. As in
[`gtrack.import`](https://tanaylab.github.io/misha/reference/gtrack.import.md),
chromosome names are matched to the database through its aliases (`chr1`
\<-\> `1`, ...), a file none of whose chromosomes the database has is an
error, and a record that runs past the end of its chromosome is clipped
to it (dropped if it starts past it) with a warning. With
`intervals = NULL`, records on chromosomes the database lacks are
skipped with a message, or a warning if one of them looks like a primary
chromosome. If no record overlaps `intervals`, the result is a data
frame of no rows with the same columns (a sparse signal often has no
record over a region); with `intervals = NULL`, a file of no records is
an error.

The result works directly as the source of a value-based virtual track
([`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)),
which computes the signal over intervals without importing the file.
`func = "sum"` adds up the values of the records an interval touches,
each counted once; `"weighted.sum"` weighs each by the number of its
bases inside the interval, which is the sum of the signal over the
interval's bases. An interval that no record overlaps gets `NaN` from
`sum`, `weighted.sum`, `avg`, `min`, `max` and `lse`, and 0 from `size`
and `exists`.

## See also

[`gtrack.import`](https://tanaylab.github.io/misha/reference/gtrack.import.md),
[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md),
[`gintervals.import_bed`](https://tanaylab.github.io/misha/reference/gintervals.import_bed.md).

## Examples

``` r

gdb.init_examples()
# dense_track.bw holds only chr1 and chr2, positions 0-5000, of the example
# database's dense_track: no record overlaps intervals outside them
f <- system.file("extdata", "dense_track.bw", package = "misha")
peaks <- gintervals(c(1, 2), c(1000, 3000), c(1500, 3600))
signal <- gintervals.import_bigwig(f, peaks)
head(signal)
#>   chrom start  end value
#> 1  chr1  1000 1050  0.10
#> 2  chr1  1050 1100  0.00
#> 3  chr1  1100 1150  0.02
#> 4  chr1  1150 1200  0.00
#> 5  chr1  1200 1250  0.00
#> 6  chr1  1250 1300  0.02

# the signal summed over each peak's bases, without importing the file
gvtrack.create("sig", src = signal, func = "weighted.sum")
gextract("sig", peaks, iterator = peaks)
#>   chrom start  end sig intervalID
#> 1  chr1  1000 1500  27          1
#> 2  chr2  3000 3600  74          2
```
