# Creates dense tracks from the columns of Parquet files

Creates one 'Dense' track per value column of one or more Parquet files.

## Usage

``` r
gtrack.import_parquet(
  description,
  path,
  binsize,
  track.prefix = NULL,
  columns = NULL,
  defval = NaN
)
```

## Arguments

- description:

  a character string description

- path:

  Parquet file path (may contain wildcards)

- binsize:

  bin size of the newly created 'Dense' tracks (\> 0)

- track.prefix:

  prefix for the track names

- columns:

  names of the value columns to import, or NULL for all

- defval:

  default track value

## Value

Names of the created tracks, invisibly.

## Details

The files must have the interval columns 'chrom' (string), 'start' and
'end' (0-based, end exclusive; integers, or floating point numbers that
are whole) and one or more value columns (integer or floating point),
such as the files
[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md)
writes with `file = "x.parquet"` or pandas' `to_parquet`. 'path' may
contain wildcards, so that per-chromosome files import in one call; all
the files must have the same value columns. Columns are matched by name,
not by position.

'columns' selects the value columns to import, in the order of the
created tracks. With 'columns = NULL' every column other than 'chrom',
'start' and 'end' is imported, and a column that is not numeric is an
error.

The track of column 'col' is named '\[track.prefix\]col'. All the track
names are checked before anything is written: a name must start with a
letter and hold only letters, digits, '\_' and '.' (columns are never
renamed; leave such a column out with 'columns'). Existing tracks are
not overwritten and no directories are created.

Values are binned as
[`gtrack.import`](https://tanaylab.github.io/misha/reference/gtrack.import.md)
bins the same rows read from a tab-delimited file: each bin holds the
mean of the values that overlap it, weighted by the overlap. Bases no
row covers, and rows whose value is null or NaN, count as 'defval'
(unless it is NaN). The rows of a chromosome may be spread over several
files and row groups in any order, but may not overlap. Chromosome names
are matched through the chromosome alias chain as in
[`gtrack.import`](https://tanaylab.github.io/misha/reference/gtrack.import.md),
with the same reporting of names the database does not have.

The import is all or nothing: if it fails or is interrupted, none of the
tracks is created.

Supported: flat schemas; PLAIN, dictionary and BYTE_STREAM_SPLIT
encodings; data pages v1 and v2; UNCOMPRESSED, SNAPPY and GZIP
compression. Other compression codecs (e.g. ZSTD) and encodings (e.g.
DELTA_BINARY_PACKED) are an error that names them. Unsigned integer
columns are read as unsigned; decimal, date, time and timestamp columns
are not numbers and are refused.

'description' is added to the created tracks as an attribute.

## See also

[`gtrack.import_set`](https://tanaylab.github.io/misha/reference/gtrack.import_set.md),
[`gtrack.import`](https://tanaylab.github.io/misha/reference/gtrack.import.md),
[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md)

## Examples

``` r

# \donttest{
gdb.init_examples()

f <- tempfile(fileext = ".parquet")
gextract("dense_track", "dense_track * 2", gintervals(1, 0, 10000), iterator = 50, file = f)
gtrack.import_parquet("from parquet", f, 50, track.prefix = "pq_", columns = "dense_track")
gextract("pq_dense_track", gintervals(1, 0, 200))
#>   chrom start end pq_dense_track intervalID
#> 1  chr1     0  50      0.1777778          1
#> 2  chr1    50 100      0.1600000          1
#> 3  chr1   100 150      0.1800000          1
#> 4  chr1   150 200      0.1600000          1
gtrack.rm("pq_dense_track", force = TRUE)
# }
```
