# Convert 1D interval set to indexed format

Converts a per-chromosome interval set to indexed format
(intervals.dat + intervals.idx) which reduces file descriptor usage.

## Usage

``` r
gintervals.convert_to_indexed(
  set.name = NULL,
  remove.old = FALSE,
  force = FALSE
)
```

## Arguments

- set.name:

  name of interval set to convert

- remove.old:

  if TRUE, removes old per-chromosome files after successful conversion

- force:

  if TRUE, re-converts even if already in indexed format

## Value

invisible NULL

## Details

The indexed format stores all chromosomes in a single intervals.dat file
with an intervals.idx index file. This reduces file descriptor usage
from N files (one per chromosome) to just 2 files.

The conversion process:

1.  Creates temporary intervals.dat.tmp and intervals.idx.tmp files

2.  Concatenates all per-chromosome files into intervals.dat.tmp

3.  Builds index with offsets and checksums

4.  Atomically renames temporary files to final names

5.  Optionally removes old per-chromosome files

The indexed format is 100

## See also

[`gintervals.save`](https://tanaylab.github.io/misha/reference/gintervals.save.md),
[`gintervals.load`](https://tanaylab.github.io/misha/reference/gintervals.load.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert an interval set
gintervals.convert_to_indexed("my_intervals")

# Convert and remove old files
gintervals.convert_to_indexed("my_intervals", remove.old = TRUE)

# Force re-conversion
gintervals.convert_to_indexed("my_intervals", force = TRUE)
} # }
```
