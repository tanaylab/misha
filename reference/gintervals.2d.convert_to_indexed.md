# Convert 2D interval set to indexed format

Converts a per-chromosome interval set to indexed format
(intervals2d.dat + intervals2d.idx) which reduces file descriptor usage.

## Usage

``` r
gintervals.2d.convert_to_indexed(
  set.name = NULL,
  remove.old = FALSE,
  force = FALSE
)
```

## Arguments

- set.name:

  name of 2D interval set to convert

- remove.old:

  if TRUE, removes old per-chromosome files after successful conversion

- force:

  if TRUE, re-converts even if already in indexed format

## Value

invisible NULL

## Details

The indexed format stores all chromosome pairs in a single
intervals2d.dat file with an intervals2d.idx index file. This
dramatically reduces file descriptor usage, especially for genomes with
many chromosomes (N\*(N-1)/2 files to just 2).

Only non-empty pairs are stored in the index, avoiding O(N^2) space
overhead.

The conversion process:

1.  Scans directory for existing per-pair files

2.  Creates temporary intervals2d.dat.tmp and intervals2d.idx.tmp files

3.  Concatenates all per-pair files into intervals2d.dat.tmp

4.  Builds index with pair offsets and checksums

5.  Atomically renames temporary files to final names

6.  Optionally removes old per-pair files

The indexed format is 100

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert a 2D interval set
gintervals.2d.convert_to_indexed("my_2d_intervals")

# Convert and remove old files
gintervals.2d.convert_to_indexed("my_2d_intervals", remove.old = TRUE)

# Force re-conversion
gintervals.2d.convert_to_indexed("my_2d_intervals", force = TRUE)
} # }
```
