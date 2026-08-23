# Change Database to Indexed Genome Format

Converts a per-chromosome database to indexed genome format with a
single consolidated genome.seq file and genome.idx index. Optionally
also converts tracks and interval sets to indexed format.

## Usage

``` r
gdb.convert_to_indexed(
  groot = NULL,
  remove_old_files = FALSE,
  force = FALSE,
  validate = TRUE,
  convert_tracks = FALSE,
  convert_intervals = FALSE,
  verbose = FALSE,
  chunk_size = 104857600,
  threads = NULL
)
```

## Arguments

- groot:

  Root directory of the database to change to indexed format. If NULL,
  uses the currently active database.

- remove_old_files:

  Logical. If TRUE, removes old per-chromosome files after successful
  conversion. Default: FALSE.

- force:

  Logical. If TRUE, forces the conversion without confirmation. Default:
  FALSE.

- validate:

  Logical. If TRUE, validates the conversion by comparing sequences.
  Default: TRUE.

- convert_tracks:

  Logical. If TRUE, also converts all eligible tracks to indexed format.
  Default: FALSE.

- convert_intervals:

  Logical. If TRUE, also converts all eligible interval sets to indexed
  format. Default: FALSE.

- verbose:

  Logical. If TRUE, prints verbose messages. Default: FALSE.

- chunk_size:

  Integer. The size of the chunk to read from the sequence files.
  Default: 104857600 (100MB). Reduce if you are running into memory
  issues.

- threads:

  Integer or NULL. Number of parallel processes used when converting
  tracks and interval sets (each worker handles one track/interval set
  via [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html)).
  If `NULL` (default), uses `min(parallel::detectCores(), 8)`. Set to
  `1` for serial execution. Falls back to serial on non-Unix platforms
  (mclapply requires fork).

## Value

Invisible NULL

## Details

This function converts a per-chromosome database (with separate .seq
files per contig) to indexed format (single genome.seq + genome.idx).
The indexed format provides better performance and scalability,
especially for genomes with many contigs.

**Important: Preserving Chromosome Order**

For exact conversion that produces bit-for-bit identical results before
and after conversion, you should load the source database first using
[`gsetroot()`](https://tanaylab.github.io/misha/reference/gdb.init.md)
or
[`gdb.init()`](https://tanaylab.github.io/misha/reference/gdb.init.md):

- If database is loaded: Uses chromosome order from ALLGENOME (exact
  preservation)

- If database is not loaded: Uses order from chrom_sizes.txt (may differ
  from ALLGENOME)

This ensures that the converted database has the exact same chromosome
ordering, which affects iteration order, interval IDs, and other
operations that depend on chromosome order.

The conversion process:

1.  Checks if database is already in indexed format

2.  Gets chromosome order from ALLGENOME (if loaded) or chrom_sizes.txt

3.  Consolidates all per-chromosome .seq files into genome.seq

4.  Creates genome.idx with CRC64 checksum

5.  Optionally validates the conversion

6.  Optionally removes old .seq files

7.  If convert_tracks=TRUE, converts all eligible tracks (1D: dense,
    sparse, array; 2D: rectangles, points)

8.  If convert_intervals=TRUE, converts all eligible interval sets (1D
    and 2D)

Tracks and intervals that cannot be converted (and are skipped):

- Tracks: virtual tracks, single-file tracks, already converted tracks

- Intervals: Single-file interval sets, already converted interval sets

## See also

[`gdb.create`](https://tanaylab.github.io/misha/reference/gdb.create.md),
[`gdb.init`](https://tanaylab.github.io/misha/reference/gdb.init.md),
[`gtrack.convert_to_indexed`](https://tanaylab.github.io/misha/reference/gtrack.convert_to_indexed.md),
[`gintervals.convert_to_indexed`](https://tanaylab.github.io/misha/reference/gintervals.convert_to_indexed.md),
[`gintervals.2d.convert_to_indexed`](https://tanaylab.github.io/misha/reference/gintervals.2d.convert_to_indexed.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Recommended: Load database first for exact conversion
gsetroot("/path/to/database")
gdb.convert_to_indexed(
    convert_tracks = TRUE,
    convert_intervals = TRUE,
    remove_old_files = TRUE,
    verbose = TRUE
)

# Convert current database to indexed format (genome only)
gdb.convert_to_indexed()

# Convert specific database without loading it first
# Note: chromosome order may differ from ALLGENOME
gdb.convert_to_indexed(groot = "/path/to/database")

# Convert genome and all tracks to indexed format
gdb.convert_to_indexed(convert_tracks = TRUE)

# Full conversion with validation and cleanup
gsetroot("/path/to/database") # Load first for exact order preservation
gdb.convert_to_indexed(
    convert_tracks = TRUE,
    convert_intervals = TRUE,
    remove_old_files = TRUE,
    validate = TRUE,
    verbose = TRUE
)
} # }
```
