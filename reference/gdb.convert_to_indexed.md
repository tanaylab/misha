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

The converted database keeps the chromosome names and order that
[`gsetroot()`](https://tanaylab.github.io/misha/reference/gdb.init.md)
gives it, whether or not it is the loaded database, so its indexed
tracks and interval sets, iteration order and interval IDs stay the
same. In a per-chromosome database whose chrom_sizes.txt names lack the
"chr" prefix of its .seq files, that order is the sorted chromosome
names, not the chrom_sizes.txt order, and chrom_sizes.txt is rewritten
in it. The conversion stops, leaving the database as it was, unless the
indexed sequence holds exactly the chromosomes of chrom_sizes.txt: a
chromosome name with a character other than a letter, digit, '\_', '-'
or '.' (which the index would replace with '\_'), or a .seq file whose
length differs from chrom_sizes.txt, stops it, and so does a chromosome
named 'genome' (its genome.seq would be overwritten). A database whose
seq/ or chrom_sizes.txt is in another database, such as a dataset saved
with `copy_seq = FALSE` or a
[`gdb.create_linked()`](https://tanaylab.github.io/misha/reference/gdb.create_linked.md)
database, is not converted: convert that database. A seq/ or
chrom_sizes.txt linked to storage of its own elsewhere is converted
there. The session loaded before the conversion is put back as it was
after it, except the converted database if it was loaded as a dataset
and the conversion changed it: it is unloaded, with a warning. If the
conversion changed the sequence it reads (its seq/ is the converted
database's, with a chrom_sizes.txt of its own), a warning says so, and
nothing is loaded when its chromosomes would read other chromosomes'
sequence. chrom_sizes.txt keeps its group when the converting user
belongs to it; otherwise it takes that user's group.

The conversion process:

1.  Checks if database is already in indexed format

2.  Gets the chromosome names and order that
    [`gsetroot()`](https://tanaylab.github.io/misha/reference/gdb.init.md)
    gives the database

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
# Convert the loaded database with its tracks and interval sets
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
gdb.convert_to_indexed(groot = "/path/to/database")

# Convert genome and all tracks to indexed format
gdb.convert_to_indexed(convert_tracks = TRUE)

# Full conversion with validation and cleanup
gsetroot("/path/to/database")
gdb.convert_to_indexed(
    convert_tracks = TRUE,
    convert_intervals = TRUE,
    remove_old_files = TRUE,
    validate = TRUE,
    verbose = TRUE
)
} # }
```
