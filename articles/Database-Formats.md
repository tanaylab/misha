# Database Formats and Multi-Contig Support

``` r

library(misha)
gdb.init_examples()
```

## Overview

Starting with misha 5.3.0, databases can be stored in two formats:

- **Indexed format** (default, recommended): Single unified files for
  sequences and tracks
- **Per-chromosome format** (legacy): Separate files for each chromosome

The indexed format provides better performance and scalability,
especially for genomes with many contigs (\>50 chromosomes).

### Key Features

- **Automatic format detection** - misha automatically detects which
  format your database uses
- **Fully backward compatible** - existing databases continue to work
  without modification
- **Transparent to users** - same API for both formats
- **Migration tools** - convert databases when convenient
- **Performance benefits** - 4-14% faster for large-scale analyses

## Database Formats

### Indexed Format (Recommended)

The indexed format uses unified files:

**Sequence data:** - `seq/genome.seq` - All chromosome sequences
concatenated - `seq/genome.idx` - Index mapping chromosome names to
positions

**Track data:** - `tracks/mytrack.track/track.dat` - All chromosome data
concatenated - `tracks/mytrack.track/track.idx` - Index with
offset/length per chromosome

**Advantages:** - Fewer file descriptors (important for genomes with
100+ contigs) - Better performance for large workloads (14% faster) -
Smaller disk footprint - Faster track creation and conversion

### Per-Chromosome Format (Legacy)

The per-chromosome format uses separate files:

**Sequence data:** - `seq/chr1.seq`, `seq/chr2.seq`, … - One file per
chromosome

**Track data:** - `tracks/mytrack.track/chr1.track`, `chr2.track`, … -
One file per chromosome

**When to use:** - Compatibility with older misha versions (\<5.3.0) -
Small genomes (\<25 chromosomes) where performance difference is
negligible

## Creating Databases

### New Databases (Indexed Format)

By default, new databases use the indexed format:

``` r

# eval = FALSE: needs a genome FASTA on disk, or a multi-GB download from UCSC.
# Create database from FASTA file
gdb.create("mydb", "/path/to/genome.fa")

# Or download pre-built genome
gdb.create_genome("hg38", path = "/path/to/install")
```

### Force Legacy Format

To create a database in legacy format (for compatibility with older
misha):

``` r

# eval = FALSE: same as above, gdb.create() needs a genome FASTA.
# Set option before creating database
options(gmulticontig.indexed_format = FALSE)
gdb.create("mydb", "/path/to/genome.fa")
```

## Checking Database Format

Use
[`gdb.info()`](https://tanaylab.github.io/misha/reference/gdb.info.md)
to check your database format:

``` r

gdb.init_examples() # the bundled examples database
info <- gdb.info()
print(info$format) # "indexed" or "per-chromosome"
#> [1] "per-chromosome"
```

The full record:

``` r

gdb.info()
#> $path
#> [1] "/tmp/RtmpNJDl0O/trackdb/test"
#> 
#> $is_db
#> [1] TRUE
#> 
#> $format
#> [1] "per-chromosome"
#> 
#> $num_chromosomes
#> [1] 3
#> 
#> $genome_size
#> [1] 1e+06
#> 
#> $chromosomes
#>   chrom  size
#> 1     1 5e+05
#> 2     2 3e+05
#> 3     X 2e+05
```

## Converting Databases

### Convert Entire Database

Convert all tracks and sequences to indexed format:

``` r

gdb.convert_to_indexed()
gdb.info()$format
#> [1] "indexed"
```

This will: 1. Convert sequence files (`chr*.seq` →
`genome.seq + genome.idx`) 2. Convert all tracks to indexed format 3.
Validate conversions 4. Remove old files after successful conversion

### Convert Individual Tracks

Convert specific tracks while keeping others in legacy format:

``` r

gdb.init_examples() # start again from a per-chromosome database

# 1D tracks
gtrack.convert_to_indexed("dense_track")

# 2D tracks (rectangles and points)
gtrack.2d.convert_to_indexed("rects_track")
```

### Convert Intervals

Convert interval sets to indexed format:

``` r

# Only "big" interval sets are stored per chromosome, so lower the threshold
# to get big sets out of the small example database (default: 1,000,000).
options(gbig.intervals.size = 10)
gintervals.save("myintervals", gscreen("dense_track > 0.3"))
gintervals.save("my2dintervals", gextract("rects_track", gintervals.2d.all())[, 1:6])
options(gbig.intervals.size = 1e6)

# 1D intervals
gintervals.convert_to_indexed("myintervals")

# 2D intervals
gintervals.2d.convert_to_indexed("my2dintervals")
```

## Migration Guide

### When to Migrate

**High priority (significant benefits):** - Genomes with many contigs
(\>50 chromosomes) - Large-scale analyses (10M+ bp regions frequently) -
2D track workflows - File descriptor limit issues

**Medium priority (moderate benefits):** - Repeated extraction
workflows - Regular analyses on medium-sized regions (1-10M bp)

**Low priority (minimal benefits):** - Small genomes (\<25
chromosomes) - One-off analyses - Simple queries on small regions

### Migration Workflow

**Step 1: Backup (optional but recommended)**

``` r

# eval = FALSE: copies a whole database outside the vignette's temp directory.
# Create backup of important database
system("cp -r /path/to/mydb /path/to/mydb.backup")
```

**Step 2: Check current format**

``` r

gdb.init_examples()
info <- gdb.info()
print(paste("Current format:", info$format))
#> [1] "Current format: per-chromosome"
```

**Step 3: Convert**

``` r

gdb.convert_to_indexed()
```

**Step 4: Verify**

``` r

# Check format changed
info <- gdb.info()
print(paste("New format:", info$format))
#> [1] "New format: indexed"

# Test a few operations
result <- gextract("dense_track", gintervals(1, 0, 1000))
print(head(result))
#>   chrom start end dense_track intervalID
#> 1  chr1     0  50   0.1777778          1
#> 2  chr1    50 100   0.1600000          1
#> 3  chr1   100 150   0.1800000          1
#> 4  chr1   150 200   0.1600000          1
#> 5  chr1   200 250   0.1600000          1
#> 6  chr1   250 300   0.2000000          1
```

**Step 5: Remove backup (after validation)**

``` r

# eval = FALSE: deletes a directory, and there is nothing here to delete.
# After thorough testing
system("rm -rf /path/to/mydb.backup")
```

## Copying Tracks Between Databases

You can freely copy tracks between databases with different formats.

### Use `gtrack.copy()`

`gtrack.copy(src, db = target_db)` is the idiom. It reads the source
track directly, reconciles per-chromosome vs. indexed storage and any
difference in the two databases’ chromosome sets, and preserves the
track’s type. `src` may be a character vector, so a batch copy is a
single call.

``` r

gdb.init_examples()
source_db <- .misha$GROOT

# A second database to copy into
target_db <- file.path(tempdir(), "target_db")
unlink(target_db, recursive = TRUE)
gdb.create_linked(target_db, parent = source_db)
#> Created linked database at /tmp/RtmpNJDl0O/target_db (linked to /tmp/RtmpNJDl0O/trackdb/test)

# One track, or many
gtrack.copy("dense_track", db = target_db)
gtrack.copy(c("sparse_track", "array_track"), db = target_db)

gsetroot(target_db)
gtrack.ls()
#> [1] "array_track"  "dense_track"  "sparse_track"
gtrack.info("dense_track")[c("type", "size.in.bytes")]
#> $type
#> [1] "dense"
#> 
#> $size.in.bytes
#> [1] 80012
```

The destination does not have to be in the same format: converting it
afterwards leaves the tracks readable.

``` r

gdb.convert_to_indexed()
gdb.info()$format
#> [1] "indexed"
head(gextract("dense_track", gintervals(1, 0, 300)))
#>   chrom start end dense_track intervalID
#> 1  chr1     0  50   0.1777778          1
#> 2  chr1    50 100   0.1600000          1
#> 3  chr1   100 150   0.1800000          1
#> 4  chr1   150 200   0.1600000          1
#> 5  chr1   200 250   0.1600000          1
#> 6  chr1   250 300   0.2000000          1
```

### Round-tripping through a text file

Exporting with `gextract(..., file = ...)` and re-importing with
[`gtrack.import()`](https://tanaylab.github.io/misha/reference/gtrack.import.md)
also works, but it is lossy: `binsize = 0` builds a **sparse** track, so
a dense source track comes back sparse and about 5x larger (20 bytes per
value against 4 bytes per bin). Reach for it only when the destination
is not a misha database, or when you actually want a different track
type.

``` r

gsetroot(source_db)
exported <- file.path(tempdir(), "dense_track.txt")
gextract("dense_track", gintervals.all(), iterator = "dense_track", file = exported)

gsetroot(target_db)
gtrack.import("dense_track_roundtrip", "Round-tripped", exported, binsize = 0)
gtrack.info("dense_track")[c("type", "size.in.bytes")]
#> $type
#> [1] "dense"
#> 
#> $size.in.bytes
#> [1] 80012
gtrack.info("dense_track_roundtrip")[c("type", "size.in.bytes")]
#> $type
#> [1] "sparse"
#> 
#> $size.in.bytes
#> [1] 400120
gtrack.rm("dense_track_roundtrip", force = TRUE)
```

## Performance Comparison

Based on comprehensive benchmarks comparing indexed vs legacy formats:

### Operations Faster with Indexed Format

- **Very large workloads (10M+ bp)**: 14% faster
- **2D track operations**: 2% faster (was 36% slower in legacy)
- **Repeated extractions**: 14% faster
- **Real workflows**: 8% faster on average

### Operations Similar Performance

- Single chromosome extraction: Within 5%
- Multi-chromosome (10-22 chr): Within 1%
- PWM operations: Within 3%
- Small workloads (\<1M bp): Within 10%

### Summary

- **64% of operations are faster** with indexed format
- **93% within 5%** of legacy format (statistically equal)
- **Average: 4% faster** across all benchmarks
- **No regressions** for common use cases

## Backward Compatibility

### Fully Compatible

- ✅ Existing databases work without modification
- ✅ Existing scripts work without changes
- ✅ Same API for both formats
- ✅ Automatic format detection
- ✅ Can mix formats in same analysis

### Example: Mixed Environment

``` r

# Work with both formats in same session
gsetroot(source_db) # per-chromosome
data1 <- gextract("dense_track", gintervals(1, 0, 1000))
head(data1)
#>   chrom start end dense_track intervalID
#> 1  chr1     0  50   0.1777778          1
#> 2  chr1    50 100   0.1600000          1
#> 3  chr1   100 150   0.1800000          1
#> 4  chr1   150 200   0.1600000          1
#> 5  chr1   200 250   0.1600000          1
#> 6  chr1   250 300   0.2000000          1

gsetroot(target_db) # indexed
data2 <- gextract("dense_track", gintervals(1, 0, 1000))
head(data2)
#>   chrom start end dense_track intervalID
#> 1  chr1     0  50   0.1777778          1
#> 2  chr1    50 100   0.1600000          1
#> 3  chr1   100 150   0.1800000          1
#> 4  chr1   150 200   0.1600000          1
#> 5  chr1   200 250   0.1600000          1
#> 6  chr1   250 300   0.2000000          1
```

## Troubleshooting

### “File descriptor limit reached”

This occurs with many-contig genomes in legacy format:

**Solution:** Convert to indexed format

``` r

gdb.init_examples()
gdb.convert_to_indexed()
```

### “Track not found after copying files”

After manually copying track directories:

**Solution:** Reload database

``` r

gdb.reload()
```

### “Conversion fails with disk space error”

Conversion needs 2x track size temporarily:

**Solution:** Free disk space or convert tracks individually

``` r

# Convert one track at a time
gdb.init_examples()
gtrack.convert_to_indexed("dense_track")
gtrack.convert_to_indexed("sparse_track")
# etc.
```

## Best Practices

### For New Projects

1.  Use indexed format (default) for new databases
2.  Use
    [`gdb.create_genome()`](https://tanaylab.github.io/misha/reference/gdb.create_genome.md)
    for standard genomes
3.  Use
    [`gdb.create()`](https://tanaylab.github.io/misha/reference/gdb.create.md)
    with multi-FASTA for custom genomes

### For Existing Projects

1.  Check format with
    [`gdb.info()`](https://tanaylab.github.io/misha/reference/gdb.info.md)
2.  Migrate if beneficial (many contigs, large analyses)
3.  No rush - legacy format remains fully supported

## Summary

- **Indexed format is the default and recommended** for new databases
- **Legacy format remains fully supported** - no forced migration
- **Automatic format detection** - users don’t need to think about it
- **Significant performance benefits** for large-scale analyses
- **Easy migration** with
  [`gdb.convert_to_indexed()`](https://tanaylab.github.io/misha/reference/gdb.convert_to_indexed.md)
- **Fully backward compatible** - existing code works unchanged
