# Creates a track from WIG / BigWig / BedGraph / BED / tab-delimited file

Creates a track from WIG / BigWig / BedGraph / BED / tab-delimited file

## Usage

``` r
gtrack.import(
  track = NULL,
  description = NULL,
  file = NULL,
  binsize = NULL,
  defval = NaN,
  attrs = NULL
)
```

## Arguments

- track:

  track name

- description:

  a character string description

- file:

  file path

- binsize:

  bin size of the newly created 'Dense' track or '0' for a 'Sparse'
  track

- defval:

  default track value

- attrs:

  a named vector or list of attributes to be set on the track after
  import

## Value

None.

## Details

This function creates a track from WIG / BigWig / BedGraph /
tab-delimited file. Zipped files are supported (file name must have
'.gz' or '.zip' suffix).

Tab-delimited files must start with a header line with the following
column names (tab-separated): 'chrom', 'start', 'end', and exactly one
value column name (e.g. 'value'). Each subsequent line provides a single
interval: - chrom: chromosome name (e.g. 'chr1') - start: 0-based start
coordinate (inclusive) - end: 0-based end coordinate (exclusive) -
value: numeric value (floating point allowed); exactly one value column
is supported

Columns must be separated by tabs. Coordinates must refer to chromosomes
existing in the current genome. Missing values can be specified as
'NaN'.

Chromosome names are matched against the genome database through the
chromosome alias chain (`chr1` \<-\> `1`, `M` \<-\> `MT`, and
`<groot>/chrom_aliases.tsv` if present). If the file names chromosomes
and none of them can be matched, 'gtrack.import' fails instead of
creating an empty track; a file with no chromosome records at all is not
an error and still imports as an empty track. If only some of the names
can be matched, the rest of the file is imported and the unmatched names
are reported: a message for contigs the database does not have (unplaced
scaffolds, patches), a warning when the unmatched name looks like a
primary chromosome, since that suggests a naming mismatch.

BED files (.bed/.bed.gz/.bed.zip) are also supported. If the BED 'score'
column (5th column) exists and is numeric, it is used as the interval
value; otherwise a constant value of 1 is used. For BED inputs,
'binsize' controls the output type: if 'binsize' is 0 the track is
'Sparse'; otherwise the track is 'Dense' with bin-averaged values based
on overlaps with BED intervals (and 'defval' for regions not covered).

If 'binsize' is 0 the resulted track is created in 'Sparse' format.
Otherwise the 'Dense' format is chosen with a bin size equal to
'binsize'. The values that were not defined in input file file are
substituted by 'defval' value.

'description' is added as a track attribute.

When multiple databases are connected via
[`gsetroot`](https://tanaylab.github.io/misha/reference/gdb.init.md),
the track is created in the current working directory (.misha\$GWD),
which defaults to the last connected database. Use
[`gdir.cd`](https://tanaylab.github.io/misha/reference/gdir.cd.md) with
an absolute path to change where new tracks are created.

## See also

[`gtrack.import_set`](https://tanaylab.github.io/misha/reference/gtrack.import_set.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md),
[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md)

## Examples

``` r

# \donttest{
gdb.init_examples()

# Create a simple WIG file for demonstration
temp_file <- tempfile(fileext = ".wig")
writeLines(c(
    "track type=wiggle_0 name=\"example track\"",
    "fixedStep chrom=chr1 start=1 step=1",
    "1.5",
    "2.0",
    "1.8",
    "3.2"
), temp_file)

# Basic import
gtrack.import("example_track", "Example track from WIG file",
    temp_file,
    binsize = 1
)
gtrack.info("example_track")
#> $type
#> [1] "dense"
#> 
#> $dimensions
#> [1] 1
#> 
#> $size.in.bytes
#> [1] 4000012
#> 
#> $format
#> [1] "per-chromosome"
#> 
#> $bin.size
#> [1] 1
#> 
gtrack.rm("example_track", force = TRUE)

# Import with custom attributes
attrs <- c("author" = "researcher", "version" = "1.0", "experiment" = "test")
gtrack.import("example_track_with_attrs", "Example track with attributes",
    temp_file,
    binsize = 1, attrs = attrs
)

# Check that attributes were set
gtrack.attr.get("example_track_with_attrs", "author")
#> [1] "researcher"
gtrack.attr.get("example_track_with_attrs", "version")
#> [1] "1.0"
gtrack.attr.get("example_track_with_attrs", "experiment")
#> [1] "test"

# Clean up
gtrack.rm("example_track_with_attrs", force = TRUE)
# }
```
