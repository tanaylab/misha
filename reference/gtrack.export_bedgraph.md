# Export a track to bedGraph format

Exports a track or track expression to a UCSC bedGraph file.

## Usage

``` r
gtrack.export_bedgraph(
  track,
  file,
  intervals = NULL,
  iterator = NULL,
  name = NULL
)
```

## Arguments

- track:

  track name or track expression (character string)

- file:

  output file path. If it ends in `.gz`, output is gzip-compressed.

- intervals:

  genomic intervals to export. If `NULL` (default), the entire genome
  (`.misha$ALLGENOME`) is used.

- iterator:

  iterator bin size. If `NULL` (default), the iterator is determined
  automatically from the track expression.

- name:

  track name for the bedGraph header line. If `NULL` (default), uses the
  `track` parameter value.

## Value

`NULL` (invisible). Called for its side effect of writing a file.

## Details

This function evaluates a track expression over the specified genomic
intervals and writes the result in standard bedGraph format (4-column,
tab-separated: chrom, start, end, value). NaN values are omitted from
the output.

The function supports physical tracks, virtual tracks, and arbitrary
track expressions (e.g. `"dense_track * 2"`). 2D tracks are not
supported.

If the output file path ends in `.gz`, the output is gzip-compressed.

## See also

[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md),
[`gtrack.export_bigwig`](https://tanaylab.github.io/misha/reference/gtrack.export_bigwig.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md)

## Examples

``` r
if (FALSE) { # \dontrun{
gdb.init_examples()

# Export a dense track
gtrack.export_bedgraph("dense_track", "/tmp/dense.bedgraph")

# Export with specific intervals
intervs <- gintervals(1, 0, 1000)
gtrack.export_bedgraph("dense_track", "/tmp/dense_chr1.bedgraph",
    intervals = intervs
)

# Export a track expression
gtrack.export_bedgraph("dense_track * 2", "/tmp/scaled.bedgraph",
    iterator = 100
)

# Export compressed
gtrack.export_bedgraph("dense_track", "/tmp/dense.bedgraph.gz")
} # }
```
