# Export a track to BigWig format

Exports a track or track expression to BigWig format by first creating a
temporary bedGraph file and then converting it using `bedGraphToBigWig`
(or `wigToBigWig` as a fallback).

## Usage

``` r
gtrack.export_bigwig(track, file, intervals = NULL, iterator = NULL)
```

## Arguments

- track:

  track name or track expression (character string)

- file:

  output file path (typically ending in `.bw` or `.bigwig`).

- intervals:

  genomic intervals to export. If `NULL` (default), the entire genome
  (`.misha$ALLGENOME`) is used.

- iterator:

  iterator bin size. If `NULL` (default), the iterator is determined
  automatically from the track expression.

## Value

`NULL` (invisible). Called for its side effect of writing a file.

## Details

This function requires the UCSC `bedGraphToBigWig` utility to be
installed and available on the system PATH, or bundled with the misha
package. If not found, the function will raise an error with
installation instructions.

## See also

[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md),
[`gtrack.export_bedgraph`](https://tanaylab.github.io/misha/reference/gtrack.export_bedgraph.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md)

## Examples

``` r
if (FALSE) { # \dontrun{
gdb.init_examples()

# Export to BigWig (requires bedGraphToBigWig)
gtrack.export_bigwig("dense_track", "/tmp/dense.bw")

# With specific region
gtrack.export_bigwig("dense_track", "/tmp/dense_chr1.bw",
    intervals = gintervals(1, 0, 1e6)
)
} # }
```
