# Export a track to BigWig format

Exports a track or track expression to a bigWig file.

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

The track expression is evaluated over the intervals as by
[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md),
one batch of chromosomes at a time, and the values are written to the
file as they are extracted, so the whole genome's values are never held
in memory. Overlapping intervals are merged first, so each base is
written once. NaN values are not written. Values are stored as 32-bit
floats, as in every bigWig file.

The file is written by misha itself; no external tool is needed. Runs of
bins of equal size (as from a dense track or a numeric iterator) are
stored as fixedStep sections, other values as bedGraph sections. The
file has zoom levels (summaries at coarser resolutions), which genome
browsers and UCSC's `bigWigSummary` use for zoomed-out views.
Chromosomes are written in the order of their names (as the bigWig
format requires), and every chromosome of the genome is listed in the
file, also those without values.

The file is written under a temporary name in the same directory and
renamed to `file` once it is complete, so an export that fails or is
interrupted leaves a file already at `file` unchanged.

## See also

[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md),
[`gtrack.export_bedgraph`](https://tanaylab.github.io/misha/reference/gtrack.export_bedgraph.md),
[`gtrack.import`](https://tanaylab.github.io/misha/reference/gtrack.import.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md)

## Examples

``` r

gdb.init_examples()

f <- tempfile(fileext = ".bw")
gtrack.export_bigwig("dense_track", f)

# A region, and a track expression in 100 bp bins
gtrack.export_bigwig("dense_track * 2", f,
    intervals = gintervals(1, 0, 1e5), iterator = 100
)

# Read it back
head(gintervals.import_bigwig(f, gintervals(1, 0, 1000)))
#>   chrom start end     value
#> 1  chr1     0 100 0.3377778
#> 2  chr1   100 200 0.3400000
#> 3  chr1   200 300 0.3600000
#> 4  chr1   300 400 0.3200000
#> 5  chr1   400 500 0.2200000
#> 6  chr1   500 600 0.0800000
unlink(f)
```
