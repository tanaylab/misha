# Fixed-width windows as arrays for machine learning

Extracts fixed-width windows around a set of loci and returns their DNA
as a one-hot array and the values of track expressions over them as a
numeric array, ready to feed a sequence model. One call replaces
[`gseq.extract`](https://tanaylab.github.io/misha/reference/gseq.extract.md),
a one-hot encoding of the strings and
[`gextract.matrix`](https://tanaylab.github.io/misha/reference/gextract.matrix.md);
the one-hot array is encoded in C++ straight from the sequence files.

## Usage

``` r
gextract.windows(
  intervals = NULL,
  tracks = NULL,
  width = NULL,
  signal_width = width,
  jitter = 0,
  binsize = 1,
  use_strand = TRUE,
  exclude = NULL,
  min_counts = NULL,
  max_counts = NULL,
  count_track = 1,
  sequence = TRUE,
  dtype = c("raw", "integer", "double")
)
```

## Arguments

- intervals:

  the loci: 1D intervals (a data frame or the name of an intervals set),
  or a list of such, e.g. `list(pos = peaks, neg = background)`, whose
  rows are concatenated and labelled by a 'set' column (the names of the
  list, or their indices)

- tracks:

  track expressions (a character vector) for 'y', or NULL for sequence
  only. Its names, if any, name the second dimension of 'y'.

- width:

  width of the sequence window

- signal_width:

  width of the signal window

- jitter:

  extra margin added to each side of both windows

- binsize:

  bin size of the signals: a positive integer that divides
  `signal_width + 2 * jitter`. Each value is the one
  [`gextract.matrix`](https://tanaylab.github.io/misha/reference/gextract.matrix.md)
  returns for the bin (e.g. the mean of a dense track; use a virtual
  track with `func = "sum"` for sums).

- use_strand:

  if TRUE, rows with strand -1 are reverse complemented (see Details)

- exclude:

  intervals whose overlapping windows are dropped, or NULL

- min_counts, max_counts:

  bounds on the sum of 'count_track' over the signal window (inclusive),
  or NULL

- count_track:

  the track expression the count filters use: an index into 'tracks', or
  one of its names

- sequence:

  if FALSE, 'x' is not computed (NULL)

- dtype:

  storage type of 'x': "raw" (one byte per value, the default),
  "integer" or "double"

## Value

A list with

- x:

  the one-hot sequences, an array of dimensions
  `c(N, 4, width + 2 * jitter)` with channels A, C, G, T (see
  [`gseq.onehot`](https://tanaylab.github.io/misha/reference/gseq.onehot.md)),
  or NULL

- y:

  the signals, a numeric array of dimensions
  `c(N, T, (signal_width + 2 * jitter) / binsize)`, or NULL without
  'tracks'

- intervals:

  the kept rows, one per window, in their original order: 'chrom',
  'start' and 'end' of the sequence window (jitter included), the
  'strand' column if any, 'signal_start' and 'signal_end' with tracks,
  the other columns of 'intervals', 'index' (the row in 'intervals',
  1-based) and, for a list, 'set'

- dropped:

  the dropped rows: 'index' and 'reason'

## Details

Each locus (a row of 'intervals') is centred at
`mid = start + (end - start) %/% 2`. Its sequence window is
`[mid - width %/% 2 - jitter, mid - width %/% 2 + width + jitter)` and
its signal window is the same with `signal_width`: an odd width reaches
one base further right than left. A positive `jitter` widens both
windows by `jitter` on each side, so that a data loader can crop
randomly shifted windows of the requested widths later without reading
the genome again (as in tangermeme's `extract_loci`).

Rows are dropped, and listed in `dropped` with the first reason that
applies, when

- "off_chrom": the sequence window or (with tracks) the signal window
  runs past either end of the chromosome;

- "excluded": either window overlaps an interval of 'exclude' (by at
  least one base; intervals of 'exclude' on chromosomes the genome does
  not have are ignored);

- "min_counts" / "max_counts": the sum of the values of track expression
  'count_track' over the signal window (as returned, i.e. after binning,
  NaN counted as 0) is below 'min_counts' or above 'max_counts'.

With `use_strand = TRUE` and a 'strand' column in 'intervals', the rows
with strand -1 (or "-") are reverse complemented: their one-hot encoding
is that of the reverse complement, and their signals are reversed along
the length axis.

The arrays have the same axes as pymisha's `gextract_windows`:
`x[n, c, p]` and `y[n, t, b]`, i.e. `(batch, channels, length)`. (R
arrays are column-major, so the memory order differs from numpy's;
indexing, and conversion to a torch tensor, are the same.)

Memory: 'x' takes `N * 4 * (width + 2 * jitter)` bytes with
`dtype = "raw"` (4 and 8 times that for "integer" and "double"), and 'y'
takes `N * T * (signal_width + 2 * jitter) / binsize * 8` bytes for T
track expressions. The signals are extracted with
[`gextract.matrix`](https://tanaylab.github.io/misha/reference/gextract.matrix.md)
in chunks of about 1e7 values, so the extraction itself adds little to
that. For 100,000 windows of 2 kb, 'x' takes 800 MB; 'y' of three tracks
at 1 bp resolution takes 4.8 GB, or 150 MB in 32 bp bins.

## See also

[`gseq.onehot`](https://tanaylab.github.io/misha/reference/gseq.onehot.md),
[`gextract.matrix`](https://tanaylab.github.io/misha/reference/gextract.matrix.md),
[`gseq.extract`](https://tanaylab.github.io/misha/reference/gseq.extract.md)

## Examples

``` r

gdb.init_examples()
loci <- data.frame(chrom = "chr1", start = c(5000, 20000, 100), end = c(5001, 20001, 101))
w <- gextract.windows(loci, tracks = c(dense = "dense_track"), width = 500, binsize = 50)
dim(w$x)
#> [1]   2   4 500
dim(w$y)
#> [1]  2  1 10
w$dropped # the third locus is too close to the chromosome start
#>   index    reason
#> 1     3 off_chrom
```
