# Returns DNA sequences

Returns DNA sequences for given intervals

## Usage

``` r
gseq.extract(intervals = NULL, pad = FALSE)
```

## Arguments

- intervals:

  intervals for which DNA sequence is returned

- pad:

  if TRUE, the part of an interval outside its chromosome is returned as
  "N" (see Details). If FALSE (the default), an interval that runs past
  a chromosome end is an error.

## Value

An array of character strings representing DNA sequence.

## Details

This function returns an array of sequence strings for each interval
from 'intervals'. If intervals contain an additional 'strand' column and
its value is '-1', the reverse-complementary sequence is returned.

By default ('pad = FALSE') each interval must lie within its chromosome:
an interval that starts before 0 or ends past the end of the chromosome
is an error.

With 'pad = TRUE' an interval may run past either end of its chromosome
('start' below 0, 'end' above the chromosome length). The positions
outside the chromosome come back as "N", so each sequence is 'end -
start' long, and an interval wholly outside its chromosome comes back as
all "N". 'start' must still be below 'end', and the chromosome must
exist. For a row with strand -1 the padded sequence is reverse
complemented as a whole, so the "N" run of the chromosome start comes
last. The bases inside the chromosome are returned as stored, lowercase
included. This gives fixed-width windows near chromosome ends, e.g. for
sequence models (see the examples).

Each interval costs one disk read, so large interval sets are dominated
by I/O latency rather than by the sequence data itself. When the
sequence files are not in the page cache, the work is split across
several processes (capped by 'gmax.processes', and in practice usually
well below it). The decision is made by timing the first few reads and
extrapolating over the number of read-ahead windows the intervals span:
cached data, and interval sets dense enough that many intervals share a
read, stay in a single process because distributing them would only add
fork overhead.

Sets of fewer than 500 intervals are never distributed. Named intervals
sets stored in the "big" format are always read sequentially.
Distributing roughly doubles peak memory while the result is being
assembled, so a call close to 'gmax.data.size' may need a lower limit or
'gmultitasking = FALSE'. Note also that this function forks: calling it
from inside your own 'mclapply()' or similar multiplies the process
count, as it does for the other multitasking misha functions.

Set 'gseq.extract.probe.usec' to override the per-interval microsecond
threshold at which distributing kicks in; 0 always distributes and also
requests the maximum number of processes, a large value never
distributes. 'gmultitasking = FALSE' disables it entirely.

## See also

[`gextract`](https://tanaylab.github.io/misha/reference/gextract.md)

## Examples

``` r

gdb.init_examples()
intervs <- gintervals(c(1, 2), 10000, 10020)
gseq.extract(intervs)
#> [1] "ccaggggccagcactgctcg" "GGGGCTCCCACCCGCCGTCC"

# Fixed-width windows centred on positions. The second window starts 30 bp
# before the start of chr1, so its first 30 bases are N.
w <- 100
centers <- c(250000, 20)
windows <- data.frame(chrom = "chr1", start = centers - w / 2, end = centers + w / 2)
seqs <- gseq.extract(windows, pad = TRUE)
nchar(seqs)
#> [1] 100 100
seqs[2]
#> [1] "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaac"
```
