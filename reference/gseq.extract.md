# Returns DNA sequences

Returns DNA sequences for given intervals

## Usage

``` r
gseq.extract(intervals = NULL)
```

## Arguments

- intervals:

  intervals for which DNA sequence is returned

## Value

An array of character strings representing DNA sequence.

## Details

This function returns an array of sequence strings for each interval
from 'intervals'. If intervals contain an additional 'strand' column and
its value is '-1', the reverse-complementary sequence is returned.

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
```
