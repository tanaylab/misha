# Promoter windows around the 5' end of stranded intervals

Returns, for each interval, 'upstream' bp before its 5' end and
'downstream' bp from the 5' base on, following the interval's strand.
The 5' base is 'start' for strand 1 and 'end - 1' for strand -1, so a
window is 'upstream + downstream' bp wide before clipping at the
chromosome ends. On a 1 bp TSS set, such as the 'tss' set that
[`gdb.install_intervals`](https://tanaylab.github.io/misha/reference/gdb.install_intervals.md)
writes, this gives promoters; on transcripts it uses their first base.

## Usage

``` r
gintervals.promoters(intervals = NULL, upstream = 500, downstream = 50)
```

## Arguments

- intervals:

  intervals set with a 'strand' column (1, -1 or 0)

- upstream:

  bp upstream of the 5' end

- downstream:

  bp downstream of the 5' end, the 5' base included

## Value

A data frame of the windows, with the other columns kept, clipped to the
chromosome boundaries.

## Details

Strand 0 marks features on opposite strands merged into one interval, as
[`gintervals.import_genes`](https://tanaylab.github.io/misha/reference/gintervals.import_genes.md)
does for a TSS that two genes share. For such an interval the window
covers both the strand 1 and the strand -1 window.

## See also

[`gintervals.expand`](https://tanaylab.github.io/misha/reference/gintervals.expand.md),
[`gintervals.import_genes`](https://tanaylab.github.io/misha/reference/gintervals.import_genes.md)

## Examples

``` r

gdb.init_examples()
tss <- data.frame(chrom = "chr1", start = c(10000, 20000), end = c(10001, 20001), strand = c(1, -1))
gintervals.promoters(tss, upstream = 500, downstream = 50)
#>   chrom start   end strand
#> 1  chr1  9500 10050      1
#> 2  chr1 19951 20501     -1
```
