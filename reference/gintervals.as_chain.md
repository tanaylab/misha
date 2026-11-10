# Transforms existing intervals to a chain format

Transforms existing intervals to a chain format by validating required
columns and adding chain attributes.

## Usage

``` r
gintervals.as_chain(
  intervals = NULL,
  src_overlap_policy = "error",
  tgt_overlap_policy = "auto",
  min_score = NULL
)
```

## Arguments

- intervals:

  a data frame with chain columns: chrom, start, end, strand, chromsrc,
  startsrc, endsrc, strandsrc, chain_id, score. Coordinates are on the
  forward strand of each genome. Strands are 1 (forward) or -1
  (reverse); a strand of 0 is treated as forward. Liftover uses only the
  relative orientation: a block is reversed when exactly one of 'strand'
  and 'strandsrc' is -1.

- src_overlap_policy:

  source overlap policy: "error", "keep", "discard", or
  "discard_intervals"

- tgt_overlap_policy:

  target overlap policy: "error", "auto", "auto_first", "auto_longer",
  "auto_score", "discard", "discard_intervals", "keep", or "agg"

- min_score:

  optional minimum alignment score threshold

## Value

A data frame in chain format with chain attributes set

## Details

This function checks that the input intervals data frame has all the
required columns for a chain format and adds the necessary attributes. A
chain format requires both target coordinates (chrom, start, end,
strand) and source coordinates (chromsrc, startsrc, endsrc, strandsrc),
as well as chain_id and score columns.

## See also

[`gintervals.load_chain`](https://tanaylab.github.io/misha/reference/gintervals.load_chain.md),
[`gintervals.liftover`](https://tanaylab.github.io/misha/reference/gintervals.liftover.md)

## Examples

``` r

gdb.init_examples()

# Create a chain from existing intervals
chain_data <- data.frame(
    chrom = "chr1",
    start = 1000,
    end = 2000,
    strand = 1,
    chromsrc = "chr1",
    startsrc = 5000,
    endsrc = 6000,
    strandsrc = 1,
    chain_id = 1L,
    score = 1000.0
)
chain <- gintervals.as_chain(chain_data)
```
