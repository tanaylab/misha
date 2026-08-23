# Generate random genome intervals

Generate random genome intervals with a specified number of regions of a
specified size. This function samples intervals uniformly across the
genome, weighted by chromosome length.

## Usage

``` r
gintervals.random(
  size,
  n,
  dist_from_edge = 3000000,
  chromosomes = NULL,
  filter = NULL
)
```

## Arguments

- size:

  The size of the intervals to generate (in base pairs)

- n:

  The number of intervals to generate

- dist_from_edge:

  The minimum distance from the edge of the chromosome for a region to
  start or end (default: 3e6)

- chromosomes:

  The chromosomes to sample from (default: all chromosomes). Can be a
  character vector of chromosome names.

- filter:

  A set of intervals to exclude from sampling (default: NULL). Generated
  intervals will not overlap with these regions.

## Value

A data.frame with columns chrom, start, and end representing genomic
intervals

## Details

The function samples intervals randomly across the genome, with
chromosomes weighted by their length. Each interval is guaranteed to:

- Be of the specified size

- Start and end at least `dist_from_edge` bases away from chromosome
  boundaries

- Fall entirely within a single chromosome

- Not overlap with any intervals in the `filter` (if provided)

When a filter is provided, the function pre-computes valid genome
segments (regions not in the filter) and samples from these segments.
Note that this can be slow when the filter contains many intervals.

The function uses R's random number generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) can be used for
reproducibility.

This function is implemented in C++ for high performance and can
generate millions of intervals quickly.

## Examples

``` r
if (FALSE) { # \dontrun{
gdb.init_examples()

# Generate 1000 random intervals of 100bp
intervals <- gintervals.random(100, 1000)
head(intervals)

# Generate intervals only on chr1 and chr2
intervals <- gintervals.random(100, 1000, chromosomes = c("chr1", "chr2"))

# Generate intervals avoiding specific regions
filter_regions <- gintervals(c("chr1", "chr2"), c(1000, 5000), c(2000, 6000))
intervals <- gintervals.random(100, 1000, filter = filter_regions)

# Verify no overlaps with filter
overlaps <- gintervals.intersect(intervals, filter_regions)
nrow(overlaps) # Should be 0

# For reproducibility
set.seed(123)
intervals1 <- gintervals.random(100, 100)
set.seed(123)
intervals2 <- gintervals.random(100, 100)
identical(intervals1, intervals2) # TRUE
} # }
```
