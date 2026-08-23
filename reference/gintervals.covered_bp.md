# Calculate total base pairs covered by intervals

Returns the total number of base pairs covered by a set of intervals.

## Usage

``` r
gintervals.covered_bp(intervals = NULL)
```

## Arguments

- intervals:

  set of one-dimensional intervals

## Value

A single numeric value representing the total number of base pairs
covered by the intervals.

## Details

This function first canonicalizes the intervals to remove overlaps and
touching intervals, then sums up the lengths of all resulting intervals.
Overlapping intervals are counted only once.

## See also

[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.canonic`](https://tanaylab.github.io/misha/reference/gintervals.canonic.md),
[`gintervals.coverage_fraction`](https://tanaylab.github.io/misha/reference/gintervals.coverage_fraction.md)

## Examples

``` r

gdb.init_examples()

# Create some intervals
intervs <- gintervals(
    c("chr1", "chr1", "chr2"),
    c(100, 150, 1000),
    c(200, 250, 2000)
)

# Calculate total bp covered
# Note: intervals [100,200) and [150,250) overlap,
# so total is (200-100) + (250-150) + (2000-1000) = 100 + 100 + 1000 = 1200
# But after canonicalization: [100,250) + [1000,2000) = 150 + 1000 = 1150
gintervals.covered_bp(intervs)
#> [1] 1150
```
