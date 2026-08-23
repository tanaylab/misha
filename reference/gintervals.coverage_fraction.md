# Calculate fraction of genomic space covered by intervals

Returns the fraction of a genomic space that is covered by a set of
intervals.

## Usage

``` r
gintervals.coverage_fraction(intervals1 = NULL, intervals2 = NULL)
```

## Arguments

- intervals1:

  set of one-dimensional intervals (the covering set)

- intervals2:

  set of one-dimensional intervals to be covered (default: NULL, meaning
  the entire genome)

## Value

A single numeric value between 0 and 1 representing the fraction of
'intervals2' (or the genome) covered by 'intervals1'.

## Details

This function calculates what fraction of 'intervals2' is covered by
'intervals1'. If 'intervals2' is NULL, it calculates the fraction of the
entire genome that is covered by 'intervals1'. Overlapping intervals in
either set are automatically unified before calculation.

## See also

[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.intersect`](https://tanaylab.github.io/misha/reference/gintervals.intersect.md),
[`gintervals.covered_bp`](https://tanaylab.github.io/misha/reference/gintervals.covered_bp.md),
[`gintervals.all`](https://tanaylab.github.io/misha/reference/gintervals.all.md)

## Examples

``` r

gdb.init_examples()

# Create some intervals
intervs1 <- gscreen("dense_track > 0.15")
intervs2 <- gintervals(c("chr1", "chr2"), 0, c(100000, 100000))

# Calculate fraction of intervs2 covered by intervs1
gintervals.coverage_fraction(intervs1, intervs2)
#> [1] 0.22325

# Calculate fraction of entire genome covered by intervs1
gintervals.coverage_fraction(intervs1)
#> [1] 0.1391
```
