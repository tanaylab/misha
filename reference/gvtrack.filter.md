# Attach or clear a genomic mask filter on a virtual track

Attaches or clears a genomic mask filter on a virtual track. When a
filter is attached, the virtual track function is evaluated only over
the unmasked regions (i.e., regions not covered by the filter
intervals).

## Usage

``` r
gvtrack.filter(vtrack = NULL, filter = NULL)
```

## Arguments

- vtrack:

  virtual track name

- filter:

  genomic mask to apply. Can be:

  - A data.frame with columns 'chrom', 'start', 'end' (intervals to
    mask)

  - A character string naming an intervals set

  - A character string naming a track (must be intervals-type track)

  - A list of any combination of the above (all will be unified)

  - NULL to clear the filter

## Value

None (invisibly).

## Details

The filter defines regions to *exclude* from virtual track evaluation.
The virtual track function will be evaluated only on the complement of
the filter. Once a filter is attached to a virtual track, it applies to
**all subsequent extractions** of that virtual track until explicitly
cleared with `filter = NULL`.

**Order of Operations:**

Filters are applied *after* iterator modifiers (sshift/eshift/dim). The
order is:

1.  Apply iterator modifiers (gvtrack.iterator with sshift/eshift)

2.  Subtract mask from the modified intervals

3.  Evaluate virtual track function over unmasked regions

**Semantics by function type:**

- **Aggregations (avg/sum/min/max/stddev/quantile):** Length-weighted
  over unmasked regions

- **coverage:** Returns (covered length in unmasked region) / (total
  unmasked length)

- **distance/distance.center:** Unaffected by mask (pure geometry)

- **PWM/kmer:** Masked bases act as hard boundaries; matches cannot span
  masked regions. **Important:** When `extend=TRUE` (the default),
  motifs at the boundaries of unmasked segments can use bases from the
  adjacent masked regions to complete the motif scoring. For example, if
  a 4bp motif starts at position 1998 in an unmasked region that ends at
  2000, and positions 2000-2002 are masked, the motif will still be
  scored using the masked bases. In other words, motif matches *starting
  positions* must be in unmasked regions, but the motif sequence itself
  can extend into masked regions when `extend=TRUE`. Set `extend=FALSE`
  to prevent any use of masked bases in scoring.

**Completely Masked Intervals:** If an entire iterator interval is
masked, the function returns `NA` (not 0).

## See also

[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md),
[`gvtrack.iterator`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.md),
[`gvtrack.info`](https://tanaylab.github.io/misha/reference/gvtrack.info.md)

## Examples

``` r

gdb.init_examples()

## Basic usage: Excluding specific regions
gvtrack.create("vtrack1", "dense_track", func = "avg")

# Create intervals to mask (e.g., repetitive regions)
repeats <- gintervals(c(1, 1), c(100, 500), c(200, 600))

# Attach filter - track will be evaluated excluding these regions
gvtrack.filter("vtrack1", filter = repeats)

# Extract values - masked regions are excluded from calculation
result_filtered <- gextract("vtrack1", gintervals(1, 0, 1000))

# Check filter info
gvtrack.info("vtrack1")
#> $src
#> [1] "dense_track"
#> 
#> $func
#> [1] "avg"
#> 
#> $filter
#> [1] "filter__tmp_RtmpCsKjkH_trackdb_test_tracks_2_7d54e8370e6eb979"
#> 
#> $filter_stats
#> $filter_stats$num_chroms
#> [1] 1
#> 
#> $filter_stats$total_bases
#> [1] 200
#> 
#> $filter_stats$empty
#> [1] FALSE
#> 
#> 

# Clear the filter and compare
gvtrack.filter("vtrack1", filter = NULL)
result_unfiltered <- gextract("vtrack1", gintervals(1, 0, 1000))

## Using multiple filter sources (combined automatically)
centromeres <- gintervals(1, 10000, 15000)
telomeres <- gintervals(1, 0, 1000)
combined_mask <- list(repeats, centromeres, telomeres)

gvtrack.filter("vtrack1", filter = combined_mask)
result_multi_filter <- gextract("vtrack1", gintervals(1, 0, 20000))

## Filters work with iterator modifiers
gvtrack.create("vtrack2", "dense_track", func = "sum")
gvtrack.filter("vtrack2", filter = repeats)
gvtrack.iterator("vtrack2", sshift = -50, eshift = 50)

# Iterator shifts applied first, then mask subtracted
result_shifted <- gextract("vtrack2", gintervals(1, 1000, 2000), iterator = 100)
```
