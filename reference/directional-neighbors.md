# Directional neighbor finding functions

These functions find neighbors using query strand directionality, where
upstream/downstream directionality is determined by the strand of the
query intervals rather than the target intervals. This is particularly
useful for TSS analysis where you want distances relative to gene
direction.

## Usage

``` r
gintervals.neighbors.upstream(
  query_intervals,
  target_intervals,
  maxneighbors = 1,
  maxdist = 1e+09,
  ...
)

gintervals.neighbors.downstream(
  query_intervals,
  target_intervals,
  maxneighbors = 1,
  maxdist = 1e+09,
  ...
)

gintervals.neighbors.directional(
  query_intervals,
  target_intervals,
  maxneighbors_upstream = 1,
  maxneighbors_downstream = 1,
  maxdist = 1e+09,
  ...
)
```

## Arguments

- query_intervals:

  intervals with strand information (query intervals)

- target_intervals:

  intervals to search for neighbors

- maxneighbors:

  maximum number of neighbors per query interval (default: 1)

- maxdist:

  maximum distance to search (default: 1e+09)

- ...:

  additional arguments passed to `gintervals.neighbors`

- maxneighbors_upstream:

  maximum upstream neighbors per query interval (default: 1)

- maxneighbors_downstream:

  maximum downstream neighbors per query interval (default: 1)

## Value

- gintervals.neighbors.upstream:

  data frame of upstream neighbors

- gintervals.neighbors.downstream:

  data frame of downstream neighbors

- gintervals.neighbors.directional:

  list with 'upstream' and 'downstream' components

## Details

\*\*Distance interpretation:\*\*

- \*\*Positive strand queries:\*\* upstream distances \< 0, downstream
  distances \> 0

- \*\*Negative strand queries:\*\* upstream distances \> 0, downstream
  distances \< 0

If no strand column is present, all intervals are treated as positive
strand.

## See also

[`gintervals.neighbors`](https://tanaylab.github.io/misha/reference/gintervals.neighbors.md)

## Examples

``` r

gdb.init_examples()

# Create TSS intervals with strand information
tss <- data.frame(
    chrom = c("chr1", "chr1", "chr1"),
    start = c(1000, 2000, 3000),
    end = c(1001, 2001, 3001),
    strand = c(1, -1, 1), # +, -, +
    gene = c("GeneA", "GeneB", "GeneC")
)

# Create regulatory features
features <- data.frame(
    chrom = "chr1",
    start = c(500, 800, 1200, 1800, 2200, 2800, 3200),
    end = c(600, 900, 1300, 1900, 2300, 2900, 3300),
    feature_id = paste0("F", 1:7)
)

# Find upstream neighbors (promoter analysis)
upstream <- gintervals.neighbors.upstream(tss, features,
    maxneighbors = 2, maxdist = 1000
)
print(upstream)
#>   chrom start  end strand  gene chrom1 start1 end1 feature_id dist
#> 1  chr1  1000 1001      1 GeneA   chr1    800  900         F2 -100
#> 2  chr1  1000 1001      1 GeneA   chr1    500  600         F1 -400
#> 3  chr1  2000 2001     -1 GeneB   chr1   2200 2300         F5 -199
#> 4  chr1  2000 2001     -1 GeneB   chr1   2800 2900         F6 -799
#> 5  chr1  3000 3001      1 GeneC   chr1   2800 2900         F6 -100
#> 6  chr1  3000 3001      1 GeneC   chr1   2200 2300         F5 -700

# Find downstream neighbors (gene body analysis)
downstream <- gintervals.neighbors.downstream(tss, features,
    maxneighbors = 2, maxdist = 5000
)
print(downstream)
#>   chrom start  end strand  gene chrom1 start1 end1 feature_id dist
#> 1  chr1  1000 1001      1 GeneA   chr1   1200 1300         F3  199
#> 2  chr1  1000 1001      1 GeneA   chr1   1800 1900         F4  799
#> 3  chr1  2000 2001     -1 GeneB   chr1   1800 1900         F4  100
#> 4  chr1  2000 2001     -1 GeneB   chr1   1200 1300         F3  700
#> 5  chr1  3000 3001      1 GeneC   chr1   3200 3300         F7  199

# Find both directions in one call
both <- gintervals.neighbors.directional(tss, features,
    maxneighbors_upstream = 1,
    maxneighbors_downstream = 1,
    maxdist = 1000
)
print(both$upstream)
#>   chrom start  end strand  gene chrom1 start1 end1 feature_id dist
#> 1  chr1  1000 1001      1 GeneA   chr1    800  900         F2 -100
#> 2  chr1  2000 2001     -1 GeneB   chr1   2200 2300         F5 -199
#> 3  chr1  3000 3001      1 GeneC   chr1   2800 2900         F6 -100
print(both$downstream)
#>   chrom start  end strand  gene chrom1 start1 end1 feature_id dist
#> 1  chr1  1000 1001      1 GeneA   chr1   1200 1300         F3  199
#> 2  chr1  2000 2001     -1 GeneB   chr1   1800 1900         F4  100
#> 3  chr1  3000 3001      1 GeneC   chr1   3200 3300         F7  199
```
