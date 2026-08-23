# Converts intervals from another assembly

Converts intervals from another assembly to the current one.

## Usage

``` r
gintervals.liftover(
  intervals = NULL,
  chain = NULL,
  src_overlap_policy = "error",
  tgt_overlap_policy = "auto",
  min_score = NULL,
  include_metadata = FALSE,
  canonic = FALSE,
  value_col = NULL,
  multi_target_agg = c("mean", "median", "sum", "min", "max", "count", "first", "last",
    "nth", "max.coverage_len", "min.coverage_len", "max.coverage_frac",
    "min.coverage_frac"),
  params = NULL,
  na.rm = TRUE,
  min_n = NULL
)
```

## Arguments

- intervals:

  intervals from another assembly

- chain:

  name of chain file or data frame as returned by
  'gintervals.load_chain'

- src_overlap_policy:

  policy for handling source overlaps: "error" (default), "keep", or
  "discard". "keep" allows one source interval to map to multiple target
  intervals, "discard" discards all source intervals that have overlaps
  and "error" throws an error if source overlaps are detected.

- tgt_overlap_policy:

  policy for handling target overlaps. One of:

  |  |  |
  |----|----|
  | Policy | Description |
  | error | Throws an error if any target overlaps are detected. |
  | auto | Default. Alias for "auto_score". |
  | auto_score | Resolves overlaps by segmenting the target region and selecting the best chain for each segment based on alignment score (highest score wins). Tie-breakers: longest span, then lowest chain_id. |
  | auto_longer | Resolves overlaps by segmenting and selecting the chain with the longest span for each segment. Tie-breakers: highest score, then lowest chain_id. |
  | auto_first | Resolves overlaps by segmenting and selecting the chain with the lowest chain_id for each segment. |
  | keep | Preserves all overlapping intervals. |
  | discard | Discards any chain interval that has a target overlap with another chain interval. |
  | agg | Segments overlaps into smaller disjoint regions where each region contains all contributing chains, allowing downstream aggregation to process multiple values per region. |
  | best_source_cluster | Best source cluster strategy based on source overlap. When multiple chains map a source interval, clusters them by source overlap: if chain source intervals overlap (indicating true duplications), all mappings are retained; if chain source intervals are disjoint (indicating conflicting/alternative mappings), only the cluster with the largest total target length is kept. |

- min_score:

  optional minimum alignment score threshold. Chains with scores below
  this value are filtered out. Useful for excluding low-quality
  alignments.

- include_metadata:

  logical; if TRUE, adds 'score' column to the output indicating the
  alignment score of the chain used for each mapping. Only applicable
  with "auto_score" or "auto" policy.

- canonic:

  logical; if TRUE, merges adjacent target intervals that originated
  from the same source interval (same intervalID) and same chain (same
  chain_id). This is useful when a source interval maps to multiple
  adjacent target blocks due to chain gaps.

- value_col:

  optional character string specifying the name of a numeric column in
  the intervals data frame to track through the liftover. When
  specified, this column's values are preserved in the output with the
  same column name. Use with multi_target_agg to aggregate values when
  multiple source intervals map to overlapping target regions.

- multi_target_agg:

  aggregation method to use when value_col is specified. One of: "mean",
  "median", "sum", "min", "max", "count", "first", "last", "nth",
  "max.coverage_len", "min.coverage_len", "max.coverage_frac",
  "min.coverage_frac". Default: "mean". Ignored when value_col is NULL.

- params:

  additional parameters for specific aggregation methods. Currently only
  used for "nth" aggregation, where it specifies which element to select
  (e.g., params = 2 for second element, or params = list(n = 2)).

- na.rm:

  logical; if TRUE (default), NA values are removed before aggregation.
  If FALSE, any NA in the values will cause the result to be NA. Only
  used when value_col is specified.

- min_n:

  optional minimum number of non-NA observations required for
  aggregation. If fewer observations are available, the result is NA.
  NULL (default) means no minimum. Only used when value_col is
  specified.

## Value

A data frame representing the converted intervals. For 1D intervals,
always includes 'intervalID' (index of original interval) and 'chain_id'
(identifier of the chain that produced the mapping) columns. The
chain_id column is essential for distinguishing results when a source
interval maps to multiple target regions via different chains
(duplications). When include_metadata=TRUE, also includes 'score'
column. When value_col is specified, includes the value column with its
original name.

## Details

This function converts 'intervals' from another assembly to the current
one. Chain file instructs how the conversion of coordinates should be
done. It can be either a name of a chain file or a data frame in the
same format as returned by 'gintervals.load_chain' function.

The converted intervals are returned. An additional column named
'intervalID' is added to the resulted data frame. For each interval in
the resulted intervals it indicates the index of the original interval.

Note: When passing a pre-loaded chain (data frame), overlap policies
cannot be specified - they are taken from the chain's attributes that
were set during loading. When passing a chain file path, policies can be
specified and will be used for loading.

## See also

[`gintervals.load_chain`](https://tanaylab.github.io/misha/reference/gintervals.load_chain.md),
[`gtrack.liftover`](https://tanaylab.github.io/misha/reference/gtrack.liftover.md),
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md)

## Examples

``` r

gdb.init_examples()
chainfile <- paste(.misha$GROOT, "data/test.chain", sep = "/")
intervs <- data.frame(
    chrom = "chr25", start = c(0, 7000),
    end = c(6000, 20000)
)
# Liftover with default policies
gintervals.liftover(intervs, chainfile)
#>   chrom start   end intervalID chain_id
#> 1  chr1 12000 12500          1        1
#> 2  chr1 12700 13500          1        1
#> 3  chr1 14100 16500          1        1
#> 4  chr1 17500 18500          2        1
#> 5  chrX  5000  7000          2        2

# Liftover keeping source overlaps (one source interval may map to multiple targets)
# gintervals.liftover(intervs, chainfile, src_overlap_policy = "keep")
```
