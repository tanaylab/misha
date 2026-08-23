# Imports a track from another assembly

Imports a track from another assembly.

## Usage

``` r
gtrack.liftover(
  track = NULL,
  description = NULL,
  src.track.dir = NULL,
  chain = NULL,
  src_overlap_policy = "error",
  tgt_overlap_policy = "auto",
  multi_target_agg = c("mean", "median", "sum", "min", "max", "count", "first", "last",
    "nth", "max.coverage_len", "min.coverage_len", "max.coverage_frac",
    "min.coverage_frac"),
  params = NULL,
  na.rm = TRUE,
  min_n = NULL,
  min_score = NULL
)
```

## Arguments

- track:

  name of a created track

- description:

  a character string description

- src.track.dir:

  path to the directory of the source track

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

- multi_target_agg:

  aggregation/selection policy for contributors that land on the same
  target locus. When multiple source intervals map to overlapping
  regions in the target genome (after applying tgt_overlap_policy),
  their values must be combined into a single value.

- params:

  additional parameters for aggregation (e.g., for "nth" aggregation)

- na.rm:

  logical indicating whether NA values should be removed before
  aggregation (default: TRUE)

- min_n:

  minimum number of non-NA values required for aggregation. If fewer
  values are available, the result will be NA.

- min_score:

  optional minimum alignment score threshold. Chains with scores below
  this value are filtered out. Useful for excluding low-quality
  alignments.

## Value

None.

## Details

This function imports a track located in 'src.track.dir' of another
assembly to the current database. Chain file instructs how the
conversion of coordinates should be done. It can be either a name of a
chain file or a data frame in the same format as returned by
'gintervals.load_chain' function. The name of the newly created track is
specified by 'track' argument and 'description' is added as a track
attribute.

Note: When passing a pre-loaded chain (data frame), overlap policies
cannot be specified - they are taken from the chain's attributes that
were set during loading. When passing a chain file path, policies can be
specified and will be used for loading. Aggregation parameters
(multi_target_agg, params, na.rm, min_n) can always be specified
regardless of chain type.

## Note

Terminology note for UCSC chain format users: In the UCSC chain format
specification, the fields prefixed with 't' (tName, tStart, tEnd, etc.)
are called "target" or "reference", while fields prefixed with 'q'
(qName, qStart, qEnd, etc.) are called "query". However, misha uses
reversed terminology: UCSC's "target/reference" corresponds to misha's
"source" (chromsrc, startsrc, endsrc), and UCSC's "query" corresponds to
misha's "target" (chrom, start, end).

## See also

[`gintervals.load_chain`](https://tanaylab.github.io/misha/reference/gintervals.load_chain.md),
[`gintervals.liftover`](https://tanaylab.github.io/misha/reference/gintervals.liftover.md)
