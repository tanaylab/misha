# Loads assembly conversion table from a chain file

Loads assembly conversion table from a chain file.

## Usage

``` r
gintervals.load_chain(
  file = NULL,
  src_overlap_policy = "error",
  tgt_overlap_policy = "auto",
  src_groot = NULL,
  min_score = NULL
)
```

## Arguments

- file:

  name of chain file

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

- src_groot:

  optional path to source genome database for validating source
  chromosomes and coordinates. If provided, the function temporarily
  switches to this database to verify that all source chromosomes exist
  and coordinates are within bounds, then restores the original
  database.

- min_score:

  optional minimum alignment score threshold. Chains with scores below
  this value are filtered out. Useful for excluding low-quality
  alignments.

## Value

A data frame representing assembly conversion table with columns: chrom,
start, end, strand, chromsrc, startsrc, endsrc, strandsrc, chain_id,
score.

## Details

This function reads a file in 'chain' format and returns assembly
conversion table that can be used in 'gtrack.liftover' and
'gintervals.liftover'.

Source overlaps occur when the same source genome position maps to
multiple target genome positions. Target overlaps occur when multiple
source positions map to overlapping regions in the target genome.

The 'src_overlap_policy' controls how source overlaps are handled:

- "error" (default): Throw an error if source overlaps are detected

- "keep": Keep all mappings, allowing one source to map to multiple
  targets

- "discard": Remove all chain intervals involved in source overlaps

The 'tgt_overlap_policy' controls how target overlaps are handled:

- "error": Throw an error if target overlaps are detected

- "auto" (default) / "auto_first": Keep the first overlapping chain
  (original file order) by trimming or discarding later overlaps while
  keeping source/target lengths consistent

- "auto_longer": Keep the longer overlapping chain and trim/drop the
  shorter ones

- "discard": Remove all chain intervals involved in target overlaps

- "keep": Allow target overlaps to remain untouched (liftover must be
  able to handle them)

## See also

[`gintervals.liftover`](https://tanaylab.github.io/misha/reference/gintervals.liftover.md),
[`gtrack.liftover`](https://tanaylab.github.io/misha/reference/gtrack.liftover.md)

## Examples

``` r

gdb.init_examples()
chainfile <- paste(.misha$GROOT, "data/test.chain", sep = "/")
# Load chain file with default policies
gintervals.load_chain(chainfile)
#>   chrom start   end strand chromsrc startsrc endsrc strandsrc chain_id score
#> 1  chr1 12000 12500      1    chr25     2000   2500         1        1 2e+05
#> 2  chr1 12700 13500      1    chr25     2500   3300         1        1 2e+05
#> 3  chr1 14100 18500      1    chr25     3600   8000         1        1 2e+05
#> 4  chrX  5000  7000      1    chr25    10000  12000         1        2 2e+05
```
