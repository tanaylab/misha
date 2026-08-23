# Install interval sets onto an existing groot

Given an existing groot and a source recipe (or registry name, or
accession), fetches the relevant annotation files and installs interval
sets - one or more of `genes / rmsk / cgi / cytoband`.

## Usage

``` r
gdb.install_intervals(
  groot,
  source,
  sets = c("genes", "rmsk", "cgi", "cytoband"),
  prefix = "",
  gene_sets = c(tss = "tss", exons = "exons", utr3 = "utr3", utr5 = "utr5"),
  gtf_priority = c("ncbiRefSeq", "bestRefSeq", "ensGene", "augustus", "xenoRefGene"),
  overwrite = FALSE,
  registry = NULL,
  target_chroms = NULL,
  target_lengths = NULL,
  min_coverage = 1,
  match_by_length = TRUE,
  force = FALSE,
  verbose = TRUE,
  prefetched_alias = NULL,
  .from_build_genome = FALSE
)
```

## Arguments

- groot:

  Path to a misha groot. `NULL` uses the active groot.

- source:

  Either a registry name, a recipe `list`, or a bare
  `GC[FA]_<digits>.<digits>` accession.

- sets:

  Subset of `c("genes", "rmsk", "cgi", "cytoband")`.

- prefix:

  Character scalar prepended verbatim to each set name. Include the
  trailing dot if you want one (e.g. `"intervs.global."`).

- gene_sets:

  Named character vector mapping `c("tss", "exons", "utr3", "utr5")` to
  the on-disk set name. `NA` value skips that role.

- gtf_priority:

  Character vector ordering GTF source preference for sources that ship
  multiple GTFs (currently `ucsc-hub`). First found wins.

- overwrite:

  If `FALSE` (default), error on existing target sets. If `TRUE`, remove
  existing sets before saving.

- registry:

  Optional path to a registry YAML; overrides the resolution chain.

- target_chroms:

  Optional character vector to pin the canonical column to (e.g. chrom
  names from `halStats --sequenceStats` for a HAL you intend to liftover
  against). When `NULL` (default) misha uses the groot's own chrom names
  (i.e. picks the alias column matching whatever is currently in the
  database). When supplied, misha picks the alias column matching
  `target_chroms` instead, and switches detection to count-weighted
  coverage (bp weighting requires lengths, which target chrom lists
  typically don't carry).

- target_lengths:

  Optional numeric vector aligned with `target_chroms`. Only honored
  when this call originates from
  [`gdb.build_genome`](https://tanaylab.github.io/misha/reference/gdb.build_genome.md)
  (the groot was just force-aligned to `target_chroms`); standalone
  calls ignore it and use the strict column gate. When honored,
  canonical is set to a synthetic `".target_chroms"` column populated by
  name match + unique-length pairing, so `chrom_aliases.tsv` writes
  `target_chroms` as canonical with all other chromAlias columns as
  aliases.

- min_coverage:

  Minimum fraction that must be covered by a chromAlias column for it to
  be picked as the canonical mapping. Default `1.0` (strict). On the
  groot side this is bp-weighted (fraction of genome basepairs
  covered) - a long-tail of small unmapped contigs (e.g. a 16 kb
  mitochondrion missing from UCSC's `genbank` column out of a 3 Gb
  genome) costs ~0.0005 (asset chroms read from a GTF/GFF) the metric is
  the count-weighted fraction of distinct names. Unmapped contigs
  receive no annotations.

- match_by_length:

  If `TRUE` (default), complement the column-based canonical detection
  with a per-row length-based fill: alias rows whose chosen column is
  empty are paired with a groot chrom of the same length, but only when
  the length is unique on both sides (ambiguous lengths are skipped,
  never guessed). Asset translation also switches to a per-row
  cross-column lookup, so a GFF in any naming scheme (or mixed schemes)
  imports cleanly. Currently honored only by the `ucsc-hub` backend
  (which ships per-contig lengths in `<acc>.chrom.sizes.txt`); other
  backends are unaffected. Set `FALSE` for the stricter
  single-column-only behavior.

- force:

  If `FALSE` (default), any requested set that the source doesn't
  provide raises an error and aborts the call before touching the groot.
  If `TRUE`, missing sets are demoted to a single summary warning and
  the available sets are installed.

- verbose:

  If `TRUE`, prints progress.

- prefetched_alias:

  Optional pre-fetched chromAlias bundle (the return value of
  `.hub_preflight_coverage`). When supplied the ucsc-hub fetcher reuses
  it instead of re-downloading. Internal; users never set this directly.

- .from_build_genome:

  Internal flag; when `TRUE`, `gdb.build_genome` signals that it has
  already rescanned the groot and we can skip the entry rescan. Users
  never set this directly.

## Value

Invisible `NULL`. Side effects: writes `.interv` files under
`<groot>/tracks/`, extends `<groot>/chrom_aliases.tsv`, appends to
`<groot>/genome_info.yaml`, and re-initializes the active groot.

The gene-derived sets (`tss`, `exons`, `utr3`, `utr5`) carry a `name`
column (transcript/RNA accession) and a `geneName` column (gene symbol
from the source annotation; blank when the source has no symbol).

## Details

Decoupled from
[`gdb.build_genome`](https://tanaylab.github.io/misha/reference/gdb.build_genome.md)
so that:

- users with a private FASTA build can layer canonical annotations onto
  it;

- failed installs can be resumed without re-fetching the FASTA;

- the same groot can host annotations from multiple sources under
  different prefixes (e.g. `intervs.global.`, `intervs.repeats.`).

## See also

[`gdb.build_genome`](https://tanaylab.github.io/misha/reference/gdb.build_genome.md),
[`gdb.install_gtf_converter`](https://tanaylab.github.io/misha/reference/gdb.install_gtf_converter.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Standalone install on an existing groot.
gdb.install_intervals(
    groot  = "/genomes/arctic_fox",
    source = "GCA_004023825.1",
    prefix = "intervs.global."
)

# Layered: private FASTA groot + intervals from a UCSC hub assembly.
gdb.install_intervals(
    groot  = "/genomes/my_private",
    source = list(source = "ucsc-hub", accession = "GCF_009806435.1"),
    sets   = c("genes", "rmsk")
)
} # }
```
