# Build a misha genome database from a name

Builds a misha genomic database for a named assembly. Resolves the name
through the registry chain (or pattern-fallback for `GC[FA]_*`
accessions), downloads the FASTA, calls
[`gdb.create`](https://tanaylab.github.io/misha/reference/gdb.create.md)
to build the seq-only groot, then dispatches to
[`gdb.install_intervals`](https://tanaylab.github.io/misha/reference/gdb.install_intervals.md)
for the requested sets.

## Usage

``` r
gdb.build_genome(
  name,
  path = name,
  registry = NULL,
  sets = c("genes", "rmsk", "cgi", "cytoband"),
  prefix = "",
  gene_sets = c(tss = "tss", exons = "exons", utr3 = "utr3", utr5 = "utr5"),
  gtf_priority = c("ncbiRefSeq", "bestRefSeq", "ensGene", "augustus", "xenoRefGene"),
  chrom_naming = NULL,
  target_chroms = NULL,
  target_lengths = NULL,
  min_coverage = 1,
  match_by_length = TRUE,
  format = NULL,
  verbose = TRUE
)
```

## Arguments

- name:

  Genome name (registry key, alias, or `GC[FA]_*` accession).

- path:

  Output directory; must not exist.

- registry:

  Optional path to an explicit registry YAML.

- sets:

  Subset of `c("genes", "rmsk", "cgi", "cytoband")`. Empty vector
  `character(0)` = sequence-only build.

- prefix:

  Character scalar prepended to set names (see
  [`gdb.install_intervals`](https://tanaylab.github.io/misha/reference/gdb.install_intervals.md)).

- gene_sets:

  Named character vector mapping the four
  [`gintervals.import_genes()`](https://tanaylab.github.io/misha/reference/gintervals.import_genes.md)
  roles to on-disk set names; `NA` skips a role.

- gtf_priority:

  Character vector ordering GTF source preference.

- chrom_naming:

  Optional override for the recipe's `chrom_naming`. Selects which name
  space the canonical chrom names should come from. For `ucsc-hub`: any
  chromAlias column (`"ucsc"`, `"genbank"`, `"refseq"`, `"ncbi"`), plus
  the friendly aliases `"sequence_name"` (= `"assembly"`) and
  `"accession"` (keep the FASTA's source column). For `ncbi`:
  `"sequence_name"` (default), `"ucsc"`, or `"accession"`. `NULL`
  (default) keeps whatever the recipe specifies.

- target_chroms:

  Optional character vector of chrom names the resulting groot should
  align to (typically the output of `halStats --sequenceStats`, the
  chrom names in a HAL file you intend to liftover against). When
  supplied, misha auto-picks the chromAlias column whose values cover
  `target_chroms` best and uses that column as the canonical naming,
  instead of `chrom_naming`. Honored only by the `ucsc-hub` backend;
  supplying it with any other source is an error (raised before any
  download).

- target_lengths:

  Optional numeric vector aligned with `target_chroms` (typically the
  second field of `halStats --sequenceStats`). When supplied alongside
  `target_chroms` and with `match_by_length = TRUE`, this is the
  strong-guarantee path: misha force-aligns the hub FASTA to
  `target_chroms`, placing every target on its chromAlias row by name
  match across columns or unique-on-both-sides length pairing. If any
  target can't be placed, the build errors (in the pre-flight, before
  the multi-GB FASTA download). On success the resulting groot's chrom
  names are exactly `target_chroms` (alias rows not in `target_chroms`
  keep their original FASTA-header accession). Honored only by the
  `ucsc-hub` backend.

- min_coverage:

  Minimum fraction of groot chroms that must appear in a chromAlias
  column for that column to be picked as canonical (forwarded to
  [`gdb.install_intervals`](https://tanaylab.github.io/misha/reference/gdb.install_intervals.md)).
  Default `1.0` (strict). Lower to e.g. `0.99` when a column has small
  gaps – typical when a target column doesn't span every contig (e.g.
  UCSC's `genbank` column has no value for the mitochondrion in many
  hubs, leaving 1 stray chrom). Honored only by the `ucsc-hub` backend;
  supplying a non-default value for any other source is an error (raised
  before any download).

- match_by_length:

  Forwarded to
  [`gdb.install_intervals`](https://tanaylab.github.io/misha/reference/gdb.install_intervals.md).
  When `TRUE` (default), complements column-based canonical detection
  with a per-row length match for alias rows the chosen column couldn't
  cover, and switches asset translation to a cross-column per-row lookup
  so GFFs in any naming scheme import cleanly. Set `FALSE` for the
  stricter single-column-only behavior.

- format:

  `"indexed"` or `"per-chromosome"`; `NULL` =\>
  `getOption("gmulticontig.indexed_format", TRUE)`.

- verbose:

  If `TRUE`, prints progress.

## Value

None (invisible `NULL`). The installed gene-derived sets (`tss`,
`exons`, `utr3`, `utr5`) carry a `name` column (transcript/RNA
accession) and a `geneName` column (gene symbol from the source
annotation; blank when the source has no symbol).

## Details

For details on resolution, sources, sets, and chromosome-alias handling,
see
[`gdb.install_intervals`](https://tanaylab.github.io/misha/reference/gdb.install_intervals.md).

## See also

[`gdb.install_intervals`](https://tanaylab.github.io/misha/reference/gdb.install_intervals.md),
[`gdb.create`](https://tanaylab.github.io/misha/reference/gdb.create.md),
[`gdb.list_genomes`](https://tanaylab.github.io/misha/reference/gdb.list_genomes.md),
[`gdb.genome_info`](https://tanaylab.github.io/misha/reference/gdb.genome_info.md).

## Examples

``` r
if (FALSE) { # \dontrun{
gdb.build_genome("hg38", path = "~/genomes/hg38")
gdb.build_genome("GCA_004023825.1",
    path   = "~/genomes/arctic_fox",
    prefix = "intervs.global."
)
# Match HAL/Cactus canonical names (GenBank accessions like JH880237.1):
gdb.build_genome("GCF_000298355.1",
    path         = "~/genomes/Bos_mutus",
    chrom_naming = "genbank",
    prefix       = "intervs.global."
)
} # }
```
