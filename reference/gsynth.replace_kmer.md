# Iteratively replace a k-mer in the genome

Performs an iterative replacement of a `target` k-mer with a
`replacement` sequence. This is useful for creating synthetic genomes
with specific motifs removed (e.g., creating a CpG-null genome by
iteratively swapping CG to GC).

## Usage

``` r
gsynth.replace_kmer(
  target,
  replacement,
  output_path = NULL,
  output_format = c("misha", "fasta", "vector"),
  intervals = NULL,
  check_composition = TRUE
)
```

## Arguments

- target:

  The k-mer sequence to remove (e.g., "CG").

- replacement:

  The replacement sequence (e.g., "GC").

- output_path:

  Path to the output file (ignored when output_format = "vector").

- output_format:

  Output format:

  - "misha": .seq binary format (default)

  - "fasta": FASTA text format

  - "vector": Return sequences as a character vector (does not write to
    file)

- intervals:

  Genomic intervals to process. If NULL, uses all chromosomes.

- check_composition:

  Logical. If TRUE (default), ensures target and replacement have the
  same nucleotide composition (preserving exact base counts).

## Value

When output_format is "misha" or "fasta", returns invisible NULL and
writes to output_path. When output_format is "vector", returns a
character vector of modified sequences.

## Details

**Bubble Sort / Iterative Logic:** The function scans the sequence and
replaces occurrences of `target` with `replacement`. If a replacement
creates a new instance of `target` (e.g., removing "CG" with "GC" in the
sequence "CCG" -\> "CGC"), the new instance is also replaced. This
continues until the sequence is free of the `target` k-mer.

When `target` and `replacement` are permutations of each other (e.g.,
"CG" and "GC"), this acts as a "bubble sort" of nucleotides, moving
bases locally without altering the total GC content or base counts of
the genome.

## Examples

``` r
if (FALSE) { # \dontrun{
# Robust removal of all CpG dinucleotides (preserving GC%)
gsynth.replace_kmer(
    target = "CG",
    replacement = "GC",
    output_path = "genome_no_cpg.seq",
    output_format = "misha"
)
} # }
```
