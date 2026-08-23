# Score DNA sequences with a k-mer over a region of interest

Counts exact matches of a k-mer in DNA sequences over a specified region
of interest (ROI). The ROI is defined by `start_pos` and `end_pos`
(1-based, inclusive), with optional extension controlled by `extend`.

## Usage

``` r
gseq.kmer(
  seqs,
  kmer,
  mode = c("count", "frac"),
  strand = 0L,
  start_pos = NULL,
  end_pos = NULL,
  extend = FALSE,
  skip_gaps = TRUE,
  gap_chars = c("-", ".")
)
```

## Arguments

- seqs:

  character vector of DNA sequences (A/C/G/T/N; case-insensitive)

- kmer:

  single character string containing the k-mer to search for (A/C/G/T
  only)

- mode:

  character; one of "count" or "frac"

- strand:

  integer; 1=forward, -1=reverse, 0=both strands (default: 0)

- start_pos:

  integer or NULL; 1-based inclusive start of ROI (default: 1)

- end_pos:

  integer or NULL; 1-based inclusive end of ROI (default: sequence
  length)

- extend:

  logical or integer; extension of allowed window starts (default:
  FALSE)

- skip_gaps:

  logical; if TRUE, treat gap characters as holes and skip them while
  scanning. Windows are k consecutive non-gap bases (default: TRUE)

- gap_chars:

  character vector; which characters count as gaps (default: c("-",
  "."))

## Value

Numeric vector with counts (for "count" mode) or fractions (for "frac"
mode). Returns 0 when sequence is too short or ROI is invalid.

## Details

This function counts k-mer occurrences in DNA sequences directly without
requiring a genomics database. For detailed documentation on k-mer
counting parameters, see
[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)
(functions "kmer.count" and "kmer.frac").

The ROI (region of interest) is defined by `start_pos` and `end_pos`.
The `extend` parameter controls whether k-mer matches can extend beyond
the ROI boundaries. For palindromic k-mers, use `strand=1` or `-1` to
avoid double counting.

When `skip_gaps=TRUE`, characters specified in `gap_chars` are treated
as gaps. Windows are defined as k consecutive non-gap bases. The `frac`
denominator counts the number of possible logical starts (non-gap
windows) in the region. `start_pos` and `end_pos` are interpreted as
physical coordinates on the full sequence.

## See also

[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)
for detailed k-mer parameter documentation

## Examples

``` r
if (FALSE) { # \dontrun{
# Example sequences
seqs <- c("CGCGCGCGCG", "ATATATATAT", "ACGTACGTACGT")

# Count CG dinucleotides on both strands
gseq.kmer(seqs, "CG", mode = "count", strand = 0)

# Count on forward strand only
gseq.kmer(seqs, "CG", mode = "count", strand = 1)

# Get CG fraction
gseq.kmer(seqs, "CG", mode = "frac", strand = 0)

# Count in a specific region
gseq.kmer(seqs, "CG", mode = "count", start_pos = 2, end_pos = 8)

# Allow k-mer to extend beyond ROI boundaries
gseq.kmer(seqs, "CG", mode = "count", start_pos = 2, end_pos = 8, extend = TRUE)

# Calculate GC content by summing G and C fractions
g_frac <- gseq.kmer(seqs, "G", mode = "frac", strand = 1)
c_frac <- gseq.kmer(seqs, "C", mode = "frac", strand = 1)
gc_content <- g_frac + c_frac
gc_content

# Compare AT counts on different strands
at_forward <- gseq.kmer(seqs, "AT", mode = "count", strand = 1)
at_reverse <- gseq.kmer(seqs, "AT", mode = "count", strand = -1)
at_both <- gseq.kmer(seqs, "AT", mode = "count", strand = 0)
data.frame(forward = at_forward, reverse = at_reverse, both = at_both)
} # }
```
