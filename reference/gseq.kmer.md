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

  integer; 1=forward, -1=reverse, 0=both strands, summed (default: 0)

- start_pos:

  integer or NULL; 1-based inclusive start of ROI (default: 1)

- end_pos:

  integer or NULL; 1-based inclusive end of ROI (default: sequence
  length)

- extend:

  logical or integer; if TRUE, counts k-mers whose anchor lies in the
  ROI even if the k-mer extends beyond its end; if FALSE (default), only
  k-mers fully contained in the ROI are considered. An integer is the
  number of bases a k-mer may run past `end_pos`: TRUE is k-1 and FALSE
  is 0. An integer above k-1 also counts k-mers that start after
  `end_pos` (and `frac` counts their starts).

- skip_gaps:

  logical; if TRUE, treat gap characters as holes and skip them while
  scanning. Windows are k consecutive non-gap bases (default: TRUE)

- gap_chars:

  character vector; which characters count as gaps (default: c("-",
  "."))

## Value

Numeric vector with counts (for "count" mode) or fractions (for "frac"
mode). Returns 0 when sequence is too short or ROI is invalid, and NA
for an `NA` sequence.

## Details

This function counts k-mer occurrences in DNA sequences directly without
requiring a genomics database. For detailed documentation on k-mer
counting parameters, see
[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)
(functions "kmer.count" and "kmer.frac").

The ROI (region of interest) is defined by `start_pos` and `end_pos`. As
in the k-mer virtual tracks, a match on either strand is anchored at its
first (leftmost) base, which must not lie before `start_pos`. `extend`
sets how many bases past `end_pos` a match may run (k-1 with TRUE, so
every anchor in the ROI counts), never past the end of the sequence.

With `strand=0` the matches on the two strands are summed, as in the
k-mer virtual tracks, so a palindromic k-mer (its own reverse
complement) counts twice per site and `strand=0` issues a warning; use
`strand=1` or `-1` to count it once.

For `mode="frac"` the denominator is the number of possible anchors in
the ROI: every anchor with `extend=TRUE` (minus the last k-1 at the end
of the sequence), the first n-k+1 of an n-base ROI with `extend=FALSE`;
with `strand=0` anchors on both strands are counted, so the fraction
stays in \[0, 1\].

When `skip_gaps=TRUE`, characters specified in `gap_chars` are treated
as gaps. Windows are k consecutive non-gap bases, anchored at their
first base, and `extend` and the `frac` denominator count non-gap bases
and windows. `start_pos` and `end_pos` are interpreted as physical
coordinates on the full sequence.

## See also

[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md)
for detailed k-mer parameter documentation

## Examples

``` r
if (FALSE) { # \dontrun{
# Example sequences
seqs <- c("CGCGCGCGCG", "ATATATATAT", "ACGTACGTACGT")

# Count CG dinucleotides on both strands (CG is its own reverse complement, so
# each site counts twice and a warning is given)
gseq.kmer(seqs, "CG", mode = "count", strand = 0)

# Count on forward strand only
gseq.kmer(seqs, "CG", mode = "count", strand = 1)

# Get CG fraction
gseq.kmer(seqs, "CG", mode = "frac", strand = 1)

# Count in a specific region
gseq.kmer(seqs, "CG", mode = "count", strand = 1, start_pos = 2, end_pos = 8)

# Also count k-mers that start in the region and end past it
gseq.kmer(seqs, "CG", mode = "count", strand = 1, start_pos = 2, end_pos = 8, extend = TRUE)

# Calculate GC content by summing G and C fractions
g_frac <- gseq.kmer(seqs, "G", mode = "frac", strand = 1)
c_frac <- gseq.kmer(seqs, "C", mode = "frac", strand = 1)
gc_content <- g_frac + c_frac
gc_content

# Compare AC counts on different strands (AC on the reverse strand is GT on the forward)
ac_forward <- gseq.kmer(seqs, "AC", mode = "count", strand = 1)
ac_reverse <- gseq.kmer(seqs, "AC", mode = "count", strand = -1)
ac_both <- gseq.kmer(seqs, "AC", mode = "count", strand = 0)
data.frame(forward = ac_forward, reverse = ac_reverse, both = ac_both)
} # }
```
