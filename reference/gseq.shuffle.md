# Shuffle sequences keeping their k-mer counts

Draws random sequences with exactly the same k-mer counts as the input
sequences - the dinucleotide counts with the default `k = 2` - for use
as a composition-matched background, e.g. in motif enrichment.

## Usage

``` r
gseq.shuffle(seqs, k = 2L, n = 1L, seed = NULL)
```

## Arguments

- seqs:

  a character vector of DNA sequences (`NA` gives `NA`), or an intervals
  data frame, whose sequences are read with
  [`gseq.extract`](https://tanaylab.github.io/misha/reference/gseq.extract.md)
  (strand-aware).

- k:

  the k-mer length whose counts are kept, 1 to 32. Default 2
  (dinucleotides).

- n:

  the number of shuffles of each sequence.

- seed:

  a whole number (\|seed\| \<= 2^53) seeding the shuffle, or `NULL` to
  draw one from R's random number generator.

## Value

With `n = 1`, a character vector like `seqs` (names kept). With `n > 1`,
a character matrix with a row per sequence and a column per shuffle.

## Details

Each shuffle is drawn uniformly at random from all the sequences that
keep the k-mer counts and the first k-1 bases. The construction is the
Eulerian-path one of Altschul & Erickson (1985) and Kandel et al.
(1996), as in uShuffle (Jiang et al. 2008): the (k-1)-mers are the
vertices of a graph whose edges are the k-mers of the sequence, a
uniformly random spanning tree of "last exits" is drawn with Wilson's
algorithm, and the remaining edges out of each vertex are put in a
random order. `k = 1` is a plain random permutation of the bases.
Because the k-mers form a path through the (k-1)-mers, the last k-1
bases are kept too.

Only A, C, G and T (in either case) are shuffled. Any other character -
N, an IUPAC code, a gap - stays at its position and splits the sequence
into segments, each shuffled on its own; the k-mers that touch such a
character are therefore unchanged, and the k-mer counts of the whole
sequence, N-containing k-mers included, are kept. Shuffling the N's into
the sequence instead would scatter assembly gaps across it. Bases are
compared ignoring case and every position keeps its case, so a
soft-masked stretch stays lowercase at the same coordinates (with new
bases in it); [`toupper()`](https://rdrr.io/r/base/chartr.html) the
input first to ignore case altogether. A segment of k bases or fewer has
only one arrangement and is returned as is.

The random numbers come from misha's own generator (xoshiro256\*\*), not
from R's: with a `seed` the result does not depend on, and does not
change, R's random number state, and is the same on every platform and
the same as PyMisha's `gseq_shuffle` for the same seed. Without a `seed`
one is drawn from R's generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) makes the call
reproducible (but the sequences then differ from PyMisha's, whose
generator differs). Each input sequence gets its own random stream,
derived from the seed and its position in `seqs`, so its shuffles do not
depend on the sequences that come after it.

## References

Altschul SF, Erickson BW (1985). Significance of nucleotide sequence
alignments: a method for random sequence permutation that preserves
dinucleotide and codon usage. *Mol Biol Evol* 2(6):526-538.

Kandel D, Matias Y, Unger R, Winkler P (1996). Shuffling biological
sequences. *Discrete Applied Mathematics* 71:171-185.

Jiang M, Anderson J, Gillespie J, Mayne M (2008). uShuffle: a useful
tool for shuffling biological sequences while preserving the k-let
counts. *BMC Bioinformatics* 9:192.

## See also

[`gseq.kmer`](https://tanaylab.github.io/misha/reference/gseq.kmer.md),
[`gseq.kmer.dist`](https://tanaylab.github.io/misha/reference/gseq.kmer.dist.md),
[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md)
for sampling from a Markov model, which keeps the expected composition
rather than the exact counts.

## Examples

``` r
seqs <- c(a = "ACGTTGCAACGTAGCTAGCTTTACG", b = "ccgNNacgtacgtTTGACG")
gseq.shuffle(seqs, seed = 1)
#>                           a                           b 
#> "AACTAGTTTGCGCACGTTAGCTACG"       "ccgNNacgacgtgTTTACG" 

# three shuffles of each sequence, keeping trinucleotide counts
gseq.shuffle(seqs, k = 3, n = 3, seed = 1)
#>   [,1]                        [,2]                       
#> a "ACGTAGCTAGCTTGCAACGTTTACG" "ACGTTTGCAACGTAGCTAGCTTACG"
#> b "ccgNNacgtacgtTTGACG"       "ccgNNacgtacgtTTGACG"      
#>   [,3]                       
#> a "ACGTTAGCTACGTAGCTTTGCAACG"
#> b "ccgNNacgtttgaCGTACG"      

# the dinucleotide counts are unchanged
sh <- gseq.shuffle(seqs[["a"]], seed = 7)
identical(
    table(substring(sh, 1:24, 2:25)),
    table(substring(seqs[["a"]], 1:24, 2:25))
)
#> [1] TRUE
```
