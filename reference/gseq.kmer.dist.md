# Compute k-mer distribution in genomic intervals

Counts the occurrence of all k-mers (of size k) within the specified
genomic intervals, optionally excluding masked regions.

## Usage

``` r
gseq.kmer.dist(intervals, k = 6L, mask = NULL)
```

## Arguments

- intervals:

  Genomic intervals to analyze

- k:

  Integer k-mer size (1-10). Default is 6.

- mask:

  Optional intervals to exclude from counting. Positions within the mask
  will not contribute to k-mer counts.

## Value

A data frame with columns:

- kmer:

  Character string representing the k-mer sequence

- count:

  Number of occurrences of this k-mer

Only k-mers with count \> 0 are included. K-mers containing N bases are
not counted.

## See also

[`gseq.extract`](https://tanaylab.github.io/misha/reference/gseq.extract.md),
[`gseq.kmer`](https://tanaylab.github.io/misha/reference/gseq.kmer.md)

## Examples

``` r
gdb.init_examples()

# Count all 6-mers in first 10kb of chr1
intervals <- data.frame(chrom = "chr1", start = 0, end = 10000)
kmer_dist <- gseq.kmer.dist(intervals, k = 6)
head(kmer_dist)
#>     kmer count
#> 1 AAAAAA     3
#> 2 AAAAAG     2
#> 3 AAAAAT     4
#> 4 AAAACA     1
#> 5 AAAACC     2
#> 6 AAAAGA     1

# Count dinucleotides
dinucs <- gseq.kmer.dist(intervals, k = 2)
dinucs
#>    kmer count
#> 1    AA   479
#> 2    AC   519
#> 3    AG   801
#> 4    AT   292
#> 5    CA   760
#> 6    CC  1087
#> 7    CG   277
#> 8    CT   881
#> 9    GA   574
#> 10   GC   835
#> 11   GG   982
#> 12   GT   475
#> 13   TA   278
#> 14   TC   565
#> 15   TG   806
#> 16   TT   388

# Count with mask
mask <- data.frame(chrom = "chr1", start = 5000, end = 6000)
kmer_dist_masked <- gseq.kmer.dist(intervals, k = 6, mask = mask)
```
