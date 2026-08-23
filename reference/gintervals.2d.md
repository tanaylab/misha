# Creates a set of 2D intervals

Creates a set of 2D intervals.

## Usage

``` r
gintervals.2d(
  chroms1 = NULL,
  starts1 = 0,
  ends1 = -1,
  chroms2 = NULL,
  starts2 = 0,
  ends2 = -1
)
```

## Arguments

- chroms1:

  chromosomes1 - an array of strings with or without "chr" prefixes or
  an array of integers (like: '1' for "chr1")

- starts1:

  an array of start1 coordinates

- ends1:

  an array of end1 coordinates. If '-1' chromosome size is assumed.

- chroms2:

  chromosomes2 - an array of strings with or without "chr" prefixes or
  an array of integers (like: '1' for "chr1"). If 'NULL', 'chroms2' is
  assumed to be equal to 'chroms1'.

- starts2:

  an array of start2 coordinates

- ends2:

  an array of end2 coordinates. If '-1' chromosome size is assumed.

## Value

A data frame representing the intervals.

## Details

This function returns a set of two-dimensional intervals. The returned
value can be used in all functions that accept 'intervals' argument.

Two-dimensional intervals is a data frame whose first six columns are
'chrom1', 'start1', 'end1', 'chrom2', 'start2' and 'end2'. Each row of
the data frame represents two genomic intervals from two chromosomes in
the range of \[start, end). Additional columns can be presented in 2D
intervals object yet these columns must be added after the six
obligatory ones.

## See also

[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.force_range`](https://tanaylab.github.io/misha/reference/gintervals.force_range.md)

## Examples

``` r

gdb.init_examples()

## the following 3 calls produce identical results
gintervals.2d(1)
#>   chrom1 start1  end1 chrom2 start2  end2
#> 1   chr1      0 5e+05   chr1      0 5e+05
gintervals.2d("1")
#>   chrom1 start1  end1 chrom2 start2  end2
#> 1   chr1      0 5e+05   chr1      0 5e+05
gintervals.2d("chrX")
#>   chrom1 start1  end1 chrom2 start2  end2
#> 1   chrX      0 2e+05   chrX      0 2e+05

gintervals.2d(1, 1000, 2000, "chrX", 400, 800)
#>   chrom1 start1 end1 chrom2 start2 end2
#> 1   chr1   1000 2000   chrX    400  800
gintervals.2d(c("chr2", "chrX"), 10, c(3000, 5000), 1)
#>   chrom1 start1 end1 chrom2 start2  end2
#> 1   chr2     10 3000   chr1      0 5e+05
#> 2   chrX     10 5000   chr1      0 5e+05
```
