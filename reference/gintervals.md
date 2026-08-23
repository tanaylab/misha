# Creates a set of 1D intervals

Creates a set of 1D intervals.

## Usage

``` r
gintervals(chroms = NULL, starts = 0, ends = -1, strands = NULL)
```

## Arguments

- chroms:

  chromosomes - an array of strings with or without "chr" prefixes or an
  array of integers (like: '1' for "chr1")

- starts:

  an array of start coordinates

- ends:

  an array of end coordinates. If '-1' chromosome size is assumed.

- strands:

  'NULL', a numeric vector of '-1', '0' or '1' values, or a
  character/factor vector with values "+", "-", ".", "\*" or ""

## Value

A data frame representing the intervals, sorted in the canonical
chromosome order (which is not the order of the arguments). Beware of
keeping per-interval data in a separate parallel vector across a
`gintervals()` call: the rows are reordered and the vector is not. Put
such data in a column of the resulting data frame instead.

## Details

This function returns a set of one-dimensional intervals. The returned
value can be used in all functions that accept 'intervals' argument.

One-dimensional intervals is a data frame whose first three columns are
'chrom', 'start' and 'end'. Each row of the data frame represents a
genomic interval of the specified chromosome in the range of \[start,
end). Additional columns can be presented in 1D intervals object yet
these columns must be added after the three obligatory ones.

If 'strands' argument is not 'NULL' an additional column "strand" is
added to the intervals. The possible values of a strand can be '1' (plus
strand), '-1' (minus strand) or '0' (unknown). Character values "+",
"-", ".", "\*" and "" (or factors with these levels) are also accepted
and converted internally to '1', '-1' and '0' respectively.

## See also

[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md),
[`gintervals.force_range`](https://tanaylab.github.io/misha/reference/gintervals.force_range.md)

## Examples

``` r

gdb.init_examples()

## the following 3 calls produce identical results
gintervals(1)
#>   chrom start   end
#> 1  chr1     0 5e+05
gintervals("1")
#>   chrom start   end
#> 1  chr1     0 5e+05
gintervals("chrX")
#>   chrom start   end
#> 1  chrX     0 2e+05

gintervals(1, 1000)
#>   chrom start   end
#> 1  chr1  1000 5e+05
gintervals(c("chr2", "chrX"), 10, c(3000, 5000))
#>   chrom start  end
#> 1  chr2    10 3000
#> 2  chrX    10 5000
```
