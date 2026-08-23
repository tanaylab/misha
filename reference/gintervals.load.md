# Loads a named intervals set

Loads a named intervals set.

## Usage

``` r
gintervals.load(
  intervals.set = NULL,
  chrom = NULL,
  chrom1 = NULL,
  chrom2 = NULL
)
```

## Arguments

- intervals.set:

  name of an intervals set

- chrom:

  chromosome for 1D intervals set

- chrom1:

  first chromosome for 2D intervals set

- chrom2:

  second chromosome for 2D intervals set

## Value

A data frame representing the intervals.

## Details

This function loads and returns intervals stored in a named intervals
set.

If intervals set contains 1D intervals and 'chrom' is not 'NULL' only
the intervals of the given chromosome are returned.

Likewise if intervals set contains 2D intervals and 'chrom1', 'chrom2'
are not 'NULL' only the intervals of the given pair of chromosomes are
returned.

For big intervals sets 'chrom' parameter (1D case) / 'chrom1', 'chrom2'
parameters (2D case) must be specified. In other words: big intervals
sets can be loaded only by chromosome or chromosome pair.

## See also

[`gintervals.save`](https://tanaylab.github.io/misha/reference/gintervals.save.md),
[`gintervals.is.bigset`](https://tanaylab.github.io/misha/reference/gintervals.is.bigset.md),
[`gintervals.exists`](https://tanaylab.github.io/misha/reference/gintervals.exists.md),
[`gintervals.ls`](https://tanaylab.github.io/misha/reference/gintervals.ls.md),
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md)

## Examples

``` r

gdb.init_examples()
gintervals.load("annotations")
#>   chrom start   end strand  remark
#> 1  chr1    20  2000      1 remarkA
#> 5  chr1  2500  2600      1 remarkE
#> 2  chr2    20  2000      1 remarkB
#> 3  chr2  3000  8000     -1 remarkC
#> 4  chr2  9000 11000      1 remarkD
#> 6  chr2 12000 12001      1 remarkF
#> 7  chr2 13000 14000     -1 remarkG
#> 8  chr2 15000 15500      1 remarkH
```
