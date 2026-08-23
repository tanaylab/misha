# Computes auto-correlation between the strands for a file of mapped sequences

Calculates auto-correlation between plus and minus strands for the given
chromosome in a file of mapped sequences.

## Usage

``` r
gcompute_strands_autocorr(
  file = NULL,
  chrom = NULL,
  binsize = NULL,
  maxread = 400,
  cols.order = c(9, 11, 13, 14),
  min.coord = 0,
  max.coord = 3e+08
)
```

## Arguments

- file:

  the name of the file containing mapped sequences

- chrom:

  chromosome for which the auto-correlation is computed

- binsize:

  calculate the auto-correlation for bins in the range of \[-maxread,
  maxread\]

- maxread:

  maximal length of the sequence used for statistics

- cols.order:

  order of sequence, chromosome, coordinate and strand columns in file

- min.coord:

  minimal coordinate used for statistics

- max.coord:

  maximal coordinate used for statistics

## Value

Statistics for each strand and auto-correlation by given bins.

## Details

This function calculates auto-correlation between plus and minus strands
for the given chromosome in a file of mapped sequences. Each line in the
file describes one read. Each column is separated by a TAB character.

The following columns must be presented in the file: sequence,
chromosome, coordinate and strand. The position of these columns are
controlled by 'cols.order' argument accordingly. The default value of
'cols.order' is a vector (9,11,13,14) meaning that sequence is expected
to be found at column number 9, chromosome - at column 11, coordinate -
at column 13 and strand - at column 14. The first column should be
referenced by 1 and not by 0.

Coordinates that are not in \[min.coord, max.coord\] range are ignored.

gcompute_strands_autocorr outputs the total statistics and the
auto-correlation given by bins. The size of the bin is indicated by
'binsize' parameter. Statistics is calculated for bins in the range of
\[-maxread, maxread\].

## Examples

``` r

gdb.init_examples()
gcompute_strands_autocorr(paste(.misha$GROOT, "reads", sep = "/"),
    "chr1", 50,
    maxread = 300
)
#> [[1]]
#>  Forward mean Forward stdev  Reverse mean Reverse stdev 
#>   0.007209372   0.193899377   0.003404426   0.109562905 
#> 
#> [[2]]
#>    bin         corr
#> 1   -6 -0.001155318
#> 2   -5  0.003557982
#> 3   -4 -0.001155318
#> 4   -3  0.012984580
#> 5   -2 -0.001155318
#> 6   -1 -0.001155318
#> 7    0  0.031837776
#> 8    1  0.045977673
#> 9    2  0.196803243
#> 10   3  0.385335205
#> 11   4  0.055404271
#> 12   5 -0.001155318
#> 
```
