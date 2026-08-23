# Divides track expression into segments

Divides the values of track expression into segments by using Wilcoxon
test.

## Usage

``` r
gsegment(
  expr = NULL,
  minsegment = NULL,
  maxpval = 0.05,
  onetailed = TRUE,
  intervals = NULL,
  iterator = NULL,
  intervals.set.out = NULL
)
```

## Arguments

- expr:

  track expression

- minsegment:

  minimal segment size

- maxpval:

  maximal P-value that separates two adjacent segments

- onetailed:

  if 'TRUE', Wilcoxon test is performed one tailed, otherwise two tailed

- intervals:

  genomic scope for which the function is applied

- iterator:

  track expression iterator of "fixed bin" type. If 'NULL' iterator is
  determined implicitly based on track expression.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a set of intervals where each interval
represents a segment.

## Details

This function divides the values of track expression into segments,
where each segment size is at least of 'minsegment' size and the P-value
of comparing the segment with the first 'minsegment' values from the
next segment is at most 'maxpval'. Comparison is done using Wilcoxon
(also known as Mann-Whitney) test.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

## See also

[`gscreen`](https://tanaylab.github.io/misha/reference/gscreen.md),
[`gwilcox`](https://tanaylab.github.io/misha/reference/gwilcox.md)

## Examples

``` r

gdb.init_examples()
gsegment("dense_track", 5000, 0.0001)
#>    chrom  start    end
#> 1   chr1      0   9200
#> 2   chr1   9200  31550
#> 3   chr1  31550  41350
#> 4   chr1  41350  65600
#> 5   chr1  65600  75800
#> 6   chr1  75800  99150
#> 7   chr1  99150 121800
#> 8   chr1 121800 129250
#> 9   chr1 129250 150650
#> 10  chr1 150650 157850
#> 11  chr1 157850 162900
#> 12  chr1 162900 224200
#> 13  chr1 224200 246200
#> 14  chr1 246200 309950
#> 15  chr1 309950 326350
#> 16  chr1 326350 346700
#> 17  chr1 346700 352200
#> 18  chr1 352200 362000
#> 19  chr1 362000 369400
#> 20  chr1 369400 379300
#> 21  chr1 379300 384800
#> 22  chr1 384800 426600
#> 23  chr1 426600 433950
#> 24  chr1 433950 500000
#> 25  chr2      0  28400
#> 26  chr2  28400  33600
#> 27  chr2  33600  38950
#> 28  chr2  38950  45800
#> 29  chr2  45800  50800
#> 30  chr2  50800  75600
#> 31  chr2  75600  85150
#> 32  chr2  85150  90600
#> 33  chr2  90600  99600
#> 34  chr2  99600 116150
#> 35  chr2 116150 122550
#> 36  chr2 122550 143400
#> 37  chr2 143400 154050
#> 38  chr2 154050 184350
#> 39  chr2 184350 194450
#> 40  chr2 194450 199600
#> 41  chr2 199600 207250
#> 42  chr2 207250 224400
#> 43  chr2 224400 229400
#> 44  chr2 229400 246000
#> 45  chr2 246000 276550
#> 46  chr2 276550 300000
#> 47  chrX      0  11300
#> 48  chrX  11300  16350
#> 49  chrX  16350  90000
#> 50  chrX  90000  95000
#> 51  chrX  95000 106200
#> 52  chrX 106200 116800
#> 53  chrX 116800 127850
#> 54  chrX 127850 135050
#> 55  chrX 135050 153450
#> 56  chrX 153450 200000
```
