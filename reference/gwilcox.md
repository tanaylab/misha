# Calculates Wilcoxon test on sliding windows over track expression

Calculates Wilcoxon test on sliding windows over the values of track
expression.

## Usage

``` r
gwilcox(
  expr = NULL,
  winsize1 = NULL,
  winsize2 = NULL,
  maxpval = 0.05,
  onetailed = TRUE,
  what2find = 1,
  intervals = NULL,
  iterator = NULL,
  intervals.set.out = NULL
)
```

## Arguments

- expr:

  track expression

- winsize1:

  number of values in the first sliding window

- winsize2:

  number of values in the second sliding window

- maxpval:

  maximal P-value

- onetailed:

  if 'TRUE', Wilcoxon test is performed one tailed, otherwise two tailed

- what2find:

  if '-1', lows are searched. If '1', peaks are searched. If '0', both
  peaks and lows are searched

- intervals:

  genomic scope for which the function is applied

- iterator:

  track expression iterator of "fixed bin" type. If 'NULL' iterator is
  determined implicitly based on track expression.

- intervals.set.out:

  intervals set name where the function result is optionally outputted

## Value

If 'intervals.set.out' is 'NULL' a data frame representing the intervals
with an additional 'pval' column where P-value is below 'maxpval'.

## Details

This function runs a Wilcoxon test (also known as a Mann-Whitney test)
over the values of track expression in the two sliding windows having an
identical center. The sizes of the windows are specified by 'winsize1'
and 'winsize2'. 'gwilcox' returns intervals where the smaller window
tested against a larger window gives a P-value below 'maxpval'. The test
can be one or two tailed.

'what2find' argument controls what should be searched: peaks, lows or
both.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
set. Use this parameter if the result size exceeds the limits of the
physical memory.

## See also

[`gscreen`](https://tanaylab.github.io/misha/reference/gscreen.md),
[`gsegment`](https://tanaylab.github.io/misha/reference/gsegment.md)

## Examples

``` r

gdb.init_examples()
gwilcox("dense_track", 100000, 1000,
    maxpval = 0.01,
    what2find = 1
)
#>    chrom  start    end         pval
#> 1   chr1      0    650 6.259068e-04
#> 2   chr1  31750  36650 3.230735e-09
#> 3   chr1  37600  39750 8.149484e-06
#> 4   chr1  60450  62100 5.452999e-05
#> 5   chr1  64600  68900 1.536396e-05
#> 6   chr1  69400  71000 5.134420e-05
#> 7   chr1  71050  77050 3.392364e-09
#> 8   chr1  89900  91700 5.214742e-03
#> 9   chr1  99150 101100 4.321116e-05
#> 10  chr1 106650 110850 7.451052e-07
#> 11  chr1 114550 117800 4.517828e-08
#> 12  chr1 119600 120900 3.503083e-03
#> 13  chr1 138050 140000 1.391299e-03
#> 14  chr1 150450 153300 4.480170e-05
#> 15  chr1 153400 155850 7.935469e-04
#> 16  chr1 166000 167800 1.995862e-03
#> 17  chr1 221450 223550 2.754194e-04
#> 18  chr1 236200 237800 6.262228e-03
#> 19  chr1 246100 247300 2.073013e-04
#> 20  chr1 253050 257200 9.573318e-07
#> 21  chr1 323900 325350 1.789731e-03
#> 22  chr1 327000 333350 3.616736e-07
#> 23  chr1 335650 338600 8.590741e-07
#> 24  chr1 339250 341950 1.306270e-04
#> 25  chr1 345600 346850 1.544731e-03
#> 26  chr1 353150 354400 1.328745e-03
#> 27  chr1 355850 357400 8.570873e-04
#> 28  chr1 358800 359900 9.412577e-03
#> 29  chr1 370500 371800 8.343761e-04
#> 30  chr1 372800 374250 2.156471e-05
#> 31  chr1 378700 385200 5.836447e-11
#> 32  chr1 398750 402200 7.187603e-04
#> 33  chr1 404600 406150 4.861715e-04
#> 34  chr1 407750 408800 5.158613e-03
#> 35  chr1 409000 411300 1.238897e-04
#> 36  chr1 413150 414600 5.000712e-04
#> 37  chr1 420050 421100 9.127494e-03
#> 38  chr1 437500 438850 7.115249e-05
#> 39  chr1 443000 444200 7.526990e-03
#> 40  chr1 450200 451450 1.502912e-03
#> 41  chr1 452800 454150 5.634265e-03
#> 42  chr1 455850 457450 4.320016e-06
#> 43  chr1 459850 461750 6.644247e-05
#> 44  chr2   2650   3950 2.982702e-03
#> 45  chr2   4100   5250 5.145261e-03
#> 46  chr2  29300  31200 7.186895e-04
#> 47  chr2  31500  32950 5.595490e-03
#> 48  chr2  37650  40800 3.808473e-05
#> 49  chr2  42200  45150 5.061067e-06
#> 50  chr2  50450  52000 1.174162e-04
#> 51  chr2  56350  59550 8.522234e-04
#> 52  chr2  76450  77500 9.085247e-03
#> 53  chr2  78500  80900 1.193580e-09
#> 54  chr2  82850  84650 8.147982e-06
#> 55  chr2  84750  89700 6.851225e-09
#> 56  chr2  99600 101550 4.784221e-04
#> 57  chr2 113300 114350 9.711064e-03
#> 58  chr2 124050 125800 4.682674e-04
#> 59  chr2 133550 135300 1.988063e-03
#> 60  chr2 138850 140100 7.080966e-04
#> 61  chr2 153350 156000 4.035042e-06
#> 62  chr2 157300 160900 3.763693e-07
#> 63  chr2 161850 163800 5.561545e-04
#> 64  chr2 164450 165800 4.125066e-03
#> 65  chr2 171000 173650 7.423077e-07
#> 66  chr2 178200 180500 2.461700e-05
#> 67  chr2 195750 196800 8.421521e-03
#> 68  chr2 197750 199700 3.419243e-05
#> 69  chr2 208000 209250 7.472793e-04
#> 70  chr2 214900 218500 5.921503e-08
#> 71  chr2 222750 224500 2.526845e-03
#> 72  chr2 230200 234050 4.451623e-10
#> 73  chr2 236150 238350 1.136839e-03
#> 74  chr2 239800 241200 2.977583e-03
#> 75  chr2 249050 250400 6.186424e-03
#> 76  chr2 251000 252500 7.463579e-04
#> 77  chr2 268750 269900 5.404431e-03
#> 78  chr2 275400 277000 2.654262e-03
#> 79  chr2 281150 282450 1.399682e-03
#> 80  chr2 288350 290950 4.984661e-09
#> 81  chr2 293750 295350 5.226153e-03
#> 82  chrX      0    900 1.768339e-04
#> 83  chrX   1000   3100 3.085641e-07
#> 84  chrX   5050   6800 1.793621e-06
#> 85  chrX  10850  16550 7.498458e-09
#> 86  chrX  84350  86200 1.268948e-06
#> 87  chrX  89800  95100 4.628085e-07
#> 88  chrX 116400 118050 3.115843e-04
#> 89  chrX 119400 122300 5.102100e-05
#> 90  chrX 124400 128100 2.769534e-08
#> 91  chrX 134900 136850 2.763108e-04
#> 92  chrX 138450 139550 4.878393e-03
#> 93  chrX 142450 144300 1.585345e-05
#> 94  chrX 149800 152200 9.806932e-05
#> 95  chrX 159300 161000 1.390129e-04
#> 96  chrX 165600 167500 9.145936e-04
```
