# Calculates quantiles of a track expression

Calculates the quantiles of a track expression for the given
percentiles.

## Usage

``` r
gquantiles(
  expr = NULL,
  percentiles = 0.5,
  intervals = get("ALLGENOME", envir = .misha),
  iterator = NULL,
  band = NULL
)
```

## Arguments

- expr:

  track expression

- percentiles:

  an array of percentiles of quantiles in \[0, 1\] range

- intervals:

  genomic scope for which the function is applied

- iterator:

  track expression iterator. If 'NULL' iterator is determined implicitly
  based on track expression.

- band:

  track expression band. If 'NULL' no band is used.

## Value

An array that represent quantiles.

## Details

This function calculates the quantiles for the given percentiles.

If data size exceeds the limit (see: 'getOption(gmax.data.size)'), the
data is randomly sampled to fit the limit. A warning message is
generated. The seed of the pseudo-random generator can be controlled
through 'grnd.seed' option.

Note: this function is capable to run in multitasking mode. Sampling may
vary according to the extent of multitasking. Since multitasking depends
on the number of available CPU cores, running the function on two
different machines might give different results. Please switch off
multitasking if you want to achieve identical results on any machine.
For more information regarding multitasking please refer "User Manual".

## See also

[`gbins.quantiles`](https://tanaylab.github.io/misha/reference/gbins.quantiles.md),
[`gintervals.quantiles`](https://tanaylab.github.io/misha/reference/gintervals.quantiles.md),
[`gdist`](https://tanaylab.github.io/misha/reference/gdist.md)

## Examples

``` r

gdb.init_examples()
gquantiles("dense_track", c(0.1, 0.6, 0.8), gintervals(c(1, 2)))
#>  0.1  0.6  0.8 
#> 0.02 0.10 0.14 
```
