# Defines modification rules for a one-dimensional iterator in a virtual track

Defines modification rules for a one-dimensional iterator in a virtual
track.

## Usage

``` r
gvtrack.iterator(vtrack = NULL, dim = NULL, sshift = 0, eshift = 0)
```

## Arguments

- vtrack:

  virtual track name

- dim:

  use 'NULL' or '0' for 1D iterators. '1' converts 2D iterator to
  (chrom1, start1, end1) , '2' converts 2D iterator to (chrom2, start2,
  end2)

- sshift:

  shift of 'start' coordinate

- eshift:

  shift of 'end' coordinate

## Value

None.

## Details

This function defines modification rules for one-dimensional iterator
intervals in a virtual track.

'dim' converts a 2D iterator interval (chrom1, start1, end1, chrom2,
start2, end2) to a 1D interval. If 'dim' is '1' the interval is
converted to (chrom1, start1, end1). If 'dim' is '2' the interval is
converted to (chrom2, start2, end2). If 1D iterator is used 'dim' must
be set to 'NULL' or '0' (meaning: no conversion is made).

Iterator interval's 'start' coordinate is modified by adding 'sshift'.
Similarly 'end' coordinate is altered by adding 'eshift'.

A shifted interval that runs past a chromosome boundary is clamped to
it. If nothing is left of it - it falls entirely outside the chromosome,
or the shifts collapse it to zero length - the virtual track returns
'NaN'.

## See also

[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md),
[`gvtrack.iterator.2d`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.2d.md)

## Examples

``` r

gdb.init_examples()

gvtrack.create("vtrack1", "dense_track")
gvtrack.iterator("vtrack1", sshift = 200, eshift = 200)
gextract("dense_track", "vtrack1", gintervals(1, 0, 500))
#>    chrom start end dense_track vtrack1 intervalID
#> 1   chr1     0  50   0.1777778    0.16          1
#> 2   chr1    50 100   0.1600000    0.20          1
#> 3   chr1   100 150   0.1800000    0.16          1
#> 4   chr1   150 200   0.1600000    0.16          1
#> 5   chr1   200 250   0.1600000    0.16          1
#> 6   chr1   250 300   0.2000000    0.06          1
#> 7   chr1   300 350   0.1600000    0.06          1
#> 8   chr1   350 400   0.1600000    0.02          1
#> 9   chr1   400 450   0.1600000    0.04          1
#> 10  chr1   450 500   0.0600000    0.00          1

gvtrack.create("vtrack2", "dense_track")
gvtrack.iterator("vtrack2", dim = 1)
gextract("vtrack2", gintervals.2d(1, 0, 1000, 1, 0, -1),
    iterator = "rects_track"
)
#> NULL
```
