# Defines modification rules for a two-dimensional iterator in a virtual track

Defines modification rules for a two-dimensional iterator in a virtual
track.

## Usage

``` r
gvtrack.iterator.2d(
  vtrack = NULL,
  sshift1 = 0,
  eshift1 = 0,
  sshift2 = 0,
  eshift2 = 0
)
```

## Arguments

- vtrack:

  virtual track name

- sshift1:

  shift of 'start1' coordinate

- eshift1:

  shift of 'end1' coordinate

- sshift2:

  shift of 'start2' coordinate

- eshift2:

  shift of 'end2' coordinate

## Value

None.

## Details

This function defines modification rules for one-dimensional iterator
intervals in a virtual track.

Iterator interval's 'start1' coordinate is modified by adding 'sshift1'.
Similarly 'end1', 'start2', 'end2' coordinates are altered by adding
'eshift1', 'sshift2' and 'eshift2' accordingly.

A shifted interval that runs past a chromosome boundary is clamped to
it. If nothing is left of it - it falls entirely outside the chromosome,
or the shifts collapse it to zero length - the virtual track returns
'NaN'.

## See also

[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md),
[`gvtrack.iterator`](https://tanaylab.github.io/misha/reference/gvtrack.iterator.md)

## Examples

``` r

gdb.init_examples()
gvtrack.create("vtrack1", "rects_track")
gvtrack.iterator.2d("vtrack1", sshift1 = 1000, eshift1 = 2000)
gextract(
    "rects_track", "vtrack1",
    gintervals.2d(1, 0, 5000, 2, 0, 5000)
)
#> NULL
```
