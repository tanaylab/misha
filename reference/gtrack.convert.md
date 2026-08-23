# Converts a track to the most current format

Converts a track (if needed) to the most current format.

## Usage

``` r
gtrack.convert(src.track = NULL, tgt.track = NULL)
```

## Arguments

- src.track:

  source track name

- tgt.track:

  target track name. If 'NULL' the source track is overwritten.

## Value

None

## Details

This function converts a track to the most current format. It should be
used if a track created by an old version of the library cannot be read
anymore by the newer version. The old track is given by 'src.track'.
After conversion a new track 'tgt.track' is created. If 'tgt.track' is
'NULL' the source track is overwritten.

## See also

[`gtrack.create`](https://tanaylab.github.io/misha/reference/gtrack.create.md),
[`gtrack.2d.create`](https://tanaylab.github.io/misha/reference/gtrack.2d.create.md),
[`gtrack.create_sparse`](https://tanaylab.github.io/misha/reference/gtrack.create_sparse.md)
