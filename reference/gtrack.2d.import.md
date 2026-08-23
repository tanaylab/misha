# Creates a 2D track from tab-delimited file

Creates a 2D track from tab-delimited file(s).

## Usage

``` r
gtrack.2d.import(track = NULL, description = NULL, file = NULL)
```

## Arguments

- track:

  track name

- description:

  a character string description

- file:

  vector of file paths

## Value

None.

## Details

This function creates a 2D track track from one or more tab-delimited
files. Each file must start with a header describing the columns. The
first 6 columns must have the following names: 'chrom1', 'start1',
'end1', 'chrom2', 'start2', 'end2'. The last column is designated for
the value and it may have an arbitrary name. The header is followed by a
list of intervals and a value for each interval. Overlapping intervals
are forbidden.

One can learn about the format of the tab-delimited file by running
'gextract' function on a 2D track with a 'file' parameter set to the
name of the file.

If all the imported intervals represent a point (i.e. end == start + 1)
a 'Points' track is created otherwise it is a 'Rectangles' track.

'description' is added as a track attribute.

Note: temporary files are created in the directory of the track during
the run of the function. A few of them need to be kept simultaneously
open. If the number of chromosomes and / or intervals is particularly
high, a few thousands files might be needed to be opened simultaneously.
Some operating systems limit the number of open files per user, in which
case the function might fail with "Too many open files" or similar
error. The workaround could be:

1\. Increase the limit of simultaneously opened files (the way varies
depending on your operating system). 2. Increase the value of
'gmax.data.size' option. Higher values of 'gmax.data.size' option will
increased memory usage of the function but create fewer temporary files.

## See also

[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md)
