# Convert a track to indexed format

Converts a per-chromosome track to indexed format (track.dat +
track.idx).

## Usage

``` r
gtrack.convert_to_indexed(track = NULL)
```

## Arguments

- track:

  track name to convert

## Value

None

## Details

This function converts a track from the per-chromosome file format to
single-file indexed format. The indexed format dramatically reduces file
descriptor usage for genomes with many contigs and provides better
performance for parallel access.

The function performs the following steps:

1.  Validates that all per-chromosome files have consistent metadata

2.  Creates track.dat by concatenating all per-chromosome files

3.  Creates track.idx with offset/length information for each chromosome

4.  Uses atomic operations (fsync + rename) to ensure data integrity

5.  Removes the old per-chromosome files after successful conversion

## See also

[`gtrack.create`](https://tanaylab.github.io/misha/reference/gtrack.create.md),
[`gtrack.create_sparse`](https://tanaylab.github.io/misha/reference/gtrack.create_sparse.md),
[`gtrack.create_dense`](https://tanaylab.github.io/misha/reference/gtrack.create_dense.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert a track to indexed format
gtrack.convert_to_indexed("my_track")
} # }
```
