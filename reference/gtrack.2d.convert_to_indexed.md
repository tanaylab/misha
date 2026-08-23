# Convert 2D track to indexed format

Converts a per-chromosome-pair 2D track (rectangles or points) to
indexed format (track.dat + track.idx). This reduces file descriptor
usage from O(N^2) to O(1), which is especially beneficial for genomes
with many contigs.

## Usage

``` r
gtrack.2d.convert_to_indexed(track = NULL, remove.old = FALSE, force = FALSE)
```

## Arguments

- track:

  track name to convert

- remove.old:

  Logical. If TRUE, removes old per-chromosome-pair files after
  successful conversion. Default: FALSE.

- force:

  Logical. If TRUE, re-converts even if already in indexed format.
  Default: FALSE.

## Value

None.

## See also

[`gtrack.2d.create`](https://tanaylab.github.io/misha/reference/gtrack.2d.create.md),
[`gtrack.2d.import`](https://tanaylab.github.io/misha/reference/gtrack.2d.import.md),
[`gtrack.2d.import_contacts`](https://tanaylab.github.io/misha/reference/gtrack.2d.import_contacts.md),
[`gtrack.convert_to_indexed`](https://tanaylab.github.io/misha/reference/gtrack.convert_to_indexed.md),
[`gdb.convert_to_indexed`](https://tanaylab.github.io/misha/reference/gdb.convert_to_indexed.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert a 2D track to indexed format
gtrack.2d.convert_to_indexed("my_2d_track")

# Convert and remove old per-pair files
gtrack.2d.convert_to_indexed("my_2d_track", remove.old = TRUE)

# Force re-conversion
gtrack.2d.convert_to_indexed("my_2d_track", force = TRUE)
} # }
```
