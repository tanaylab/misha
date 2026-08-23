# Choose the multitasking strategy for gextract

Read `getOption("gmultitasking.strategy", "auto")` and resolve "auto" to
either "tiles" (the historical misha multitask: each kid handles a tile
range across all tracks) or "tracks" (R-side mclapply: each worker
handles a track subset across all tiles). Track-parallel is a major win
on cold-NFS many-track gextract workloads (~5× measured on Tamar's
30-track × 2.19M-bin bench) because each worker only mmap-faults its own
files, avoiding the working-set thrashing that 24 kids × 500 tracks
creates on the kernel page cache.

## Usage

``` r
.gmultitasking_strategy(
  tracks,
  intervals,
  iterator = NULL,
  file = NULL,
  intervals.set.out = NULL,
  band = NULL
)
```

## Arguments

- tracks:

  character vector of track expressions

- intervals:

  data.frame of intervals (or NULL for big intervals sets)

- file:

  optional output file path (NULL = data.frame return)

- intervals.set.out:

  optional bigset output name

- band:

  optional 2D band parameter

## Value

one of "tiles" or "tracks"

## Details

Falls back to "tiles" when track-parallel doesn't apply (single track,
file/bigset output, 2D band iteration, or workload too small to amortize
fork overhead).
