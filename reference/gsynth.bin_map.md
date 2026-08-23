# Create a bin mapping from value-based merge specifications

Converts value-based bin merge specifications into a bin_map named
vector that can be used with
[`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md).
This allows you to specify merges using actual track values rather than
bin indices.

## Usage

``` r
gsynth.bin_map(breaks, merge_ranges = NULL)
```

## Arguments

- breaks:

  Numeric vector of bin boundaries (same as used in
  [`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md))

- merge_ranges:

  List of merge specifications. Each specification is a named list with:

  from

  :   Numeric vector of length 2 `c(min, max)` defining the source value
      range to merge. Use `-Inf` or `Inf` for open-ended ranges. Can
      also be a single number (shorthand for `c(value, Inf)`).

  to

  :   Numeric vector of length 2 `c(min, max)` defining the target bin
      that source bins should map to. Must match an existing bin defined
      by `breaks`.

## Value

A named vector (bin_map) compatible with `bin_map` parameter in
[`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md).
The names are source bin indices (1-based), and values are target bin
indices (1-based).

## See also

[`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md),
[`gsynth.cell_merge`](https://tanaylab.github.io/misha/reference/gsynth.cell_merge.md),
[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md)

## Examples

``` r
# Define breaks for GC content [0, 1] in 0.025 increments
breaks <- seq(0, 1, 0.025)

# Merge all GC content above 70% (0.7) into the bin (0.675, 0.7]
bin_map <- gsynth.bin_map(
    breaks = breaks,
    merge_ranges = list(
        list(from = 0.7, to = c(0.675, 0.7))
    )
)

# Multiple merges: merge low GC (< 0.3) and high GC (> 0.7) into middle bins
bin_map2 <- gsynth.bin_map(
    breaks = breaks,
    merge_ranges = list(
        list(from = c(-Inf, 0.3), to = c(0.4, 0.425)), # low GC -> (0.4, 0.425]
        list(from = 0.7, to = c(0.675, 0.7)) # high GC -> (0.675, 0.7]
    )
)
```
