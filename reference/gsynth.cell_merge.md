# Resolve a cell-level merge specification into flat bin indices

Unlike
[`gsynth.bin_map`](https://tanaylab.github.io/misha/reference/gsynth.bin_map.md),
which merges bins independently along each dimension,
`gsynth.cell_merge` resolves per-*joint-cell* redirects: each entry
redirects one specific training cell (identified by its per-dimension
values) to another specific training cell. This lets callers redirect
arbitrary cells whose Cartesian position cannot be expressed via
per-axis merges (e.g., “cell (GC=0.725, CG=0.05) –\> cell (GC=0.70,
CG=0.08)”).

## Usage

``` r
gsynth.cell_merge(model, cell_merge, bin_merge = NULL)
```

## Arguments

- model:

  A `gsynth.model` object from
  [`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md).

- cell_merge:

  A list of redirect specifications. Each entry is a named list with:

  from

  :   Numeric vector of length `n_dims`; one representative value per
      dimension that identifies the *source* cell.

  to

  :   Numeric vector of length `n_dims`; one representative value per
      dimension that identifies the *target* cell.

  A single data frame with columns
  `from_1, from_2, ..., to_1, to_2, ...` is also accepted and converted
  internally.

- bin_merge:

  Optional sampling-time bin merge specification (same format as in
  [`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md)).
  When supplied, source and target per-dimension bin indices are
  remapped through the resulting per-axis maps before being combined
  into flat indices, so that cell values reference post-bin_merge cells.

## Value

A data frame with one row per `cell_merge` entry and columns:

- `from_<d>`, `to_<d>` for each dimension `d`: 1-based bin index after
  any bin_merge remapping.

- `source_flat`, `target_flat`: 1-based flat bin index into
  `model$model_data$cdf`.

## Details

This function is primarily a utility for inspecting / debugging what
[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md)
will do when invoked with the `cell_merge` argument. It returns a data
frame describing, for every entry, the resolved source and target cells
in both per-dimension-bin and flat-bin space.

## See also

[`gsynth.bin_map`](https://tanaylab.github.io/misha/reference/gsynth.bin_map.md),
[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Resolve a redirect table before handing it to gsynth.sample:
redirects <- list(
    list(from = c(0.725, 0.05), to = c(0.70, 0.08)),
    list(from = c(0.75, 0.06), to = c(0.70, 0.08))
)
resolved <- gsynth.cell_merge(model, redirects)
print(resolved)

gsynth.sample(model, "out.fa",
    output_format = "fasta",
    cell_merge = redirects
)
} # }
```
