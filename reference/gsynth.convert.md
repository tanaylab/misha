# Convert a legacy RDS gsynth model to .gsm format

Reads a gsynth.model from a legacy RDS file and saves it in the
cross-platform .gsm format.

## Usage

``` r
gsynth.convert(input_file, output_file, compress = FALSE)
```

## Arguments

- input_file:

  Path to the legacy RDS model file

- output_file:

  Path for the output .gsm model (directory or zip)

- compress:

  Logical. If `TRUE`, save as a ZIP archive. If `FALSE` (default), save
  as a directory.

## Value

Invisibly returns the output file path.

## See also

[`gsynth.save`](https://tanaylab.github.io/misha/reference/gsynth.save.md),
[`gsynth.load`](https://tanaylab.github.io/misha/reference/gsynth.load.md)
