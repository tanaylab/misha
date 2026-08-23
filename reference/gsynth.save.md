# Save a gsynth.model to disk in .gsm format

Saves a trained Markov model in the cross-platform .gsm format, which
consists of a metadata YAML file and raw binary arrays for counts and
CDFs. The .gsm format can be stored as a directory (default) or a ZIP
archive.

## Usage

``` r
gsynth.save(model, file, compress = FALSE)
```

## Arguments

- model:

  A gsynth.model object from
  [`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md)

- file:

  Path to save the model (directory or .zip file)

- compress:

  Logical. If `TRUE`, save as a ZIP archive. If `FALSE` (default), save
  as a directory.

## Value

Invisibly returns the file path.

## See also

[`gsynth.load`](https://tanaylab.github.io/misha/reference/gsynth.load.md),
[`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md),
[`gsynth.convert`](https://tanaylab.github.io/misha/reference/gsynth.convert.md)
