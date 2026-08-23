# Load a gsynth.model from disk

Loads a previously saved Markov model. Auto-detects the format:

- If `file` is a directory, reads the .gsm directory format

- If `file` is a file, tries ZIP .gsm format first, then falls back to
  legacy RDS format

## Usage

``` r
gsynth.load(file)
```

## Arguments

- file:

  Path to the saved model (directory, .gsm zip, or legacy .rds)

## Value

A gsynth.model object

## See also

[`gsynth.save`](https://tanaylab.github.io/misha/reference/gsynth.save.md),
[`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md),
[`gsynth.convert`](https://tanaylab.github.io/misha/reference/gsynth.convert.md)
