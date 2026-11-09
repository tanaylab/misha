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

## Details

A .gsm model saved by a misha version whose `gsynth.save` did not store
the iterator has none: before
[`gsynth.sample`](https://tanaylab.github.io/misha/reference/gsynth.sample.md)
or
[`gsynth.score`](https://tanaylab.github.io/misha/reference/gsynth.score.md)
can use a stratified model, set `model$iterator` to the iterator used in
[`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md).
An unstratified model samples without it.

## See also

[`gsynth.save`](https://tanaylab.github.io/misha/reference/gsynth.save.md),
[`gsynth.train`](https://tanaylab.github.io/misha/reference/gsynth.train.md),
[`gsynth.convert`](https://tanaylab.github.io/misha/reference/gsynth.convert.md)
