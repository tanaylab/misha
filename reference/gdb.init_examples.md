# Initialise the example Genomic Database

Extracts the bundled `testdb` tarball under `dir` and points the active
groot at `<dir>/trackdb/test`. Used by examples and tests.

## Usage

``` r
gdb.init_examples(dir = NULL)
```

## Arguments

- dir:

  Directory under which to extract the example trackdb. Created if it
  does not exist. Defaults to `Sys.getenv("MISHA_EXAMPLES_DIR")` if set,
  otherwise [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

## Value

None.

## Details

If `dir` is `NULL` (default), uses the value of environment variable
`MISHA_EXAMPLES_DIR`, falling back to
[`tempdir()`](https://rdrr.io/r/base/tempfile.html). Set
`MISHA_EXAMPLES_DIR` to a roomy path when `/tmp` is tight.

## See also

[`gdb.init`](https://tanaylab.github.io/misha/reference/gdb.init.md)
