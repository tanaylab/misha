# Export a database genome as FASTA

Writes all contigs from a misha database to a multi-FASTA file.

## Usage

``` r
gdb.export_fasta(
  file = NULL,
  groot = NULL,
  line_width = 80L,
  chunk_size = 1000000L,
  overwrite = FALSE,
  verbose = FALSE
)
```

## Arguments

- file:

  Output FASTA file path

- groot:

  Optional database root path. If NULL, uses current database.

- line_width:

  Number of bases per FASTA line. Default: 80.

- chunk_size:

  Number of bases to extract per chunk while writing. Default: 1000000.

- overwrite:

  Logical. If TRUE, overwrite existing output file. Default: FALSE.

- verbose:

  Logical. If TRUE, prints progress messages. Default: FALSE.

## Value

Invisibly returns `file`.

## Details

By default, the currently active database is used. You can also provide
`groot` to export another database without changing the caller's active
database.

## See also

[`gdb.init`](https://tanaylab.github.io/misha/reference/gdb.init.md),
[`gseq.extract`](https://tanaylab.github.io/misha/reference/gseq.extract.md)

## Examples

``` r
if (FALSE) { # \dontrun{
gdb.init_examples()
out <- tempfile(fileext = ".fa")
gdb.export_fasta(out)
head(readLines(out))
} # }
```
