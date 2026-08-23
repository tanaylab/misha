# Create and Load a Genome Database

This function downloads, extracts, and loads a misha genome database for
the specified genome.

## Usage

``` r
gdb.create_genome(genome, path = getwd(), tmpdir = tempdir())
```

## Arguments

- genome:

  A character string specifying the genome to download. Supported
  genomes are "mm9", "mm10", "mm39", "hg19", and "hg38".

- path:

  A character string specifying the directory where the genome will be
  extracted. Defaults to genome name (e.g. "mm10") in the current
  working directory.

- tmpdir:

  A character string specifying the directory for storing temporary
  files. This is used for storing the downloaded genome file.

## Value

None.

## Details

The function checks if the specified genome is available. If tmpdir, it
constructs the download URL, downloads the genome file, extracts it to
the specified directory, and loads the genome database using `gsetroot`.
The function also calls `gdb.reload` to reload the genome database.

## Examples

``` r
if (FALSE) { # \dontrun{
mm10_dir <- tempdir()
gdb.create_genome("mm10", path = mm10_dir)
list.files(file.path(mm10_dir, "mm10"))
gsetroot(file.path(mm10_dir, "mm10"))
gintervals.ls()
} # }
```
