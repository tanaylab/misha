# Pre-install UCSC's gff3ToGenePred binary

Downloads UCSC's `gff3ToGenePred` static binary (~25 MB) into
`tools::R_user_dir("misha", "cache")/bin/`, verifies its SHA256, and
makes it executable. Used by
[`gdb.build_genome`](https://tanaylab.github.io/misha/reference/gdb.build_genome.md)
when the `ncbi` backend (or the `manual` backend with
`genes_format: gff3`) is invoked. Calling it directly is useful in CI or
in non-interactive scripts where the consent prompt would otherwise
fail.

## Usage

``` r
gdb.install_gff3_converter(force = FALSE)
```

## Arguments

- force:

  If `TRUE`, skip the consent prompt and re-download even if the binary
  is already cached.

## Value

The cache path of the installed binary (invisibly).

## Details

Override the binary location by setting environment variable
`MISHA_GFF3_TO_GENEPRED` to a binary you provide (for example, one
installed via `conda install -c bioconda ucsc-gff3togenepred`). This is
the recommended workaround on systems whose glibc is older than the one
UCSC's prebuilt binary requires.

## Examples

``` r
if (FALSE) { # \dontrun{
gdb.install_gff3_converter()
Sys.setenv(MISHA_GFF3_TO_GENEPRED = "/path/to/your/gff3ToGenePred")
} # }
```
