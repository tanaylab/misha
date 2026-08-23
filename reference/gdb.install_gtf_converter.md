# Pre-install UCSC's gtfToGenePred binary

Mirrors
[`gdb.install_gff3_converter`](https://tanaylab.github.io/misha/reference/gdb.install_gff3_converter.md).
Required for the `ucsc-hub` backend's `genes` set (UCSC mammal hubs ship
GTFs).

## Usage

``` r
gdb.install_gtf_converter(force = FALSE)
```

## Arguments

- force:

  If `TRUE`, skip consent prompt and re-download even if cached.

## Value

The cache path (invisibly).

## Details

Override the binary location by setting environment variable
`MISHA_GTF_TO_GENEPRED`.

## Examples

``` r
if (FALSE) { # \dontrun{
gdb.install_gtf_converter()
Sys.setenv(MISHA_GTF_TO_GENEPRED = "/path/to/your/gtfToGenePred")
} # }
```
