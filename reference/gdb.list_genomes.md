# List resolvable genome names

Returns a data frame describing every genome resolvable from the active
registry chain (see
[`gdb.build_genome`](https://tanaylab.github.io/misha/reference/gdb.build_genome.md)
for the chain order).

## Usage

``` r
gdb.list_genomes(registry = NULL)
```

## Arguments

- registry:

  Optional path to an explicit registry YAML, overriding the resolution
  chain.

## Value

A data frame with columns:

- `name` - registry key.

- `source` - backend (`ucsc`, `ncbi`, `s3`, `local`, `manual`).

- `detail` - assembly / accession / path.

- `resolved_from` - which registry the entry came from.

## See also

[`gdb.build_genome`](https://tanaylab.github.io/misha/reference/gdb.build_genome.md),
[`gdb.genome_info`](https://tanaylab.github.io/misha/reference/gdb.genome_info.md).

## Examples

``` r
gdb.list_genomes()
#>        name source   detail resolved_from
#> 1      hg19   ucsc     hg19      built-in
#> 2      hg38   ucsc     hg38      built-in
#> 3       mm9   ucsc      mm9      built-in
#> 4      mm10   ucsc     mm10      built-in
#> 5      mm39   ucsc     mm39      built-in
#> 6       rn6   ucsc      rn6      built-in
#> 7       rn7   ucsc      rn7      built-in
#> 8       dm6   ucsc      dm6      built-in
#> 9      ce11   ucsc     ce11      built-in
#> 10  sacCer3   ucsc  sacCer3      built-in
#> 11 danRer11   ucsc danRer11      built-in
```
