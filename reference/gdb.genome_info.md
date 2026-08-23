# Inspect a resolved genome recipe without building

Resolves `name` through the registry chain and returns the recipe (a
list) along with the source it was resolved from. Useful for previewing
what
[`gdb.build_genome`](https://tanaylab.github.io/misha/reference/gdb.build_genome.md)
would do.

## Usage

``` r
gdb.genome_info(name, registry = NULL)
```

## Arguments

- name:

  Genome name.

- registry:

  Optional path to an explicit registry YAML.

## Value

A list with components `recipe` (the resolved recipe) and
`resolved_from` (the registry source).

## See also

[`gdb.build_genome`](https://tanaylab.github.io/misha/reference/gdb.build_genome.md),
[`gdb.list_genomes`](https://tanaylab.github.io/misha/reference/gdb.list_genomes.md).

## Examples

``` r
gdb.genome_info("hg38")
#> $recipe
#> $recipe$source
#> [1] "ucsc"
#> 
#> $recipe$assembly
#> [1] "hg38"
#> 
#> 
#> $resolved_from
#> [1] "built-in (inst/genomes.yaml)"
#> 
gdb.genome_info("GCF_009806435.1")
#> $recipe
#> $recipe$source
#> [1] "ucsc-hub"
#> 
#> $recipe$accession
#> [1] "GCF_009806435.1"
#> 
#> 
#> $resolved_from
#> [1] "pattern fallback (UCSC mammal hub accession)"
#> 
```
