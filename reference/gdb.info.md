# Get Database Information

Returns information about a misha genome database including format,
number of chromosomes, total genome size, and whether it uses the
indexed format.

## Usage

``` r
gdb.info(groot = NULL)
```

## Arguments

- groot:

  Root directory of the database. If NULL, uses the currently active
  database.

## Value

A list with database information:

- `path` - Full path to the database

- `is_db` - TRUE if this is a valid misha database

- `format` - "indexed" or "per-chromosome"

- `num_chromosomes` - Number of chromosomes/contigs

- `genome_size` - Total length of genome in bases

- `chromosomes` - Data frame with chromosome names and sizes

## Examples

``` r
if (FALSE) { # \dontrun{
# Get info about currently active database
info <- gdb.info()
cat("Database format:", info$format, "\n")
cat("Genome size:", info$genome_size / 1e6, "Mb\n")

# Get info about specific database
info <- gdb.info("/path/to/database")
} # }
```
