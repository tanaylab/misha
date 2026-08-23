# Transplant sequences from one genome into another

Sugar function that extracts sequences from a source genome at the given
intervals and implants them into a target genome. Equivalent to calling
`ggenome.implant` with `donor = source_genome` and
`genome_fasta = target_genome`.

## Usage

``` r
ggenome.transplant(
  intervals,
  source_genome,
  target_genome = NULL,
  output,
  create_trackdb = TRUE,
  trackdb_path = NULL,
  line_width = 80L,
  overwrite = FALSE
)
```

## Arguments

- intervals:

  A data.frame with `chrom`, `start`, `end` columns specifying the
  regions to transplant.

- source_genome:

  Path to the misha database root containing the donor sequences.

- target_genome:

  Path to the target reference FASTA file. If `NULL`, the current misha
  database is used.

- output:

  Path for the output FASTA file.

- create_trackdb:

  Logical. If `TRUE`, creates a misha trackdb.

- trackdb_path:

  Path for the new trackdb.

- line_width:

  Integer. Number of bases per FASTA line.

- overwrite:

  Logical. If `TRUE`, overwrite existing output.

## Value

Invisibly returns the output FASTA path.

## See also

[`ggenome.implant`](https://tanaylab.github.io/misha/reference/ggenome.implant.md),
[`gdb.export_fasta`](https://tanaylab.github.io/misha/reference/gdb.export_fasta.md),
[`gseq.extract`](https://tanaylab.github.io/misha/reference/gseq.extract.md)

## Examples

``` r
gdb.init_examples()

# Create a "donor" DB with different sequence (all T's)
donor_fasta <- tempfile(fileext = ".fa")
cat(">chr1\n", paste(rep("T", 500000), collapse = ""), "\n",
    ">chr2\n", paste(rep("T", 300000), collapse = ""), "\n",
    file = donor_fasta, sep = ""
)
donor_db <- tempfile()
gdb.create(donor_db, fasta = donor_fasta, verbose = FALSE)

# Export the current DB as the target FASTA
gdb.init_examples()
ref_fasta <- tempfile(fileext = ".fa")
gdb.export_fasta(ref_fasta)

# Transplant donor sequence into positions 100-200 of chr1
intervals <- data.frame(chrom = "chr1", start = 100, end = 200)
out <- tempfile(fileext = ".fa")
trackdb <- tempfile()
ggenome.transplant(intervals,
    source_genome = donor_db,
    target_genome = ref_fasta,
    output = out,
    create_trackdb = TRUE,
    trackdb_path = trackdb
)

# Verify: positions 100-200 should now be all T's
gdb.init(trackdb)
gseq.extract(data.frame(chrom = "chr1", start = 100, end = 200))
#> [1] "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

# Clean up
unlink(
    c(
        donor_fasta, donor_db, ref_fasta, out,
        paste0(out, ".fai"), trackdb
    ),
    recursive = TRUE
)
```
