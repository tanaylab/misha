# Implant donor sequences into a reference genome

Replaces specified intervals in a reference genome with donor DNA
sequences and writes the result as a new FASTA file. Optionally creates
a misha trackdb from the output.

## Usage

``` r
ggenome.implant(
  intervals,
  donor,
  output,
  genome_fasta = NULL,
  create_trackdb = TRUE,
  trackdb_path = NULL,
  line_width = 80L,
  overwrite = FALSE
)
```

## Arguments

- intervals:

  A data.frame with `chrom`, `start`, `end` columns specifying the
  regions to replace. Coordinates are 0-based, half-open (standard misha
  convention).

- donor:

  Either a character vector of DNA sequences (one per interval row), or
  a single string path to a misha database root from which sequences
  will be extracted at the same intervals.

- output:

  Path for the output FASTA file.

- genome_fasta:

  Path to the reference FASTA file to edit. If `NULL`, the current misha
  database is exported via `gdb.export_fasta`.

- create_trackdb:

  Logical. If `TRUE`, creates a misha trackdb from the output FASTA
  using `gdb.create`.

- trackdb_path:

  Path for the new trackdb. Defaults to `<dirname(output)>/trackdb`.

- line_width:

  Integer. Number of bases per FASTA line. Default: 80.

- overwrite:

  Logical. If `TRUE`, overwrite existing output file. Default: `FALSE`.

## Value

Invisibly returns the output FASTA path.

## Details

This function fills the gap between manual sequence extraction and
genome perturbation by providing a single-call interface for editing
genomes.

The `donor` parameter controls where replacement sequences come from:

- If `donor` is a character vector with one sequence per interval row,
  those literal sequences are used directly.

- If `donor` is a single string pointing to an existing directory, it is
  treated as a misha database root. Sequences are extracted from that
  database at the same coordinates as `intervals`.

Perturbations are applied in reverse coordinate order within each
chromosome so that earlier coordinates remain valid when later ones are
replaced.

## See also

[`ggenome.transplant`](https://tanaylab.github.io/misha/reference/ggenome.transplant.md),
[`gdb.export_fasta`](https://tanaylab.github.io/misha/reference/gdb.export_fasta.md),
[`gseq.extract`](https://tanaylab.github.io/misha/reference/gseq.extract.md),
[`gdb.create`](https://tanaylab.github.io/misha/reference/gdb.create.md)

## Examples

``` r
gdb.init_examples()

# Export the example DB to a reference FASTA
ref_fasta <- tempfile(fileext = ".fa")
gdb.export_fasta(ref_fasta)

# Replace two regions with literal sequences
intervals <- data.frame(
    chrom = c("chr1", "chr1"),
    start = c(100, 200),
    end = c(110, 210)
)
donors <- c("AAAAAAAAAA", "CCCCCCCCCC")
out <- tempfile(fileext = ".fa")
trackdb <- tempfile()
ggenome.implant(intervals, donors,
    output = out,
    genome_fasta = ref_fasta,
    create_trackdb = TRUE,
    trackdb_path = trackdb
)

# Verify the implanted sequences via the new trackdb
gdb.init(trackdb)
gseq.extract(data.frame(chrom = "chr1", start = 100, end = 110))
#> [1] "AAAAAAAAAA"

# Clean up
unlink(c(ref_fasta, out, paste0(out, ".fai"), trackdb), recursive = TRUE)
```
