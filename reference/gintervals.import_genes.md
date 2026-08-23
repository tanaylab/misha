# Imports genes and annotations from files

Imports genes and annotations from files.

## Usage

``` r
gintervals.import_genes(
  genes.file = NULL,
  annots.file = NULL,
  annots.names = NULL
)
```

## Arguments

- genes.file:

  name or URL of file that contains genes

- annots.file:

  name of URL file that contains annotations. If 'NULL' no annotations
  are imported

- annots.names:

  annotations names

## Value

A list of four intervals sets named 'tss', 'exons', 'utr3' and 'utr5'.
'strand' column and annotations are attached to the intevals.

## Details

This function reads a definition of genes from 'genes.file' and returns
four sets of intervals: TSS, exons, 3utr and 5utr. In addition to the
regular intervals columns 'strand' column is added. It contains '1'
values for '+' strands and '-1' values for '-' strands.

If annotation file 'annots.file' is given then annotations are attached
too to the intervals. The names of the annotations as they would appear
in the return value must be specified in 'annots.names' argument.

Both 'genes.file' and 'annots.file' can be either a file path or URL in
a form of 'ftp://\[address\]/\[file\]'. Files that these arguments point
to can be zipped or unzipped.

Examples of 'genes.file' and 'annots.file' can be found here:

`ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz`
`ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz`

If a few intervals overlap (for example: two TSS regions) they are all
unified to an interval that covers the whole overlapping region.
'strand' value is set to '0' if two or more of the overlapping intervals
have different strands. The annotations of the overlapping intervals are
concatenated to a single character string separated by semicolons.
Identical values of overlapping intervals' annotation are eliminated.

## See also

[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gdb.create`](https://tanaylab.github.io/misha/reference/gdb.create.md)
