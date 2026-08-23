# Creates a track from a file of mapped sequences

Creates a track from a file of mapped sequences.

## Usage

``` r
gtrack.import_mappedseq(
  track = NULL,
  description = NULL,
  file = NULL,
  pileup = 0,
  binsize = -1,
  cols.order = c(9, 11, 13, 14),
  remove.dups = TRUE
)
```

## Arguments

- track:

  track name

- description:

  a character string description

- file:

  path to the mapped sequences file (see Description for accepted
  formats)

- pileup:

  interval expansion

- binsize:

  bin size of a dense track

- cols.order:

  order of sequence, chromosome, coordinate and strand columns in mapped
  sequences file or NULL if SAM file is used. For BAM input the only
  accepted values are \`NULL\` or omission; an explicit non-NULL
  \`cols.order\` with BAM input is an error.

- remove.dups:

  if 'TRUE' the duplicated coordinates are counted only once.

## Value

A list of conversion process statistics.

## Details

This function creates a track from a file of mapped sequences. The file
can be in SAM format, in a general TAB delimited text format where each
line describes a single read, in gzipped variants of either
(\`.sam.gz\`, \`.tsv.gz\`), or in BAM format (auto-detected by bgzip
magic; requires \`samtools\` on \`PATH\`).

For a SAM file 'cols.order' must be set to 'NULL'. For BAM input the
default \`cols.order = c(9, 11, 13, 14)\` is treated as SAM mode because
\`samtools view\` emits SAM-format payload; passing a non-NULL
\`cols.order\` explicitly with BAM input is an error.

For a general TAB delimited text format the following columns must be
presented in the file: sequence, chromosome, coordinate and strand. The
position of these columns should be specified in 'cols.order' argument.
The default value of 'cols.order' is an array of (9, 11, 13, 14) meaning
that sequence is expected to be found at column number 9, chromosome -
at column 11, coordinate - at column 13 and strand - at column 14. The
column indices are 1-based, i.e. the first column is referenced by 1.
Chromosome needs a prefix 'chr' e.g. 'chr1'. Valid strand values are '+'
or 'F' for forward strand and '-' or 'R' for the reverse strand.

Each read at given coordinate can be "expanded" to cover an interval
rather than a single point. The length of the interval is controlled by
'pileup' argument. The direction of expansion depends on the strand
value. If 'pileup' is '0', no expansion is performed and the read is
converted to a single point. The track is created in sparse format. If
'pileup' is greater than zero, the output track is in dense format.
'binsize' controls the bin size of the dense track.

If 'remove.dups' is 'TRUE' the duplicated coordinates are counted only
once.

'description' is added as a track attribute.

'gtrack.import_mappedseq' returns the statistics of the conversion
process.

## See also

[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md)
