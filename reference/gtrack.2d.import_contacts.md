# Creates a track from a file of inter-genomic contacts

Creates a track from a file of inter-genomic contacts.

## Usage

``` r
gtrack.2d.import_contacts(
  track = NULL,
  description = NULL,
  contacts = NULL,
  fends = NULL,
  allow.duplicates = TRUE
)
```

## Arguments

- track:

  track name

- description:

  a character string description

- contacts:

  vector of contacts files

- fends:

  name of fragment ends file

- allow.duplicates:

  if 'TRUE' duplicated contacts are allowed

## Value

None.

## Details

This function creates a 'Points' (two-dimensional) track from contacts
files. If 'allow.duplicates' is 'TRUE' duplicated contacts are allowed
and summed up, otherwise an error is reported.

Contacts (coord1, coord2) within the same chromosome are automatically
doubled to include also '(coord2, coord1)' unless 'coord1' equals to
'coord2'.

Contacts may come in one or more files.

If 'gmultitasking' is on (the default), the input files are read by
several processes in parallel, each process reading a contiguous subset
of the files, and the chromosome pairs are then written in parallel, one
pair per process at a time. A single input file is read by one process,
so split a large input into several files to read it in parallel. The
number of processes is bounded by 'gmax.processes'. The resulting track
is identical to the one created with 'gmultitasking' off. A pair is
started only when its estimated memory, added to that of the pairs being
written, stays within 'gmax.mem.usage'. Larger pairs are taken first,
but a smaller pair that fits may start while a larger one waits, and a
pair estimated above 'gmax.mem.usage' is written alone. A pair with at
least 4 x 'gmax.data.size' contacts (2 x for a chromosome with itself,
whose contacts are stored mirrored) is written in parts and estimated by
its largest part, so lowering 'gmax.data.size' lets more pairs be
written at once. Lower 'gmax.mem.usage' to use less memory. Reading the
input is not covered by the estimate.

If 'fends' is 'NULL' contacts file is expected to be in
"intervals-value" tab-separated format. The file starts with a header
defining the column names. The first 6 columns must have the following
names: 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'. The last
column is designated for the value and it may have an arbitrary name.
The header is followed by a list of intervals and a value for each
interval. An interval of form (chrom1, start1, end1, chrom2, start2,
end2) is added as a point (X, Y) to the resulted track where X =
(start1 + end1) / 2 and Y = (start2 + end2) / 2.

One can see an example of "intervals-value" format by running 'gextract'
function on a 2D track with a 'file' parameter set to the name of the
file.

If 'fends' is not 'NULL' contacts file is expected to be in
"fends-value" tab-separated format. It should start with a header
containing at least 3 column names 'fend1', 'fend2' and 'count' in
arbitrary order followed by lines each defining a contact between two
fragment ends.

|        |         |                                   |
|--------|---------|-----------------------------------|
| COLUMN | VALUE   | DESCRIPTION                       |
| fend1  | Integer | ID of the first fragment end      |
| fend2  | Integer | ID of the second fragment end     |
| count  | Numeric | Value associated with the contact |

A fragment ends file is also in tab-separated format. It should start
with a header containing at least 3 column names 'fend', 'chr' and
'coord' in arbitrary order followed by lines each defining a single
fragment end.

|  |  |  |
|----|----|----|
| COLUMN | VALUE | DESCRIPTION |
| fend | Unique integer | ID of the fragment end |
| chr | Chromosome name | Can be specified with or without "chr" prefix, like: "X" or "chrX" |
| coord | Integer | Coordinate |

'description' is added as a track attribute.

Note: temporary files are created in the directory of the track during
the run of the function. A few of them need to be kept simultaneously
open. If the number of chromosomes and / or contacts is particularly
high, a few thousands files might be needed to be opened simultaneously.
Some operating systems limit the number of open files per user, in which
case the function might fail with "Too many open files" or similar
error. The workaround could be:

1\. Increase the limit of simultaneously opened files (the way varies
depending on your operating system). 2. Increase the value of
'gmax.data.size' option. Higher values of 'gmax.data.size' option will
increased memory usage of the function but create fewer temporary files.

## See also

[`gtrack.2d.import`](https://tanaylab.github.io/misha/reference/gtrack.2d.import.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md)
