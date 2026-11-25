# Creates tracks from the tracks of a figwig file

Creates one misha track per selected track of a figwig file ('.fw', the
multi-track format of the figwig package,
<https://github.com/jmschrei/figwig>). The file is read by misha itself;
neither Python nor figwig is needed.

## Usage

``` r
gtrack.import_figwig(
  file,
  binsize = 0,
  track.prefix = NULL,
  tracks = NULL,
  select = NULL,
  names = NULL,
  defval = NaN,
  description = "",
  attrs = NULL
)
```

## Arguments

- file:

  path of the figwig file

- binsize:

  bin size of 'Dense' tracks, or 0 (or NULL) for 'Sparse' tracks

- track.prefix:

  prefix of the track names

- tracks:

  figwig track names or 1-based indices to import, or NULL for all

- select:

  a named list of attribute values to choose tracks by

- names:

  names of the created tracks (before the prefix), instead of the
  sanitized figwig names

- defval:

  value of the bases with no value in 'Dense' tracks

- description:

  description of every created track; when empty, a figwig 'description'
  attribute is kept

- attrs:

  a named vector or list of attributes set on every created track, over
  those of the file

## Value

Names of the created tracks, invisibly.

## Details

With 'binsize' 0 (or NULL) each track is 'Sparse': one interval per
maximal run of bases with the same value, leaving out the bases whose
value is 0 or NaN. That suits tracks that are mostly zero, such as
counts or cut sites; a dense float track makes one interval per base.
With a positive 'binsize' each track is 'Dense', its bins computed as
[`gtrack.import`](https://tanaylab.github.io/misha/reference/gtrack.import.md)
computes the bins of a bigWig: the mean of the values of the bin's
bases, each run of equal values weighted by its length. A figwig file
has no "no data": a base that was never written is 0, as is one written
as 0 (figwig's 'bw2fw' writes a bigWig's uncovered bases as 0 unless
given '–missing nan'). NaN bases, bases past the end of the file's
chromosome and chromosomes the file does not have have no value: they
count as 'defval' unless it is NaN, and a bin with no value is NaN.

Tracks are chosen with 'tracks' (figwig track names, or 1-based
indices), in the order given, or all of them in the file's order, and
then with 'select', a named list of attribute values (a vector of
several values matches any of them): a track is imported when every
attribute given matches. Numbers are compared by value, whether given as
numbers or as text, so `list(rep = 1)`, `list(rep = "1")` and
`list(rep = "1.0")` all match a `"rep": 1` or `"rep": 1.0` record; other
values are compared as text. 'name', 'dtype' and 'group' can be selected
on too.

The track made from figwig track 'n' is named '\[track.prefix\]n2',
where 'n2' is the entry of 'names' for it if given (a character vector,
one per imported track, or one named by figwig track names) and
otherwise 'n' sanitized: a trailing '+' or '-' (a strand) is spelled
'plus' or 'minus', each run of characters other than ASCII letters,
digits and '\_' becomes one '\_', and without a prefix a name that does
not start with a letter gets an 'X' in front ('atac.+' becomes
'atac_plus', 'H3K27ac-rep1' 'H3K27ac_rep1'). Dots in 'track.prefix' are
directories, which must exist. All the names are checked before anything
is written, and existing tracks are not overwritten.

Each figwig track's attributes become attributes of its track, except
'name', 'dtype' and 'group'; 'created.by', 'created.date' and
'created.user', which the import sets; 'type' and 'binsize'; and
'description', unless the 'description' given is not empty. A string
attribute is kept as it is, a null one is left out, and a number,
boolean, list or object becomes its JSON text. 'attrs' win over them.

Chromosome names are matched through the chromosome alias chain as in
[`gtrack.import`](https://tanaylab.github.io/misha/reference/gtrack.import.md),
with the same reporting of names the database does not have; a
chromosome whose length differs from the database's is warned about. The
import is all or nothing: if it fails or is interrupted, none of the
tracks is created.

Every layout of the figwig format version 1 is read: position and track
order, stored or zlib-compressed chunks, dense, sparse, runs and matrix
chunks and every element type. PyMisha's `gtrack_import_figwig` takes
the same arguments and makes the same tracks.

## See also

[`gtrack.import`](https://tanaylab.github.io/misha/reference/gtrack.import.md),
[`gtrack.import_parquet`](https://tanaylab.github.io/misha/reference/gtrack.import_parquet.md),
[`gtrack.attr.get`](https://tanaylab.github.io/misha/reference/gtrack.attr.get.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# every track of the file, as sparse tracks named after the figwig tracks
gtrack.import_figwig("tracks.fw")

# the ATAC-seq tracks, as dense tracks of 20 bp bins
gtrack.import_figwig("tracks.fw",
    binsize = 20, track.prefix = "fw_",
    select = list(assay = "ATAC-seq")
)
} # }
```
