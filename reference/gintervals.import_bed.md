# Import intervals from a BED file

Reads a BED file, plain or gzipped, and returns a misha 1D intervals
data frame. Track/browser/comment header lines are skipped
automatically. Chromosome names are normalized through the active
database's `CHROM_ALIAS` mechanism (so `chr1` \<-\> `1` works without
explicit configuration).

## Usage

``` r
gintervals.import_bed(
  file = NULL,
  name = TRUE,
  score = TRUE,
  strand = TRUE,
  keep_order = FALSE,
  format = NULL
)
```

## Arguments

- file:

  path to a BED file (`.bed` or `.bed.gz`) or a bigBed file. Zip
  archives are not supported - unzip them first.

- name:

  if `TRUE` and a 4th column exists, include it as `name`.

- score:

  if `TRUE` and a 5th (numeric) column exists, include it as `score`.

- strand:

  if `TRUE` and a 6th column exists, include it as `strand` (mapped to
  `1`/`-1`/`0`).

- keep_order:

  if `TRUE`, return the intervals in the file's order (header and
  comment lines skipped) instead of sorting them.

- format:

  one of `"bed"`, `"narrowPeak"` or `"broadPeak"`. `NULL` guesses from
  the file name: a name ending in `.narrowPeak` or `.narrowPeak.gz` (any
  case) is narrowPeak, `.broadPeak` or `.broadPeak.gz` is broadPeak,
  anything else is BED.

## Value

A 1D intervals data frame, sorted by chrom and start unless `keep_order`
is `TRUE`.

## Details

BED is already 0-based half-open, so coordinates are taken as-is.

A bigBed file (recognized by its first bytes, whatever its name) is read
natively. Its columns are chrom, start, end and then its other fields in
the file's order, named by the file's autoSql: fields of a numeric
autoSql type (int, uint, short, ushort, byte, ubyte, float, double,
bigint) are numeric, the others (string, lstring, char\[n\], enum, set,
arrays) are kept as the text stored, and so is itemRgb (BED field 9,
"R,G,B" text although autoSql declares it a uint). `name`, `score` and
`strand` apply to the fields of those names. The fields of a bigBed
without autoSql get the names UCSC's tools give them (name, score,
strand, thickStart, thickEnd, reserved, blockCount, blockSizes,
chromStarts, ...). As in
[`gtrack.import`](https://tanaylab.github.io/misha/reference/gtrack.import.md)
of a bigWig, records on chromosomes the database lacks are skipped with
a message (a warning if one of them looks like a primary chromosome),
and records past the end of their chromosome are clipped to it (dropped
if they start past it) with a warning. `format` must be `NULL` for a
bigBed.

ENCODE narrowPeak (BED6+4) and broadPeak (BED6+3) files carry extra
columns, which are added after the BED columns under their UCSC names:
`signalValue`, `pValue` and `qValue` (numeric) for both, and `peak`
(integer) for narrowPeak. `peak` is the summit's offset from `start`; -1
means no summit was called, and is kept as-is.

By default the intervals are sorted by chrom and start. With
`keep_order = TRUE` they stay in the file's order, so row `i` is the
`i`-th data line, as peak IDs that are row numbers expect. The functions
below accept unsorted intervals. `gextract`, `gscreen`,
`gintervals.canonic` and the set operations return their rows in genomic
order (`gextract`'s `intervalID` column gives the input row), while
`gseq.extract`, `gintervals.annotate`, `gintervals.normalize` and
`gintervals.summary` return one row per input row in input order.
[`gintervals()`](https://tanaylab.github.io/misha/reference/gintervals.md)
sorts its input.

## See also

[`gintervals.import_gff`](https://tanaylab.github.io/misha/reference/gintervals.import_gff.md),
[`gintervals.import_vcf`](https://tanaylab.github.io/misha/reference/gintervals.import_vcf.md).

## Examples

``` r

gdb.init_examples()
f <- tempfile(fileext = ".narrowPeak")
writeLines(c(
    "chr1\t20000\t20400\tpeak_1\t500\t.\t8.1\t12.3\t10.2\t150",
    "chr1\t10000\t10300\tpeak_2\t800\t.\t9.4\t15.0\t12.8\t120"
), f)
peaks <- gintervals.import_bed(f, keep_order = TRUE)
peaks
#>   chrom start   end strand   name score signalValue pValue qValue peak
#> 1  chr1 20000 20400      0 peak_1   500         8.1   12.3   10.2  150
#> 2  chr1 10000 10300      0 peak_2   800         9.4   15.0   12.8  120

# 200 bp windows centred on the summits, in the file's order: the summit
# is at (start + end) %/% 2. Drop rows with peak = -1 (no summit) first.
summit <- peaks$start + peaks$peak
windows <- transform(peaks, start = summit - 100, end = summit + 100, peak = 100L)
windows
#>   chrom start   end strand   name score signalValue pValue qValue peak
#> 1  chr1 20050 20250      0 peak_1   500         8.1   12.3   10.2  100
#> 2  chr1 10020 10220      0 peak_2   800         9.4   15.0   12.8  100
```
