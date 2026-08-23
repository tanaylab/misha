# Creates 1D intervals from coordinate strings

Creates a set of 1D intervals by parsing UCSC-style coordinate strings.

## Usage

``` r
gintervals.from_strings(regions = NULL)
```

## Arguments

- regions:

  a character vector of region strings. Accepted formats per element:
  `"chrom"`, `"chrom:start-end"`, `"chrom:start..end"`, optionally
  followed by `":+"` or `":-"`.

## Value

A data frame representing the intervals, sorted as by
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md).
A "strand" column is added only when at least one region specifies a
strand.

## Details

Parses strings of the form `"chrom:start-end"` (`"chrom:start..end"` is
also accepted) into a 1D intervals data frame. A chromosome-only string
such as `"chr1"` expands to the full chromosome extent. An optional
trailing strand suffix `":+"` or `":-"` adds a "strand" column.

Coordinates are zero-based and half-open, matching the convention used
by
[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md).

## See also

[`gintervals`](https://tanaylab.github.io/misha/reference/gintervals.md),
[`gintervals.2d`](https://tanaylab.github.io/misha/reference/gintervals.2d.md)

## Examples

``` r

gdb.init_examples()
gintervals.from_strings("chr1:1000-2000")
#>   chrom start  end
#> 1  chr1  1000 2000
gintervals.from_strings(c("chr1:1000-2000:+", "chr2:50-60:-"))
#>   chrom start  end strand
#> 1  chr1  1000 2000      1
#> 2  chr2    50   60     -1
```
