# Misha Basics (Short Guide)

This page gives a compact mental model for misha. Use it as the first
quick read before the full `Manual` vignette.

## The Core Idea

Most analyses follow the same pattern:

1.  Choose **where** to look (intervals / scope).
2.  Choose **how** to walk through it (iterator).
3.  Evaluate a **track expression** over those iterator intervals.

In misha this is usually one call to `gextract`, `gscreen`, or
`gsummary`.

You are not limited to raw track names. You can pass full expressions,
for example `log(dense_track + 1)`, `dense_track / (chip.sum + 1e-6)`,
or `pmin(dense_track, 2)`.

All examples below assume the bundled examples database:

``` r

library(misha)
gdb.init_examples()
```

## Four Concepts You Need First

### 1) Track

A **track** is genomic signal organized over coordinates.

- Dense track: value for each fixed-size bin (for example `dense_track`
  in the examples DB).
- Sparse track: values on intervals (for example peaks).
- 2D track: values on genomic rectangles (for example contact matrices).

Useful starter commands:

``` r

gtrack.ls() # list tracks in the examples DB
#> [1] "array_track"         "dense_track"         "rects_track"        
#> [4] "sparse_track"        "subdir.dense_track2"
gtrack.info("dense_track") # inspect type/metadata
#> $type
#> [1] "dense"
#> 
#> $dimensions
#> [1] 1
#> 
#> $size.in.bytes
#> [1] 80012
#> 
#> $format
#> [1] "per-chromosome"
#> 
#> $bin.size
#> [1] 50
gtrack.info("sparse_track")
#> $type
#> [1] "sparse"
#> 
#> $dimensions
#> [1] 1
#> 
#> $size.in.bytes
#> [1] 32012
#> 
#> $format
#> [1] "per-chromosome"
```

For intuition, you can think of `dense_track` as a ChIP-seq-like
coverage signal.

### 2) Intervals

An **interval set** defines genomic regions (`chrom`, `start`, `end`)
where you want to work.

- Intervals can come from files, annotations, peak calls, or be
  generated in code.
- Intervals often act as a **scope**: “analyze only here.”

``` r

regions <- gintervals(1, c(0, 250000), c(100000, 260000))
regions
#>   chrom  start    end
#> 1  chr1      0 100000
#> 2  chr1 250000 260000
```

### 3) Iterator

The **iterator** is the stepping policy inside the scope.

- `iterator = 100` -\> fixed 100 bp bins
- `iterator = "some_sparse_track"` -\> iterate over that track’s
  intervals
- `iterator = some_intervals_df` -\> iterate over explicit regions
- `iterator = "my_intervals_set"` -\> iterate directly over an intervals
  set

Think of it as: scope says *where*, iterator says *in what chunks*.

``` r

out <- gextract("dense_track", regions, iterator = 100)
head(out)
#>   chrom start end dense_track intervalID
#> 1  chr1     0 100   0.1688889          1
#> 2  chr1   100 200   0.1700000          1
#> 3  chr1   200 300   0.1800000          1
#> 4  chr1   300 400   0.1600000          1
#> 5  chr1   400 500   0.1100000          1
#> 6  chr1   500 600   0.0400000          1
log_out <- gextract("log(dense_track + 1)", regions, iterator = 100)
head(log_out)
#>   chrom start end log(dense_track + 1) intervalID
#> 1  chr1     0 100           0.15605364          1
#> 2  chr1   100 200           0.15700375          1
#> 3  chr1   200 300           0.16551444          1
#> 4  chr1   300 400           0.14842000          1
#> 5  chr1   400 500           0.10436001          1
#> 6  chr1   500 600           0.03922071          1
```

Create and use an intervals set as an iterator:

``` r

gintervals.save("my_intervals_set", regions) # name first, then the intervals
out2 <- gextract("dense_track", gintervals.all(), iterator = "my_intervals_set")
out2
#>   chrom  start    end dense_track intervalID
#> 1  chr1      0 100000   0.1011889          1
#> 2  chr1 250000 260000   0.1263158          1
```

### 4) Virtual Track

A **virtual track** is a named on-the-fly transformation, not stored as
a physical track file.

Examples:

- Local sum of a source track
- Distance to nearest annotation interval
- Quantile-like or nearest-neighbor summaries

``` r

gvtrack.create("chip.sum", "dense_track", "sum")
out <- gextract("chip.sum", regions, iterator = 200)
head(out)
#>   chrom start  end   chip.sum intervalID
#> 1  chr1     0  200 0.67777777          1
#> 2  chr1   200  400 0.68000001          1
#> 3  chr1   400  600 0.29999998          1
#> 4  chr1   600  800 0.04000000          1
#> 5  chr1   800 1000 0.09999999          1
#> 6  chr1  1000 1200 0.12000000          1
```

You can also shift the iterator window used by the virtual track:

``` r

gvtrack.create("chip.shifted", "dense_track", "sum")
gvtrack.iterator("chip.shifted", sshift = -100, eshift = 100)
out <- gextract("chip.shifted", regions, iterator = 200)
head(out)
#>   chrom start  end chip.shifted intervalID
#> 1  chr1     0  200     1.037778          1
#> 2  chr1   200  400     1.240000          1
#> 3  chr1   400  600     0.660000          1
#> 4  chr1   600  800     0.140000          1
#> 5  chr1   800 1000     0.200000          1
#> 6  chr1  1000 1200     0.220000          1
```

Here, each iterator interval is expanded by 100 bp on both sides before
evaluating `dense_track`.

Virtual tracks are session objects (easy to list with `gvtrack.ls` and
delete with `gvtrack.rm`).

## Minimal Workflow

``` r

library(misha)
gdb.init_examples()

# 1) pick scope
regions <- gintervals(1, 0, 50000)

# 2) inspect available tracks
gtrack.ls()
#> [1] "array_track"         "dense_track"         "rects_track"        
#> [4] "sparse_track"        "subdir.dense_track2"

# 3) extract signal with a chosen iterator
chip <- gextract("dense_track", regions, iterator = 100)
head(chip)
#>   chrom start end dense_track intervalID
#> 1  chr1     0 100   0.1688889          1
#> 2  chr1   100 200   0.1700000          1
#> 3  chr1   200 300   0.1800000          1
#> 4  chr1   300 400   0.1600000          1
#> 5  chr1   400 500   0.1100000          1
#> 6  chr1   500 600   0.0400000          1

# 4) screen high-signal bins (as a simple peak-like filter).
#    Pick the threshold from the data: dense_track only reaches 0.34 in this
#    scope, so a threshold above that would silently screen nothing.
hi_chip <- gscreen("dense_track > 0.2", regions, iterator = 100)
head(hi_chip)
#>   chrom start   end
#> 1  chr1 17200 17300
#> 2  chr1 20000 20100
#> 3  chr1 23300 23400
#> 4  chr1 26200 26300
#> 5  chr1 32600 32800
#> 6  chr1 32900 33000

# 5) summarize distribution/coverage
gsummary("dense_track", regions, iterator = 100)
#> Total intervals   NaN intervals             Min             Max             Sum 
#>    500.00000000      0.00000000      0.00000000      0.34000000     43.68888862 
#>            Mean         Std dev 
#>      0.08737778      0.06208607
```

## PWM in One Minute

A PWM/PSSM is a motif model over A/C/G/T. In misha, a common pattern is:

1.  Extract sequence from intervals.
2.  Score those sequences with a PWM.

``` r

regions <- gintervals(1, c(1000, 2000), c(1020, 2020))
seqs <- gseq.extract(regions)
seqs
#> [1] "cctcagtaatccgaaaagcc" "CTGCATGTAACTTAATACCA"

pssm <- matrix(c(
    0.80, 0.05, 0.10, 0.05,
    0.10, 0.10, 0.70, 0.10,
    0.05, 0.80, 0.05, 0.10,
    0.10, 0.10, 0.10, 0.70
), ncol = 4, byrow = TRUE)
colnames(pssm) <- c("A", "C", "G", "T")

scores <- gseq.pwm(seqs, pssm, mode = "lse")
scores
#> [1] -2.198000 -2.119781
```

If your database has motif files under `pssms/`, you can create a
genome-wide PWM-energy track with `gtrack.create_pwm_energy(...)`.
