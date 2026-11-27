# misha - core operations

Compact reference for the everyday misha workflow. Concepts first, then chooser tables, then one recipe per common task. Each recipe ends with short `Avoid:` callouts; the full anti-pattern catalogue is in [`misha-anti-patterns.md`](misha-anti-patterns.md). Advanced topics (2D / Hi-C, PWM, import/export, new genomes) live in [`misha-advanced.md`](misha-advanced.md).

## Contents

- [1. Concepts](#1-concepts)
- [2. Session bootstrap](#2-session-bootstrap)
- [3. Function chooser tables](#3-function-chooser-tables)
- [4. Recipes](#4-recipes)
  - [4.1 Build and manipulate intervals](#41-build-and-manipulate-intervals)
  - [4.2 Annotate intervals with nearby features](#42-annotate-intervals-with-nearby-features)
  - [4.3 Distance flavors - picking the right primitive](#43-distance-flavors)
  - [4.4 Extract per-region values](#44-extract-per-region-values)
  - [4.5 Aggregate signal in a flanking window](#45-aggregate-signal-in-a-flanking-window)
  - [4.6 Call regions of interest (gscreen + thresholding)](#46-call-regions-of-interest)
  - [4.7 Joint distributions, summaries, correlations](#47-joint-distributions-summaries-correlations)
  - [4.8 Materialize derived tracks](#48-materialize-derived-tracks)

## 1. Concepts

**Mental model.** Think of every misha query as a *walk across the genome*: an iterator visits one position (or one interval, or one pair of positions) at a time, in genomic order, gathering values from any tracks or virtual tracks you ask for at each stop. You don't load tracks into R and slice them; you describe *where to walk* (`intervals=`), *what step size* (`iterator=`), and *what to collect at each step* (the expression / vtracks). The engine streams the genome once, top to bottom, in C++. Almost everything else in this guide is a variation on that idea.

### 1.1 Track database

A *track DB* is a directory tree on disk under `<GROOT>`. Each subdirectory is a namespace; each leaf is a **track**. Track names use dots, mirroring the layout: `epi.k27me3.es` lives at `<GROOT>/epi/k27me3/es/`. Attach with `gsetroot(<path>)`. Tracks are not loaded into R; every query streams off disk in C++. The common project-portability convention is a relative symlink at the project root (`db -> /full/path/to/groot`) so `gsetroot("db")` works from any check-out. Project-local intervals can live in a separate **dataset** attached with `gdataset.load(<path>)`.

### 1.2 Tracks

Three storage shapes:

- **Dense** - value per bin (typical 20bp). `gtrack.create`, `gtrack.create_dense`. Continuous signal.
- **Sparse** - values at named intervals, NA elsewhere. `gtrack.create_sparse`. Peaks, per-CpG methylation.
- **2D** - value per rectangle in `(chrom1×coord1) × (chrom2×coord2)`. Hi-C / capture-C. Covered in [misha-advanced.md](misha-advanced.md).

`gtrack.ls("pattern", perl=TRUE)` lists tracks; `gtrack.exists(name)` is the idempotency guard. If a track you know exists on disk doesn't appear in `gtrack.ls` (for example, someone else added it after your session attached the DB, or you just dropped a file in by hand), run `gdb.reload()` to rescan and refresh the cached track / intervals index.

**On-disk format.** Older trackdbs store one file per chromosome inside each track's directory. Newer trackdbs use an **indexed** format - `track.dat` + `track.idx` per track - which scales to thousands of contigs without dragging the filesystem. The R API hides this entirely: `gextract`, `gscreen`, etc. work on both formats. **Do not reason from per-chromosome files on disk; always go through the API.** Convert an existing trackdb (or a single track) to the indexed form with `gdb.convert_to_indexed(convert_tracks = TRUE, convert_intervals = TRUE)` (both flags default to `FALSE`, so the bare call is a no-op), `gtrack.convert_to_indexed(track)`, or `gtrack.2d.convert_to_indexed(track, ...)` - important for fragmented assemblies (Zoonomia-style multi-contig genomes) where the per-chromosome layout is unworkable.

### 1.3 Intervals

A 1D intervals frame has columns `chrom, start, end` (+ optional `strand` and user columns); 2D has `chrom1, start1, end1, chrom2, start2, end2`. Build with `gintervals(...)` / `gintervals.2d(...)`. Genome-wide scopes: `gintervals.all()` (1D), `gintervals.2d.all()` (2D); whole-chrom 2D shorthand is `gintervals.2d(chrom, 0, -1, chrom, 0, -1)`.

Canonical sets shipped per genome:
- **Annotation**: `intervs.global.tss`, `intervs.global.exon`, `intervs.global.tad_names`, `intervs.global.rmsk_<class>` (LINE, SINE, LTR, ...).
- **Sequence content** (extract via vtracks for GC / CpG features): `seq.GC_500_mean`, `seq.CG_500_mean`, `seq.GC500_bin20`, `seq.G_or_C`, `seq.CG`. The `_<W>_mean` form is a windowed average; `_<W>_bin20` is binned for fast threshold queries.

### 1.4 Virtual tracks (vtracks)

A vtrack is a parameterized view onto a track. One call composes aggregator and window:

```r
gvtrack.create("k27_flank", src="epi.k27me3.es", func="sum",
               sshift = -2000, eshift = 2000)
```

Common aggregators (`func=`): `sum`, `area`, `avg`, `max`, `min`, `lse` (log-sum-exp), `stdev`, `quantile`, `global.percentile.max`, `weighted.sum`. Non-track sources: `distance` / `distance.center` / `distance.edge` (three distance-to-nearest-interval flavors - see [§4.3](#43-distance-flavors)), `pwm` / `pwm.max` / `pwm.edit_distance` (PSSM energies - see [misha-advanced.md](misha-advanced.md)), `kmer.frac` (sequence content).

For 2D: window args become `sshift1=, eshift1=, sshift2=, eshift2=`. To bind a 1D vtrack to one axis of a 2D iteration: `dim = 1` or `dim = 2`.

Re-creating a vtrack name silently overwrites - no `gvtrack.rm` needed first.

### 1.5 Iterators

Every query takes `iterator =` defining what one row of the answer is:

- `iterator = 20` - tile the scope at 20bp. One row per bin. *Chromosome-relative*, not interval-relative.
- `iterator = "some.track"` - use the track's stored bin grid (one row per stored position; common with CpG tracks).
- `iterator = <intervals frame>` - one row per interval. With `intervals = <same frame>`, this is the universal "one row per peak" recipe.
- `iterator = giterator.cartesian_grid(...)` - 2D pair grid (see [misha-advanced.md](misha-advanced.md)).

`intervals=` is the *scope*; `iterator=` is the *resolution*. They are independent.

**`intervals =` is mandatory for most query functions.** `gextract`, `gscreen`, `gsummary`, `gdist`, `gintervals.neighbors`, `gintervals.annotate`, etc. don't default to "whole genome" - if you want the walk to cover everything, pass `intervals = gintervals.all()` (1D) or `intervals = gintervals.2d.all()` (2D) explicitly. The same applies to `iterator =` when there's no natural per-row resolution.

### 1.6 Two rules that apply everywhere

**Rule 1 - prefer a virtual track over a hand-rolled computation.** Before round-tripping through R for any per-bin or per-region computation (sums, percentiles, distances, kmer counts, PWM scores, log-sum-exp aggregates), check whether a vtrack `func =` already does it. Virtual tracks evaluate in C++ at the engine level and stay inside the single-genome-pass model that makes misha fast.

**Rule 2 - unify into one `gextract` call.** When you need several values per region (or per bin), pass *all* the expressions to a single `gextract` as a `c("expr1", "expr2", ...)` vector. The engine evaluates them jointly in one genome pass; N separate calls are N× slower and require a manual join on chrom/start/end (or `intervalID`) afterwards.

## 2. Session bootstrap

```r
library(misha)                       # for installed misha
# devtools::load_all(export_all = FALSE)  # when working in the misha source tree

gsetroot("db")                       # or absolute path; symlink "db" from project root
gdataset.load("db_extra")            # optional, project-local dataset

options(gmax.data.size = 1e10)       # raise from default for large gextract returns
options(gmultitasking  = TRUE)       # parallel across chroms (default; explicit)
options(scipen         = 1e3)        # non-scientific notation when pasting coordinates
```

Attach the genome you need once at the top and stay on it. Switching `gsetroot` mid-script while intervals from the previous genome are still in scope is the most reliable way to produce silent coordinate drift.

**`gmultitasking` and external parallelism.** Misha's `gmultitasking = TRUE` forks worker processes internally. When you then wrap misha calls in `mclapply` (or other forking parallelism), the nested forks can deadlock or oversubscribe cores. Set `options(gmultitasking = FALSE)` inside the `mclapply` body (or under `gcluster.run`) and let the outer layer do the parallelism.

**Multitasked values.** A multitasked call returns the same values as a call in one process, except the sums of `gsummary` and Pearson's `gcor` (merged across the processes, last-digit differences) and the sampled results of `gquantiles`, `gintervals.quantiles` and `gcor(method = "spearman")`: each of their processes samples at most `gmax.data.size / processes` values, so past that their results are approximate even where one process would give them exactly (the two quantile functions warn, Spearman does not). For exact results raise `gmax.data.size`, set `options(gmultitasking = FALSE)` for the call, or use `method = "spearman.exact"`.

**Memory balance.** `gmax.data.size` caps the *result size* (rows) of any single `gextract`, all its processes together. It does not change how much memory a result under it takes; lowering it only makes a larger call fail early. With `gmultitasking = TRUE` (the "tiles" strategy) each process collects its rows in its own memory and then copies them to shared memory, which the parent turns into the result, freeing it as it goes. While the processes copy, a large result can take up to about twice its size, which can surface as `cannot allocate memory`. If you hit it: write the result to a file (`gextract(..., file = "x.tsv")` or `file = "x.parquet"`), drop `gmultitasking` for the call, or partition by chromosome.

**Many-track strategy.** `options(gmultitasking.strategy = "auto" | "tracks" | "tiles")` controls how the engine parallelizes a multi-track `gextract`. `"tiles"` (the historical default) parallelizes across genome tiles within one expression at a time - good for a few heavy tracks. `"tracks"` parallelizes across the *expressions* in `c("...", "...")` instead - dramatically faster for thousands of motif / feature tracks where each is cheap individually. `"auto"` picks per call. For motif-scan workloads (hundreds of PSSMs over the genome) set `"tracks"` explicitly.

**Avoid:**
- `options(gmutitasking = FALSE)` - typo of `gmultitasking`. R's `options()` silently accepts unknown names; this is a no-op.
- `options(gparam.type = 'string')` - vestigial legacy option; current misha accepts both quoted-string and bare-symbol track names.

## 3. Function chooser tables

Pick the verb by *what's returned*, not by alphabetical proximity.

### Workhorse verbs

| Verb | Returns | When |
|---|---|---|
| `gextract` | data.frame | Values to look at, plot, or join. |
| `gextract.matrix` | N x L matrix (N x L x E array) | Per-bin values over equal-width windows, one row per window in input order. |
| `gscreen`  | intervals frame | Filter to a region set (peaks above threshold). |
| `gquantiles` | quantile cutoffs | Pick a threshold, normalize. Argument is `percentiles=` (0..1), **not** `probs=`. |
| `gdist`    | N-d count array | Binned joint distribution (signal × distance, etc.). |
| `gsummary` | per-scope summary | One-shot min/max/mean/sum/Nbin (often before deciding `gdist` breaks). |
| `gintervals.summary` | per-*interval* min/max/mean/sum/sd/nbin | One row per input interval; one-call alternative to multi-expression `gextract`. |
| `gcor`     | scalar correlation (or matrix) | Pearson / Spearman correlation between two+ expressions. |
| `gintervals.neighbors` | intervals + nearest-feature cols | Nearest-feature annotation; supports `maxdist`/`mindist`. |
| `gintervals.annotate` | input rows + selected annot cols + `dist` | Column-attach without changing row count or order. |
| `gintervals.canonic` | merged intervals frame | After rbind / union when you want overlaps merged. |
| `gintervals.to_mat` / `gintervals.from_mat` | `intervs_mat` (matrix subclass) / data.frame | Round-trip intervals + values to a numeric matrix for clustering / heatmaps. |
| `gtrack.create` / `gtrack.smooth` | (side effect) | Materialize a derived track. |

### Distance - which primitive?

| Use | Returns | When |
|---|---|---|
| `gvtrack.create(_, src, "distance")` / `"distance.center"` / `"distance.edge"` | per-bin numeric column, signed only when the source has a strand | You're sweeping the genome at fixed resolution and want distance as a feature alongside other tracks. |
| `gintervals.neighbors(a, b)` | paired-row frame (a + b + `dist`, signed only by a strand) | k-NN, distance bands via `mindist`/`maxdist`, asymmetric upstream/downstream. |
| `gintervals.annotate(intervs, annot, annotation_columns)` | input rows + selected annot cols + `dist` | Attach gene symbol + distance to a fixed query frame without changing row count. |

### Aggregators (`func =` on `gvtrack.create`)

| `func` | Use when |
|---|---|
| `sum`  | Count signals (ATAC, ChIP, contacts). Adds across bins in the window. |
| `area` | Width-aware sum for sparse / 2D sources where bin widths vary. |
| `weighted.sum` | Sparse / 2D contact aggregation against rectangles. |
| `avg`  | Continuous signal where you want a per-bin mean. **Do not use for count tracks.** |
| `max`  | Peak detection (PWM energies, percentile maxima). |
| `min`  | Local-min insulation border calling, diagnostic windows. |
| `lse`  | log-sum-exp - for summing motif energies (log-space). Composes correctly under nested aggregation. |
| `global.percentile.max` | Genome-wide percentile of max-in-window. Comparable across tracks. |

### Argument-name gotchas (silent footguns)

| Wrong | Right | Why it bites |
|---|---|---|
| Bare unquoted track names: `gextract(epi.k27me3.es, ...)` | `gextract("epi.k27me3.es", ...)` | The character-vector form is unambiguous and works with sprintf'd names. |

## 4. Recipes

### 4.1 Build and manipulate intervals

**Build.**

```r
# 1D - vectorized; chrom can be a scalar or per-row.
peaks <- gintervals(chrom = c("chr1", "chr2"),
                    start = c(100,    200),
                    end   = c(500,    700))

# 2D - three canonical forms.
gintervals.2d(chrom1, start1, end1, chrom2, start2, end2)   # explicit rectangles
gintervals.2d(chroms1 = chrs, chroms2 = chrs)               # whole-chrom cartesian
gintervals.2d(chrom, 0, -1, chrom, 0, -1)                   # one whole chrom (cis)
```

**Load from file.** Dedicated readers for the three common annotation formats - output is a validated misha intervals frame:

```r
peaks <- gintervals.import_bed("peaks.bed", name = TRUE, score = TRUE, strand = TRUE)
genes <- gintervals.import_gff("refseq.gff", feature = "exon", attrs = TRUE)
snps  <- gintervals.import_vcf("snps.vcf",  info = TRUE)
```

`gintervals.import_bed` sorts by default; `keep_order = TRUE` keeps the file's row order (header lines skipped), which peak IDs that are row numbers need. `format = "narrowPeak"` / `"broadPeak"` (guessed from a `.narrowPeak` / `.broadPeak` name, `.gz` too) adds `signalValue`, `pValue`, `qValue` and, for narrowPeak, `peak`, the summit's offset from `start` (`-1` = no summit). Summit-centred windows, with the summit at `(start + end) %/% 2L`:

```r
peaks <- gintervals.import_bed("ENCFF000ABC.narrowPeak.gz", keep_order = TRUE)
peaks <- peaks[peaks$peak >= 0, ]
win <- transform(peaks, start = start + peak - 250L, end = start + peak + 250L, peak = 250L)
```

`gintervals.import_bed` also reads a bigBed (detected by its magic number, not its name): columns are chrom, start, end, then the file's fields in order, named and typed by its autoSql (numeric autoSql types -> numeric; strings, arrays, enums and itemRgb stay text). `format` must be `NULL` for a bigBed.

**bigWig signal over intervals, without importing a track.** `gintervals.import_bigwig(file, intervals)` returns the records (chrom, start, end, value) overlapping the intervals, whole and each once, reading only the data blocks under them. Feed it to a value-based vtrack; `weighted.sum` is the per-base sum (what `bbi.stackup(..., summary = "sum")` gives), `sum` adds each touched record's value once:

```r
gvtrack.create("sig", src = gintervals.import_bigwig("signal.bigWig", peaks), func = "weighted.sum")
gextract("sig", peaks, iterator = peaks)   # NaN where no record overlaps
```

**Sort and merge overlaps (`gintervals.canonic`).** When you want a non-overlapping, sorted set - typical after rbind'ing per-condition `gscreen` results. To fold per-row metadata into the merged set, use `gintervals.mark_overlaps` (tags each source row with its merge-group ID) and aggregate per group:

```r
merged <- rbind(set_a, set_b) |>
    gintervals.mark_overlaps() |>
    group_by(overlap_group) |>
    summarise(chrom  = first(chrom),
              start  = min(start),
              end    = max(end),
              strand = if (n_distinct(strand) == 1L) first(strand) else 0L,
              .groups = "drop") |>
    select(-overlap_group)
```

`gintervals.mark_overlaps` is `gintervals.canonic` plus a join-friendly group column; reach for it when you need custom per-merged-region aggregation. Bare `gintervals.canonic(...)` is fine when you only want the merged intervals frame and don't need to carry source metadata.

Both also merge intervals that only touch (`end` of one == `start` of the next). To merge overlapping sites but keep adjacent ones apart, pass `unify_touching_intervals = FALSE` to either. `attr(gintervals.canonic(x), "mapping")` gives, for each row of `x`, the row of the merged set it went into; that vector is the `overlap_group` column.

If you want to *keep* duplicate / overlapping rows (paired matches, parallel rows per sample), skip canonic.

**Set operations.**

```r
gintervals.union(a, b)         # union
gintervals.intersect(a, b)     # intersection
gintervals.diff(a, b)          # a minus b
```

**Fixed-width peaks (center ± W).**

```r
peaks_2k <- gintervals.normalize(peaks, 2000)   # 2kb-wide, centered on the original midpoint
```

For asymmetric expansion, do it manually with integer arithmetic - *never* float division:

```r
peaks$mid <- (peaks$start + peaks$end) %/% 2L
peaks <- gintervals(chrom = peaks$chrom,
                    start = peaks$mid - W_up,
                    end   = peaks$mid + W_down + 1L)
peaks <- gintervals.canonic(peaks)
```

**Clip to genome bounds.** `gintervals.force_range(intervs)` clamps `start >= 0` and `end <= chrom_len`. Use after any expansion that may push past chromosome ends. It shortens windows silently, so for fixed-width windows (model inputs, `gextract.matrix`) drop the rows whose width changed: `win <- gintervals.force_range(win); win <- win[win$end - win$start == width, ]`.

**Symmetric expansion shortcut.** `gintervals.expand(intervs, 100)` adds N bp to each side and clamps with `gintervals.force_range` in one call. Use it in place of the explicit arithmetic when expansion is symmetric.

**Persist.**

```r
gintervals.save("intervs.my.peaks", peaks)
my_peaks <- gintervals.load("intervs.my.peaks")
```

**Avoid:**
- `(start + end) / 2` - float division yields non-integer centers. Always `(start + end) %/% 2L`.
- `gintervals.canonic` reflexively after every `rbind` - only run it when you actually want overlaps merged.
- Base `data.frame(chrom = ..., start = ..., end = ...)` - checked against the genome only when a query reads it. Pass through `gintervals(...)` for validation, unless the row order matters (A9).
- Hand-rolling fixed-width peaks via mutate when `gintervals.normalize(intervs, size)` does the same in one call.

### 4.2 Annotate intervals with nearby features

Two primitives, different return shapes:

- **`gintervals.neighbors(query, target, ...)`** - paired-row output, one input row × matched neighbor; `dist` column (sign rule below). Good for k-NN, distance-banded matching, asymmetric upstream/downstream.
- **`gintervals.annotate(intervals, annotation_intervals, annotation_columns=, ...)`** - column-attach output, preserves row count and order. Good for "add gene symbol + distance".

```r
tss <- gintervals.load("intervs.global.tss")

# Nearest neighbor, dropping unmatched rows:
near <- gintervals.neighbors(peaks, tss)

# Distance-banded - keep all peaks (na.if.notfound is the key flag):
near <- gintervals.neighbors(peaks, tss, maxdist = 50e3, na.if.notfound = TRUE)

# k nearest neighbors:
knn <- gintervals.neighbors(peaks, tss, maxneighbors = 5)

# Every overlapping pair (e.g. motif hits inside peaks). maxdist = 0 also returns
# intervals that only touch (end == start), so keep the true overlaps:
ov <- gintervals.neighbors(peaks, hits, maxneighbors = 1e6, mindist = 0, maxdist = 0)
ov <- ov[ov$start1 < ov$end & ov$end1 > ov$start, ]

# Strand-aware variants:
gintervals.neighbors.upstream(peaks, tss, maxdist = 100e3)
gintervals.neighbors.downstream(peaks, tss, maxdist = 100e3)
gintervals.neighbors.directional(peaks, tss,
    maxneighbors_upstream = 1, maxneighbors_downstream = 1)

# Column attach - keeps row count and order:
peaks <- gintervals.annotate(peaks, tss,
                             annotation_columns = "geneSymbol",
                             dist_column = "tss_dist",
                             max_dist    = 100e3,
                             na_value    = NA)
```

Sign of `dist` (the same for `gintervals.annotate`'s distance column). It is never signed by genome coordinates alone:

- **Default:** it follows the target's (`intervals2`) strand: positive when the query lies downstream of the target in the target's direction, negative upstream. A peak 1 kb downstream of a + strand TSS has `dist = 1000`.
- **Target without a `strand` column:** `dist` is unsigned (never negative), so a negative `mindist` or `dist < 0` matches nothing.
- **`use_intervals1_strand = TRUE`** (and `gintervals.neighbors.upstream/downstream/directional`): it follows the query's strand: positive when the target lies downstream of the query, negative upstream, on either strand. A query without a `strand` column counts as +, which gives the sign in genome coordinates (positive = target to the right).

**Promoters from TSS in one call.** Strand-aware promoter windows come up constantly:

```r
promoters <- gintervals.promoters("intervs.global.tss", upstream = 500, downstream = 50)
```

returns one window per TSS, `upstream` bp before the TSS and `downstream` bp from it on, oriented by the `strand` column (a set built by `gdb.install_intervals` names it `"tss"`). Strand 0, a TSS shared by genes on both strands, gets a window covering both directions. Use it instead of rolling the strand math by hand.

**Avoid:**
- `gintervals.neighbors` without `na.if.notfound = TRUE` when annotating a fixed set - rows with no neighbor in range are dropped silently.
- `gintervals.neighbors` when you only want to *attach* columns to an existing frame - `gintervals.annotate` is the right tool.
- Assuming `dist` is signed by genome position. It is signed only by a strand (the target's, by default), so check which one before filtering on its sign. Filter on size with `abs(dist) < W`; with a strand, bare `dist < W` accepts arbitrary upstream distances.

### 4.3 Distance flavors

`gvtrack.create` accepts three `func` values for distance, with different semantics:

| `func` | Measures | Returns 0 when | Notes |
|---|---|---|---|
| `distance` | Iterator-bin *center* → nearest source *edge* (outside); *normalized* fractional position when inside. | Bin center exactly on a source edge | Mixed semantics - fine for genome-wide profiles, surprising if you assumed pure edge-to-edge. |
| `distance.center` | Iterator-bin *center* → *center* of the source interval containing it; NaN when the bin center is outside every source interval. | Bin center coincides with a source center | Use for offsets inside source intervals (e.g. position within a peak); for distance to anchors outside them use `distance` or `distance.edge`. |
| `distance.edge` | *Edge-to-edge*, same as `gintervals.neighbors`. | Iterator interval overlaps the source | Use when you specifically want `gintervals.neighbors` semantics inside a vtrack. |

All three return *signed* distance when the source has a `strand` column (sign = direction relative to strand); unsigned otherwise. All three return NA when the chromosome has no source intervals, and `distance.center` also when the bin center is outside every source interval.

```r
gvtrack.create("d_tss", src = "intervs.global.tss", func = "distance")

# Per-bin distance - one column alongside other tracks:
gextract(c("epi.atac.es", "d_tss"), intervals = peaks, iterator = 20)

# Distance-banded screen:
gscreen("epi.atac.es > 5 & abs(d_tss) > 5000",
        intervals = gintervals.all(), iterator = 20)
```

For 2D, bind a 1D `distance` vtrack to one axis of a 2D iteration with `dim = 1` or `dim = 2`.

**Avoid:**
- Treating distance vtrack values as unsigned - `abs(d_tss) < W` is almost always what you want.
- Forgetting NA on chromosomes with no source intervals. Wrap with `ifelse(is.na(d_x), Inf, d_x)` if you want "no neighbor = infinitely far".

### 4.4 Extract per-region values

Signature: `gextract(..., intervals, colnames, iterator, band, file, intervals.set.out, intervals_join)`. The expressions (`...`) are a character vector - bare track names, vtrack names, or arithmetic involving them.

**Many tracks, one call - the central concept.** A single `gextract` over a `c("expr1", "expr2", ...)` vector makes ONE pass over the genome:

```r
profs <- gextract(c("epi.k27me3.es", "epi.k4me3.es", "atac.es"),
                  intervals = peaks,
                  iterator  = peaks,
                  colnames  = c("k27", "k4", "atac"))
```

**`colnames =` is how you name the output columns.** Without it, each column is named by the literal expression string - fine for bare track names, ugly for anything with arithmetic. Always pass `colnames` when expressions are not bare track names.

**One row per region.** Pass the same intervals frame as both `intervals=` and `iterator=`:

```r
gvtrack.create("k27", "epi.k27me3.es", "sum", sshift = -2000, eshift = 2000)
gvtrack.create("atac", "atac.es",       "sum", sshift = -250,  eshift = 250)

per_peak <- gextract(c("k27", "atac"),
                     intervals = peaks, iterator = peaks,
                     colnames  = c("k27_flank", "atac_summit")) |>
            arrange(intervalID) |>
            select(-intervalID)
```

`arrange(intervalID) |> select(-intervalID)` recovers input row order and drops the helper column.

**Attach the input intervals' columns directly (no R join).** When you want the input `peaks` columns - coords, gene names, scores - on every output row, pass `intervals_join = "intervals"`:

```r
per_peak <- gextract(c("k27", "atac"),
                     intervals = peaks, iterator = peaks,
                     colnames  = c("k27_flank", "atac_summit"),
                     intervals_join = "intervals")
# per_peak has chrom/start/end + k27_flank + atac_summit + chrom1/start1/end1
# + every other column of `peaks` (e.g. gene_id, score, ...)
```

`intervalID` is dropped; conflicting names from `peaks` get a `"1"` suffix (`chrom` -> `chrom1`). Output rows are in genomic-sort order, not original `peaks` order. The C++ side copies the input columns by `intervalID` index, replacing the older `misha.ext::gextract.left_join` R-side wrapper at ~2x the speed on wide-window workloads (e.g. TSS +/- 5kb at 50bp iterator) and dropping the misha.ext dependency. Use `intervals_join = "none"` to just drop `intervalID` without attaching anything. Only available when returning to memory; combining with `file=` or `intervals.set.out=` raises an error.

**NA handling - crucial gotcha.** Sparse tracks return NA at unmeasured bins; many dense tracks also have NAs. Fill at extract time, inside the expression:

```r
gextract(c("ifelse(is.na(epi.atac.es), 0, epi.atac.es)",
           "ifelse(is.na(cpgs.meth),  0, cpgs.meth)"),
         intervals = peaks, iterator = peaks,
         colnames  = c("atac", "meth"))
```

Treat this as a default, not an edge case.

**Expression functions must be vectorized.** The C++ engine passes track *vectors* to your expression and expects a same-length vector back. Standard vectorized R ops all work (basic arithmetic, `ifelse`, `pmin`/`pmax`, log/exp). A function that returns one value, such as `mean(track)` or `min(a, b)`, is not an error: misha sees that the result is shorter than its input and evaluates the expression one bin at a time instead. `min(a, b)` then gives the per-bin minimum, but the call is much slower (about 40x slower than `pmin(a, b)` over 1M bins), and `mean(track)` is just the bin's own value, not a window mean - use an `avg` vtrack for that. Base R has no `pmean`; for per-element averaging across two tracks, use `tgutil::pmean(a, b, na.rm = TRUE)` (`(a + b) / 2` is wrong: it propagates NA whenever *either* input is NA, so any position where one track has data and the other doesn't comes out as NA instead of the value you have).

**Tiled bins inside intervals.** For sub-region resolution (per-bin profile), pass an integer iterator:

```r
prof <- gextract("epi.k27me3.es",
                 intervals = regions,
                 iterator  = 20)        # one row per 20bp bin within regions
```

**Iterator coordinate gotcha.** `iterator = 20` tiles bins relative to **chromosome start**, not to each interval's start - so two peaks at different chromosome offsets get *misaligned* bin grids. For interval-relative bin indices use `giterator.intervals(intervals = peaks, iterator = 20, interval_relative = TRUE)`.

**Equal-width windows as an N x L matrix (`gextract.matrix`).** For per-base (or per-bin) profiles over many windows of the same width, e.g. model inputs or motif-hit masks, skip the long table:

```r
# windows around the summits, in the peaks' row order (gintervals() would sort them, A14b)
peaks <- gintervals.import_bed("peaks.narrowPeak.gz", keep_order = TRUE)
peaks <- peaks[peaks$peak >= 0, ]
win <- transform(peaks, start = start + peak - 1057L, end = start + peak + 1057L)    # 2114bp
sig  <- gextract.matrix("atac.es", intervals = win)                      # nrow(win) x 2114
gvtrack.create("hit_cov", motif_hits, func = "coverage")
mask <- gextract.matrix("hit_cov", intervals = win) > 0                  # logical mask
win32 <- transform(peaks, start = start + peak - 1056L, end = start + peak + 1056L)  # 2112 = 66 x 32
arr  <- gextract.matrix(c("atac.es", "k27"), intervals = win32, iterator = 32)  # N x 66 x 2
```

Row `i` is input row `i` (unsorted, duplicated and overlapping rows are fine); columns are in genomic order whatever the strand. Several expressions give an N x L x E array named by `colnames`. All windows must have the same width and `iterator` must divide it. Unlike gextract's numeric iterator, the bins start at each window's start (as `giterator.intervals(..., interval_relative = TRUE)`), so any start works. Values are what gextract gives for each bin as an iterator interval, NaN included; a dense-track bin that overlaps several track bins gets their plain mean, not weighted by overlap. A window past a chromosome end is an error, as in gextract.

**One row per stored position (`iterator = <track>`).** For sparse-aligned data - per-CpG methylation, per-fragment coverage:

```r
m <- gextract("cpgs.meth", intervals = gintervals.all(), iterator = "cpgs.cov")
```

**Per-interval descriptive stats in one call.** When you want `nbins / n_nan / min / max / sum / mean / sd` per input region, `gintervals.summary(expr, intervals, iterator)` returns one row per interval with all of them appended:

```r
stats <- gintervals.summary("epi.k27me3.es", intervals = peaks, iterator = 20)
# stats has chrom/start/end + nbins, nbins.nan, min, max, sum, mean, stdev
```

This subsumes a multi-expression `gextract` for the standard descriptive-stats case. Reach for the multi-expression form only when you need non-stat aggregates (custom expressions, NA-fills, ratios) at the same time.

**Persist as a named intervals set.**

```r
gextract("score > 5", intervals = gintervals.all(), iterator = 20,
         intervals.set.out = "intervs.high_score")
```

**Reshape to a matrix for downstream analysis.** After `gextract`, the typical next step for clustering / heatmaps is a numeric matrix indexed by interval. `gintervals.to_mat` and its inverse `gintervals.from_mat` round-trip cleanly:

```r
mat <- gintervals.to_mat(profs)                # N peaks x K tracks; rownames = "chr:start-end"
pcs <- prcomp(mat, center = TRUE)$x[, 1:5]
out <- gintervals.from_mat(mat[order(pcs[, 1]), ])  # back to intervals + values
```

The intervals are carried in `attr(mat, "intervals")` - rownames are display-only and never parsed back, so chrom names containing underscores (`chrUn_KI270442v1`, scaffolds) round-trip safely and `start`/`end` keep integer type. `[`, `head`, `tail`, and `rbind` keep the attribute aligned. Pass `labels = FALSE` to skip rowname construction (~25x faster `to_mat` on large inputs when you don't need to look at them). Pass `id_col = "gene"` (or any column name) to use that column for rownames instead of coordinates. Pass `value_cols = c("t1", "t2")` to subset which numeric columns go into the matrix.

**Avoid:**
- Looping `gextract` over a vector of tracks. Always one call with `c("track1", "track2", ...)`.
- Bare unquoted expression arguments. Always use the character-vector form.
- Post-extract `prof[is.na(prof)] <- 0` when you can write `ifelse(is.na(track), 0, track)` inside the expression.
- Scalar-result functions (`mean(track)`, `min(a, b)`, `sum(track)`) inside an expression - misha falls back to evaluating them one bin at a time, which is slow, and `mean(track)` is just the bin's value. Use `pmin` / `pmax` / `ifelse` or a vtrack.
- `iterator = peaks` when you want sub-region resolution. Use `iterator = 20`.
- Forgetting that `iterator = 20` is chromosome-relative, not interval-relative.
- Skipping `arrange(intervalID)` when downstream code assumes input row order.
- `misha.ext::intervs_to_mat` / `mat_to_intervs` - deprecated. Use `gintervals.to_mat` / `gintervals.from_mat` (in misha >= 5.7.0). The legacy versions corrupt chrom names containing underscores and lose integer typing on round-trip.
- `misha.ext::gextract.left_join(tracks, intervals, ...)` - deprecated. Use `gextract(tracks, intervals = ..., intervals_join = "intervals", ...)` (in misha >= 5.7.1). The built-in attaches the input columns at the C++ layer (~2x faster on wide-window workloads) and drops the misha.ext dependency.
- Any other `misha.ext` helper - the package is retired. `get_promoters` -> `gintervals.promoters`, `gintervals.expand` / `gintervals.remove_overlaps` -> the misha functions of the same name, `gintervals.centers(x)` -> `gintervals.normalize(x, 1)`, `gintervals.neighbors1(...)` -> `gintervals.neighbors(..., na.if.notfound = TRUE)`, `grandom_genome` -> `gintervals.random`. misha.ext's own `gintervals.normalize` gives `size - 1` bp for odd sizes; misha's gives `size`.

### 4.5 Aggregate signal in a flanking window

One `gvtrack.create` per aggregator + window combination, then `gextract`:

```r
# Window-summed ChIP signal in ±2kb:
gvtrack.create("k27_flank", src = "epi.k27me3.es", func = "sum",
               sshift = -2000, eshift = 2000)

# Sum ATAC counts in the summit (-250 / +250bp):
gvtrack.create("atac_summit", src = "atac.es", func = "sum",
               sshift = -250, eshift = 250)

# Max PWM energy across ±100bp:
gvtrack.create("ctcf_max", src = "motifs.ctcf", func = "max",
               sshift = -100, eshift = 100)

per_peak <- gextract(c("k27_flank", "atac_summit", "ctcf_max"),
                     intervals = peaks, iterator = peaks)
```

**Window placement.** Symmetric `sshift = -W, eshift = W` is the default. For asymmetric (upstream-only DI, strand-bias diagnostics) use `sshift = -W, eshift = 0` or `sshift = 0, eshift = W`. Point-sample at the iterator bin: omit shift args. 2D rectangle: `sshift1=, eshift1=, sshift2=, eshift2=`.

**Avoid:**
- `gvtrack.create(...); gvtrack.iterator(name, sshift=, eshift=)` as two calls - current misha takes `sshift / eshift` directly in `gvtrack.create`. The two-call form is legacy.
- Defensive `if (length(gvtrack.ls(n)) == 1) gvtrack.rm(n)` before every `gvtrack.create` - current misha silently overwrites on re-create.
- Picking `avg` when you mean `sum` - for count tracks, `avg` divides by bin count and discards magnitude.
- Forgetting that `sshift` / `eshift` are in *base pairs*, not bin counts.

### 4.6 Call regions of interest

`gscreen(expr, intervals, iterator, intervals.set.out)`. Bins where the expression evaluates to non-zero / non-NA become intervals (adjacent bins auto-merged).

**Pick threshold from a quantile, then screen.** Almost never hardcode a threshold - derive it from the empirical distribution:

```r
thr <- gquantiles("epi.k27me3.es",
                  percentiles = 0.99,
                  intervals   = gintervals.all())

peaks <- gscreen(sprintf("epi.k27me3.es > %g", thr),
                 intervals = gintervals.all(),
                 iterator  = 20)
```

`gquantiles` argument is **`percentiles =`** (0..1), not base R's `probs =`. For *per-region* quantiles (one row per interval), reach for `gintervals.quantiles` instead - same arg name, dedicated to the per-interval case.

**Multi-track OR.** Build a `" | "`-joined expression so one `gscreen` returns the union of per-track hits:

```r
expr <- paste(sprintf("(%s > %g)", track_names, thrs), collapse = " | ")
peaks <- gscreen(expr, intervals = gintervals.all(), iterator = 20)
peaks <- gintervals.canonic(peaks)
```

**Refine peak centers by argmax.**

```r
summits <- gextract("epi.k27me3.es", intervals = peaks, iterator = 20) |>
    group_by(intervalID) |>
    slice_max(epi.k27me3.es, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(chrom, start, end) |>
    gintervals.normalize(280)
```

**Persist as a named set.**

```r
gscreen("epi.k27me3.es > 5",
        intervals          = gintervals.all(),
        iterator           = 20,
        intervals.set.out  = "intervs.k27.peaks")
```

**Distance-band screen.**

```r
gvtrack.create("d_tss",  "intervs.global.tss",         "distance")
gvtrack.create("d_rmsk", "intervs.global.rmsk_LINE",   "distance")
gscreen("abs(d_tss) > 5000 & abs(d_rmsk) > 100",
        intervals = peaks, iterator = 20,
        intervals.set.out = "intervs.peaks.intergenic")
```

**Avoid:**
- Hardcoded literal thresholds (`gscreen("track > 1870", ...)`) - uninterpretable across datasets; derive from `gquantiles`.
- `gquantiles(x, probs = 0.9)` - wrong argument name. Always `percentiles =`.
- `gquantiles` on a track with many NAs assumed to include zeros. NA positions are silently dropped. Wrap with `ifelse(is.na(track), 0, track)` if zeros should count.

### 4.7 Joint distributions, summaries, correlations

Four primitives:

- **`gdist(track1, breaks1, track2, breaks2, ..., intervals, iterator)`** - N-d count histogram. Pass `dataframe = TRUE` for a long-format frame with one column per axis + count column `n`.
- **`gbins.summary(strat_track, breaks, expr = value_track, ...)`** - per-bin n/mean/sum/var/min/max of `value_track` stratified by `strat_track` bins.
- **`gintervals.quantiles(expr, percentiles, intervals, iterator)`** - per-region quantile vector.
- **`gcor(expr1, expr2, ..., method)`** - Pearson / Spearman correlation between two+ expressions.

```r
# Joint signal × distance histogram:
gvtrack.create("d_tss", "intervs.global.tss", "distance")
h <- gdist("epi.k27me3.es", c(-Inf, seq(0, 10, 0.5), Inf),
           "d_tss",          c(-Inf, -5e4, -1e4, -1e3, 0, 1e3, 1e4, 5e4, Inf),
           intervals = gintervals.all(), iterator = 20,
           dataframe = TRUE, names = c("signal", "dist_tss"))

# Stratified mean signal per distance bin:
summ <- gbins.summary("d_tss", seq(-1e6, 1e6, 1e4),
                      expr = "epi.k27me3.es",
                      intervals = gintervals.all(), iterator = 20)

# Per-region quantiles (one row per peak):
q <- gintervals.quantiles("epi.atac.es",
                          percentiles = c(0.5, 0.9, 0.99),
                          intervals   = peaks, iterator = 20)

# Correlation between two tracks (Spearman, genome-wide):
gcor("epi.k27me3.es", "epi.atac.es",
     intervals = gintervals.all(), iterator = 20,
     method    = "spearman")
```

**Picking `gdist` breaks.** Continuous tracks: `seq(0, 1, length.out = 21)`. Count tracks: zero-vs-positive split + log spacing. Always `gsummary` first to pick range; never hardcode upper bound.

**Avoid:**
- Reaching for `gdist` when one axis is continuous and you want a per-bin mean - `gbins.summary` is the right tool.
- Hardcoded `breaks` without `gsummary` first - produces overflow / underflow bins silently.
- Mistaking `gdist`'s output for a probability - it's *counts*, not a density.

### 4.8 Materialize derived tracks

Four creation primitives:

- **`gtrack.create(track, description, expr, iterator)`** - evaluate an expression genome-wide and write the result as a dense track.
- **`gtrack.create_sparse(track, description, intervals, values)`** - values at irregular positions; NA elsewhere.
- **`gtrack.create_dense(track, description, intervals, values, binsize, defval, func)`** - fully dense fixed-bin track from interval/value pairs. `func` ∈ {`"weighted.mean"` (default), `"weighted.sum"`, `"coverage"`, `"max"`, `"min"`, `"median"`, `"count"`}.
- **`gtrack.smooth(track, description, expr, winsize, alg)`** - windowed-mean track. `alg = "LINEAR_RAMP"` (triangular, default) or `"MEAN"` (boxcar).

```r
gtrack.create("epi.k27me3.log",
              description = "log2(K27me3 + 1)",
              expr        = "log2(epi.k27me3.es + 1)",
              iterator    = 20)
```

**Indicator from an intervals frame - no materialization.** Add a `value` column and pass the frame as `src` of a vtrack:

```r
peaks$value <- 1
gvtrack.create("peak_ind", src = peaks, func = "max")
gextract("ifelse(is.na(peak_ind), 0, peak_ind)",
         intervals = gintervals.all(), iterator = 20)
```

Materialize an on-disk indicator track only when the indicator is reused across sessions or a downstream consumer specifically needs a stored track.

**Smooth an existing track.**

```r
gtrack.smooth("epi.k27me3.smooth",
              description = "K27me3 ±20kb LINEAR_RAMP",
              expr        = "epi.k27me3.es",
              winsize     = 20000)
```

**Namespace via directories.** Long names with dots map to directories on disk:

```r
gtrack.create_dirs("epi.derived")
gtrack.create("epi.derived.k27_over_atac",
              description = "K27me3 over ATAC ratio",
              expr        = "epi.k27me3.es / (atac.es + 1)",
              iterator    = 20)
```

**Idempotency guard.**

```r
if (gtrack.exists("derived")) gtrack.rm("derived", force = TRUE)
gtrack.create("derived", "...", expr = "...", iterator = 20)
```

**Persist track metadata.** `gtrack.attr.set(track, key, value)` / `gtrack.attr.get(track, key)` stores per-track scalar metadata (training parameters, source paths, dates). The intervals-set parallel is `gintervals.attr.set(set, key, value)` / `gintervals.attr.get(set, key)` for annotation metadata on named intervals sets (source build, filter version, etc.) - same convention.

**Avoid:**
- Materializing a "track" used once or twice - for ephemeral signals, a vtrack (in-memory) is the right tool.
- `gtrack.create(... expr = "a / b")` without NA-safe divisor - wrap as `"a / (b + 1e-6)"` at creation time.
- `gtrack.rm("name")` without `force = TRUE` in a script - the interactive confirmation prompt hangs unattended runs.
- Aggregating a stored `pwm` (LSE) track with `func = "sum"` - LSE doesn't compose under summation. Use `func = "lse"`.
