---
name: importing-tracks
description: Use when importing external data (WIG / bigWig / BED / TSV / BAM / 2D contacts / array matrices) into a misha track DB, when picking which `gtrack.*import*` / `gtrack.create_*` / `gtrack.liftover` / `gtrack.copy` variant to call, or when a track import fails opaquely or produces an all-zero / all-NaN track.
---

# Importing tracks into misha

## Overview

A "track import" is any operation that materializes a persistent on-disk track from an external source — a single file, a directory or URL of files, a BAM, an in-R data frame, a track in another misha DB, or a track in another genome assembly. There are eight+ entry points; picking the wrong one is the most common ingestion problem after silent chrom-name mismatches.

Two failure modes worth pre-empting:

1. **Picking the wrong variant** for the source shape.
2. **Silent partial imports.** WIG/bigWig importers run with the C++-internal flag `ignore_unknown_chroms` hard-wired to true (`src/GenomeTrackImportWig.cpp` instantiates `Wig` with that argument; there is no R parameter to flip it) - chroms not in `gintervals.all()$chrom` are dropped without error. Mismatched naming (`1` vs `chr1`) yields an "all NaN" track, not a failure.

## Precondition: chrom-name compatibility (mandatory)

Before any file-based import, verify the source chroms resolve in the active gdb:

```r
gdb.init("/path/to/trackdb")
setdiff(unique(input$chrom), gintervals.all()$chrom)   # MUST be empty for chroms you care about
```

For huge WIG/bigWig inputs you can't load into R, extract chroms from the source:

```bash
grep -oE 'chrom=[^ ]+' file.wig | sort -u | sed 's/chrom=//'   # WIG
bigWigInfo -chroms file.bw                                     # bigWig
```

Non-empty diff → wrong gdb, naming-convention mismatch, or a gdb missing scaffolds. Fix before importing. Silent partial imports are the #1 cause of "the track is all zeros".

## Format chooser

Convention used throughout: `binsize > 0` produces a dense fixed-bin track (one value per `binsize`-wide genomic window); `binsize = 0` produces a sparse track (one value per source interval, no inter-record bookkeeping). Picking the wrong shape for the source is the second-most-common ingestion mistake after chrom-name mismatch.

| Source | Function | Key args / caveats |
|---|---|---|
| In-memory R intervals + values | `gtrack.create_dense(name, desc, intervals, values, binsize, func = "weighted.mean"/"coverage"/"max"/"min"/"sum"/"median"/"count")` (5.6.31+) or `gtrack.create_sparse(name, desc, intervs, values)` | preferred whenever data is already in R — no TSV round-trip and no `scipen` traps |
| Single WIG / bigWig file | `gtrack.import(name, desc, file, binsize, defval = NaN)` | auto-dispatch on `.wig` / `.bw` / `.bigWig` / magic bytes; pass `defval = NaN` (not `0`) for continuous signal so uncovered bins read as NA |
| Single BED file (or 4-col TSV) | same `gtrack.import` — pass `binsize` for a dense track, omit it for a sparse track | misha ≥ 5.7.1 required if you pass `binsize` (5.7.0 errors with `func argument is not a string`) |
| R-computed values via TSV (chrom/start/end/value) | `gtrack.import(name, desc, file, binsize)` | run `options(scipen = 20)` BEFORE `write.table` — scientific-notation coords parse as malformed |
| Directory or URL-with-wildcards of WIG/bigWig files | `gtrack.import_set(description, path, binsize, track.prefix = "", defval = NaN)` | bulk import; `path` may be a glob (`/data/*.bw`) or an `ftp://` URL; continues on per-file errors and returns successes/failures |
| ChIP/ATAC pileup from read intervals (data.frame in R) | `gtrack.create_dense(name, desc, intervals = reads, values = rep(1, nrow(reads)), binsize = 20, defval = 0, func = "coverage")` | one-call replacement for the old `gtrack.import_mappedseq` route when reads are already in R |
| Mapped reads (SAM, tab-delimited, or gzipped `.sam.gz` / `.tsv.gz`) | `gtrack.import_mappedseq(name, desc, file, pileup, binsize, cols.order = c(9, 11, 13, 14), remove.dups = TRUE)` | for SAM pass `cols.order = NULL`; gzip is auto-detected by magic bytes (misha ≥ 5.8.0); the default `remove.dups = TRUE` is correct for ChIP/ATAC/CUT&RUN |
| Single BAM file, no pre-filter (misha ≥ 5.8.0) | `gtrack.import_mappedseq(name, desc, "<file>.bam", pileup, binsize, remove.dups = TRUE)` | bgzip magic auto-detected; misha spawns `samtools view` internally so `samtools` must be on `PATH`. `cols.order` is forced to SAM mode; passing a non-NULL `cols.order` is an error |
| Multi-BAM concat or MAPQ-filtered BAM | `misha.ext::gtrack.import_mappedseq_bam(bam_files, track, min_mapq = NULL, ...)` | misha.ext wrapper: stages `samtools view` (with `samtools cat` for multi-BAM) through a named fifo into `gtrack.import_mappedseq`. Use this when you need MAPQ filtering or multiple BAMs concatenated; for a single unfiltered BAM the misha 5.8.0+ native path above is simpler |
| 2D Hi-C / capture-C contacts (per-contact adj + fends) | `gtrack.2d.import_contacts(name, desc, contacts, fends = "<redb>/<RE>.fends", allow.duplicates = TRUE)` | `contacts =` accepts a vector (shaman per-rect scores, scHi-C per-core stagings - see "scHi-C pooling" below). `allow.duplicates = TRUE` (the default) sums repeated fend pairs; use `allow.duplicates = FALSE` only for pre-scored shaman per-rectangle outputs where any duplicate is an upstream bug |
| Restriction-enzyme fragment + fend prerequisites (any new 4C/Hi-C bootstrap) | `gtrack.import('redb.<RE>_flen', flen_file, 0)` + `gtrack.import('redb.<RE>_gc', gc_file, 0)` + `gtrack.create('redb.<RE>_map', mapab_track, iterator = 'redb.<RE>_flen')` | per-fragment iterator promotes a per-base mapability track to per-fragment mean. See "RE fragment-track bootstrap" below |
| Pooled multi-cell scHi-C track | `gtrack.2d.import_contacts(pool_name, desc, contacts = c(stage_1, ..., stage_N), fends = "<redb>/<RE>.fends")` | parallel-stage per-cell `gextract`s to N TSVs, then a single pool call. See "scHi-C pooling" below |
| Per-coordinate 2D matrix (already binned, no fends) | `gtrack.2d.import(name, desc, file)` | lower-level than `_contacts`; takes a coord-resolution TSV directly. Used for 4DN-style 1kb contacts. Skip when you have adj+fends |
| 2D track from a per-rectangle R data.frame | `gtrack.2d.create(name, desc, intervals, values)` | rectangles must be non-overlapping |
| Deep-learning model predictions (Borzoi / IceQream / equivalent) | `gtrack.import(track_name, desc, per_bin_tsv, binsize = res)` where the TSV is the model's per-bin output | lab convention: predicted tracks live under `seq.IQ.<model>.<genome>.<config>_<target>` (e.g. `seq.IQ.pcg.flashzoi.mm10.rf524k_EB4_cnt`). See "Deep-learning bridge" below |
| Per-CpG methylation from bismark `.cov.gz` (per sample) | Read TSV in R, build the 4-track quartet (`.cov` / `.meth` / `.unmeth` / `.avg`) via `gtrack.create_sparse`. See "Methylation tracks" below | bismark is 1-based, misha is 0-based half-open — convert. Always CG-validate against the genome before importing |
| 450k / EPIC array methylation matrix (TCGA-style wide TSV) | `gtrack.array.import(name, desc, ...)` + `gtrack.var.set` for per-sample metadata | varargs: each `...` entry is either a wide TSV (chrom/start/end + N value columns, one per sample) or another array track; columns from multiple sources are merged by interval. Use `gtrack.array.extract` to query later |
| Wide feature matrix (non-methylation, e.g. signature scores) | `gtrack.array.import(name, desc, ...)` + `gtrack.var.set` for per-feature metadata | same varargs shape as above - one wide TSV per source |
| Track from another misha DB (same assembly) | `gtrack.copy(src, dest, db = "/path/to/other/trackdb", overwrite = FALSE)` | format conversion (per-chrom ↔ indexed) and chrom-order remap handled automatically; vector `src` supported |
| Intervals from another assembly (liftover via chain) | `gintervals.liftover(intervs, chain = "<srcToTgt>.over.chain.fixed1")` after `gsetroot(target_assembly)` | lifted intervals may extend past chrom ends and may overlap each other — always follow with `gintervals.force_range()` and (if you need disjoint output) `gintervals.canonic()`. See "Liftover" below |
| Track from another assembly (liftover via chain) | `gtrack.liftover(track, desc, src.track.dir, chain, multi_target_agg = "mean")` | requires a chain (`gintervals.load_chain`); aggregation policy controls how multiple source bins folded into one target bin combine — pick deliberately. Much less used in the lab than the intervals-side liftover above |
| Intervals (BED / GFF / GTF / VCF) as an interval set, not a track | `gintervals.import_bed` / `gintervals.import_genes` / `gintervals.import_gff` / `gintervals.import_vcf` | for a sparse-track route instead, build the data frame in R and pass to `gtrack.create_sparse` |

## Pileup tracks: prefer the in-R path

For ChIP/ATAC/CUT&Tag pileup tracks, prefer `gtrack.create_dense(..., func = "coverage")` over the legacy `gtrack.import_mappedseq` round-trip. With `values = rep(1, nrow(reads))`, bin value = `sum(overlap_i) / binsize` = average per-base read count — exactly the ChIP-seq pileup definition, in one C++ pass over the data frame. Use `gtrack.import_mappedseq` only when reads are too large to load into R or when you specifically want the `pileup =` read-extension feature.

`gtrack.create_dense` `func` choices (5.6.31+): `"weighted.mean"` (default), `"weighted.sum"`, `"coverage"`, `"max"`, `"min"`, `"median"`, `"count"`. Note that 5.6.32 fixed a bug where overlapping intervals of mixed lengths in the same bin gave plausible-looking but wrong means under the old default — re-import affected tracks if they predate 5.6.32.

## Methylation tracks: the 4-track quartet convention

Lab-wide convention (high-frequency, cross-user — 95+ files / 23+ projects in the corpus): every per-sample methylation dataset materialises as **four sparse tracks** sharing a common stem, not one:

| Track | Value per CpG |
|---|---|
| `<sample>.cov`    | total reads at the CpG (`meth + unmeth`) |
| `<sample>.meth`   | methylated read count |
| `<sample>.unmeth` | unmethylated read count |
| `<sample>.avg`    | methylation level `meth / cov` (often materialised; sometimes recomputed as an expression-derived vtrack at the CpG iterator) |

Downstream tooling (lab helpers like `combine.libs` / `copy.lib`, methylation-aware `misha.ext` wrappers) walks the four suffixes by convention. Materialising only `.avg` and skipping the count pair breaks coverage-aware aggregations and re-binning.

### Importing from a bismark `.cov.gz` (per sample)

The lab's current methylation entry point is bismark's `.cov` output piped into R. (The older `gpatterns`-based BAM-import path is deprecated and not used in recent projects; new methylation cohorts go through this `.cov` route.)

Bismark coverage files are 1-based, six columns: `chrom, start, end, avg, meth_count, unmeth_count` (start == end). Misha is 0-based half-open. **Always CG-validate against the reference** before importing — bismark coordinates can sit on either strand and ambiguous calls land on non-CG positions:

```r
library(data.table); library(dplyr); library(glue)

df <- fread(filename,
            col.names = c("chrom", "start", "end", "avg", "meth", "unmeth")) %>%
    select(-avg) %>%
    filter(chrom %in% gintervals.all()$chrom) %>%
    # 1-based -> 0-based half-open: written as two mutates so the second
    # `start + 1` uses the post-shift `start`, not the original bismark value.
    mutate(start = start - 1) %>%
    mutate(end = start + 1)

# Anchor to the C of the CpG: if the reference base at start is 'C' keep it; if it's
# 'G' (= bismark called the minus strand) shift back by one so start points at the C.
df <- df %>%
    mutate(s = toupper(gseq.extract(.))) %>%
    mutate(start = ifelse(s == "C", start, start - 1)) %>%
    mutate(end = start + 1)

# Validate that the 2bp dinucleotide at the anchored position really is CG;
# drop everything else (strand-ambiguous calls, soft-masked, sequencing errors).
df1 <- df %>% mutate(end = start + 2) %>% mutate(s = toupper(gseq.extract(.)))
df  <- df %>% filter(df1$s == "CG") %>%
    group_by(chrom, start, end) %>%
    summarise(meth = sum(meth), unmeth = sum(unmeth), .groups = "drop") %>%
    mutate(cov = meth + unmeth, avg = meth / cov) %>%
    filter(cov > 0)

# Materialise the quartet
for (suffix in c("cov", "meth", "unmeth", "avg")) {
    tr <- glue("{track_stem}.{suffix}")
    gtrack.create_sparse(track       = tr,
                         description = sprintf("%s %s", description, suffix),
                         intervals   = df[, c("chrom", "start", "end")],
                         values      = df[[suffix]])
}
```

Three points the corpus emphasises and that bite if skipped:

1. **`start = bismark_pos - 1`.** Bismark uses 1-based coordinates; misha is 0-based. Off-by-one silently shifts every CpG one base.
2. **The `gseq.extract` CG validation pair.** First call anchors to the C of the CpG; second call confirms the 2bp dinucleotide really is CG. Drops strand-ambiguous and miscalled rows. Skipping this leaks non-CG noise into the track.
3. **Group-summarise before creating.** The same CpG can appear twice in a bismark file (plus/minus calls both pointing to the same C after anchoring); deduplicate by summing counts.

### Parallelising over many samples

For a cohort, wrap the per-sample block in a function and dispatch via `misha.ext::gcluster.run2`:

```r
import_one <- function(filename, track_stem, description) {
    # ... the per-sample body above ...
}

cmds <- samples %>%
    mutate(cmd = glue("import_one('{cov_path}', '{track_stem}', '{description}')")) %>%
    pull(cmd)
misha.ext::gcluster.run2(command_list = cmds, max.jobs = 25, threads = 5, io_saturation = 1)
gdb.reload()
```

`gdb.reload()` afterwards is mandatory - `gcluster.run2` writes the new tracks from worker processes, and the parent gdb won't see them until reload.

### Population matrices: 450k / EPIC / PBAT cohorts

When you have a per-position × per-sample matrix (TCGA 450k arrays, PBAT cohort summaries, EWAS outputs), don't create N quartets — use the array form:

```r
# Wide TSV: chrom, start, end, sample1, sample2, ..., sampleN
gtrack.array.import("meth.tcga_brca", "TCGA BRCA 450k beta values", "/path/to/wide.tsv")
gtrack.var.set("meth.tcga_brca", "donors", donors_df)     # per-sample metadata
```

Query with `gtrack.array.extract(track, names_or_indices, intervals)` — pulls a sample-subset slice in one call. Build the wide TSV from N per-sample sparse `.avg` tracks via `gextract(c(samples), gintervals.all(), file = tempfile())` + `gtrack.array.import(... file = that_tempfile)` (`gextract` writes the right shape directly).

## RE fragment-track bootstrap (4C/Hi-C prerequisites)

Before any 2D contact import, the destination DB needs the restriction-enzyme fragment family for the RE used in library prep. Standard layout: `redb.<RE>_{flen, gc, map}` per-fragment 1D tracks plus a `fe<RE>_*` family produced from them. Lab corpus: 90 files across 25 projects — every new 4C/Hi-C bootstrap touches this.

```r
# Per-fragment tracks (binsize = 0 -> sparse, one value per RE fragment).
# flen_file / gc_file are produced by an external perl pipeline
# (e.g. map3c/TG3C/gen_re_frags.pl); not part of misha proper.
gtrack.import("redb.DpnII_flen", flen_file, binsize = 0)
gtrack.import("redb.DpnII_gc",   gc_file,   binsize = 0)

# Per-fragment mapability: aggregate a per-base mapability track using
# the flen track as iterator — each output bin is one RE fragment.
gtrack.create("redb.DpnII_map", mapab_per_base_track,
              iterator = "redb.DpnII_flen")

# Export back to text so re_frags_to_fends.pl can build the fends file
# consumed by gtrack.2d.import_contacts later.
gextract("redb.DpnII_map", gintervals.all(),
         iterator = "redb.DpnII_map", file = fragmap_file)
```

The `<redb>/<RE>.fends` file produced downstream (`GATC.fends` for DpnII, etc.) is the `fends =` argument of every `gtrack.2d.import_contacts` call that follows.

## scHi-C: pool many single-cell tracks into one 2D track

Universal pattern (166 files, 19 projects): partition the cell list across cores, each core `gextract`s its cells' contacts into one TSV stage, then a single `gtrack.2d.import_contacts` call with a **vector** of stage files pools everything. Avoids both N separate import calls and concatenating gigabytes through R.

```r
library(parallel); library(glue)

tmp_dir <- tempfile(); dir.create(tmp_dir)
num.cores <- 16
cells.partition <- split(cells, rep(1:num.cores, length.out = length(cells)))

mclapply(seq_along(cells.partition), function(cur.core) {
    # do.call(rbind, lapply(...)) - NOT a NULL-init growing accumulator.
    # An rbind-in-a-for-loop here is O(N^2) and is the standard scHi-C
    # anti-pattern flagged across the lab corpus.
    all.cell.conts <- do.call(rbind, lapply(cells.partition[[cur.core]], function(cell) {
        gextract(cell, gintervals.2d.all())[, 1:6]
    }))
    all.cell.conts$value <- 1
    write.table(format(all.cell.conts, scientific = FALSE),
                file.path(tmp_dir, cur.core),
                sep = "\t", quote = FALSE, row.names = FALSE)
}, mc.cores = num.cores)

gtrack.2d.import_contacts(pool_track, "Pooled scHi-C contacts",
                          contacts = file.path(tmp_dir, seq_len(num.cores)),
                          fends    = file.path(redb_dir, "GATC.fends"))
                          # allow.duplicates defaults to TRUE - the right
                          # value for pooling; setting FALSE here would
                          # error on the first cross-cell repeat.
```

Key points: `value = 1` per-contact (counts get aggregated by the importer); `format(..., scientific = FALSE)` is the 2D analogue of `options(scipen = 20)` for the TSV form - large coordinates serialised in scientific notation reject silently or partially. `allow.duplicates = TRUE` is correct (and is the function default - source: `R/track-2d.R:295`): when two cells contribute the same fend pair, the importer sums their counts, which is what pooling means. Setting `FALSE` makes the importer error on the first repeated fend pair, which kills any non-trivial pool. `FALSE` is appropriate only for pre-scored shaman per-rectangle imports, where each rectangle should be unique. The lab's canonical pooling idiom (e.g. `combine_adjs.pl` + import) simply omits the argument and relies on the default.

## Deep-learning bridge: model predictions → misha tracks

Borzoi / IceQream / similar deep-learning genome models round-trip through misha: the training half uses `gextract(..., intervals_join = "intervals")` (added in misha 5.7.1) to assemble per-bin training data with the input intervals' columns attached at the C++ layer; the inference half writes per-bin predictions back as a track. By lab convention, predicted tracks live under `seq.IQ.<model>.<genome>.<config>_<target>` (e.g. `seq.IQ.pcg.flashzoi.mm10.rf524k_EB4_cnt`).

Training-side extraction:

```r
borzoi_data <- gextract(borzoi_vt, intervals = regs, iterator = 32,
                        intervals_join = "intervals") %>%
    mutate(across(all_of(borzoi_vt), ~ ifelse(is.na(.), 0, .))) %>%
    group_by(chrom) %>%
    mutate(across(all_of(borzoi_vt), ~ iceqream::norm01(log2(1 + .))))
fwrite(borzoi_data, here("output/data-for-borzoi/borzoi_data.tsv"), sep = "\t")
```

`misha.ext::gextract.left_join` is the older R-side wrapper that does the same join via `dplyr::left_join`; deprecated in misha 5.7.1 (the built-in is ~2x faster on wide-window workloads and drops the misha.ext dependency). Don't introduce new calls.

Prediction-side import (after the external model run produces per-bin TSVs):

```r
options(scipen = 20)
gtrack.import("seq.IQ.pcg.flashzoi.mm10.rf524k_EB4_cnt",
              description = "Flashzoi mm10 rf524k EB4_cnt predictions",
              file        = pred_tsv,
              binsize     = 32)         # bin width must match the model's output resolution
stopifnot(gtrack.exists("seq.IQ.pcg.flashzoi.mm10.rf524k_EB4_cnt"))
```

The DL-prediction track behaves like any other dense track downstream — `gcor`, `gscreen`, peak overlap, vtrack composition all work. The `binsize` must match the model's per-bin resolution exactly (typically 1, 20, 32, or 128); a binsize mismatch silently re-bins predictions onto the wrong grid.

## Liftover: intervals vs tracks

Two related operations. By corpus frequency the intervals form dominates (~76 files across 17+ projects: ronisto, nettam, atanay, nimrodra, ofirr, aviezerl, effi, evghenic). `gintervals.liftover` was substantially extended in the last year — picking the right policy combo is now the load-bearing decision.

### Canonical idiom: lift then recenter

```r
gsetroot("/path/to/mm9")   # target assembly
lifted <- gintervals.liftover(intervs_mm10, chain = "mm10ToMm9.over.chain.fixed1")

# Recenter to a fixed 300bp window so stretched/shrunk intervals don't confound downstream extracts:
mid          <- round((lifted$start + lifted$end) / 2)
lifted$start <- mid - 150
lifted$end   <- mid + 150
lifted       <- gintervals.force_range(lifted)        # clamp to chrom bounds
# If overlap-free output is required (e.g. for sparse-track creation):
# lifted    <- gintervals.canonic(lifted)
```

Two non-negotiable post-liftover steps: `gintervals.force_range` (lifted intervals can extend past chrom ends - `gintervals()` then rejects them with `end coordinate exceeds chromosome boundaries`) and, if disjoint output is required, `gintervals.canonic` (lifted intervals may overlap). `gintervals.canonic` defaults to `unify_touching_intervals = TRUE`, which merges intervals that share an endpoint as well as overlapping ones - pass `unify_touching_intervals = FALSE` if you need to preserve adjacency as separate rows.

Output always carries `intervalID` (index into the input) and `chain_id` (which chain produced the mapping). The `chain_id` is essential when a source interval lifted via multiple chains (duplications): two rows with the same `intervalID` but different `chain_id` are distinct mappings.

### Policy chooser

**`src_overlap_policy` — source position mapped by multiple chains (paralogs / duplications):**

| Value | When to use |
|---|---|
| `"error"` (default) | Strict: fail if input contains source duplications. Right for curated chains and 1:1 lift between close assemblies (mm9↔mm10, hg19↔hg38) |
| `"keep"` | Cross-species / paralog-rich (e.g. cross-primate). One source → many output rows; distinguish by `chain_id` |
| `"discard"` | Drop any source that has duplications. Use when downstream tooling can't tolerate ambiguity |

**`tgt_overlap_policy` - multiple chains converging on overlapping target loci:**

| Value | When to use |
|---|---|
| `"error"` | Strict: fail if any target overlap is detected. Use to catch unexpected chain-file structure before downstream work assumes there are none |
| `"auto"` / `"auto_score"` (default) | Segment target overlaps and pick the highest-scoring chain per segment. Right default for almost everything |
| `"auto_longer"` | Prefer the longer chain per segment. Useful when alignment scores are noisy or uniform across the chain file |
| `"auto_first"` | Prefer the lowest `chain_id` per segment (= original chain-file order). Use when the chain file encodes ranked preference |
| `"keep"` | Leave target overlaps untouched (downstream must handle them) |
| `"discard"` | Drop any chain involved in a target overlap |
| `"agg"` | Segment into disjoint regions and retain all contributors per region. **Pair with `value_col` + `multi_target_agg`** for proper value aggregation |
| `"best_source_cluster"` | Cluster chains by source-side overlap: keep every "true duplication" cluster (chains whose source intervals overlap) but pick only the largest "conflicting alternative" cluster (disjoint source intervals). For paralog-aware lifts where you want duplications but not alternative competing mappings |
| `"best_cluster_union"` / `"best_cluster_sum"` / `"best_cluster_max"` | Cluster-best variants for specialised cases (cluster scoring by union / sum / max). Most users want `"best_source_cluster"` or `"auto_score"` |

**`min_score = N`** — filter out chains with alignment score below N. Useful with public chain files that include low-confidence alignments (e.g. liftOver's chains from older UCSC builds).

### Lifting a value column

When source intervals carry a value (peak score, methylation level, deep-learning prediction, …) and you need that value in the output, pass `value_col`. Pair with `tgt_overlap_policy = "agg"` so target-side conflicts are segmented and contributors aggregated:

```r
gintervals.liftover(intervs_with_score, chain,
                    value_col          = "score",
                    tgt_overlap_policy = "agg",
                    multi_target_agg   = "mean",      # or sum / max / max.coverage_len / ...
                    na.rm              = TRUE,
                    min_n              = 1)            # require >=1 non-NA contributor per row
```

`multi_target_agg` choices: `mean` (default) / `median` / `sum` / `min` / `max` / `count` / `first` / `last` / `nth` (set `params = N`) / `max.coverage_len` / `min.coverage_len` / `max.coverage_frac` / `min.coverage_frac`. Pick deliberately — the default is `"mean"`, which is wrong for count-like or coverage-weighted values. The `na.rm` and `min_n` knobs are only consulted when `value_col` is set.

### Other flags

- `include_metadata = TRUE` adds a `score` column to the output (from the winning chain per row). Only meaningful with `tgt_overlap_policy = "auto"` / `"auto_score"`. Use when you want to threshold or weight lifted intervals by alignment quality post-lift.
- `canonic = TRUE` merges adjacent target intervals that share the same `(intervalID, chain_id)` — collapses chain-gap splits into one row per source/chain. Cheaper than running `gintervals.canonic` after the fact when that's the only canonicalisation you need.

### Pre-loaded chains + hand-built chains

`gintervals.load_chain(file, src_overlap_policy, tgt_overlap_policy, src_groot, min_score)` parses a chain file once and stamps the policies as attributes on the returned data.frame. Reusing it across many `gintervals.liftover` calls is much faster than re-parsing per call. Once policies are baked in this way, `gintervals.liftover` rejects attempts to override them — set them at load time.

`src_groot` is a path to the source-assembly groot; when provided, `gintervals.load_chain` validates that chain source chroms and coordinates exist there. Cheap safety net before a long batch.

`gintervals.as_chain(intervals, src_overlap_policy, tgt_overlap_policy, min_score)` promotes an arbitrary data.frame with the chain columns (`chrom, start, end, strand, chromsrc, startsrc, endsrc, strandsrc, chain_id, score`) into a chain — for synthesised or hand-edited chains.

### gtrack.liftover

Whole-track form. Shares the `multi_target_agg` enum with `gintervals.liftover` (no `value_col` argument — the track's value IS the value column). Same caveat: the default `"mean"` can hide signal when many source bins fold to one target locus; for coverage / count tracks prefer `"sum"` or one of `max.coverage_*` / `min.coverage_*`. On indexed destination DBs the writer skips per-chrom placeholder files (5.6.30+).

## Validating a multi-file concatenated source before import

If the input is a `cat`-style concatenation of pipeline chunks, validate before launching `gtrack.import` (on large inputs it can take 10+ minutes before erroring):

**Trailing newlines on every chunk.** A missing `\n` glues a value line to the next chunk's `fixedStep` header → invalid WIG. Detect and normalize:

```bash
# detect chunks missing terminal newline
for f in chunk_*.wig; do
  [ "$(tail -c1 "$f" | od -An -c | tr -d ' ')" != '\n' ] && echo "MISSING NL: $f"
done

# concatenate with numeric chunk order + newline-normalization in one pass
ls chunk_*.wig | sort -V | xargs awk 1 > combined.wig
```

`awk 1` re-emits every line with a terminating `\n`, so missing newlines are fixed in-flight; `sort -V` keeps `chunk_0, chunk_1, ..., chunk_10` in numeric order (not the lexical `0, 1, 10, 11, ..., 2` that plain `ls` or shell globbing gives).

**No malformed lines.**

```bash
LC_ALL=C grep -cvE '^(fixedStep|variableStep|track|#|-?[0-9.eE+-]+$|$)' combined.wig   # should print 0
```

**No chrom revisits.** Each known chrom must appear contiguously across the concat — the WIG parser rejects with `not sorted by chromosomes` otherwise. Sanity-check by walking the headers (`grep -nE '^(fixed|variable)Step' | awk ...`).

**No chunk overlaps within a chrom.** For consecutive same-chrom chunks, require `next_start_1based - 1 >= prev_start_1based - 1 + prev_value_count`. Overlaps surface at `gextract` time (not import), with `overlapping intervals` from the WIG `get_data` path.

## Post-import sanity

Sample multiple chromosomes — a single read at `chr1:0-1000` can be legitimately zero in continuous-signal tracks and hides chrom-name mismatches:

```r
stopifnot(gtrack.exists(name))
gtrack.info(name)            # check bin.size / dimensions / size.in.bytes vs expectation

# Sample ~1000 random positions on each of the first few chroms in the gdb.
# Driving chrom selection from gintervals.all() keeps the check portable to any genome.
all_chr <- gintervals.all()
for (i in seq_len(min(5, nrow(all_chr)))) {
    c <- all_chr$chrom[i]; cend <- all_chr$end[i]
    starts <- sort(sample.int(max(cend - 100, 1L), min(1000, cend - 100)))
    pts <- gintervals(c, starts, starts + 100)
    v <- gextract(name, pts)[[name]]
    cat(sprintf("%-12s nz=%d/%d mean=%.4f range=[%.3f, %.3f]\n",
                c, sum(v != 0, na.rm = TRUE), length(v),
                mean(v, na.rm = TRUE), min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
}
```

All-zero across every chrom you expect signal on -> chrom-name mismatch (revisit the precondition). Non-zero with a sane range -> import succeeded.

## Post-import metadata (provenance)

Every imported track should carry enough metadata for someone (you, in six months; a collaborator; a paper companion repo) to figure out *what it is* without rereading the import script. The lab convention is `gtrack.attr.set` with a small, consistent attribute vocabulary:

```r
gtrack.attr.set(track, "source",      "ENCODE")               # or GEO, in-house, ...
gtrack.attr.set(track, "experiment",  "chip-seq")             # chip-seq, atac, wgbs, hic, ...
gtrack.attr.set(track, "protein",     "Oct4")                 # for chip-seq / cut&run
gtrack.attr.set(track, "PMID",        "28212747")
gtrack.attr.set(track, "data_link",   "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1910646")
gtrack.attr.set(track, "liftover",    "mm9")                  # source assembly if lifted
gtrack.attr.set(track, "pileup",      200)                    # for BAM/mappedseq imports
```

`created.by` and `created.date` are populated automatically. List everything attached to a track with `gtrack.attr.export(track)`; query one with `gtrack.attr.get`.

Use `gtrack.var.set(track, var, value)` (not `attr.set`) when you need to attach a structured R object - a donors data frame for an array track, a per-sample metadata table, an arbitrary list. `attr` is for short character/numeric labels; `var` is for R objects.

## Failure modes

| Symptom / error | Cause | Fix |
|---|---|---|
| `Unrecognized format of file <X>` (misha < 5.7.2) | Real WIG/CSV parser error is being swallowed by a two-step fallback | Upgrade to misha ≥ 5.7.2; the real diagnostic (file + line) is then surfaced |
| `Invalid format of WIG file ..., line N` | Malformed line N — often a value glued to a `fixedStep` header from a cat-without-newline concat | Re-concat with `awk 1` (see Validation §1) |
| `not sorted by chromosomes` | Concat order interleaves chroms across files | `sort -V` chunk files before `cat`; verify with a header pass |
| `overlapping intervals` (at `gextract`, not import) | Same-chrom chunks overlap | Validation §4 |
| `func argument is not a string` (BED + `binsize`, misha < 5.7.1) | Internal BED-to-dense path missed a required arg | Upgrade to misha ≥ 5.7.1 |
| Import succeeds, `gextract` returns all 0 or all NaN | Chrom names didn't match gdb (silent skip); or `defval = 0` with uncovered bases on a continuous-signal source | Precondition; use `defval = NaN` for continuous signal |
| `gtrack.import_mappedseq` track gives 5–20× expected counts at peaks | Legacy call site explicitly passed `remove.dups = FALSE` (current default is `TRUE`) | Drop the explicit `FALSE`; rely on the safe default |
| TSV import errors with non-numeric / parse error on coords | R serialized large coordinates in scientific notation | `options(scipen = 20)` BEFORE `write.table` |
| Re-imports under the same name behave inconsistently | Importer doesn't fully overwrite | `gtrack.rm(name, force = TRUE)` before re-import |
| `gtrack.create_dense` `func = "coverage"` returns implausibly large values (misha < 5.6.32) | Overlapping intervals of mixed lengths inflated the per-bin coverage | Upgrade to misha ≥ 5.6.32 |
| `gtrack.liftover` produces NA-heavy track | Many source bins mapped to overlapping target loci and `multi_target_agg` discarded most | Try `multi_target_agg = "mean"` / `"max.coverage_len"` deliberately; inspect the chain coverage |
| `gtrack.copy(..., db = X)` errors on 2D | Chrom orders differ between source and destination | Reorder the destination DB, or stage via R (`gextract` + `gtrack.2d.create`) |
| `BAD_CHROM` (paths that don't ignore unknown chroms) | Source has chroms the gdb doesn't know | Extend the gdb (see `gdb.install_intervals` / `gdb.build_genome`), add aliases, or filter input |

## Common mistakes

- Mistaking `gtrack.import_set` (current bulk importer) for the legacy `gtrack.import.wigs` name. The dot-form is not an exported function in current misha — calls fail at the R level. Use `gtrack.import_set(description, path, binsize, track.prefix)`.
- Round-tripping data already in R through a temp TSV. `gtrack.create_dense` / `gtrack.create_sparse` take the data.frame directly — faster and avoids `scipen` traps.
- Using `gtrack.import_mappedseq` for reads already loaded in R, instead of `gtrack.create_dense(..., func = "coverage")`. The latter is the modern one-call path.
- Using `defval = 0` for continuous signal (phyloP, normalized ratios). Downstream `is.na` filters become useless; everything looks like "data".
- Skipping post-import sampling because the function returned without error. The most common silent failure is chrom-name mismatch → all-NaN track.
- Treating `gtrack.copy` as a same-DB-only operation. It supports cross-DB copy via `db = "/other/groot"` (5.6.28+), including format conversion.

## Cross-references

- `agent-guides/misha-advanced.md` §5 — short recipe view of the same task. This skill is the full reference.
- `agent-guides/misha-core.md` §1.1–1.2 — what a trackdb / track is on disk.
- `R/track-import.R` / `R/track-management.R` / `R/track-liftover.R` / `R/track-2d.R` / `src/GenomeTrackImportWig.cpp` — implementation.
