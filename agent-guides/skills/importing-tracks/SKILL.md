---
name: importing-tracks
description: Use when importing external data (WIG / bigWig / BED / TSV / BAM / 2D contacts / array matrices) into a misha track DB, when picking which `gtrack.*import*` / `gtrack.create_*` / `gtrack.liftover` / `gtrack.copy` variant to call, or when a track import fails opaquely or produces an all-zero / all-NaN track.
---

# Importing tracks into misha

## Overview

A "track import" is any operation that materializes a persistent on-disk track from an external source — a single file, a directory or URL of files, a BAM, an in-R data frame, a track in another misha DB, or a track in another genome assembly. There are eight+ entry points; picking the wrong one is the most common ingestion problem after silent chrom-name mismatches.

Two failure modes worth pre-empting:

1. **Picking the wrong variant** for the source shape.
2. **Silent partial imports.** WIG/bigWig importers run with `ignore_unknown_chroms = TRUE` — chroms not in `gintervals.all()$chrom` are dropped without error. Mismatched naming (`1` vs `chr1`) yields an "all NaN" track, not a failure.

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

| Source | Function | Key args / caveats |
|---|---|---|
| In-memory R intervals + values | `gtrack.create_dense(name, desc, intervals, values, binsize, func = "weighted.mean"/"coverage"/"max"/"min"/"sum"/"median"/"count")` (5.6.31+) or `gtrack.create_sparse(name, desc, intervs, values)` | preferred whenever data is already in R — no TSV round-trip and no `scipen` traps |
| Single WIG / bigWig file | `gtrack.import(name, desc, file, binsize, defval = NaN)` | auto-dispatch on `.wig` / `.bw` / `.bigWig` / magic bytes; pass `defval = NaN` (not `0`) for continuous signal so uncovered bins read as NA |
| Single BED file (or 4-col TSV) | same `gtrack.import` — pass `binsize` for a dense track, omit it for a sparse track | misha ≥ 5.7.1 required if you pass `binsize` (5.7.0 errors with `func argument is not a string`) |
| R-computed values via TSV (chrom/start/end/value) | `gtrack.import(name, desc, file, binsize)` | run `options(scipen = 20)` BEFORE `write.table` — scientific-notation coords parse as malformed |
| Directory or URL-with-wildcards of WIG/bigWig files | `gtrack.import_set(description, path, binsize, track.prefix = "", defval = NaN)` | bulk import; `path` may be a glob (`/data/*.bw`) or an `ftp://` URL; continues on per-file errors and returns successes/failures |
| ChIP/ATAC pileup from read intervals (data.frame in R) | `gtrack.create_dense(name, desc, intervals = reads, values = rep(1, nrow(reads)), binsize = 20, defval = 0, func = "coverage")` | one-call replacement for the old `gtrack.import_mappedseq` route when reads are already in R |
| Mapped reads (TSV from `samtools view`, or pre-staged file) | `gtrack.import_mappedseq(name, desc, file, pileup, binsize, cols.order = c(9, 11, 13, 14), remove.dups = TRUE)` | for BAM input, pipe through `samtools view` first; the default `remove.dups = TRUE` is correct for ChIP/ATAC/CUT&RUN |
| 2D Hi-C / capture-C contacts | `gtrack.2d.import_contacts(name, desc, contacts, fends, allow.duplicates = FALSE)` | `contacts =` accepts a vector (e.g. shaman per-rect scores are concatenated by the importer) |
| 2D track from a per-rectangle R data.frame | `gtrack.2d.create(name, desc, intervals, values)` | rectangles must be non-overlapping |
| Wide feature matrix (e.g. 450k array) | `gtrack.array.import(name, desc, file)` + `gtrack.var.set` for per-feature metadata | one wide TSV: chrom/start/end + N value columns |
| Track from another misha DB (same assembly) | `gtrack.copy(src, dest, db = "/path/to/other/trackdb", overwrite = FALSE)` | format conversion (per-chrom ↔ indexed) and chrom-order remap handled automatically; vector `src` supported |
| Intervals from another assembly (liftover via chain) | `gintervals.liftover(intervs, chain = "<srcToTgt>.over.chain.fixed1")` after `gsetroot(target_assembly)` | lifted intervals may extend past chrom ends and may overlap each other — always follow with `gintervals.force_range()` and (if you need disjoint output) `gintervals.canonic()`. See "Liftover" below |
| Track from another assembly (liftover via chain) | `gtrack.liftover(track, desc, src.track.dir, chain, multi_target_agg = "mean")` | requires a chain (`gintervals.load_chain`); aggregation policy controls how multiple source bins folded into one target bin combine — pick deliberately. Much less used in the lab than the intervals-side liftover above |
| Intervals (BED / GFF / GTF / VCF) as an interval set, not a track | `gintervals.import_bed` / `gintervals.import_genes` / `gintervals.import_gff` / `gintervals.import_vcf` | for a sparse-track route instead, build the data frame in R and pass to `gtrack.create_sparse` |

## Pileup tracks: prefer the in-R path

For ChIP/ATAC/CUT&Tag pileup tracks, prefer `gtrack.create_dense(..., func = "coverage")` over the legacy `gtrack.import_mappedseq` round-trip. With `values = rep(1, nrow(reads))`, bin value = `sum(overlap_i) / binsize` = average per-base read count — exactly the ChIP-seq pileup definition, in one C++ pass over the data frame. Use `gtrack.import_mappedseq` only when reads are too large to load into R or when you specifically want the `pileup =` read-extension feature.

`gtrack.create_dense` `func` choices (5.6.31+): `"weighted.mean"` (default), `"weighted.sum"`, `"coverage"`, `"max"`, `"min"`, `"median"`, `"count"`. Note that 5.6.32 fixed a bug where overlapping intervals of mixed lengths in the same bin gave plausible-looking but wrong means under the old default — re-import affected tracks if they predate 5.6.32.

## Liftover: intervals vs tracks

Two related but distinct operations. By corpus frequency, the intervals form is dominant (~76 files, 17+ projects, used by ronisto, nettam, atanay, nimrodra, ofirr, aviezerl); the track form is much rarer.

**Intervals across assemblies** — `gintervals.liftover(intervs, chain = "<srcToTgt>.over.chain.fixed1")`. Standard lab idiom: switch to the target assembly first, lift, then recenter to a fixed width and clamp:

```r
gsetroot("/path/to/mm9")   # target assembly
lifted <- gintervals.liftover(intervs_mm10, chain = "mm10ToMm9.over.chain.fixed1")

# Recenter to a fixed 300bp window so stretched / shrunk intervals don't
# confound downstream extracts:
mid          <- round((lifted$start + lifted$end) / 2)
lifted$start <- mid - 150
lifted$end   <- mid + 150
lifted       <- gintervals.force_range(lifted)   # clamp to chrom bounds
# If overlap-free output is required (e.g. for sparse-track creation):
# lifted    <- gintervals.canonic(lifted)
```

Two non-negotiable post-liftover steps: `gintervals.force_range` (lifted intervals can extend past chrom ends — `gintervals()` then rejects them with `end coordinate exceeds chromosome boundaries`) and, if you need a non-overlapping set, `gintervals.canonic` (lifted intervals may overlap). Chain files are basename-passed and resolved via misha's chain-search path; `gintervals.load_chain(file)` lets you preload + inspect a chain.

**Whole tracks across assemblies** — `gtrack.liftover(track, desc, src.track.dir, chain, multi_target_agg = "mean")`. Reads every bin in the source track, projects via the chain, and aggregates contributors landing on the same target locus per `multi_target_agg` (`mean` / `median` / `sum` / `min` / `max` / `max.coverage_len` / etc.). The aggregation policy is load-bearing: defaults can drop most of the signal if many bins fold into one. On indexed destination DBs the writer skips the per-chrom placeholder files (5.6.30+).

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

All-zero across every chrom you expect signal on → chrom-name mismatch (revisit the precondition). Non-zero with a sane range → import succeeded.

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
