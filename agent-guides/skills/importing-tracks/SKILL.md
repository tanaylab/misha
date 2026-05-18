---
name: importing-tracks
description: Use when importing external data (WIG / bigWig / BED / TSV / BAM / 2D contacts / array matrices) into a misha track DB, when picking which `gtrack.*import*` variant to call, or when a track import fails opaquely or produces an all-zero / all-NaN track.
---

# Importing tracks into misha

## Overview

`gtrack.import` and its siblings turn external files into persistent on-disk tracks queryable via `gextract` and vtracks. The two failure modes worth pre-empting:

1. **Picking the wrong variant** for the source shape — there are 6+ entry points, each tuned to a different input.
2. **Silent partial imports** — WIG/bigWig importers run with `ignore_unknown_chroms = TRUE`, so chrom-name mismatches yield an "all NaN" track with no error.

## Precondition: chrom-name compatibility (mandatory)

Before any import, verify input chroms resolve in the active gdb:

```r
gdb.init("/path/to/trackdb")
setdiff(unique(input$chrom), gintervals.all()$chrom)   # MUST be empty for chroms you care about
```

For huge WIG/bigWig inputs you can't load into R, extract chroms from headers:

```bash
grep -oE 'chrom=[^ ]+' file.wig | sort -u | sed 's/chrom=//'   # WIG
bigWigInfo -chroms file.bw                                     # bigWig
```

Non-empty diff → wrong gdb, naming-convention mismatch (`1` vs `chr1`), or missing scaffolds. Fix before importing. Silent partial imports are the #1 cause of "the track is all zeros".

## Format chooser

Ordered by lab-corpus frequency (most used first).

| Source | Function | Key args / caveats |
|---|---|---|
| URL list of bigWigs (ENCODE / public data) | `gtrack.import.wigs(urls, track.prefix=, binsize=)` | one track per URL; lab-internal `import.tracks()` wrappers add `mkdir`+`gsetroot`-by-genome boilerplate |
| R-computed values via TSV (chrom/start/end/value) | `gtrack.import(name, desc, file, binsize=res)` | run `options(scipen = 20)` BEFORE `write.table` — scientific-notation coords produce opaque parser errors |
| In-memory R intervals + values | `gtrack.create_dense(..., func = "coverage")` for pileups; `gtrack.create_sparse(intervs, values)` for sparse | no TSV round-trip; preferred whenever data is already in R |
| Single WIG / bigWig file | `gtrack.import(name, desc, file, binsize=, defval = NaN)` | auto-dispatch on `.wig` / `.bw` / `.bigWig` / magic bytes; pass `defval = NaN` (not `0`) for continuous signal so uncovered bins read as NA |
| Mapped reads (TSV from `samtools view`) | `gtrack.import_mappedseq(..., pileup =, binsize =, cols.order =, remove.dups = TRUE)` | `remove.dups = TRUE` is mandatory for ChIP/ATAC/CUT&RUN — duplicates inflate every threshold call |
| BAM directly | `gtrack.import_mappedseq_bam(...)` | same caveats; saves the `samtools view` step |
| Sparse track from a 4-col BED | `gintervals.canonic()` + `tapply(values, key, mean)` → `gtrack.create_sparse` | canonicalize first if intervals can overlap |
| 2D Hi-C / capture-C contacts | `gtrack.2d.import_contacts(name, desc, contacts =, fends =, allow.duplicates = FALSE)` | `contacts =` accepts a vector (e.g. shaman per-rect scores get concatenated by the importer) |
| Wide feature matrix (e.g. 450k array) | `gtrack.array.import(...)` + `gtrack.var.set` for metadata | one wide TSV: chrom/start/end + N value columns |

## Validating a multi-file concatenated source before import

If the input is a `cat`-style concatenation of pipeline chunks, validate before launching `gtrack.import` (which on large inputs can take 10+ minutes before erroring):

**Trailing newlines on every chunk.** A missing `\n` glues a value line to the next chunk's `fixedStep` header → invalid WIG. Detect and normalize:

```bash
# detect chunks missing terminal newline
for f in chunk_*.wig; do
  [ "$(tail -c1 "$f" | od -An -c | tr -d ' ')" != '\n' ] && echo "MISSING NL: $f"
done

# concatenate with numeric chunk order + newline-normalization
ls chunk_*.wig | sort -V | xargs awk 1 > combined.wig
```

`awk 1` re-emits every line with a terminating `\n`, so missing newlines are fixed in-flight; `sort -V` keeps `chunk_0, chunk_1, chunk_2, ..., chunk_10` in numeric (not lexical) order.

**No malformed lines.**

```bash
LC_ALL=C grep -cvE '^(fixedStep|variableStep|track|#|-?[0-9.eE+-]+$|$)' combined.wig   # should print 0
```

**No chrom revisits.** Each known chrom must appear contiguously across the concat. WIG parser rejects with `not sorted by chromosomes` otherwise. Verify with an awk pass over the headers (`grep -nE '^(fixed|variable)Step'`).

**No chunk overlaps within a chrom.** For consecutive same-chrom chunks, require `next_start - 1 >= prev_start - 1 + prev_value_count`. Overlaps surface at `gextract` time (not import), with `overlapping intervals` from the WIG `get_data` path.

## Post-import sanity

```r
stopifnot(gtrack.exists(name))
gtrack.info(name)            # check bin.size / dimensions / size.in.bytes vs expectation

# Sample MULTIPLE chroms — a first-position read can be legitimately zero in continuous-signal tracks
for (c in c("chr1", "chr7", "chrX")) {
    pts <- gintervals(c, seq(0, 1e8, length.out = 1000), seq(1, 1e8 + 1, length.out = 1000))
    pts <- pts[pts$start < pts$end, ]
    v <- gextract(name, pts)[[name]]
    cat(sprintf("%-6s nz=%d/%d mean=%.4f range=[%.3f, %.3f]\n",
                c, sum(v != 0, na.rm = TRUE), length(v),
                mean(v, na.rm = TRUE), min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
}
```

All-zero across every chrom you expect signal on → chrom-name mismatch (revisit the precondition). Non-zero with a sane range → import succeeded.

## Failure modes

| Symptom / error | Cause | Fix |
|---|---|---|
| `Unrecognized format of file <X>` (misha < 5.7.2) | Real WIG/CSV parser error is being swallowed by a two-step fallback | Upgrade to misha ≥ 5.7.2; the real diagnostic (file + line) is then surfaced |
| `Invalid format of WIG file ..., line N` | Malformed line N — typically a value glued to a `fixedStep` header from a cat-without-newline concat | Re-concat with `awk 1` (Validation §1) |
| `not sorted by chromosomes` | Concat order interleaves chroms across files | `sort -V` chunk files before `cat`; verify with header pass |
| `overlapping intervals` (at `gextract`, not import) | Same-chrom chunks overlap | Validation §4 |
| `func argument is not a string` (BED + `binsize`, misha < 5.7.1) | Internal BED-to-dense path missed a required arg | Upgrade to ≥ 5.7.1 |
| Import succeeds, `gextract` returns all 0 or all NaN | Chrom names didn't match gdb (silent skip); or `defval = 0` with uncovered bases | Precondition; for continuous signal use `defval = NaN` |
| `gtrack.import_mappedseq` track gives 5–20× expected counts at peaks | `remove.dups = FALSE` (the default) | Re-import with `remove.dups = TRUE` |
| Re-imports under the same name behave inconsistently | Importer doesn't fully overwrite | `gtrack.rm(name, force = TRUE)` before re-import |
| TSV import errors with non-numeric / parse error on coords | R serialized large coordinates in scientific notation | `options(scipen = 20)` BEFORE `write.table` |
| `BAD_CHROM` (when `ignore_unknown_chroms = FALSE` paths) | Source has chroms the gdb doesn't know | Either extend the gdb (`gdb.create` with full chrom set), add aliases, or filter input |

## Common mistakes

- Treating `gtrack.import` and `gtrack.import.wigs` as synonyms. They're different functions: the first imports one file, the second iterates a URL/file vector under a `track.prefix`. Picking the wrong one wastes a loop or kills a batch.
- Writing R-computed TSVs without `options(scipen = 20)`. Scientific-notation coordinates parse as malformed.
- Round-tripping data already in R through a temp TSV. `gtrack.create_dense` / `gtrack.create_sparse` take the data.frame directly — faster and avoids `scipen` traps.
- Using `defval = 0` for continuous signal (phyloP, coverage normalised in ratios). Downstream `is.na` filters become useless; everything is "data".
- Skipping post-import sampling because the function returned without error. The most common silent failure is chrom-name mismatch → all-NaN track.

## Cross-references

- `agent-guides/misha-advanced.md` §3.13 — short recipe view of the same task. This skill is the full reference.
- `agent-guides/misha-core.md` §1.1–1.2 — what a trackdb / track is on disk.
- `R/track-import.R` / `R/track-2d.R` / `src/GenomeTrackImportWig.cpp` — implementation.
