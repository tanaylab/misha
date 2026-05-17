# misha - anti-patterns

The silent footguns - code runs, output looks plausible, downstream conclusions are subtly wrong. Inline `Avoid:` callouts in [`misha-core.md`](misha-core.md) and [`misha-advanced.md`](misha-advanced.md) reference entries here.

## Contents

- [Configuration](#configuration)
- [Virtual track lifecycle](#virtual-track-lifecycle)
- [Intervals](#intervals)
- [Extraction & expressions](#extraction--expressions)
- [Quantiles / thresholds](#quantiles--thresholds)
- [Aggregator choice](#aggregator-choice)

## Configuration

**A1. `options(gparam.type = 'string')` set in modern scripts.** A real switch in misha circa 2014 - bare track names required quoting only when this was set. Current misha accepts both quoted and bare forms transparently; setting the option is dead weight. Drop it from new code.

**A2. Mid-script `gsetroot` swap.** Switching genomes while in-memory intervals reference the previous genome silently corrupts coordinates - every subsequent `gintervals.load`, vtrack lookup, or sequence extract resolves against the new reference. If you must swap, wrap with `on.exit(gsetroot(.misha$GROOT))` so a partial run restores cleanly (see misha-advanced.md §6).

**A3. `gmultitasking = TRUE` inside `mclapply` / forked workers.** Misha's multitasking forks worker processes to parallelize across chromosomes. When the misha call is *itself* inside a forked context (`mclapply`, `parallel::mcMap`, `gcluster.run` body), the nested forks deadlock or oversubscribe cores. Fix: set `options(gmultitasking = FALSE)` inside the inner body and let the outer layer own the parallelism. Symptom: hung workers that never make progress and don't error.

**A4. High `gmax.data.size` combined with `gmultitasking = TRUE` → `cannot allocate memory`.** `gmax.data.size` caps the result size of *each* `gextract`; with multitasking, several extractions hold that much simultaneously. On a heavy genome-wide pull (`gmax.data.size = 1e10` × multiple worker copies) the resident set exceeds available RAM. Three fixes, pick by what's tightest: lower `gmax.data.size` to what one extract actually needs, set `gmultitasking = FALSE` for the heavy call, or partition the query by chromosome and process the partitions serially.

## Virtual track lifecycle

**A5. Defensive `gvtrack.rm` before every `gvtrack.create`.** Boilerplate of the form `if (length(gvtrack.ls(n)) == 1) gvtrack.rm(n)` (or `try(gvtrack.rm(n), TRUE)`) before every `gvtrack.create` is dead weight - `gvtrack.create` silently overwrites on name clash with another vtrack. The error only fires when the name clashes with a *regular track* or *intervals set*, which the defensive `gvtrack.rm` doesn't fix anyway.

## Intervals

**A6. Iterating `gintervals.normalize` + `gintervals.canonic` in a `while` loop until "stable".** Old code does `repeat { x <- gintervals.canonic(gintervals.normalize(x, W)); if (...) break }` to "settle" overlapping peaks. A single pass suffices by construction: `gintervals.normalize` writes width-`W` intervals centered on the input midpoints, and `gintervals.canonic` then merges overlaps in one sort+merge. The second iteration changes nothing because the input is already fixed-width and sorted.

**A7. `(start + end) / 2` for the midpoint.** Float division yields a non-integer center; downstream operations round inconsistently (sometimes via `as.integer`, sometimes via `floor`, sometimes via `round`). Always integer-divide: `(start + end) %/% 2L`.

**A8. Reflexive `gintervals.canonic` after every `rbind`.** Canonicalizing destroys per-row identity - paired matched-control rows, per-sample replicate identity, parallel rows for downstream joins all collapse into the merged set. Run `gintervals.canonic` *only* when you actually want overlaps merged (typical: union-of-peak-sets across conditions). When the rows mean "distinct things at the same location", skip canonic.

**A9. Base `data.frame(chrom, start, end)` instead of `gintervals(...)`.** `data.frame` produces unvalidated coordinates and turns `chrom` into a factor with stale levels (especially after subsetting). Equality filters then silently miss rows. Always pass through `gintervals(chrom, start, end)` which validates against the active genome and gives integer-typed coordinates.

## Extraction & expressions

**A10. Looping `gextract` over track names.** Calling `gextract` once per track and joining on `chrom/start/end` is N× slower than one `gextract(c("track1", "track2", ...), ...)`. The C++ engine processes all expressions jointly in one genome pass; the merge step is unnecessary work and a source of join bugs (factor levels, NA-vs-missing-row).

**A11. Post-extract `prof[is.na(prof)] <- 0` instead of inside the expression.** Sparse tracks (and many dense tracks) return NA at unmeasured bins. Filling at extract time via `gextract("ifelse(is.na(track), 0, track)", ...)` keeps the NA handling inside the C++ engine; the post-hoc base-R fill works but is slower and easier to forget.

**A12. `iterator = peaks` when you meant per-bin resolution.** With `iterator = <intervals>`, gextract collapses to one row per interval (with vtrack aggregation). For a per-bin profile inside the same intervals, use `iterator = 20` (or another bin size). For interval-relative bin indices use `giterator.intervals(... interval_relative = TRUE)`.

**A13. `iterator = 20` and assuming bin offsets are interval-relative.** Integer iterators are *chromosome-relative* - two peaks at different chromosome offsets get *misaligned* bin grids. If you need bin indices anchored to each peak's start, use `giterator.intervals(... interval_relative = TRUE)` (with the caveat that dense source tracks still align to their stored bin grid).

**A14. Skipping `arrange(intervalID)` (or `order(.$intervalID)`) on gextract output.** With `gmultitasking = TRUE`, `gextract` returns rows in chunked-chromosome order, not input order. If the downstream code indexes by row position into the original peaks frame, the result is silently shuffled. Always `arrange(intervalID) |> select(-intervalID)` or equivalent.

## Quantiles / thresholds

**A15. Hardcoded threshold values.** `gscreen("track > 1870", ...)` makes the analysis uninterpretable across datasets and silently breaks on retraining or track replacement. Derive the threshold from `gquantiles(track, percentiles = p, ...)` so the rule is "top p-percentile of the empirical distribution" - re-evaluable on any new track.

**A16. `gquantiles(x, probs = 0.9, ...)` instead of `percentiles =`.** R does not partial-match across function definitions; `probs` is silently dropped and `gquantiles` returns the default `percentiles = 0.5` (median) regardless of what you passed. The threshold is wrong; the downstream screen call returns the wrong intervals; nothing errors.

**A17. `gquantiles` on a track with many NAs assumed to include zeros.** `gquantiles` silently drops NA positions when computing percentiles. For "top 1% of covered positions" that's right; for "top 1% of all positions including zero-coverage", wrap with `ifelse(is.na(track), 0, track)` first.

## Aggregator choice

**A18. `func = "avg"` for count-like tracks.** ATAC reads, ChIP coverage, contact counts have meaningful magnitudes - `avg` divides by bin count and discards the count. For enrichment analyses you almost always want `sum` (1D) or `weighted.sum` (2D / sparse).

**A19. Aggregating a stored `pwm` (LSE) track with `func = "sum"`.** LSE does not compose under summation. Re-aggregating per-bin LSE values needs `func = "lse"` so the result is the LSE over the region's windows. `sum` over LSE values produces a number that has no probabilistic interpretation.

**A20. `bidirect = TRUE` when motif orientation matters.** Bidirectional PWM scoring returns the better-strand score per position, discarding the strand. For oriented-motif analyses (CTCF orientation, stranded peak calling), pair two vtracks with `strand = +1` and `strand = -1`, `bidirect = FALSE`, so the orientation survives.
