# misha - advanced topics

Specialized recipes that aren't part of the everyday flow. Read [`misha-core.md`](misha-core.md) first for concepts and the common workhorse recipes.

## Contents

- [1. Meta-profile around peak centers](#1-meta-profile-around-peak-centers)
- [2. 2D contact pile-ups](#2-2d-contact-pile-ups)
- [3. Insulation, directionality index, domain borders](#3-insulation-directionality-index-domain-borders)
- [4. Sequence and PWM (motif) tracks](#4-sequence-and-pwm-motif-tracks)
- [5. Bulk import and export](#5-bulk-import-and-export)
- [6. Side topics (pointers)](#6-side-topics-pointers)

## 1. Meta-profile around peak centers

**Goal.** For a peak set, get a `(peak × offset)` signal matrix you can heatmap, row-cluster, or average into a smooth meta-profile.

```r
# 1. Anchor: per-peak center as a 1bp interval frame.
# With misha.ext available, this is one call: center_1bp <- misha.ext::gintervals.centers(peaks)
peaks$center <- (peaks$start + peaks$end) %/% 2L
center_1bp   <- gintervals(chrom = peaks$chrom,
                            start = peaks$center,
                            end   = peaks$center + 1L)

# 2. Distance vtrack -> signed offset at each iterator bin.
gvtrack.create("d_center", src = center_1bp, func = "distance.center")

# 3. Expand region symmetrically and extract one or more tracks at fine bins.
region <- transform(peaks,
                    start = center - W,
                    end   = center + W + 20L)
spat <- gextract(c("d_center", track_names),
                 intervals = region,
                 iterator  = 20,
                 colnames  = c("d_center", track_names))

# 4. Reshape long -> wide (peak x offset bin) for one track at a time.
mat <- spat |>
    mutate(bin = floor(d_center / 20)) |>
    filter(abs(bin) <= W / 20) |>
    pivot_wider(id_cols     = intervalID,
                names_from  = bin,
                values_from = !!sym(track_names[1]),
                values_fn   = mean,
                values_fill = 0) |>
    arrange(intervalID) |>
    select(-intervalID) |>
    as.matrix()
```

`distance.center` is the right flavor here - anchor → anchor.

**Variants.**
- **Multi-track in one pass.** Pass the whole vector to `gextract`; reshape once per track from the same long-format frame.
- **Pre-aggregated sum (no offset axis).** If you don't need the offset axis, skip `d_center` and build a window vtrack: `gvtrack.create("v", track, "sum", sshift = -W, eshift = W)`, then `gextract("v", intervals = peaks, iterator = peaks)`.
- **Argmax-recenter first.** When peak centers are coarse, refine via the `slice_max` recipe in misha-core §4.6 before building `center_1bp`.
- **Row-clustering.** `tglkmeans::TGL_kmeans_tidy(mat, k = 8)` is the lab-standard meta-profile clustering pass.

**Avoid:**
- `(peaks$start + peaks$end) / 2` - float division yields non-integer centers. Always `%/% 2L`.
- Forgetting to wrap the signal column with `ifelse(is.na(...), 0, ...)` inside the extract expression.
- Using `func = "distance"` instead of `func = "distance.center"` - the default `distance` mixes edge and normalized-in-interval semantics, which is *not* what an anchor-relative offset matrix needs.

## 2. 2D contact pile-ups

**Goal.** From a 2D contact track (Hi-C / capture-HiC), build a feature-pair contact density matrix - e.g. CTCF × CTCF, enhancer × promoter - stratified by genomic distance.

**Build a 2D pair iterator.** `giterator.cartesian_grid(intervals1, expansion1, intervals2, expansion2, min.band.idx, max.band.idx)` produces a 2D iterator over the cartesian product of two 1D intervals sets, each row expanded by `expansionN`:

```r
sites_x <- gintervals.load("intervs.ctcf_peaks")
sites_y <- gintervals.load("intervs.ctcf_peaks")

it <- giterator.cartesian_grid(
    intervals1   = sites_x, expansion1 = 1000,
    intervals2   = sites_y, expansion2 = 1000,
    min.band.idx = 1, max.band.idx = Inf)         # cis-only, off-diagonal
```

`min.band.idx` / `max.band.idx` restrict to a diagonal band - `1, Inf` is "off-diagonal cis only", `1, 1` is "near-diagonal", `0, 0` is "diagonal only".

**Aggregate contacts in the 2D window per pair.**

```r
gvtrack.create("obs", src = "hic.es.score", func = "area")

pair_mat <- gextract("obs", intervals = it, iterator = it,
                     band = c(-1e8, -1e3))   # exclude near-diagonal trivial contacts
```

`band = c(-1e8, -1e3)` keeps only contacts at least 1kb off-diagonal. The `band =` arg is a 2D-specific feature of `gextract` / `gscreen` / `gdist`.

**Per-axis distance binning** (asymmetric feature pile-up):

```r
gvtrack.create("dx", src = sites_x, func = "distance.center"); gvtrack.iterator("dx", dim = 1)
gvtrack.create("dy", src = sites_y, func = "distance.center"); gvtrack.iterator("dy", dim = 2)

h <- gdist("dx", dist_breaks,
           "dy", dist_breaks,
           "obs", obs_breaks,
           intervals = gintervals.2d.all(), iterator = it, include.lowest = TRUE)
# Collapse the obs-value axis to (dx, dy) mean obs, weighting each obs-bin by
# its midpoint (not its lower edge) and dividing by the total count per cell:
mids   <- (head(obs_breaks, -1) + tail(obs_breaks, -1)) / 2
mat    <- apply(h, c(1, 2), function(w) sum(w * mids) / max(sum(w), 1))
```

**Materialize a 2D track from a contacts text file.**

```r
gtrack.2d.import_contacts("hic.cell_a", "Cell A Hi-C contacts",
                          contacts = "/path/to/contacts.txt",
                          fends    = "/path/to/redb.fends",
                          allow.duplicates = FALSE)
```

For pre-scored rectangle outputs (shaman `.score` files), `contacts =` accepts a vector of file paths.

**`gtrack.2d.create` for sparse rectangles from R.**

```r
rects$value <- rects$obs
gtrack.2d.create("hic.derived", "Per-rectangle scores",
                 intervals = rects, values = rects$value)
```

**Avoid:**
- `giterator.cartesian.grid` (dot-form) - legacy name, no longer exported. Always underscore.
- 2D iterator without `min.band.idx` / `max.band.idx` when you only care about cis - the implicit "everything" scope is much slower.
- `gtrack.2d.create` with overlapping rectangles - 2D tracks expect non-overlapping cells.
- Forgetting `band = c(-X, -Y)` on cis extracts - the near-diagonal dominates the signal numerically.

## 3. Insulation, directionality index, domain borders

**Goal.** From a 2D contact track, build 1D per-bin scores capturing local TAD structure - *insulation* (low across a TAD border, high inside) and *directionality index* (asymmetry between upstream and downstream contact counts) - then call borders from local minima.

**Insulation via paired 2D vtracks.** Define a square 2D window on the diagonal:

```r
W <- 1e5
gvtrack.create("obs_ins",
               src    = "hic.es.score",
               func   = "weighted.sum",
               sshift1 = -W, eshift1 = W,
               sshift2 = -W, eshift2 = W)

ins <- gextract("obs_ins",
                intervals = gintervals.all(),
                iterator  = 2e4)            # 20kb diagonal iterator
```

The engine sweeps the diagonal when you pass an integer `iterator` against a 2D track.

**Directionality index** - same shape, two asymmetric windows:

```r
gvtrack.create("obs_up",
               src = "hic.es.score", func = "weighted.sum",
               sshift1 = -W, eshift1 = 0,
               sshift2 =  0, eshift2 = W)
gvtrack.create("obs_dn",
               src = "hic.es.score", func = "weighted.sum",
               sshift1 =  0, eshift1 = W,
               sshift2 = -W, eshift2 = 0)

di <- gextract("(obs_up - obs_dn) / (obs_up + obs_dn + 1)",
               intervals = gintervals.all(), iterator = 2e4)
```

**Persist as a track.**

```r
gtrack.create("hic.es.ins_1e5",
              description = "Insulation, 100kb window, 20kb diagonal",
              expr        = "obs_ins",
              iterator    = 2e4)
```

**Call borders from local-min of insulation.**

```r
gvtrack.create("ins_min",
               src    = "hic.es.ins_1e5",
               func   = "min",
               sshift = -W, eshift = W)

borders <- gscreen("hic.es.ins_1e5 == ins_min & hic.es.ins_1e5 < threshold",
                   intervals = gintervals.all(),
                   iterator  = 2e4)
```

Pick `threshold` with `gquantiles("hic.es.ins_1e5", percentiles = 0.1, ...)`.

**Multi-scale insulation.** Loop the same recipe over a vector of window sizes:

```r
for (W in c(5e4, 1e5, 2e5, 5e5)) {
    track_nm <- sprintf("hic.es.ins_%g", W)
    if (gtrack.exists(track_nm)) next
    gvtrack.create("obs_ins",
                   src = "hic.es.score", func = "weighted.sum",
                   sshift1 = -W, eshift1 = W,
                   sshift2 = -W, eshift2 = W)
    gtrack.create(track_nm, sprintf("Insulation, %gbp window", W),
                  expr = "obs_ins", iterator = 2e4)
}
```

**Avoid:**
- A single hardcoded `W = 1e5` window - TAD sizes vary; multi-scale is the lab norm.
- Picking `func = "sum"` over `weighted.sum` - for 2D contact tracks where rectangles vary in width, `weighted.sum` is correct.
- Computing DI without normalizing by `(up + dn + 1)` - raw difference scales with chromosome-arm coverage.
- Calling borders directly from `ins < threshold` without the local-min constraint - every bin in a deep valley qualifies, inflating border counts.

## 4. Sequence and PWM (motif) tracks

**Goal.** Extract genome sequence under intervals, score PSSMs on the fly via vtracks, and persist motif energies as tracks when reuse warrants it.

**Before writing string-level sequence code, check `gseq.*`.** Misha ships a sequence-manipulation family - `gseq.extract`, `gseq.rev`, `gseq.comp` (reverse / complement), `gseq.kmer` / `gseq.kmer.dist` (kmer counting / distance), `gseq.pwm` / `gseq.pwm_edits` (direct PWM scoring on R character strings), `gseq.read_homer` / `.read_jaspar` / `.read_meme` (motif file readers). Reach for these before hand-rolling reverse complement / kmer scans / PWM scoring in base R.

**Extract sequence.**

```r
seqs <- gseq.extract(intervals)
# Character vector, one DNA string per row. Reverse-complements when
# intervals$strand == -1; positive strand if `strand` is absent.
```

For 2D intervals frames, returns paired sequences (one per axis).

**Sequence-content vtracks.** Built-in `seq.GC500_bin20`, `seq.CG_500_mean`, etc. For custom windows, build `kmer.frac` / `kmer.count` vtracks. The `kmer` argument is a **single** k-mer string (e.g. `"G"`, `"CG"`, `"GATC"`); to count two or more, build one vtrack per k-mer and combine in the expression:

```r
gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G", sshift = -250, eshift = 250)
gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C", sshift = -250, eshift = 250)
gextract("g_frac + c_frac", intervals = peaks, iterator = peaks, colnames = "gc_w")

# Palindromic / strand-explicit k-mer: pass strand = 1 (or -1) to avoid double-counting.
gvtrack.create("cg_frac", NULL, "kmer.frac", kmer = "CG", strand = 1)
```

`src = NULL` selects the genome sequence as the source. Optional named args via `...`: `strand` (0 = both, default; 1 / -1 = stranded), `extend` (TRUE by default, scan past iterator bin edges so a k-mer straddling the boundary still gets counted).

**Materialize a k-mer track in one call.** When the multi-k-mer-sum-over-window pattern is something you want to persist as a track (not just an ephemeral vtrack), `misha.ext::gtrack.create_kmer` does it in one call - and accepts a *vector* of k-mers, so GC content is direct:

```r
# 200bp GC fraction at 20bp resolution, written as a dense track:
misha.ext::gtrack.create_kmer("seq.GC",
                              description = "GC fraction in 200bp window",
                              kmer        = c("G", "C"),
                              iterator    = 20,
                              window      = 200,
                              mode        = "frac",
                              strand      = 1)
```

For ephemeral analyses, the two-vtrack-plus-expression form above is fine; for tracks you'll reuse, prefer `gtrack.create_kmer`.

**Where motifs come from.** `prego` is an optional companion package (not in misha's `Imports` / `Suggests`); install it separately when you need curated PSSM libraries:

```r
prego::all_motif_datasets()       # union of HOMER + JASPAR + JOLMA + HOCOMOCO
prego::HOMER_motifs               # HOMER only
prego::JASPAR_motifs              # JASPAR only
ctcf_pssm <- prego::get_motif_pssm("CTCF", dataset = prego::JASPAR_motifs)
```

For de-novo motifs, train a PSSM from labelled sequences with `prego::regress_pwm(...)`.

**Dependency-free path.** When `prego` isn't installed, read PSSMs directly from motif files with misha's own `gseq.read_homer(file)`, `gseq.read_jaspar(file)`, or `gseq.read_meme(file)`. The resulting matrix plugs straight into `params = list(pssm = ...)`.

**PWM vtrack semantics.** A PWM vtrack scans the source sequence with sliding windows of width `nrow(pssm)`. Per-window score = log-likelihood under the PSSM relative to a uniform background (modulated by `prior`). Different funcs aggregate the per-window scores differently:

| `func` | What it returns per iterator bin | Aggregation across windows |
|---|---|---|
| `pwm`     | log-sum-exp of all window scores in the iterator interval | LSE (soft sum) |
| `pwm.max` | the maximum window score in the iterator interval | best-window |
| `pwm.max.pos` | the position (bp from interval start) of the best-scoring window | argmax |
| `pwm.count` | number of windows with score ≥ `score.thresh` | count above threshold |
| `pwm.edit_distance` | min #edits to raise (`direction = "above"`) or disrupt (`direction = "below"`) the best-window score across `score.thresh` | search over edit budget |
| `pwm.edit_distance.lse` | same but LSE of windows, not the max | LSE-based |
| `pwm.edit_distance.pos` / `pwm.edit_distance.lse.pos` | 1-based position of the most impactful single edit (signed by strand when `bidirect = TRUE`) | argmax over edit positions |

**Edit-distance knobs.** `direction = "above"` (default) finds the minimum edits to *raise* the score across `score.thresh` - "how close is this site to becoming a hit". `direction = "below"` finds the minimum edits to *disrupt* an existing hit - "how robust is this site". `max_edits` caps the total budget (exhaustive search for `max_edits <= 2`, greedy heuristic for `max_edits >= 3`). `max_indels` allows insertions / deletions in addition to substitutions. `score.min` / `score.max` clip the per-window score range before the edit search to keep the optimization bounded. With `bidirect = TRUE`, an edit's effect is evaluated against both strands and the better orientation is kept.

The iterator-bin width and the **PSSM width** interact. With `iterator = 20` and a 12bp PSSM, every 20bp bin contains 9 candidate windows (offsets 0..8 within the bin). With `iterator = 1`, each bin contains a single window.

**`extend`** (default `0`): extends the window scan past the iterator-bin edges by N bp, so a motif straddling a bin boundary still gets scored. Set to `nrow(pssm) - 1` to ensure no candidate window is missed at boundaries.

**`bidirect` and `strand`.**

```r
# bidirect = TRUE (default): both strands at every position, keep the better.
gvtrack.create("ctcf_max",
               func   = "pwm.max",
               params = list(pssm = ctcf_pssm, prior = 0.01,
                             bidirect = TRUE))

# Strand-specific: pair two vtracks, strand = +1 / -1, bidirect = FALSE.
gvtrack.create("ctcf_fwd", func = "pwm.max",
               params = list(pssm = ctcf_pssm, prior = 0.01,
                             strand = +1, bidirect = FALSE))
gvtrack.create("ctcf_rev", func = "pwm.max",
               params = list(pssm = ctcf_pssm, prior = 0.01,
                             strand = -1, bidirect = FALSE))
```

`prior` controls the strength of the uniform-background regularizer (smaller = sharper PWM; typical range 0.001-0.05). Tune by AUC on labelled data.

**Materialize a PWM energy track when reuse is warranted.** Preferred path: feed a `pwm` vtrack into `gtrack.create` - composes cleanly with any aggregator and lets you control `iterator`, `prior`, `bidirect`, etc. in one place:

```r
gvtrack.create("ctcf_lse",
               func   = "pwm",
               params = list(pssm = ctcf_pssm, prior = 0.01))

gtrack.create("motifs.ctcf",
              description = "CTCF JASPAR PWM energy (mm10)",
              expr        = "ctcf_lse",
              iterator    = 20)
```

Standard recipe: **materialize a dense PWM track at 20bp with the `pwm` (LSE) func**, then aggregate over arbitrary regions with `func = "lse"` vtracks on top - log-sum-exp composes correctly across bins. There's also a one-shot `gtrack.create_pwm_energy(track, description, pssmset, pssmid, prior, iterator)` for PSSMs stored under `<GROOT>/pssms/`, but the vtrack-then-`gtrack.create` route is more flexible (custom PSSM in `params`, control of `bidirect` / `strand` / `prior`, no `pssms/` directory required) and is the recommended form.

```r
# Re-aggregate a stored LSE track over arbitrary windows:
gvtrack.create("ctcf_region",
               src    = "motifs.ctcf",
               func   = "lse",
               sshift = -250, eshift = 250)
```

**Background calibration.** `gintervals.random(size, n, chromosomes, filter)` draws `n` fixed-width random intervals - feed to `gseq.extract` + `prego::regress_pwm` for de-novo training, or score with a PWM vtrack for empirical-quantile calibration.

**Avoid:**
- Aggregating a stored `pwm` (LSE) track with `func = "sum"` - LSE doesn't compose under summation. Use `func = "lse"`.
- `bidirect = TRUE` when downstream code cares about motif orientation - use paired `strand = +1` / `strand = -1` vtracks.
- Manually reverse-complementing strings from `gseq.extract` to score the opposite strand - the PWM engine handles strand internally.
- Picking `iterator = 1` "to be safe" - for most motifs `iterator = 10` or `20` with `extend = nrow(pssm) - 1` gives the same hits with 10-20× less storage.

## 5. Bulk import and export

**Goal.** Get external data into the misha track DB and round-trip tracks back out to standard formats for browser viewing, deposition, sharing.

**bigWig / WIG / TSV → dense fixed-bin track.**

```r
gtrack.import("epi.k27me3.es",
              description = "K27me3 ES bigWig (ENCODE ENCSR...)",
              file        = "k27me3_es.bw",
              binsize     = 20,
              defval      = NaN)
```

`gtrack.import` auto-dispatches on extension (`.bw`, `.wig`, `.bed`, `.txt`/`.tsv`).

**BED-style read intervals (R data.frame) → pileup track in one call.** `gtrack.create_dense(..., func = "coverage")`: bin value = `sum(value_i * overlap_i / binsize)` = average per-base signal in the bin. With `values = rep(1, n)`, this is a ChIP-seq-style pileup in a single call:

```r
gtrack.create_dense("atac.es",
                    description = "ATAC ES pileup",
                    intervals   = read_intervals,
                    values      = rep(1, nrow(read_intervals)),
                    binsize     = 20,
                    defval      = 0,
                    func        = "coverage")
```

Preferred path for any read-intervals-in-R source.

**Mapped reads from a SAM-style text file.**

```r
gtrack.import_mappedseq("atac.es",
                        description = "ATAC ES read coverage",
                        file        = "atac_es.txt.gz",
                        pileup      = 200,
                        binsize     = 20,
                        cols.order  = c(9, 11, 13, 14),
                        remove.dups = TRUE)
```

**Per-fragment 2D contacts → 2D track.**

```r
gtrack.2d.import_contacts("hic.cell_a",
                          description = "Cell A Hi-C",
                          contacts    = "contacts.txt",
                          fends       = "redb.fends",
                          allow.duplicates = FALSE)
```

**Export a track to bigWig / bedGraph.** For UCSC browser uploads, GEO depositions, sharing with collaborators on non-misha stacks:

```r
gtrack.export_bigwig("epi.k27me3.es",
                     file      = "k27me3_es.bw",
                     intervals = gintervals.all(),
                     iterator  = 20)

gtrack.export_bedgraph("epi.k27me3.es",
                       file      = "k27me3_es.bg",
                       intervals = peaks,                       # subset to a region set
                       iterator  = 20,
                       name      = "K27me3 ES (peaks only)")    # track header line
```

**Persist provenance at import.** Pass `attrs = list(source = ..., date = ...)` to `gtrack.import` (or `gtrack.attr.set` post-import) so the track survives someone trying to reproduce your analysis a year later.

**Copy a track between databases.** `gtrack.copy(src, dest = NULL, db = NULL, overwrite = FALSE)` copies a track from the *active* trackdb to `db` (another trackdb on disk), preserving indexed-format files, attributes, and metadata. `dest = NULL` keeps the same name in the destination; pass a vector of `src` names with `dest` as a namespace prefix to copy many at once. Useful for staging tracks from a personal trackdb to a shared one, or for promoting a derived track from a project-local trackdb (attached via `gdataset.load`) into the main `<GROOT>`. Distinct from `gsetroot` (which only switches which DB is *active*) and `gdataset.load` (which *attaches* a secondary DB read-only).

**Avoid:**
- `gtrack.import` without `binsize` for a continuous-signal source - resulting bin grid may not align with the rest of your trackdb.
- `gtrack.import_mappedseq` without `remove.dups = TRUE` for ChIP / ATAC - PCR duplicates inflate per-bin counts.
- Per-chromosome `mclapply` + `rbind` + `gtrack.create_sparse` to build a coverage track from a BAM - `gtrack.create_dense(..., func = "coverage")` does it in one call.

## 6. Side topics (pointers)

These have full conventions outside this guide.

**Methylation tracks (WGBS / PBAT / RRBS).** The lab convention is a four-suffix family per sample - `<sample>.cov`, `.meth`, `.unmeth`, `.avg`. For per-region methylation, the rule is `sum(meth) / sum(cov)` from the *count* tracks - **never** `mean(.avg)`, which is biased by low-coverage CpGs and undefined where `cov = 0`. The `.avg` track is for per-CpG views (heatmaps, scatter plots), not region-level aggregation.

With `misha.ext` installed, `gextract_meth` packages this convention in one call:

```r
meth <- misha.ext::gextract_meth(
    tracks    = c("wgbs.sample_a", "wgbs.sample_b"),
    intervals = peaks,
    iterator  = "intervs.global.seq_CG",   # default: per-CpG resolution
    d_expand  = NULL,                       # smooth ±d_expand bp before the ratio
    min_cov   = 5                           # mask samples with cov < 5 to NA
)
# Columns: chrom/start/end + one column per sample (the ratio) + ".cov" suffix per sample.
```

Without `misha.ext`, the hand-rolled equivalent is one `func = "sum"` vtrack per `.meth` and `.cov`, then `gextract(c("meth_v / cov_v", "cov_v"), ...)` and an explicit NA-fill on `.cov`.

**New genome bootstrap.** Three layers, from low-level to convenience:

- `gdb.create(groot, fasta, genes.file, annots.file, annots.names)` - write directory structure, index chromosomes from a FASTA, optionally seed annotation tracks from files you provide.
- `gdb.install_intervals(groot, source, sets, prefix, ...)` - install / refresh `intervs.global.*` annotation sets into an existing groot from an upstream source (UCSC, UCSC track hub, NCHI, local files, S3 backend). `sets =` is a whitelist (default `c("genes", "rmsk", "cgi", "cytoband")`); `prefix` lets you namespace alternative versions.
- `gdb.build_genome(name, path, registry, sets, ...)` - top-level convenience: builds an assembly end-to-end from a registry entry. Looks up the FASTA + annotation sources via `registry` (defaults to the lab `misha.yaml` registry) and runs `gdb.create` + `gdb.install_intervals` in one call.

Registry / discovery helpers: `gdb.list_genomes(registry)` shows what's available; `gdb.genome_info("hg38")` returns the registry entry (FASTA URL, annotation sources, alias rules). The default `misha.yaml` chain points at canonical UCSC / NCBI sources; override with a project-local `misha.yaml` for custom assemblies.

Naming note: in current misha the annotation-set whitelist arg is `sets =` (renamed from `annotations =`); the canonical CpG-island set is `cgi` (renamed from `cpgIsland`). Older scripts may still reference the old names.

For multi-contig / fragmented assemblies (Zoonomia primates, draft genomes), use chrom-alias matching to map source-track chromosome names onto the new genome's contigs: `target_chroms`, `match_by_length`, `min_coverage` arguments on `gdb.install_intervals` control the heuristic. Pair with `gdb.convert_to_indexed()` after build so the resulting trackdb scales.

Parallelize over many species with `parallel::mclapply` (set `gmultitasking = FALSE` inside the worker - see §2 in misha-core).

**Cross-assembly liftover.** `gintervals.load_chain(file)` + `gintervals.liftover(intervs, chain)`; always follow with `gintervals.force_range()` (lifted intervals can extend past chromosome ends) and `gintervals.canonic()` (lifted intervals may overlap).

**Reentrant genome swap.** When a script must temporarily switch genomes:

```r
cur_db <- .misha$GROOT
on.exit(gsetroot(cur_db), add = TRUE)
gsetroot("/path/to/other_db")
# ... analysis ...
```

Without the `on.exit` restore, a partially-completed run leaves the session pointing at the wrong DB and subsequent `gintervals.load` calls silently load against the wrong reference.

**Synthetic genomes via `gsynth`.** Train a k-th-order Markov sequence model on a real genome (optionally stratified by tracks such as GC content / CpG fraction), then sample synthetic sequence with matched per-base composition. Useful for ML background-distribution work.

- `gsynth.train(..., mask = NULL, intervals = NULL, iterator = NULL, k = 5, prior = "marginal")` - fit a k-th-order Markov model. The `...` slot takes zero or more stratification specs (each `list(expr = ..., breaks = ..., bin_merge = ...)` over an existing vtrack) so the model has per-stratum transition tables. With no `...`, fits a single global model. `mask =` is an intervals frame whose positions are skipped (e.g. repeats).
- `gsynth.save(model, file, compress = FALSE)` / `gsynth.load(file)` - canonical persistence (writes a `.gsm` directory, or a zip with `compress = TRUE`). `gsynth.load` also reads legacy `.rds` files.
- `gsynth.score(model, track, description, intervals, resolution = NULL)` - materialize the model's per-bin log-probability as a dense misha track (`track =` is the new track name; `resolution =` defaults to the model's training iterator).
- `gsynth.cell_merge(model, cell_merge, bin_merge = NULL)` - merge specified stratification cells *within* a single stratified model (not a multi-model combiner). `cell_merge` is a list of `list(from = c(...), to = c(...))` mappings.
- `gsynth.forbid_kmer(model, pattern)` - zero out transitions that would produce `pattern` (e.g. `"CG"` for a CpG-depleted background). `pattern` must be at most `k + 1` bases.
- `gsynth.sample(model, output_path, output_format = c("misha", "fasta", "vector"), intervals = NULL, n_samples = 1, preserve_n = TRUE, seed = NULL)` - emit synthetic sequence. Length is determined by `intervals` (defaults to the whole genome). `output_format = "fasta"` writes FASTA + `.fai`, `"misha"` writes a sequence track into the current trackdb, `"vector"` returns an R character vector.

**Synthetic perturbation genomes (`ggenome.*`).** Build a derived FASTA by editing the source genome at chosen intervals, then optionally bootstrap a new trackdb on top of it for downstream extraction.

- `ggenome.implant(intervals, donor, output, create_trackdb = TRUE, trackdb_path = ...)` - replace `intervals` in the source genome with the corresponding donor sequences (one per row). Used for CRE-destroy / motif-shuffle ablations.
- `ggenome.transplant(intervals, source_genome, target_genome, output, ...)` - splice intervals from one assembly into another at matched positions (cross-species transplants).

Typical workflow: load the source genome, write the derived FASTA via `ggenome.implant`, point `gdb.create` or `gdb.build_genome` at it to make a new trackdb, then run your usual extract / vtrack analysis against the perturbed genome.
