#include "GlmFeatureExtractor.h"
#include "GenomeTrack1D.h"  // for lse_accumulate
#include "rdbutils.h"
#include "GenomeTrack.h"
#include "GenomeTrackFixedBin.h"

#include <Rinternals.h>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <atomic>
#include <thread>
#include <stdexcept>
#include <mutex>

using namespace rdb;
using namespace std;

// ---------------------------------------------------------------------------
// Track opening — handles both FixedBin and Sparse tracks
// ---------------------------------------------------------------------------
void GlmFeatureExtractor::open_track_static(TrackHandle &handle, int chromid,
                                             const GenomeChromKey &chromkey)
{
    string resolved = GenomeTrack::find_existing_1d_filename(chromkey, handle.track_dir, chromid);
    string filename = handle.track_dir + "/" + resolved;

    handle.fixedbin = nullptr;
    handle.sparse = nullptr;
    handle.bin_size = 0;

    if (handle.type == GenomeTrack::FIXED_BIN) {
        auto t = make_shared<GenomeTrackFixedBin>();
        t->init_read(filename.c_str(), chromid);
        handle.track = t;
        handle.fixedbin = t.get();
        handle.bin_size = t->get_bin_size();
    } else if (handle.type == GenomeTrack::SPARSE) {
        auto t = make_shared<GenomeTrackSparse>();
        t->init_read(filename.c_str(), chromid);
        handle.track = t;
        handle.sparse = t.get();
    } else {
        // Worker context: throw a plain runtime_error so we don't longjmp
        // out of a non-main thread via verror().
        throw std::runtime_error(
            std::string("Track ") + handle.track_dir +
            " has unsupported type (expected dense or sparse)");
    }
}

void GlmFeatureExtractor::open_track(TrackHandle &handle, int chromid)
{
    open_track_static(handle, chromid, m_iu.get_chromkey());
}

// ---------------------------------------------------------------------------
// Aggregation: log-sum-exp over a window
// Uses the same float-precision lse_accumulate as misha's vtrack func="lse"
// (GenomeTrack1D.h) to produce bit-identical results.
// ---------------------------------------------------------------------------
double GlmFeatureExtractor::aggregate_lse(const TrackHandle &handle,
                                           int64_t start, int64_t end)
{
    if (start < 0) start = 0;
    if (start >= end) return numeric_limits<double>::quiet_NaN();

    if (handle.fixedbin) {
        unsigned bin_size = handle.bin_size;
        if (bin_size == 0) return numeric_limits<double>::quiet_NaN();

        int64_t sbin = start / (int64_t)bin_size;
        int64_t ebin = (int64_t)ceil(end / (double)bin_size);
        int64_t num_bins = ebin - sbin;
        if (num_bins <= 0) return numeric_limits<double>::quiet_NaN();

        int64_t out_count = 0;
        const float *ptr = handle.fixedbin->get_mmap_bins_ptr(sbin, num_bins, out_count);
        if (!ptr || out_count <= 0) return numeric_limits<double>::quiet_NaN();

        // Match misha's vtrack LSE: sequential float-precision accumulation
        // (GenomeTrackFixedBin.cpp line 497-506, uses lse_accumulate(float&, float))
        float lse = -numeric_limits<float>::infinity();
        uint64_t num_vs = 0;
        for (int64_t i = 0; i < out_count; i++) {
            if (!isnan(ptr[i])) {
                lse_accumulate(lse, ptr[i]);
                num_vs++;
            }
        }
        if (num_vs == 0) return numeric_limits<double>::quiet_NaN();
        return (double)lse;
    }

    if (handle.sparse) {
        const GIntervals &intervals = handle.sparse->get_intervals();
        const vector<float> &vals = handle.sparse->get_vals();

        size_t idx = sparse_lower_bound(intervals, start);

        float lse = -numeric_limits<float>::infinity();
        uint64_t num_vs = 0;
        for (size_t i = idx; i < intervals.size(); i++) {
            if (intervals[i].start >= end) break;
            if (!isnan(vals[i])) {
                lse_accumulate(lse, vals[i]);
                num_vs++;
            }
        }
        if (num_vs == 0) return numeric_limits<double>::quiet_NaN();
        return (double)lse;
    }

    return numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------------------------
// Aggregation: sum over a window (supports both FixedBin and Sparse)
// ---------------------------------------------------------------------------
double GlmFeatureExtractor::aggregate_sum(const TrackHandle &handle,
                                           int64_t start, int64_t end)
{
    if (start < 0) start = 0;
    if (start >= end) return numeric_limits<double>::quiet_NaN();

    if (handle.fixedbin) {
        unsigned bin_size = handle.bin_size;
        if (bin_size == 0) return numeric_limits<double>::quiet_NaN();

        int64_t sbin = start / (int64_t)bin_size;
        int64_t ebin = (int64_t)ceil(end / (double)bin_size);
        int64_t num_bins = ebin - sbin;
        if (num_bins <= 0) return numeric_limits<double>::quiet_NaN();

        int64_t out_count = 0;
        const float *ptr = handle.fixedbin->get_mmap_bins_ptr(sbin, num_bins, out_count);
        if (!ptr || out_count <= 0) return numeric_limits<double>::quiet_NaN();

        double acc = 0.0;
        bool has_val = false;
        for (int64_t i = 0; i < out_count; i++) {
            float v = ptr[i];
            if (!isnan(v) && !isinf(v)) {
                acc += (double)v;
                has_val = true;
            }
        }
        return has_val ? acc : numeric_limits<double>::quiet_NaN();
    }

    if (handle.sparse) {
        const GIntervals &intervals = handle.sparse->get_intervals();
        const vector<float> &vals = handle.sparse->get_vals();

        size_t idx = sparse_lower_bound(intervals, start);

        double acc = 0.0;
        bool has_val = false;
        for (size_t i = idx; i < intervals.size(); i++) {
            if (intervals[i].start >= end) break;
            float v = vals[i];
            if (!isnan(v) && !isinf(v)) {
                int64_t overlap_start = max((int64_t)intervals[i].start, start);
                int64_t overlap_end = min((int64_t)intervals[i].end, end);
                int64_t overlap_len = overlap_end - overlap_start;
                acc += (double)v * overlap_len;
                has_val = true;
            }
        }
        return has_val ? acc : numeric_limits<double>::quiet_NaN();
    }

    return numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------------------------
// Binary search: find first interval where intervals[i].end > pos
// (i.e., the first interval that could overlap a query starting at pos)
// ---------------------------------------------------------------------------
size_t GlmFeatureExtractor::sparse_lower_bound(const GIntervals &intervals, int64_t pos)
{
    size_t lo = 0, hi = intervals.size();
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (intervals[mid].end <= pos)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

// ---------------------------------------------------------------------------
// Main extraction
// ---------------------------------------------------------------------------
void GlmFeatureExtractor::extract(
    const vector<string> &track_names,
    const int *peak_chromids,
    const int64_t *peak_starts,
    const int64_t *peak_ends,
    int n_peaks,
    int tile_size,
    int flank_size,
    const vector<LmScalingConfig> &scaling,
    const vector<LmTransformConfig> &transforms,
    const string &gc_track_name,
    double gc_scale_factor,
    double *output,
    int n_cols,
    int n_threads)
{
    int n_motifs = (int)track_names.size();
    int n_transforms = (int)transforms.size();

    // Tile geometry: peaks of uniform width get n_tiles symmetric tiles of
    // tile_size, with flank_size of outer flank on each side.
    int first_peak_size = (int)(peak_ends[0] - peak_starts[0]);
    int extended = first_peak_size + 2 * flank_size;
    int n_tiles = extended / tile_size;

    int n_gc_inter = n_tiles * (n_tiles - 1) / 2;
    int expected_cols = n_motifs * n_tiles * n_transforms + n_tiles + n_gc_inter;
    if (n_cols != expected_cols) {
        verror("Column count mismatch: expected %d, got %d", expected_cols, n_cols);
    }

    // Resolve all track paths up-front on the main thread (R-env access is
    // not safe to call from worker threads — track2path / get_type touch the
    // R interpreter via the misha environment).
    SEXP envir = m_iu.get_env();
    const GenomeChromKey &chromkey = m_iu.get_chromkey();

    vector<string> motif_dirs(n_motifs);
    vector<GenomeTrack::Type> motif_types(n_motifs);
    for (int m = 0; m < n_motifs; m++) {
        motif_dirs[m] = track2path(envir, track_names[m]);
        motif_types[m] = GenomeTrack::get_type(motif_dirs[m].c_str(), chromkey, false);
    }
    string gc_dir = track2path(envir, gc_track_name);
    GenomeTrack::Type gc_type =
        GenomeTrack::get_type(gc_dir.c_str(), chromkey, false);

    // Sort peaks by (chrom, start). Each contiguous chromosome run lets a
    // worker open all motif/GC tracks once and process many peaks before
    // re-opening for the next chromosome. The atomic counter below also
    // hands out chunks of consecutive sorted indices, so within a chunk a
    // worker typically stays on a single chromosome.
    vector<int> peak_order(n_peaks);
    iota(peak_order.begin(), peak_order.end(), 0);
    sort(peak_order.begin(), peak_order.end(), [&](int a, int b) {
        if (peak_chromids[a] != peak_chromids[b])
            return peak_chromids[a] < peak_chromids[b];
        return peak_starts[a] < peak_starts[b];
    });

    memset(output, 0, sizeof(double) * (int64_t)n_peaks * n_cols);

    int motif_block_size = n_motifs * n_tiles * n_transforms;
    int gc_block_start = motif_block_size;
    int gc_inter_start = gc_block_start + n_tiles;

    if (n_threads < 1) n_threads = 1;
    if (n_threads > n_peaks) n_threads = n_peaks;

    // Group peak_order indices by chromosome run. peak_order is already
    // sorted by (chromid, start), so each chromosome is a contiguous slice.
    struct ChromRun { int chromid; int oi_lo; int oi_hi; };
    std::vector<ChromRun> chrom_runs;
    {
        int run_lo = 0;
        while (run_lo < n_peaks) {
            int chromid = peak_chromids[peak_order[run_lo]];
            int run_hi = run_lo + 1;
            while (run_hi < n_peaks &&
                   peak_chromids[peak_order[run_hi]] == chromid)
                ++run_hi;
            chrom_runs.push_back({chromid, run_lo, run_hi});
            run_lo = run_hi;
        }
    }

    std::mutex err_mtx;
    std::string first_error;

    // Per-chromosome processing: open all motif + GC handles ONCE (single
    // thread), pre-materialize sparse tracks, then fan out the peak loop
    // across worker threads that share the read-only handles. This avoids
    // the previous per-thread re-open cost (191 mmap setups × N_threads ×
    // N_chroms) which made >4-thread runs scale negatively.
    //
    // Thread safety:
    //   - GenomeTrackFixedBin::get_mmap_bins_ptr is `const` and only reads
    //     the mmap region — safe for concurrent calls.
    //   - GenomeTrackSparse::get_intervals/get_vals are lazy on first call;
    //     we trigger them from the main thread before fan-out so the
    //     workers see fully-loaded vectors.
    //   - Output rows are disjoint per peak (no write contention).
    //
    // Per-chromosome handle release: shared_ptr<GenomeTrack> goes out of
    // scope at the end of each iteration, freeing the previous chrom's
    // mmap mappings before the next chrom is opened.
    constexpr int CHUNK = 256;

    for (const auto &run : chrom_runs) {
        int chromid = run.chromid;

        std::vector<TrackHandle> motif_handles(n_motifs);
        TrackHandle gc_handle;
        try {
            for (int m = 0; m < n_motifs; m++) {
                motif_handles[m].track_dir = motif_dirs[m];
                motif_handles[m].type = motif_types[m];
                open_track_static(motif_handles[m], chromid, chromkey);
                if (motif_handles[m].sparse) {
                    motif_handles[m].sparse->get_intervals();
                    motif_handles[m].sparse->get_vals();
                }
            }
            gc_handle.track_dir = gc_dir;
            gc_handle.type = gc_type;
            open_track_static(gc_handle, chromid, chromkey);
            if (gc_handle.sparse) {
                gc_handle.sparse->get_intervals();
                gc_handle.sparse->get_vals();
            }
        } catch (TGLException &e) {
            verror("%s", e.msg());
        } catch (const std::exception &e) {
            verror("%s", e.what());
        }

        std::atomic<int> next_chunk{0};
        const int n_run_peaks = run.oi_hi - run.oi_lo;

        auto worker = [&]() {
            try {
                std::vector<double> gc_vals(n_tiles, 0.0);
                while (true) {
                    int chunk_idx = next_chunk.fetch_add(1, std::memory_order_relaxed);
                    int rel_lo = chunk_idx * CHUNK;
                    if (rel_lo >= n_run_peaks) break;
                    int rel_hi = std::min(rel_lo + CHUNK, n_run_peaks);

                    for (int rel = rel_lo; rel < rel_hi; rel++) {
                        int pi = peak_order[run.oi_lo + rel];

                        int64_t peak_center = (peak_starts[pi] + peak_ends[pi]) / 2;
                        int peak_size = (int)(peak_ends[pi] - peak_starts[pi]);
                        int half = peak_size / 2;

                        for (int ti = 0; ti < n_tiles; ti++) {
                            int64_t tile_start = peak_center - half - flank_size +
                                                 (int64_t)ti * tile_size;
                            int64_t tile_end = tile_start + tile_size;

                            double gc_raw = aggregate_sum(gc_handle, tile_start, tile_end);
                            double gc_scaled = isnan(gc_raw) ? 0.0
                                : (gc_raw / tile_size) * gc_scale_factor;
                            gc_vals[ti] = gc_scaled;

                            int gc_col = gc_block_start + ti;
                            output[pi + (int64_t)gc_col * n_peaks] = gc_scaled;

                            for (int m = 0; m < n_motifs; m++) {
                                double raw = aggregate_lse(motif_handles[m], tile_start, tile_end);
                                double scaled = (isnan(raw) || !isfinite(raw))
                                    ? 0.0 : apply_scaling(raw, scaling[m]);

                                for (int t = 0; t < n_transforms; t++) {
                                    double transformed = apply_transform(scaled, transforms[t]);
                                    if (!isfinite(transformed)) transformed = 0.0;
                                    int col = (ti * n_motifs + m) * n_transforms + t;
                                    output[pi + (int64_t)col * n_peaks] = transformed;
                                }
                            }
                        }

                        int inter_idx = 0;
                        for (int a = 0; a < n_tiles; a++) {
                            for (int b = a + 1; b < n_tiles; b++) {
                                int col = gc_inter_start + inter_idx;
                                output[pi + (int64_t)col * n_peaks] =
                                    gc_vals[a] * gc_vals[b] / gc_scale_factor;
                                inter_idx++;
                            }
                        }
                    }
                }
            } catch (TGLException &e) {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (first_error.empty()) first_error = e.msg();
            } catch (const std::exception &e) {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (first_error.empty()) first_error = e.what();
            } catch (...) {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (first_error.empty()) first_error = "unknown error in worker thread";
            }
        };

        int run_threads = std::min(n_threads, n_run_peaks);
        if (run_threads <= 1) {
            worker();
        } else {
            std::vector<std::thread> threads;
            threads.reserve(run_threads);
            for (int i = 0; i < run_threads; i++) threads.emplace_back(worker);
            for (auto &t : threads) t.join();
        }
        if (!first_error.empty()) break;
    }

    if (!first_error.empty()) verror("%s", first_error.c_str());
}

// ---------------------------------------------------------------------------
// .Call entry point
// ---------------------------------------------------------------------------
extern "C" SEXP C_glm_extract_features(
    SEXP _track_names,     // character vector
    SEXP _chroms,          // character vector (chromosome names)
    SEXP _starts,          // numeric vector (int64 as double)
    SEXP _ends,            // numeric vector (int64 as double)
    SEXP _tile_size,       // integer scalar
    SEXP _flank_size,      // integer scalar
    SEXP _max_caps,        // numeric vector (per-motif, same order as track_names)
    SEXP _dis_from_cap,    // numeric scalar
    SEXP _scale_factor,    // numeric scalar
    SEXP _transforms,      // numeric matrix: n_transforms x 5 (L, k, x_0, pre_shift, post_shift)
    SEXP _gc_track,        // character scalar
    SEXP _gc_scale_factor, // numeric scalar
    SEXP _n_threads,       // integer scalar (0 = auto)
    SEXP _envir)           // R environment
{
    try {
        RdbInitializer rdb_init;
        IntervUtils iu(_envir);

        // Parse track names
        int n_motifs = Rf_length(_track_names);
        vector<string> track_names(n_motifs);
        for (int i = 0; i < n_motifs; i++) {
            track_names[i] = CHAR(STRING_ELT(_track_names, i));
        }

        // Parse intervals: resolve chromosome names to chromkey IDs here so
        // callers don't have to know misha's internal chromid layout (which
        // is the chromkey insertion order, not any subset ordering).
        if (TYPEOF(_chroms) != STRSXP) {
            verror("'chroms' must be a character vector of chromosome names");
        }
        int n_peaks = Rf_length(_chroms);
        const GenomeChromKey &chromkey_for_ids = iu.get_chromkey();
        vector<int> chromids(n_peaks);
        for (int i = 0; i < n_peaks; i++) {
            SEXP s = STRING_ELT(_chroms, i);
            if (s == NA_STRING) {
                verror("chroms[%d] is NA", i + 1);
            }
            chromids[i] = chromkey_for_ids.chrom2id(CHAR(s));
        }
        const double *starts_d = REAL(_starts);
        const double *ends_d = REAL(_ends);

        // Convert doubles to int64
        vector<int64_t> starts(n_peaks), ends(n_peaks);
        for (int i = 0; i < n_peaks; i++) {
            starts[i] = (int64_t)starts_d[i];
            ends[i] = (int64_t)ends_d[i];
        }

        int tile_size = INTEGER(_tile_size)[0];
        int flank_size = INTEGER(_flank_size)[0];

        // Parse scaling config
        const double *max_caps = REAL(_max_caps);
        double dis_from_cap = REAL(_dis_from_cap)[0];
        double scale_factor = REAL(_scale_factor)[0];

        vector<LmScalingConfig> scaling(n_motifs);
        for (int i = 0; i < n_motifs; i++) {
            scaling[i].max_cap = max_caps[i];
            scaling[i].dis_from_cap = dis_from_cap;
            scaling[i].scale_factor = scale_factor;
        }

        // Parse transforms (matrix: n_transforms x 5, column-major)
        int n_transforms = Rf_nrows(_transforms);
        const double *tdata = REAL(_transforms);
        vector<LmTransformConfig> transforms(n_transforms);
        for (int i = 0; i < n_transforms; i++) {
            transforms[i].L          = tdata[i + 0 * n_transforms];
            transforms[i].k          = tdata[i + 1 * n_transforms];
            transforms[i].x_0        = tdata[i + 2 * n_transforms];
            transforms[i].pre_shift  = tdata[i + 3 * n_transforms];
            transforms[i].post_shift = tdata[i + 4 * n_transforms];
        }

        string gc_track = CHAR(STRING_ELT(_gc_track, 0));
        double gc_scale_factor = REAL(_gc_scale_factor)[0];

        // Thread count. 0 = auto: min(n_peaks, hardware_concurrency, 40).
        // Negative or NA values fall back to 1 (silent, matching glm_batch_quantiles).
        int n_threads_req = INTEGER(_n_threads)[0];
        int n_threads;
        if (n_threads_req == NA_INTEGER || n_threads_req < 0) {
            n_threads = 1;
        } else if (n_threads_req == 0) {
            unsigned hw = std::thread::hardware_concurrency();
            if (hw == 0) hw = 4;
            n_threads = (int)std::min((unsigned)n_peaks, std::min(hw, 40u));
            if (n_threads < 1) n_threads = 1;
        } else {
            n_threads = std::min(n_threads_req, n_peaks);
            if (n_threads < 1) n_threads = 1;
        }

        // Compute output dimensions
        int first_peak_size = (int)(ends[0] - starts[0]);
        int extended = first_peak_size + 2 * flank_size;
        int n_tiles = extended / tile_size;
        int n_gc_inter = n_tiles * (n_tiles - 1) / 2;
        int n_cols = n_motifs * n_tiles * n_transforms + n_tiles + n_gc_inter;

        // Allocate output matrix
        SEXP result;
        rprotect(result = Rf_allocMatrix(REALSXP, n_peaks, n_cols));
        double *output = REAL(result);

        // Run extraction
        GlmFeatureExtractor extractor(iu);
        extractor.extract(
            track_names,
            chromids.data(),
            starts.data(),
            ends.data(),
            n_peaks,
            tile_size,
            flank_size,
            scaling,
            transforms,
            gc_track,
            gc_scale_factor,
            output,
            n_cols,
            n_threads
        );

        return result;

    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    } catch (const exception &e) {
        rerror("%s", e.what());
    }
    return R_NilValue;
}
