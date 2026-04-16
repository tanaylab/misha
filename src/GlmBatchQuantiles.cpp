#include "GenomeTrack1D.h"  // for lse_accumulate
#include "rdbutils.h"
#include "GenomeTrack.h"
#include "GenomeTrackFixedBin.h"
#include "GenomeTrackSparse.h"

#include <Rinternals.h>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <limits>
#include <algorithm>
#include <thread>

using namespace rdb;
using namespace std;

// ---------------------------------------------------------------------------
// Track handle — same pattern as GlmFeatureExtractor
// ---------------------------------------------------------------------------
struct BqTrackHandle {
    string track_dir;
    GenomeTrack::Type type;
    shared_ptr<GenomeTrack> track;
    GenomeTrackFixedBin *fixedbin = nullptr;
    GenomeTrackSparse *sparse = nullptr;
    unsigned bin_size = 0;
};

static void bq_open_track(BqTrackHandle &handle, int chromid,
                           const GenomeChromKey &chromkey)
{
    string resolved = GenomeTrack::find_existing_1d_filename(
        chromkey, handle.track_dir, chromid);
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
        throw runtime_error(string("Track ") + handle.track_dir +
                            " has unsupported type (expected dense or sparse)");
    }
}

// ---------------------------------------------------------------------------
// LSE aggregation over [start, end) for sparse tracks
// ---------------------------------------------------------------------------
static double bq_aggregate_lse_sparse(const BqTrackHandle &handle,
                                       int64_t start, int64_t end)
{
    if (start < 0) start = 0;
    if (start >= end) return numeric_limits<double>::quiet_NaN();

    const GIntervals &intervals = handle.sparse->get_intervals();
    const vector<float> &vals = handle.sparse->get_vals();

    size_t lo = 0, hi = intervals.size();
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (intervals[mid].end <= start)
            lo = mid + 1;
        else
            hi = mid;
    }

    float lse = -numeric_limits<float>::infinity();
    uint64_t num_vs = 0;
    for (size_t i = lo; i < intervals.size(); i++) {
        if (intervals[i].start >= end) break;
        if (!isnan(vals[i])) {
            lse_accumulate(lse, vals[i]);
            num_vs++;
        }
    }
    if (num_vs == 0) return numeric_limits<double>::quiet_NaN();
    return (double)lse;
}

// ---------------------------------------------------------------------------
// Fast FixedBin scan: mmap entire chromosome, iterate with tight inner loop
// ---------------------------------------------------------------------------
static void bq_scan_fixedbin_fast(
    GenomeTrackFixedBin *fb,
    unsigned bin_size,
    uint64_t chrom_size,
    int iterator_step,
    int sshift,
    int eshift,
    vector<float> &values)
{
    int64_t total_bins = (int64_t)((chrom_size + bin_size - 1) / bin_size);
    if (total_bins <= 0) return;

    int64_t out_count = 0;
    const float *all_bins = fb->get_mmap_bins_ptr(0, total_bins, out_count);
    if (!all_bins || out_count <= 0) return;

    for (int64_t center = 0; center < (int64_t)chrom_size;
         center += iterator_step) {
        int64_t win_start = center + sshift;
        int64_t win_end = center + eshift;

        if (win_start < 0) win_start = 0;
        if (win_end <= win_start) continue;

        int64_t sbin = win_start / (int64_t)bin_size;
        int64_t ebin = (win_end + (int64_t)bin_size - 1) / (int64_t)bin_size;
        if (sbin < 0) sbin = 0;
        if (ebin > out_count) ebin = out_count;
        if (sbin >= ebin) continue;

        float lse = -numeric_limits<float>::infinity();
        uint64_t num_vs = 0;
        for (int64_t b = sbin; b < ebin; b++) {
            float v = all_bins[b];
            if (!isnan(v)) {
                lse_accumulate(lse, v);
                num_vs++;
            }
        }
        if (num_vs > 0) {
            values.push_back(lse);
        }
    }
}

// ---------------------------------------------------------------------------
// Per-track worker: scan genome, collect values, compute quantiles
// ---------------------------------------------------------------------------
struct TrackQuantileTask {
    // Input (read-only, shared across threads)
    string track_dir;
    GenomeTrack::Type type;
    int n_chroms;
    const GenomeChromKey *chromkey;
    int iterator_step;
    int sshift;
    int eshift;
    int n_pctiles;
    const double *percentiles;
    uint64_t estimated_positions;

    // Output (written by thread, read after join)
    vector<double> results;
    string error_msg;
};

static void worker_process_track(TrackQuantileTask &task) {
    try {
        BqTrackHandle handle;
        handle.track_dir = task.track_dir;
        handle.type = task.type;

        vector<float> values;
        values.reserve(task.estimated_positions);

        for (int chrom = 0; chrom < task.n_chroms; chrom++) {
            uint64_t chrom_size = task.chromkey->get_chrom_size(chrom);
            if (chrom_size == 0) continue;

            string resolved = GenomeTrack::find_existing_1d_filename(
                *task.chromkey, handle.track_dir, chrom);
            string filename = handle.track_dir + "/" + resolved;
            if (access(filename.c_str(), F_OK) != 0)
                continue;

            bq_open_track(handle, chrom, *task.chromkey);

            if (handle.fixedbin) {
                bq_scan_fixedbin_fast(
                    handle.fixedbin, handle.bin_size, chrom_size,
                    task.iterator_step, task.sshift, task.eshift, values);
            } else if (handle.sparse) {
                for (int64_t center = 0; center < (int64_t)chrom_size;
                     center += task.iterator_step) {
                    int64_t win_start = center + task.sshift;
                    int64_t win_end = center + task.eshift;
                    double val = bq_aggregate_lse_sparse(handle,
                                                          win_start, win_end);
                    if (!isnan(val) && isfinite(val)) {
                        values.push_back((float)val);
                    }
                }
            }

            handle.track.reset();
            handle.fixedbin = nullptr;
            handle.sparse = nullptr;
        }

        // Compute quantiles using nth_element (O(N) average, exact)
        int64_t N = (int64_t)values.size();
        task.results.resize(task.n_pctiles);

        for (int p = 0; p < task.n_pctiles; p++) {
            if (N == 0) {
                task.results[p] = numeric_limits<double>::quiet_NaN();
            } else {
                int64_t idx = (int64_t)floor(task.percentiles[p] * (double)(N - 1));
                if (idx < 0) idx = 0;
                if (idx >= N) idx = N - 1;

                nth_element(values.begin(), values.begin() + idx, values.end());
                task.results[p] = (double)values[idx];
            }
        }
    } catch (const exception &e) {
        task.error_msg = e.what();
    } catch (...) {
        task.error_msg = "Unknown error in worker thread";
    }
}

// ---------------------------------------------------------------------------
// .Call entry point: C_glm_batch_quantiles
//
// Computes genome-wide quantiles for multiple tracks using parallel threads.
// Each thread processes one track at a time, scanning the genome and computing
// LSE over the specified window at each iterator position.
//
// Memory: ~540 MB per concurrent thread (one track's values at 20bp).
// Default threads = min(n_tracks, hardware_concurrency, 40).
// ---------------------------------------------------------------------------
extern "C" SEXP C_glm_batch_quantiles(
    SEXP _track_names,   // character vector of track names
    SEXP _percentiles,   // numeric vector of percentiles (e.g., 0.9999)
    SEXP _iterator,      // integer scalar: iterator step size in bp
    SEXP _sshift,        // integer scalar: start shift for LSE window
    SEXP _eshift,        // integer scalar: end shift for LSE window
    SEXP _n_threads,     // integer scalar: number of parallel threads
    SEXP _envir)         // R environment
{
    try {
        RdbInitializer rdb_init;
        IntervUtils iu(_envir);

        // Parse inputs
        int n_tracks = Rf_length(_track_names);
        vector<string> track_names(n_tracks);
        for (int i = 0; i < n_tracks; i++) {
            track_names[i] = CHAR(STRING_ELT(_track_names, i));
        }

        int n_pctiles = Rf_length(_percentiles);
        const double *percentiles = REAL(_percentiles);

        int iterator_step = INTEGER(_iterator)[0];
        int sshift = INTEGER(_sshift)[0];
        int eshift = INTEGER(_eshift)[0];
        int n_threads_req = INTEGER(_n_threads)[0];

        if (iterator_step <= 0)
            verror("iterator must be a positive integer");

        // Get chromosome info
        const GenomeChromKey &chromkey = iu.get_chromkey();
        int n_chroms = (int)chromkey.get_num_chroms();

        // Resolve all track paths and types upfront (main thread only)
        SEXP envir = iu.get_env();
        vector<string> track_dirs(n_tracks);
        vector<GenomeTrack::Type> track_types(n_tracks);
        for (int m = 0; m < n_tracks; m++) {
            track_dirs[m] = track2path(envir, track_names[m]);
            track_types[m] = GenomeTrack::get_type(
                track_dirs[m].c_str(), chromkey, false);
        }

        // Pre-compute total genome size
        uint64_t total_genome_bp = 0;
        for (int chrom = 0; chrom < n_chroms; chrom++) {
            total_genome_bp += chromkey.get_chrom_size(chrom);
        }
        uint64_t estimated_positions = total_genome_bp / (uint64_t)iterator_step;

        // Determine number of threads
        int n_threads;
        if (n_threads_req > 0) {
            n_threads = min(n_threads_req, n_tracks);
        } else {
            // Auto-detect
            unsigned hw = thread::hardware_concurrency();
            if (hw == 0) hw = 4;
            n_threads = (int)min((unsigned)n_tracks, min(hw, 40u));
        }
        // Ensure at least 1 thread
        if (n_threads < 1) n_threads = 1;

        // Allocate output
        SEXP result;
        if (n_pctiles == 1) {
            rprotect(result = Rf_allocVector(REALSXP, n_tracks));
        } else {
            rprotect(result = Rf_allocMatrix(REALSXP, n_tracks, n_pctiles));
        }
        double *out = REAL(result);

        // Process tracks in parallel batches
        for (int batch_start = 0; batch_start < n_tracks;
             batch_start += n_threads) {
            check_interrupt();

            int batch_end = min(batch_start + n_threads, n_tracks);
            int batch_size = batch_end - batch_start;

            vector<TrackQuantileTask> tasks(batch_size);
            for (int i = 0; i < batch_size; i++) {
                int m = batch_start + i;
                tasks[i].track_dir = track_dirs[m];
                tasks[i].type = track_types[m];
                tasks[i].n_chroms = n_chroms;
                tasks[i].chromkey = &chromkey;
                tasks[i].iterator_step = iterator_step;
                tasks[i].sshift = sshift;
                tasks[i].eshift = eshift;
                tasks[i].n_pctiles = n_pctiles;
                tasks[i].percentiles = percentiles;
                tasks[i].estimated_positions = estimated_positions;
            }

            vector<thread> threads;
            threads.reserve(batch_size);
            for (int i = 0; i < batch_size; i++) {
                threads.emplace_back(worker_process_track, ref(tasks[i]));
            }

            for (auto &t : threads) {
                t.join();
            }

            for (int i = 0; i < batch_size; i++) {
                if (!tasks[i].error_msg.empty()) {
                    verror("Error processing track %s: %s",
                           track_names[batch_start + i].c_str(),
                           tasks[i].error_msg.c_str());
                }

                int m = batch_start + i;
                for (int p = 0; p < n_pctiles; p++) {
                    if (n_pctiles == 1) {
                        out[m] = tasks[i].results[p];
                    } else {
                        out[m + (int64_t)p * n_tracks] = tasks[i].results[p];
                    }
                }
            }
        }

        // Set names
        if (n_pctiles == 1) {
            SEXP names;
            rprotect(names = Rf_allocVector(STRSXP, n_tracks));
            for (int i = 0; i < n_tracks; i++) {
                SET_STRING_ELT(names, i, STRING_ELT(_track_names, i));
            }
            Rf_setAttrib(result, R_NamesSymbol, names);
        } else {
            SEXP rownames;
            rprotect(rownames = Rf_allocVector(STRSXP, n_tracks));
            for (int i = 0; i < n_tracks; i++) {
                SET_STRING_ELT(rownames, i, STRING_ELT(_track_names, i));
            }

            SEXP colnames;
            rprotect(colnames = Rf_allocVector(STRSXP, n_pctiles));
            for (int i = 0; i < n_pctiles; i++) {
                char buf[64];
                snprintf(buf, sizeof(buf), "%.6g", percentiles[i]);
                SET_STRING_ELT(colnames, i, Rf_mkChar(buf));
            }

            SEXP dimnames;
            rprotect(dimnames = Rf_allocVector(VECSXP, 2));
            SET_VECTOR_ELT(dimnames, 0, rownames);
            SET_VECTOR_ELT(dimnames, 1, colnames);
            Rf_setAttrib(result, R_DimNamesSymbol, dimnames);
        }

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
