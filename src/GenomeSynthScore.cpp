/*
 * GenomeSynthScore.cpp
 *
 * C++ implementation for scoring a reference sequence under a trained
 * stratified Markov-k model. Writes a misha fixed-bin dense track whose
 * value at each output bin is the summed natural-log conditional
 * probability of the reference sequence under the model.
 */

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <vector>

#include "GenomeChromKey.h"
#include "GenomeSeqFetch.h"
#include "GenomeTrackFixedBin.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"
#include "StratifiedMarkovModel.h"

using namespace std;
using namespace rdb;

extern "C" {

/**
 * C_gsynth_score: Score reference sequence under a trained Markov model.
 *
 * @param _log_p_list   List of log-p matrices (one per bin); each is
 *                      num_kmers x 4 in column-major R layout.
 * @param _bin_indices  Integer vector of bin indices for each
 *                      iter-position (0-based; -1 = NA).
 * @param _iter_starts  Integer vector of iter-position start coords.
 * @param _iter_chroms  Integer vector of iter-position chrom IDs.
 * @param _intervals    R intervals object (regions to score).
 * @param _track_dir    String: misha track NAME (kernel calls
 *                      create_track_dir(envir, name) to mkdir).
 * @param _binsize      Integer: output bin size in bp (resolution).
 * @param _k            Integer: Markov order k.
 * @param _iter_size    Integer: stratum iterator bin size in bp.
 * @param _n_policy     Integer: 0 = NA, 1 = uniform.
 * @param _sparse_policy Integer: 0 = NA, 1 = uniform.
 * @param _envir        R environment.
 *
 * @return R_NilValue on success.
 */
SEXP C_gsynth_score(SEXP _log_p_list, SEXP _bin_indices,
                    SEXP _iter_starts, SEXP _iter_chroms,
                    SEXP _intervals, SEXP _track_dir,
                    SEXP _binsize, SEXP _k, SEXP _iter_size,
                    SEXP _n_policy, SEXP _sparse_policy,
                    SEXP _envir) {
    try {
        RdbInitializer rdb_init;
        IntervUtils iu(_envir);

        // ---- Parse scalars ----
        int k = Rf_asInteger(_k);
        if (k < 1 || k > StratifiedMarkovModel::MAX_K) {
            verror("k must be between 1 and %d",
                   StratifiedMarkovModel::MAX_K);
        }
        int num_kmers = 1;
        for (int i = 0; i < k; ++i) num_kmers *= NUM_BASES;

        unsigned binsize = (unsigned)Rf_asInteger(_binsize);
        if (binsize <= 0) verror("binsize must be a positive integer");

        int iter_size = Rf_asInteger(_iter_size);
        if (iter_size <= 0)
            verror("iter_size must be a positive integer; got %d",
                   iter_size);

        int n_policy = Rf_asInteger(_n_policy);            // 0=NA,1=unif
        int sparse_policy = Rf_asInteger(_sparse_policy);  // 0=NA,1=unif
        const float UNIFORM_LOGP = (float)std::log(0.25);

        // ---- Parse log_p list (num_bins entries, each num_kmers x 4) ----
        int num_bins = Rf_length(_log_p_list);
        if (num_bins <= 0) verror("log_p_list is empty");

        // Flat layout: log_p[bin][ctx * NUM_BASES + base].
        // Track per-bin "is sparse" flag: bin sparse if any cell is NaN.
        vector<vector<float>> log_p(num_bins);
        vector<bool> bin_is_sparse(num_bins, false);
        for (int b = 0; b < num_bins; ++b) {
            SEXP m = VECTOR_ELT(_log_p_list, b);
            if (!Rf_isReal(m)) verror("log_p_list[%d] is not numeric", b + 1);
            if (Rf_length(m) != num_kmers * NUM_BASES) {
                verror("log_p_list[%d] has wrong length (%d, expected %d)",
                       b + 1, Rf_length(m), num_kmers * NUM_BASES);
            }
            double *p = REAL(m);
            log_p[b].resize(num_kmers * NUM_BASES);
            for (int ctx = 0; ctx < num_kmers; ++ctx) {
                for (int base = 0; base < NUM_BASES; ++base) {
                    // R matrices are column-major: [row + col * nrow]
                    double v = p[ctx + base * num_kmers];
                    if (std::isnan(v)) bin_is_sparse[b] = true;
                    log_p[b][ctx * NUM_BASES + base] = (float)v;
                }
            }
        }

        // ---- Parse bin indices ----
        int num_iter_positions = Rf_length(_bin_indices);
        int *bin_indices = INTEGER(_bin_indices);
        int *iter_starts = INTEGER(_iter_starts);
        int *iter_chroms = INTEGER(_iter_chroms);

        // ---- Parse intervals ----
        const GenomeChromKey &chromkey = iu.get_chromkey();
        int num_chroms = chromkey.get_num_chroms();
        vector<vector<GInterval>> sample_per_chrom(num_chroms);
        if (!Rf_isNull(_intervals)) {
            GIntervalsFetcher1D *si = NULL;
            iu.convert_rintervs(_intervals, &si, NULL);
            unique_ptr<GIntervalsFetcher1D> guard(si);
            si->sort();
            for (si->begin_iter(); !si->isend(); si->next()) {
                const GInterval &iv = si->cur_interval();
                if (iv.chromid >= 0 && iv.chromid < num_chroms)
                    sample_per_chrom[iv.chromid].push_back(iv);
            }
        }

        // ---- Track name -> create dir (matches gtrack_create_dense) ----
        const char *track_name = CHAR(STRING_ELT(_track_dir, 0));
        string track_dir_s = create_track_dir(_envir, track_name);

        // ---- Per-chromosome 200bp bin lookup (sorted by start) ----
        vector<vector<pair<int64_t, int>>> chrom_bins(num_chroms);
        for (int i = 0; i < num_iter_positions; ++i) {
            int cid = iter_chroms[i];
            if (cid >= 0 && cid < num_chroms) {
                chrom_bins[cid].push_back({iter_starts[i], bin_indices[i]});
            }
        }
        for (int c = 0; c < num_chroms; ++c) {
            sort(chrom_bins[c].begin(), chrom_bins[c].end());
        }

        GenomeSeqFetch seqfetch;
        seqfetch.set_seqdir(string(rdb::get_groot(_envir)) + "/seq");

        // Total bp for progress reporting.
        uint64_t total_bp = 0;
        for (int c = 0; c < num_chroms; ++c) {
            for (const auto &iv : sample_per_chrom[c])
                total_bp += iv.end - iv.start;
        }
        Progress_reporter progress;
        progress.init(total_bp, 1000000);

        const float NaN_FLOAT = numeric_limits<float>::quiet_NaN();

        // Iterate ALL chromosomes in chromkey order (matching
        // GenomeTrackCreateDense): every chrom gets a track file, even
        // if no input intervals fall on it (writes all-NaN). Required
        // so gextract on out-of-input chroms doesn't fail.
        for (int chromid = 0; chromid < num_chroms; ++chromid) {
            int64_t chrom_size = chromkey.get_chrom_size(chromid);
            if (chrom_size <= 0) continue;

            int64_t num_out_bins =
                ((int64_t)chrom_size + binsize - 1) / binsize;

            // Per-bin running state. We can't simply sum log_p as we go
            // because if any bp is NA the bin must be NA (strict
            // policy preserves T(a) = log P(seq[a..a+r-1] | ...)).
            //   running_sum[i]  -- accumulated log_p sum for bin i
            //   any_na[i]       -- true if any bp in bin i is NA
            //   covered_bp[i]   -- bp count in bin i inside any input
            //                      interval (0 => bin is outside data)
            vector<double> running_sum(num_out_bins, 0.0);
            vector<bool> any_na(num_out_bins, false);
            vector<int64_t> covered_bp(num_out_bins, 0);

            const vector<GInterval> &ivs = sample_per_chrom[chromid];
            const vector<pair<int64_t, int>> &bins = chrom_bins[chromid];

            for (size_t iv_idx = 0; iv_idx < ivs.size(); ++iv_idx) {
                const GInterval &iv = ivs[iv_idx];
                int64_t start = max<int64_t>(0, iv.start);
                int64_t end = min<int64_t>(chrom_size, iv.end);
                if (end <= start) continue;

                // Read sequence with k bp upstream context where
                // available. seq[0] = read_start; seq[up_pad] = start.
                int64_t read_start = max<int64_t>(0, start - k);
                GInterval read_iv(chromid, read_start, end, 0);
                vector<char> seq;
                seqfetch.read_interval(read_iv, chromkey, seq);

                // Forward cursor over the 200bp bin lookup.
                size_t bin_cursor = 0;
                if (!bins.empty()) {
                    while (bin_cursor + 1 < bins.size() &&
                           start >= bins[bin_cursor + 1].first) {
                        ++bin_cursor;
                    }
                }

                for (int64_t pos = start; pos < end; ++pos) {
                    int64_t out_bin = pos / binsize;
                    if (out_bin < 0 || out_bin >= num_out_bins) continue;

                    covered_bp[out_bin]++;

                    if (any_na[out_bin]) {
                        // already poisoned — no need to score this bp
                        continue;
                    }

                    // k upstream bases available?
                    int64_t rel = pos - read_start;
                    if (rel < k) {
                        // Hit chromosome start with < k context.
                        // Boundary NA is unconditional (not gated by
                        // n_policy or sparse_policy).
                        any_na[out_bin] = true;
                        continue;
                    }

                    // Encode k-mer context (k bases ending at pos-1).
                    int ctx_idx = StratifiedMarkovModel::encode_kmer(
                        &seq[rel - k], k);
                    int base_idx = StratifiedMarkovModel::encode_base(
                        seq[rel]);

                    // Stratum bin lookup at this 200bp window.
                    int stratum_bin = -1;
                    if (!bins.empty()) {
                        while (bin_cursor + 1 < bins.size() &&
                               pos >= bins[bin_cursor + 1].first) {
                            ++bin_cursor;
                        }
                        if (pos >= bins[bin_cursor].first &&
                            pos < bins[bin_cursor].first + iter_size) {
                            stratum_bin = bins[bin_cursor].second;
                        }
                    }

                    bool n_invalid = (ctx_idx < 0 || base_idx < 0);
                    bool stratum_invalid =
                        (stratum_bin < 0 || stratum_bin >= num_bins);
                    bool sparse = (!stratum_invalid &&
                                   bin_is_sparse[stratum_bin]);

                    float contrib;
                    if (n_invalid) {
                        if (n_policy == 1) {
                            contrib = UNIFORM_LOGP;
                        } else {
                            any_na[out_bin] = true;
                            continue;
                        }
                    } else if (stratum_invalid) {
                        // Stratum NA is unconditional — there's no
                        // sensible default to score against.
                        any_na[out_bin] = true;
                        continue;
                    } else if (sparse) {
                        if (sparse_policy == 1) {
                            contrib = UNIFORM_LOGP;
                        } else {
                            any_na[out_bin] = true;
                            continue;
                        }
                    } else {
                        contrib = log_p[stratum_bin]
                                       [ctx_idx * NUM_BASES + base_idx];
                        // Defensive: cell could be NaN (p=0 even after
                        // pseudocount, very rare).
                        if (std::isnan(contrib)) {
                            any_na[out_bin] = true;
                            continue;
                        }
                    }
                    running_sum[out_bin] += contrib;
                }

                progress.report(end - start);
                check_interrupt();
            }

            // Stream per-bin output values through GenomeTrackFixedBin.
            char filename[FILENAME_MAX];
            snprintf(filename, sizeof(filename), "%s/%s",
                     track_dir_s.c_str(),
                     GenomeTrack::get_1d_filename(chromkey,
                                                   chromid).c_str());

            GenomeTrackFixedBin gtrack;
            gtrack.init_write(filename, binsize, chromid);

            vector<float> batch;
            batch.reserve(10000);
            for (int64_t i = 0; i < num_out_bins; ++i) {
                float v;
                if (covered_bp[i] == 0)
                    v = NaN_FLOAT;
                else if (any_na[i])
                    v = NaN_FLOAT;
                else
                    v = (float)running_sum[i];
                batch.push_back(v);
                if (batch.size() >= 10000) {
                    gtrack.write_next_bins(batch.data(), batch.size());
                    batch.clear();
                    check_interrupt();
                }
            }
            if (!batch.empty()) {
                gtrack.write_next_bins(batch.data(), batch.size());
            }
        }

        progress.report_last();
        return R_NilValue;

    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

}  // extern "C"
