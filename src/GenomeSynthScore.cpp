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
        const char *track_dir = track_dir_s.c_str();

        // Stub: nothing more to do until Tasks 4-5. Avoid unused-var warnings.
        (void)num_iter_positions;
        (void)bin_indices;
        (void)iter_starts;
        (void)iter_chroms;
        (void)track_dir;
        (void)UNIFORM_LOGP;
        (void)bin_is_sparse;
        (void)n_policy;
        (void)sparse_policy;
        (void)iter_size;

        return R_NilValue;

    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

}  // extern "C"
