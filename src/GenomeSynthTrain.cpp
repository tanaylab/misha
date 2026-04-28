/*
 * GenomeSynthTrain.cpp
 *
 * C++ implementation for training a stratified Markov-k model
 * from genomic sequences.
 */

#include <algorithm>
#include <climits>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "Filter.h"
#include "GenomeChromKey.h"
#include "GenomeSeqFetch.h"
#include "GIntervalsBigSet1D.h"
#include "MaskUtils.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"
#include "StratifiedMarkovModel.h"

using namespace std;
using namespace rdb;

extern "C" {

/**
 * C_gsynth_train: Train a stratified Markov-k model from genome sequences.
 *
 * @param _chrom_ids Integer vector of chromosome IDs to process
 * @param _chrom_starts Integer vector of start positions for each chromosome
 * @param _chrom_ends Integer vector of end positions for each chromosome
 * @param _bin_indices Integer vector of bin indices for each position
 *                     (aligned with the iterator positions)
 * @param _iter_starts Integer vector of iterator interval start positions
 * @param _iter_chroms Integer vector of iterator interval chromosome IDs
 * @param _breaks Numeric vector of bin boundaries
 * @param _bin_map Integer vector mapping source bins to target bins (-1 = keep)
 * @param _mask R intervals object for mask regions (NULL if no mask)
 * @param _pseudocount Numeric pseudocount for normalization
 * @param _k Markov order
 * @param _prior_mode Character: "uniform", "marginal", "global", or "explicit"
 * @param _prior_matrix Numeric matrix (n_bins x 4) for "explicit" mode, else NULL
 * @param _envir R environment
 *
 * @return A list containing the trained model data
 */
SEXP C_gsynth_train(SEXP _chrom_ids, SEXP _chrom_starts, SEXP _chrom_ends,
                     SEXP _bin_indices, SEXP _iter_starts, SEXP _iter_chroms,
                     SEXP _breaks, SEXP _bin_map, SEXP _mask,
                     SEXP _pseudocount, SEXP _k,
                     SEXP _prior_mode, SEXP _prior_matrix, SEXP _envir) {
    try {
        RdbInitializer rdb_init;
        IntervUtils iu(_envir);

        // Parse Markov order k (default 5 for backward compatibility)
        int k = Rf_isNull(_k) ? 5 : INTEGER(_k)[0];
        if (k < 1 || k > StratifiedMarkovModel::MAX_K) {
            verror("k must be between 1 and %d", StratifiedMarkovModel::MAX_K);
        }
        int kmer_len = k + 1;  // sliding window length: k context + 1 next base

        // Extract chromosome ranges
        int num_chroms = Rf_length(_chrom_ids);
        if (num_chroms == 0) {
            verror("No chromosomes to process");
        }

        int* chrom_ids = INTEGER(_chrom_ids);
        int* chrom_starts = INTEGER(_chrom_starts);
        int* chrom_ends = INTEGER(_chrom_ends);

        // Extract bin indices from track extraction
        int num_iter_positions = Rf_length(_bin_indices);
        int* bin_indices = INTEGER(_bin_indices);
        int* iter_starts = INTEGER(_iter_starts);
        int* iter_chroms = INTEGER(_iter_chroms);

        // Extract breaks
        int num_breaks = Rf_length(_breaks);
        int num_bins = num_breaks - 1;
        if (num_bins <= 0) {
            verror("breaks must have at least 2 elements");
        }

        double* breaks = REAL(_breaks);
        vector<double> breaks_vec(breaks, breaks + num_breaks);

        // Extract bin_map (or create identity mapping)
        vector<int> bin_map_vec(num_bins);
        if (!Rf_isNull(_bin_map) && Rf_length(_bin_map) == num_bins) {
            int* bin_map = INTEGER(_bin_map);
            for (int i = 0; i < num_bins; ++i) {
                bin_map_vec[i] = bin_map[i];
            }
        } else {
            // Identity mapping
            for (int i = 0; i < num_bins; ++i) {
                bin_map_vec[i] = i;
            }
        }

        double pseudocount = Rf_asReal(_pseudocount);

        // Parse prior mode (default uniform if NULL/missing)
        std::string prior_mode = "uniform";
        if (!Rf_isNull(_prior_mode) && Rf_length(_prior_mode) > 0) {
            prior_mode = CHAR(STRING_ELT(_prior_mode, 0));
        }

        // Parse mask intervals if provided
        vector<vector<GInterval>> mask_per_chrom;
        mask_per_chrom.resize(iu.get_chromkey().get_num_chroms());

        if (!Rf_isNull(_mask)) {
            GIntervalsFetcher1D* mask_intervals = NULL;
            iu.convert_rintervs(_mask, &mask_intervals, NULL);
            unique_ptr<GIntervalsFetcher1D> mask_guard(mask_intervals);
            mask_intervals->sort();

            for (mask_intervals->begin_iter(); !mask_intervals->isend();
                 mask_intervals->next()) {
                const GInterval& iv = mask_intervals->cur_interval();
                if (iv.chromid >= 0 &&
                    iv.chromid < (int)mask_per_chrom.size()) {
                    mask_per_chrom[iv.chromid].push_back(iv);
                }
            }
        }

        // Initialize the Markov model
        StratifiedMarkovModel model;
        model.init(num_bins, breaks_vec, k);

        int iter_size = 0;
        if (num_iter_positions > 0) {
            // Compute iterator bin size from first two positions on same chrom
            for (int i = 1; i < num_iter_positions; ++i) {
                if (iter_chroms[i] == iter_chroms[i - 1]) {
                    iter_size = iter_starts[i] - iter_starts[i - 1];
                    break;
                }
            }
            if (iter_size <= 0) iter_size = 1;
        }

        // Set up sequence fetcher
        GenomeSeqFetch seqfetch;
        seqfetch.set_seqdir(string(rdb::get_groot(_envir)) + "/seq");

        // Compute total range for progress reporting
        uint64_t total_range = 0;
        for (int c = 0; c < num_chroms; ++c) {
            total_range += chrom_ends[c] - chrom_starts[c];
        }

        Progress_reporter progress;
        progress.init(total_range, 1000000);

        // Statistics
        uint64_t total_valid = 0;
        uint64_t total_masked = 0;
        uint64_t total_n = 0;

        // Build per-chromosome bin lookup
        int num_chroms_key = iu.get_chromkey().get_num_chroms();
        vector<vector<pair<int64_t, int>>> chrom_bins(num_chroms_key);
        for (int i = 0; i < num_iter_positions; ++i) {
            int chromid = iter_chroms[i];
            if (chromid >= 0 && chromid < num_chroms_key) {
                chrom_bins[chromid].push_back({iter_starts[i], bin_indices[i]});
            }
        }
        for (int c = 0; c < num_chroms_key; ++c) {
            sort(chrom_bins[c].begin(), chrom_bins[c].end());
        }

        // Process each chromosome
        for (int c = 0; c < num_chroms; ++c) {
            int chromid = chrom_ids[c];
            int64_t start = chrom_starts[c];
            int64_t end = chrom_ends[c];

            // Load chromosome sequence
            GInterval chrom_interval(chromid, start, end, 0);
            vector<char> seq;
            seqfetch.read_interval(chrom_interval, iu.get_chromkey(), seq);

            if ((int)seq.size() < kmer_len) continue;  // Need at least (k+1) bp

            // Get mask intervals for this chromosome
            const vector<GInterval>& mask_ivs = mask_per_chrom[chromid];
            size_t mask_cursor = 0;

            const vector<pair<int64_t, int>>& bins = chrom_bins[chromid];
            size_t bin_cursor = 0;

            // Scan with (k+1)-mer sliding window
            int64_t seq_size = static_cast<int64_t>(seq.size());
            for (int64_t pos = 0; pos <= seq_size - kmer_len; ++pos) {
                int64_t genome_pos = start + pos;

                // Check if masked
                if (is_position_masked(genome_pos, mask_ivs, mask_cursor)) {
                    ++total_masked;
                    continue;
                }

                // Check for N's in the (k+1)-mer
                bool has_n = false;
                for (int i = 0; i < kmer_len; ++i) {
                    char base = seq[pos + i];
                    if (StratifiedMarkovModel::encode_base(base) < 0) {
                        has_n = true;
                        break;
                    }
                }
                if (has_n) {
                    ++total_n;
                    continue;
                }

                // Find bin for this position using a forward cursor
                int bin_idx = -1;
                if (!bins.empty()) {
                    while (bin_cursor + 1 < bins.size() &&
                           genome_pos >= bins[bin_cursor + 1].first) {
                        ++bin_cursor;
                    }
                    if (genome_pos >= bins[bin_cursor].first &&
                        genome_pos < bins[bin_cursor].first + iter_size) {
                        bin_idx = bins[bin_cursor].second;
                    }
                }

                if (bin_idx < 0 || bin_idx >= num_bins) {
                    continue;  // Position not covered by iterator
                }

                // Encode k-mer context and next base
                int context_idx = StratifiedMarkovModel::encode_kmer(&seq[pos], k);
                int next_base_idx =
                    StratifiedMarkovModel::encode_base(seq[pos + k]);

                if (context_idx >= 0 && next_base_idx >= 0) {
                    // Add forward strand count
                    model.increment_count(bin_idx, context_idx, next_base_idx);

                    // Add reverse complement count for strand symmetry
                    // This ensures the model learns symmetric transition probabilities
                    int revcomp_context_idx, revcomp_next_idx;
                    StratifiedMarkovModel::revcomp_kmer(
                        context_idx, next_base_idx, k,
                        revcomp_context_idx, revcomp_next_idx);
                    model.increment_count(bin_idx, revcomp_context_idx, revcomp_next_idx);

                    // Count as 2 k-mers (forward + reverse complement)
                    total_valid += 2;
                }
            }

            progress.report(end - start);
            check_interrupt();
        }

        progress.report_last();

        // Apply bin mapping
        model.apply_bin_mapping(bin_map_vec);

        // Resolve prior pi(b) from m_counts (post-merge), explicit matrix,
        // or as a uniform/global broadcast — depending on prior_mode.
        int marginal_fallbacks = 0;
        if (prior_mode == "uniform") {
            model.set_prior_uniform();
        } else if (prior_mode == "marginal") {
            marginal_fallbacks = model.set_prior_from_marginal();
        } else if (prior_mode == "global") {
            model.set_prior_from_global_marginal();
        } else if (prior_mode == "explicit") {
            if (Rf_isNull(_prior_matrix)) {
                verror("prior_mode='explicit' requires a non-NULL prior_matrix");
            }
            int nrow = Rf_nrows(_prior_matrix);
            int ncol = Rf_ncols(_prior_matrix);
            int n_bins_out = model.get_num_bins();
            if (nrow != n_bins_out || ncol != NUM_BASES) {
                verror("prior_matrix must be %d x %d (got %d x %d)",
                       n_bins_out, NUM_BASES, nrow, ncol);
            }
            std::vector<std::array<double, NUM_BASES>> pi_rows(n_bins_out);
            double* mat = REAL(_prior_matrix);
            for (int b = 0; b < n_bins_out; ++b) {
                for (int a = 0; a < NUM_BASES; ++a) {
                    // R column-major: [b + a * nrow]
                    pi_rows[b][a] = mat[b + a * nrow];
                }
            }
            model.set_prior_explicit(pi_rows);
        } else {
            verror("Unknown prior_mode: %s", prior_mode.c_str());
        }

        // Normalize and build CDFs (uses m_prior set above)
        model.normalize_and_build_cdf(pseudocount);

        // Build return list
        SEXP answer, names;
        rprotect(answer = Rf_allocVector(VECSXP, 9));
        rprotect(names = Rf_allocVector(STRSXP, 9));

        // num_bins
        SET_VECTOR_ELT(answer, 0, Rf_ScalarInteger(num_bins));
        SET_STRING_ELT(names, 0, Rf_mkChar("num_bins"));

        // breaks
        SEXP r_breaks;
        rprotect(r_breaks = Rf_allocVector(REALSXP, num_breaks));
        memcpy(REAL(r_breaks), breaks_vec.data(), num_breaks * sizeof(double));
        SET_VECTOR_ELT(answer, 1, r_breaks);
        SET_STRING_ELT(names, 1, Rf_mkChar("breaks"));

        // total_kmers
        SET_VECTOR_ELT(answer, 2, Rf_ScalarReal(static_cast<double>(model.get_total_kmers())));
        SET_STRING_ELT(names, 2, Rf_mkChar("total_kmers"));

        // per_bin_kmers
        SEXP r_per_bin;
        rprotect(r_per_bin = Rf_allocVector(REALSXP, num_bins));
        for (int i = 0; i < num_bins; ++i) {
            REAL(r_per_bin)[i] = static_cast<double>(model.get_bin_kmers(i));
        }
        SET_VECTOR_ELT(answer, 3, r_per_bin);
        SET_STRING_ELT(names, 3, Rf_mkChar("per_bin_kmers"));

        // total_masked
        SET_VECTOR_ELT(answer, 4, Rf_ScalarReal(static_cast<double>(total_masked)));
        SET_STRING_ELT(names, 4, Rf_mkChar("total_masked"));

        // total_n
        SET_VECTOR_ELT(answer, 5, Rf_ScalarReal(static_cast<double>(total_n)));
        SET_STRING_ELT(names, 5, Rf_mkChar("total_n"));

        // Store the model data as a raw vector for later use
        // We'll serialize it to a temporary buffer
        // For now, we store counts and CDFs as nested lists

        // counts: list of matrices (num_bins x (num_kmers * 4))
        SEXP r_counts, r_cdf;
        rprotect(r_counts = Rf_allocVector(VECSXP, num_bins));
        rprotect(r_cdf = Rf_allocVector(VECSXP, num_bins));

        const auto& model_counts = model.get_counts();
        const auto& model_cdf = model.get_cdf();
        int num_kmers = model.get_num_kmers();

        for (int b = 0; b < num_bins; ++b) {
            // Counts matrix: num_kmers rows x 4 cols (column-major for R)
            SEXP count_mat;
            rprotect(count_mat = Rf_allocMatrix(REALSXP, num_kmers, NUM_BASES));
            double* count_data = REAL(count_mat);

            SEXP cdf_mat;
            rprotect(cdf_mat = Rf_allocMatrix(REALSXP, num_kmers, NUM_BASES));
            double* cdf_data = REAL(cdf_mat);

            for (int ctx = 0; ctx < num_kmers; ++ctx) {
                for (int base = 0; base < NUM_BASES; ++base) {
                    // R matrices are column-major: [row + col * nrow]
                    count_data[ctx + base * num_kmers] =
                        static_cast<double>(model_counts[b][ctx * NUM_BASES + base]);
                    cdf_data[ctx + base * num_kmers] =
                        static_cast<double>(model_cdf[b][ctx * NUM_BASES + base]);
                }
            }

            SET_VECTOR_ELT(r_counts, b, count_mat);
            SET_VECTOR_ELT(r_cdf, b, cdf_mat);
        }

        // Create a combined data structure with counts and CDFs
        SEXP r_model_data;
        rprotect(r_model_data = Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(r_model_data, 0, r_counts);
        SET_VECTOR_ELT(r_model_data, 1, r_cdf);

        SEXP model_names;
        rprotect(model_names = Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(model_names, 0, Rf_mkChar("counts"));
        SET_STRING_ELT(model_names, 1, Rf_mkChar("cdf"));
        Rf_setAttrib(r_model_data, R_NamesSymbol, model_names);

        SET_VECTOR_ELT(answer, 6, r_model_data);
        SET_STRING_ELT(names, 6, Rf_mkChar("model_data"));

        // Resolved per-bin prior matrix (n_bins x 4, column-major for R)
        SEXP r_prior;
        rprotect(r_prior = Rf_allocMatrix(REALSXP, num_bins, NUM_BASES));
        {
            const auto& prior = model.get_prior();
            double* prior_data = REAL(r_prior);
            for (int b = 0; b < num_bins; ++b) {
                for (int a = 0; a < NUM_BASES; ++a) {
                    prior_data[b + a * num_bins] = prior[b][a];
                }
            }
        }
        SET_VECTOR_ELT(answer, 7, r_prior);
        SET_STRING_ELT(names, 7, Rf_mkChar("prior"));

        SET_VECTOR_ELT(answer, 8, Rf_ScalarInteger(marginal_fallbacks));
        SET_STRING_ELT(names, 8, Rf_mkChar("marginal_fallbacks"));

        Rf_setAttrib(answer, R_NamesSymbol, names);

        return answer;

    } catch (TGLException& e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc& e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

}  // extern "C"
