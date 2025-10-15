#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif

#include <cstdint>
#include "port.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>

#include "DnaPSSM.h"
#include "PWMScorer.h"
#include "GenomeUtils.h"

// Helper struct for score results
struct ScoreResult {
    double value;
    bool has_pos;
    int pos_1based;
    int strand;  // +1 or -1

    ScoreResult() : value(0.0), has_pos(false), pos_1based(NA_INTEGER), strand(1) {}
};

// Gap character set for fast lookup
struct GapCharset {
    bool is_gap[256];

    GapCharset() {
        for (int i = 0; i < 256; ++i) {
            is_gap[i] = false;
        }
    }

    void add_gap_char(char c) {
        is_gap[static_cast<unsigned char>(c)] = true;
    }
};

// Gap projector: creates a compacted view of a sequence with gaps removed
struct GapProjector {
    std::string comp;                // compacted sequence (non-gap bases, uppercased)
    std::vector<int> log_to_phys;    // mapping from logical (compacted) to physical indices (0-based)

    // Build projector from a physical range (inclusive, 0-based)
    static GapProjector build(const std::string& full, const GapCharset& gaps,
                             int phys_lo0, int phys_hi0) {
        GapProjector gp;
        gp.comp.reserve(phys_hi0 - phys_lo0 + 1);
        gp.log_to_phys.reserve(phys_hi0 - phys_lo0 + 1);

        for (int i = phys_lo0; i <= phys_hi0; ++i) {
            char c = full[i];
            if (!gaps.is_gap[static_cast<unsigned char>(c)]) {
                gp.comp.push_back(::toupper(c));
                gp.log_to_phys.push_back(i);
            }
        }

        return gp;
    }

    // Find first logical index j where log_to_phys[j] >= phys0
    // Returns -1 if no such index exists
    int first_log_ge_phys(int phys0) const {
        auto it = std::lower_bound(log_to_phys.begin(), log_to_phys.end(), phys0);
        if (it == log_to_phys.end()) return -1;
        return static_cast<int>(it - log_to_phys.begin());
    }

    // Find last logical index j such that log_to_phys[j+w-1] <= phys0
    // Returns -1 if no such index exists
    int last_log_window_end_le_phys(int phys0, int w) const {
        if (static_cast<int>(comp.size()) < w) return -1;

        // We need log_to_phys[j+w-1] <= phys0
        // So we need j+w-1 <= index where log_to_phys[index] <= phys0
        // Which means j <= (index - w + 1)
        auto it = std::upper_bound(log_to_phys.begin(), log_to_phys.end(), phys0);
        if (it == log_to_phys.begin()) return -1;
        --it;  // Now points to last element <= phys0

        int max_end_idx = static_cast<int>(it - log_to_phys.begin());
        int j_max = max_end_idx - w + 1;

        if (j_max < 0) return -1;
        return j_max;
    }
};

// Helper to convert R character to upper case std::string
static std::string r_string_to_upper(const char* cstr) {
    std::string s(cstr);
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

// Helper to get integer from SEXP with bounds checking
static int get_int_elem(SEXP r_vec, int idx, int default_val) {
    if (Rf_isNull(r_vec) || Rf_length(r_vec) == 0) {
        return default_val;
    }
    int len = Rf_length(r_vec);
    int actual_idx = idx % len;  // recycle
    int val = INTEGER(r_vec)[actual_idx];
    return (val == NA_INTEGER) ? default_val : val;
}

// Compute ROI bounds and allowed window starts
static void compute_bounds(int L, int w, int roi_start, int roi_end, int extend_val,
                          int& start_min0, int& start_max0) {
    // Convert 1-based ROI to 0-based for internal use
    int roi_start0 = roi_start - 1;
    int roi_end0 = roi_end - 1;
    
    // Compute allowed window starts (0-based)
    start_min0 = std::max(0, roi_start0 - extend_val);
    start_max0 = std::min(L - w, roi_end0 - w + 1 + extend_val);
}

// Score PWM over a range (with optional gap support)
static ScoreResult score_pwm_over_range(
    const std::string& seq,
    int start_min0, int start_max0,
    const DnaPSSM& pssm,
    const std::string& mode,
    bool bidirect,
    int strand_mode,  // 0=both, +1=fwd, -1=rev
    double score_thresh,
    int roi_start1,
    const std::vector<float>& spat_log_factors,
    int spat_bin_size,
    bool skip_gaps = false,
    const GapCharset* gaps = nullptr
) {
    ScoreResult result;
    int w = pssm.size();
    bool use_spat = !spat_log_factors.empty();

    // Helper lambda for log-sum-exp
    auto log_sum_exp_add = [](double a, double b) -> double {
        if (a == R_NegInf) return b;
        if (b == R_NegInf) return a;
        if (a > b) {
            return a + log1p(exp(b - a));
        } else {
            return b + log1p(exp(a - b));
        }
    };

    // Branch on gap-skipping mode
    if (!skip_gaps || gaps == nullptr) {
        // ========== Original contiguous path ==========
        // Handle empty range
        if (start_max0 < start_min0) {
            if (mode == "count") {
                result.value = 0.0;
            } else {
                result.value = R_NaReal;
            }
            return result;
        }

        double best_score = R_NegInf;
        int best_start0 = -1;
        int best_strand = 1;
        double total = R_NegInf;
        int count = 0;

        // Scan all allowed starts
        for (int s0 = start_min0; s0 <= start_max0; ++s0) {
            // Compute spatial weight
            int offset_from_roi = (s0 + 1) - roi_start1;
            float spat_log = 0.0f;
            if (use_spat) {
                int spat_bin = offset_from_roi / spat_bin_size;
                if (spat_bin < 0) spat_bin = 0;
                if (spat_bin >= (int)spat_log_factors.size()) {
                    spat_bin = spat_log_factors.size() - 1;
                }
                spat_log = spat_log_factors[spat_bin];
            }

            // Score forward and/or reverse strand
            float fwd_score = R_NegInf;
            float rev_score = R_NegInf;

            // Extract substring for this window
            std::string window = seq.substr(s0, w);

            if (strand_mode >= 0) {  // forward or both
                float logp = 0;
                pssm.calc_like(window, logp);
                fwd_score = logp;
            }

            if ((bidirect || strand_mode <= 0) && strand_mode != 1) {  // reverse or both
                float logp = 0;
                pssm.calc_like_rc(window, logp);
                rev_score = logp;
            }

            // Apply spatial weighting
            if (use_spat) {
                if (fwd_score > R_NegInf) fwd_score += spat_log;
                if (rev_score > R_NegInf) rev_score += spat_log;
            }

            // Process based on mode
            if (mode == "lse") {
                if (bidirect) {
                    total = log_sum_exp_add(total, fwd_score);
                    total = log_sum_exp_add(total, rev_score);
                } else {
                    double this_score = (strand_mode >= 0) ? fwd_score : rev_score;
                    total = log_sum_exp_add(total, this_score);
                }
            } else if (mode == "max") {
                double this_max = std::max(fwd_score, rev_score);
                if (this_max > best_score) {
                    best_score = this_max;
                }
            } else if (mode == "count") {
                double this_max = std::max(fwd_score, rev_score);
                if (this_max >= score_thresh) {
                    count++;
                }
            } else if (mode == "pos") {
                // Tie-breaking: prefer smaller start, then prefer forward strand
                if (fwd_score > best_score ||
                    (fwd_score == best_score && s0 < best_start0) ||
                    (fwd_score == best_score && s0 == best_start0 && best_strand == -1)) {
                    best_score = fwd_score;
                    best_start0 = s0;
                    best_strand = 1;
                }
                if (rev_score > best_score ||
                    (rev_score == best_score && s0 < best_start0)) {
                    best_score = rev_score;
                    best_start0 = s0;
                    best_strand = -1;
                }
            }
        }

        // Set result based on mode
        if (mode == "lse") {
            result.value = total;
        } else if (mode == "max") {
            result.value = best_score;
        } else if (mode == "count") {
            result.value = count;
        } else if (mode == "pos") {
            if (best_start0 >= 0) {
                result.has_pos = true;
                result.pos_1based = best_start0 + 1;
                result.strand = best_strand;
            }
            result.value = best_score;
        }

        return result;

    } else {
        // ========== Gap-skipping path ==========
        // Compute the full physical range that windows can span
        int end_lim_phys0 = start_max0 + w - 1;

        // Build gap projector over the full physical range
        GapProjector gp = GapProjector::build(seq, *gaps, start_min0, end_lim_phys0);

        // Compute logical bounds
        int j_min = gp.first_log_ge_phys(start_min0);
        int j_max = gp.last_log_window_end_le_phys(end_lim_phys0, w);

        // Check if we have enough non-gap bases
        if (j_min < 0 || j_max < j_min || static_cast<int>(gp.comp.size()) < w) {
            if (mode == "count") {
                result.value = 0.0;
            } else {
                result.value = R_NaReal;
            }
            return result;
        }

        double best_score = R_NegInf;
        int best_j = -1;
        int best_strand = 1;
        double total = R_NegInf;
        int count = 0;

        // Scan logical starts
        for (int j = j_min; j <= j_max; ++j) {
            int start_phys0 = gp.log_to_phys[j];
            int start_phys1 = start_phys0 + 1;

            // Compute spatial weight based on physical position
            int offset_from_roi = start_phys1 - roi_start1;
            float spat_log = 0.0f;
            if (use_spat) {
                int spat_bin = offset_from_roi / spat_bin_size;
                if (spat_bin < 0) spat_bin = 0;
                if (spat_bin >= (int)spat_log_factors.size()) {
                    spat_bin = spat_log_factors.size() - 1;
                }
                spat_log = spat_log_factors[spat_bin];
            }

            // Score using compacted sequence
            float fwd_score = R_NegInf;
            float rev_score = R_NegInf;

            // Extract window from compacted sequence
            std::string window = gp.comp.substr(j, w);

            if (strand_mode >= 0) {
                float logp = 0;
                pssm.calc_like(window, logp);
                fwd_score = logp;
            }

            if ((bidirect || strand_mode <= 0) && strand_mode != 1) {
                float logp = 0;
                pssm.calc_like_rc(window, logp);
                rev_score = logp;
            }

            // Apply spatial weighting
            if (use_spat) {
                if (fwd_score > R_NegInf) fwd_score += spat_log;
                if (rev_score > R_NegInf) rev_score += spat_log;
            }

            // Process based on mode
            if (mode == "lse") {
                if (bidirect) {
                    total = log_sum_exp_add(total, fwd_score);
                    total = log_sum_exp_add(total, rev_score);
                } else {
                    double this_score = (strand_mode >= 0) ? fwd_score : rev_score;
                    total = log_sum_exp_add(total, this_score);
                }
            } else if (mode == "max") {
                double this_max = std::max(fwd_score, rev_score);
                if (this_max > best_score) {
                    best_score = this_max;
                }
            } else if (mode == "count") {
                double this_max = std::max(fwd_score, rev_score);
                if (this_max >= score_thresh) {
                    count++;
                }
            } else if (mode == "pos") {
                // Tie-breaking: prefer smaller physical start, then prefer forward strand
                if (fwd_score > best_score ||
                    (fwd_score == best_score && start_phys0 < (best_j >= 0 ? gp.log_to_phys[best_j] : INT_MAX)) ||
                    (fwd_score == best_score && start_phys0 == (best_j >= 0 ? gp.log_to_phys[best_j] : 0) && best_strand == -1)) {
                    best_score = fwd_score;
                    best_j = j;
                    best_strand = 1;
                }
                if (rev_score > best_score ||
                    (rev_score == best_score && start_phys0 < (best_j >= 0 ? gp.log_to_phys[best_j] : INT_MAX))) {
                    best_score = rev_score;
                    best_j = j;
                    best_strand = -1;
                }
            }
        }

        // Set result based on mode
        if (mode == "lse") {
            result.value = total;
        } else if (mode == "max") {
            result.value = best_score;
        } else if (mode == "count") {
            result.value = count;
        } else if (mode == "pos") {
            if (best_j >= 0) {
                result.has_pos = true;
                result.pos_1based = gp.log_to_phys[best_j] + 1;  // physical position, 1-based
                result.strand = best_strand;
            }
            result.value = best_score;
        }

        return result;
    }
}

// Score k-mer over a range (with optional gap support)
// Returns count or fraction depending on return_frac
static ScoreResult score_kmer_over_range(
    const std::string& seq,
    int start_min0, int start_max0,
    const std::string& kmer,
    int strand_mode,
    bool return_frac,
    bool skip_gaps = false,
    const GapCharset* gaps = nullptr
) {
    ScoreResult result;
    int w = kmer.length();

    std::string kmer_upper = kmer;
    std::transform(kmer_upper.begin(), kmer_upper.end(), kmer_upper.begin(), ::toupper);

    // Compute reverse complement if needed
    std::string kmer_rc;
    if (strand_mode <= 0) {  // reverse or both
        kmer_rc = seq2reverse_complementary(kmer_upper);
    }

    // Branch on gap-skipping mode
    if (!skip_gaps || gaps == nullptr) {
        // ========== Original contiguous path ==========
        // Handle empty range
        if (start_max0 < start_min0) {
            result.value = 0.0;
            return result;
        }

        int count = 0;
        int nstarts = start_max0 - start_min0 + 1;

        // Scan all allowed starts
        for (int s0 = start_min0; s0 <= start_max0; ++s0) {
            bool match_fwd = false;
            bool match_rev = false;

            // Check forward strand
            if (strand_mode >= 0) {
                match_fwd = (seq.compare(s0, w, kmer_upper) == 0);
            }

            // Check reverse strand
            if (strand_mode <= 0) {
                match_rev = (seq.compare(s0, w, kmer_rc) == 0);
            }

            if (match_fwd || match_rev) {
                count++;
            }
        }

        if (return_frac && nstarts > 0) {
            result.value = static_cast<double>(count) / static_cast<double>(nstarts);
        } else {
            result.value = count;
        }
        return result;

    } else {
        // ========== Gap-skipping path ==========
        // Compute the full physical range that windows can span
        int end_lim_phys0 = start_max0 + w - 1;

        // Build gap projector over the full physical range
        GapProjector gp = GapProjector::build(seq, *gaps, start_min0, end_lim_phys0);

        // Compute logical bounds
        int j_min = gp.first_log_ge_phys(start_min0);
        int j_max = gp.last_log_window_end_le_phys(end_lim_phys0, w);

        // Check if we have enough non-gap bases
        if (j_min < 0 || j_max < j_min || static_cast<int>(gp.comp.size()) < w) {
            result.value = 0.0;
            return result;
        }

        int count = 0;
        int nstarts = j_max - j_min + 1;

        // Scan logical starts
        for (int j = j_min; j <= j_max; ++j) {
            bool match_fwd = false;
            bool match_rev = false;

            // Check forward strand
            if (strand_mode >= 0) {
                match_fwd = (gp.comp.compare(j, w, kmer_upper) == 0);
            }

            // Check reverse strand
            if (strand_mode <= 0) {
                match_rev = (gp.comp.compare(j, w, kmer_rc) == 0);
            }

            if (match_fwd || match_rev) {
                count++;
            }
        }

        if (return_frac && nstarts > 0) {
            result.value = static_cast<double>(count) / static_cast<double>(nstarts);
        } else {
            result.value = count;
        }
        return result;
    }
}

extern "C" {

// Main PWM scoring function
SEXP C_gseq_pwm(SEXP r_seqs, SEXP r_pssm, SEXP r_mode, SEXP r_bidirect,
                SEXP r_strand_mode, SEXP r_score_thresh,
                SEXP r_roi_start, SEXP r_roi_end, SEXP r_extend,
                SEXP r_spat_params, SEXP r_return_strand,
                SEXP r_skip_gaps, SEXP r_gap_chars) {
    try {
        // Validate and extract inputs
        if (!Rf_isString(r_seqs)) {
            Rf_error("seqs must be a character vector");
        }
        if (!Rf_isMatrix(r_pssm) || !Rf_isReal(r_pssm)) {
            Rf_error("pssm must be a numeric matrix");
        }

        int n_seqs = Rf_length(r_seqs);
        if (n_seqs == 0) {
            return Rf_allocVector(REALSXP, 0);
        }

        // Extract PSSM
        DnaPSSM pssm = PWMScorer::create_pssm_from_matrix(r_pssm);
        int w = pssm.size();

        // Extract parameters
        std::string mode = CHAR(STRING_ELT(r_mode, 0));
        bool bidirect = Rf_asLogical(r_bidirect);
        int strand_mode = Rf_asInteger(r_strand_mode);
        double score_thresh = Rf_asReal(r_score_thresh);
        bool return_strand = Rf_asLogical(r_return_strand);

        // Handle extend parameter
        int extend_val;
        if (Rf_isLogical(r_extend)) {
            extend_val = Rf_asLogical(r_extend) ? (w - 1) : 0;
        } else {
            extend_val = Rf_asInteger(r_extend);
            if (extend_val < 0) extend_val = 0;
        }

        // Extract gap parameters
        bool skip_gaps = Rf_asLogical(r_skip_gaps);
        GapCharset gap_charset;
        if (skip_gaps && !Rf_isNull(r_gap_chars) && Rf_isString(r_gap_chars)) {
            int n_gaps = Rf_length(r_gap_chars);
            for (int i = 0; i < n_gaps; ++i) {
                const char* gap_str = CHAR(STRING_ELT(r_gap_chars, i));
                if (gap_str && gap_str[0] != '\0') {
                    gap_charset.add_gap_char(gap_str[0]);
                }
            }
        }

        // Extract spatial parameters
        std::vector<float> spat_log_factors;
        int spat_bin_size = 1;
        if (!Rf_isNull(r_spat_params) && Rf_isNewList(r_spat_params)) {
            SEXP spat_factor = VECTOR_ELT(r_spat_params, 0);  // spat.factor
            SEXP spat_bin = VECTOR_ELT(r_spat_params, 1);      // spat.bin
            // spat.min and spat.max are at indices 2, 3 but we don't use them here

            if (!Rf_isNull(spat_factor) && Rf_length(spat_factor) > 0) {
                int n_spat = Rf_length(spat_factor);
                spat_log_factors.resize(n_spat);
                double* spat_vals = REAL(spat_factor);
                for (int i = 0; i < n_spat; ++i) {
                    spat_log_factors[i] = log(std::max(1e-30, spat_vals[i]));
                }
            }

            if (!Rf_isNull(spat_bin) && Rf_length(spat_bin) > 0) {
                spat_bin_size = std::max(1, Rf_asInteger(spat_bin));
            }
        }
        
        // Prepare result vectors
        SEXP r_result;
        SEXP r_pos = R_NilValue;
        SEXP r_strand_out = R_NilValue;
        
        if (mode == "pos") {
            PROTECT(r_pos = Rf_allocVector(INTSXP, n_seqs));
            if (return_strand) {
                PROTECT(r_strand_out = Rf_allocVector(INTSXP, n_seqs));
            }
        } else {
            PROTECT(r_result = Rf_allocVector(REALSXP, n_seqs));
        }
        
        // Process each sequence
        for (int i = 0; i < n_seqs; ++i) {
            std::string seq = r_string_to_upper(CHAR(STRING_ELT(r_seqs, i)));
            int L = seq.length();
            
            // Get ROI bounds
            int roi_start = get_int_elem(r_roi_start, i, 1);
            int roi_end = get_int_elem(r_roi_end, i, L);
            
            // Validate bounds
            if (roi_start < 1) roi_start = 1;
            if (roi_end > L) roi_end = L;
            if (roi_start > roi_end) {
                // Invalid ROI
                if (mode == "pos") {
                    INTEGER(r_pos)[i] = NA_INTEGER;
                    if (return_strand) INTEGER(r_strand_out)[i] = NA_INTEGER;
                } else if (mode == "count") {
                    REAL(r_result)[i] = 0.0;
                } else {
                    REAL(r_result)[i] = R_NaReal;
                }
                continue;
            }
            
            // Compute allowed window starts
            int start_min0, start_max0;
            compute_bounds(L, w, roi_start, roi_end, extend_val, start_min0, start_max0);
            
            // Score the sequence
            ScoreResult result = score_pwm_over_range(
                seq, start_min0, start_max0, pssm, mode, bidirect, strand_mode,
                score_thresh, roi_start, spat_log_factors, spat_bin_size,
                skip_gaps, skip_gaps ? &gap_charset : nullptr
            );
            
            // Store result
            if (mode == "pos") {
                INTEGER(r_pos)[i] = result.pos_1based;
                if (return_strand) {
                    INTEGER(r_strand_out)[i] = result.has_pos ? result.strand : NA_INTEGER;
                }
            } else {
                REAL(r_result)[i] = result.value;
            }
        }
        
        // Return appropriate structure
        if (mode == "pos") {
            if (return_strand) {
                // Return data.frame with pos and strand
                SEXP df;
                PROTECT(df = Rf_allocVector(VECSXP, 2));
                SET_VECTOR_ELT(df, 0, r_pos);
                SET_VECTOR_ELT(df, 1, r_strand_out);
                
                SEXP names;
                PROTECT(names = Rf_allocVector(STRSXP, 2));
                SET_STRING_ELT(names, 0, Rf_mkChar("pos"));
                SET_STRING_ELT(names, 1, Rf_mkChar("strand"));
                Rf_setAttrib(df, R_NamesSymbol, names);
                
                UNPROTECT(4);
                return df;
            } else {
                UNPROTECT(1);
                return r_pos;
            }
        } else {
            UNPROTECT(1);
            return r_result;
        }
        
    } catch (std::exception& e) {
        Rf_error("Error in C_gseq_pwm: %s", e.what());
    } catch (...) {
        Rf_error("Unknown error in C_gseq_pwm");
    }
    return R_NilValue;
}

// Main k-mer scoring function
SEXP C_gseq_kmer(SEXP r_seqs, SEXP r_kmer, SEXP r_mode, SEXP r_strand_mode,
                 SEXP r_roi_start, SEXP r_roi_end, SEXP r_extend,
                 SEXP r_skip_gaps, SEXP r_gap_chars) {
    try {
        // Validate and extract inputs
        if (!Rf_isString(r_seqs)) {
            Rf_error("seqs must be a character vector");
        }
        if (!Rf_isString(r_kmer) || Rf_length(r_kmer) != 1) {
            Rf_error("kmer must be a single character string");
        }

        int n_seqs = Rf_length(r_seqs);
        if (n_seqs == 0) {
            return Rf_allocVector(REALSXP, 0);
        }

        // Extract parameters
        std::string kmer = CHAR(STRING_ELT(r_kmer, 0));
        std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
        int w = kmer.length();

        std::string mode = CHAR(STRING_ELT(r_mode, 0));
        int strand_mode = Rf_asInteger(r_strand_mode);

        // Handle extend parameter
        int extend_val;
        if (Rf_isLogical(r_extend)) {
            extend_val = Rf_asLogical(r_extend) ? (w - 1) : 0;
        } else {
            extend_val = Rf_asInteger(r_extend);
            if (extend_val < 0) extend_val = 0;
        }

        // Extract gap parameters
        bool skip_gaps = Rf_asLogical(r_skip_gaps);
        GapCharset gap_charset;
        if (skip_gaps && !Rf_isNull(r_gap_chars) && Rf_isString(r_gap_chars)) {
            int n_gaps = Rf_length(r_gap_chars);
            for (int i = 0; i < n_gaps; ++i) {
                const char* gap_str = CHAR(STRING_ELT(r_gap_chars, i));
                if (gap_str && gap_str[0] != '\0') {
                    gap_charset.add_gap_char(gap_str[0]);
                }
            }
        }
        
        // Prepare result vector
        SEXP r_result;
        PROTECT(r_result = Rf_allocVector(REALSXP, n_seqs));
        
        // Process each sequence
        for (int i = 0; i < n_seqs; ++i) {
            std::string seq = r_string_to_upper(CHAR(STRING_ELT(r_seqs, i)));
            int L = seq.length();
            
            // Get ROI bounds
            int roi_start = get_int_elem(r_roi_start, i, 1);
            int roi_end = get_int_elem(r_roi_end, i, L);
            
            // Validate bounds
            if (roi_start < 1) roi_start = 1;
            if (roi_end > L) roi_end = L;
            if (roi_start > roi_end || L < w) {
                REAL(r_result)[i] = 0.0;
                continue;
            }
            
            // Compute allowed window starts
            int start_min0, start_max0;
            compute_bounds(L, w, roi_start, roi_end, extend_val, start_min0, start_max0);
            
            // Score the sequence
            bool return_frac = (mode == "frac");
            ScoreResult result = score_kmer_over_range(
                seq, start_min0, start_max0, kmer, strand_mode, return_frac,
                skip_gaps, skip_gaps ? &gap_charset : nullptr
            );
            
            REAL(r_result)[i] = result.value;
        }
        
        UNPROTECT(1);
        return r_result;
        
    } catch (std::exception& e) {
        Rf_error("Error in C_gseq_kmer: %s", e.what());
    } catch (...) {
        Rf_error("Unknown error in C_gseq_kmer");
    }
    return R_NilValue;
}

}  // extern "C"

