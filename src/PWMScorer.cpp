#include "PWMScorer.h"
#include <algorithm>
#include <limits>
#include <cstring> // For strcmp
#include <cmath>   // For log

PWMScorer::PWMScorer(const DnaPSSM& pssm, const std::string& genome_root, bool extend,
                     ScoringMode mode, char strand,
                     const std::vector<float>& spat_factor, int spat_bin_size)
    : GenomeSeqScorer(genome_root, extend, strand), m_pssm(pssm), m_mode(mode)
{
    if (!spat_factor.empty()) {
        m_use_spat = true;
        m_spat_bin_size = std::max(1, spat_bin_size);

        // Precompute log of spatial factors
        m_spat_log_factors.resize(spat_factor.size());
        for (size_t i = 0; i < spat_factor.size(); ++i) {
            m_spat_log_factors[i] = std::log(std::max(1e-30f, spat_factor[i]));
        }
    }
}

float PWMScorer::score_interval(const GInterval& interval, const GenomeChromKey& chromkey)
{
    // Calculate expanded interval to include full motif coverage
    int64_t motif_length = m_pssm.size();
    GInterval expanded_interval = calculate_expanded_interval(interval, chromkey, motif_length);
    expanded_interval.strand = m_strand;

    // Check if interval is too small for motif
    if (!m_extend && (expanded_interval.end - expanded_interval.start) < motif_length) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    std::vector<char> seq;
    try {
        m_seqfetch.read_interval(expanded_interval, chromkey, seq);
        std::string target(seq.begin(), seq.end());

        if (!m_use_spat) {
            // Original behavior (no spatial weighting)
            if (m_mode == TOTAL_LIKELIHOOD) {
                float energy;
                m_pssm.integrate_like(target, energy);
                return energy;
            } else {
                float best_logp;
                int best_dir;
                bool combine_strands = (m_mode == MAX_LIKELIHOOD);
                std::string::const_iterator best_pos = m_pssm.max_like_match(target, best_logp, best_dir, combine_strands);

                if (m_mode == MAX_LIKELIHOOD){
                    return best_logp;
                } else { // MAX_LIKELIHOOD_POS
                    float pos = best_pos - target.begin();
                    pos = pos + 1; // return a 1-based position

                    if (m_strand == -1){
                        // The position is now according to the reverse complement sequence, change it to be according to the original sequence (the plus strand)
                        pos = target.length() - pos - motif_length + 1;
                    }

                    if (m_pssm.is_bidirect()) {
                        // Return signed position - negative for reverse strand match
                        pos = pos * best_dir;
                    }

                    return pos;
                }
            }
        }

        // Spatial path - use integrate_energy with log(spat_factor) or max variant
        if (m_spat_log_factors.empty()) {
            // Safety check - should never happen if m_use_spat is true
            return std::numeric_limits<float>::quiet_NaN();
        }

        if (m_mode == TOTAL_LIKELIHOOD) {
            // Use integrate_energy_logspat for efficiency (avoids recomputing logs)
            float energy;
            m_pssm.integrate_energy_logspat(target, energy, m_spat_log_factors, m_spat_bin_size);
            return energy;
        }

        if (m_mode == MAX_LIKELIHOOD) {
            float energy;
            m_pssm.integrate_energy_max_logspat(target, energy, m_spat_log_factors, m_spat_bin_size);
            return energy;
        }

        // m_mode == MAX_LIKELIHOOD_POS with spatial weighting:
        // Find argmax over positions of (log-likelihood + log(spatial_weight))
        // Use the same logic as integrate_energy_max_logspat but track position
        float best_val = -std::numeric_limits<float>::infinity();
        size_t best_index = 0;
        int best_dir = 1;

        if (target.length() < (size_t)motif_length) {
            return std::numeric_limits<float>::quiet_NaN();
        }

        size_t max_i_idx = std::min<size_t>(m_pssm.get_max_range(), target.size() - motif_length);
        size_t min_i_idx = std::min<size_t>(std::max(0, m_pssm.get_min_range()), max_i_idx);

        int pos = 0;
        for (size_t i = min_i_idx; i <= max_i_idx; ++i) {
            int spat_bin = int(pos / m_spat_bin_size);
            if (spat_bin >= (int)m_spat_log_factors.size()) {
                spat_bin = m_spat_log_factors.size() - 1;
            }
            pos++;

            float spat_log = m_spat_log_factors[spat_bin];

            // Forward strand
            if (m_strand != -1) {
                float logp = 0;
                std::string::const_iterator it = target.begin() + i;
                m_pssm.calc_like(it, logp);
                float val = logp + spat_log;
                if (val > best_val) {
                    best_val = val;
                    best_index = i;
                    best_dir = 1;
                }
            }

            // Reverse strand (if bidirectional)
            if (m_pssm.is_bidirect() && m_strand != 1) {
                float logp_rc = 0;
                std::string::const_iterator it2 = target.begin() + i;
                m_pssm.calc_like_rc(it2, logp_rc);
                float val_rc = logp_rc + spat_log;
                if (val_rc > best_val) {
                    best_val = val_rc;
                    best_index = i;
                    best_dir = -1;
                }
            }
        }

        float pos_result = float(best_index) + 1.0f; // 1-based
        if (m_strand == -1) {
            pos_result = target.length() - pos_result - motif_length + 1;
        }
        if (m_pssm.is_bidirect()) {
            pos_result = pos_result * best_dir;
        }
        return pos_result;

    } catch (TGLException &e) {
        return std::numeric_limits<float>::quiet_NaN();
    }
}

DnaPSSM PWMScorer::create_pssm_from_matrix(SEXP matrix)
{
    DnaPSSM pssm;
    int row_num = Rf_nrows(matrix);  // Number of positions in the PWM
    pssm.resize(row_num);

    // Get colnames properly using R_DimNames
    SEXP dimnames = Rf_getAttrib(matrix, R_DimNamesSymbol);
    if (Rf_isNull(dimnames) || VECTOR_ELT(dimnames, 1) == R_NilValue) {
        rdb::verror("PWM matrix must have column names 'A', 'C', 'G', and 'T'");
    }

    SEXP colnames = VECTOR_ELT(dimnames, 1);
    if (!Rf_isString(colnames) || Rf_length(colnames) != 4) {
        rdb::verror("PWM matrix must have exactly 4 columns labeled 'A', 'C', 'G', and 'T'");
    }

    // Get column indices for A,C,G,T
    int a_idx = -1, c_idx = -1, g_idx = -1, t_idx = -1;
    for (int i = 0; i < 4; i++) {
        const char *name = CHAR(STRING_ELT(colnames, i));
        if (strcmp(name, "A") == 0) a_idx = i;
        else if (strcmp(name, "C") == 0) c_idx = i;
        else if (strcmp(name, "G") == 0) g_idx = i;
        else if (strcmp(name, "T") == 0) t_idx = i;
    }

    // Verify we found all required indices
    if (a_idx == -1 || c_idx == -1 || g_idx == -1 || t_idx == -1) {
        rdb::verror("PWM matrix must have columns labeled 'A', 'C', 'G', and 'T'");
    }

    double *matrix_ptr = REAL(matrix);

    // Initialize PSSM probabilities for each position
    for (int i = 0; i < row_num; i++) {
        float pa = matrix_ptr[i + a_idx * row_num];  // Get A probability for position i
        float pc = matrix_ptr[i + c_idx * row_num];  // Get C probability for position i
        float pg = matrix_ptr[i + g_idx * row_num];  // Get G probability for position i
        float pt = matrix_ptr[i + t_idx * row_num];  // Get T probability for position i
        
        pssm[i] = DnaProbVec(pa, pc, pg, pt);
    }

    return pssm;
}