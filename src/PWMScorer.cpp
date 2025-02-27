#include "PWMScorer.h"
#include <algorithm>
#include <limits>
#include <cstring> // For strcmp

PWMScorer::PWMScorer(const DnaPSSM& pssm, const std::string& genome_root, bool extend, 
                     ScoringMode mode, char strand)
    : GenomeSeqScorer(genome_root, extend, strand), m_pssm(pssm), m_mode(mode)
{
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