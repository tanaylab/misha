#ifndef PWM_SCORER_H_
#define PWM_SCORER_H_

#include <string>
#include <memory>
#include <R.h>
#include <Rinternals.h>
#include "DnaPSSM.h"
#include "rdbutils.h"
#include "GenomeSeqFetch.h"
#include "GInterval.h"
#include "GenomeChromKey.h"

class PWMScorer {
public:
    enum ScoringMode {
        TOTAL_LIKELIHOOD,  // For PWM function
        MAX_LIKELIHOOD,    // For PWM_MAX function
        MAX_LIKELIHOOD_POS // For PWM_MAX_POS function - returns position
    };

    PWMScorer(const DnaPSSM& pssm, const std::string& genome_root, const bool extend = true,ScoringMode mode = TOTAL_LIKELIHOOD) 
        : m_pssm(pssm), m_extend(extend), m_mode(mode) {
        m_seqfetch.set_seqdir(genome_root + "/seq");
    }

    // Calculate PWM score for a given interval
    float score_interval(const GInterval& interval, const GenomeChromKey& chromkey) {
        // Calculate expanded interval to include full motif coverage
        int64_t motif_length = m_pssm.size();
        GInterval expanded_interval = interval;
        
        if (m_extend){
            expanded_interval.end = std::min(expanded_interval.end + (motif_length - 1), 
                                       (int64_t)chromkey.get_chrom_size(interval.chromid));
        } else {
            if ((expanded_interval.end - expanded_interval.start) < motif_length) {
                return std::numeric_limits<float>::quiet_NaN();
            }
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
                string::const_iterator best_pos = m_pssm.max_like_match(target, best_logp, best_dir, combine_strands);
                if (m_mode == MAX_LIKELIHOOD){
                    return best_logp;
                } else { // MAX_LIKELIHOOD_POS
                    // Return signed position - negative for reverse strand match                    
                    float pos = best_pos - target.begin();
                    pos = pos + 1; // return a 1-based position
                    return pos;
                }
            }
        } catch (TGLException &e) {
            return std::numeric_limits<float>::quiet_NaN();
        }
    }

    // Static helper to create PWM from R matrix
    static DnaPSSM create_pssm_from_matrix(SEXP matrix) {
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

private:
    DnaPSSM m_pssm;
    bool m_extend = true;
    GenomeSeqFetch m_seqfetch;
    ScoringMode m_mode;
};

#endif // PWM_SCORER_H_