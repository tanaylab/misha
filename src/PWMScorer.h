// PWMScorer.h
#ifndef PWM_SCORER_H_
#define PWM_SCORER_H_

#include <string>
#include <memory>
#include <R.h>
#include <Rinternals.h>
#include "DnaPSSM.h"
#include "GenomeSeqFetch.h"
#include "GInterval.h"
#include "GenomeChromKey.h"

class PWMScorer {
public:
    enum ScoringMode {
        TOTAL_LIKELIHOOD,  // For PWM function
        MAX_LIKELIHOOD     // For PWM_MAX function
    };

    PWMScorer(const DnaPSSM& pssm, const std::string& genome_root, ScoringMode mode = TOTAL_LIKELIHOOD) 
        : m_pssm(pssm), m_mode(mode) {
        m_seqfetch.set_seqdir(genome_root + "/seq");
    }

    // Calculate PWM score for a given interval
    double score_interval(const GInterval& interval, const GenomeChromKey& chromkey) {
        std::vector<char> seq;
        try {
            m_seqfetch.read_interval(interval, chromkey, seq);

            if (seq.size() < m_pssm.size()) {
                return std::numeric_limits<double>::quiet_NaN();
            }

            std::string target(seq.begin(), seq.end());

            if (m_mode == TOTAL_LIKELIHOOD) {
                float energy;
                m_pssm.integrate_like(target, energy);
                return exp(energy);
            } else {
                float best_logp;
                int best_dir;
                m_pssm.max_like_match(target, best_logp, best_dir);
                return exp(best_logp);
            }
        } catch (TGLException &e) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Static helper to create PWM from R matrix
    static DnaPSSM create_pssm_from_matrix(SEXP matrix) {
        DnaPSSM pssm;
        int col_num = Rf_ncols(matrix);
        pssm.resize(col_num);

        // Get row indices for A,C,G,T
        SEXP rownames = Rf_getAttrib(matrix, R_RowNamesSymbol);
        int a_idx = -1, c_idx = -1, g_idx = -1, t_idx = -1;
        for (int i = 0; i < Rf_length(rownames); i++) {
            const char *name = CHAR(STRING_ELT(rownames, i));
            if (strcmp(name, "A") == 0) a_idx = i;
            else if (strcmp(name, "C") == 0) c_idx = i;
            else if (strcmp(name, "G") == 0) g_idx = i;
            else if (strcmp(name, "T") == 0) t_idx = i;
        }

        double *matrix_ptr = REAL(matrix);

        // Initialize PSSM probabilities for each position
        for (int i = 0; i < col_num; i++) {
            std::vector<float> probs(4);
            probs[0] = matrix_ptr[a_idx + i * 4]; // A
            probs[1] = matrix_ptr[c_idx + i * 4]; // C
            probs[2] = matrix_ptr[g_idx + i * 4]; // G
            probs[3] = matrix_ptr[t_idx + i * 4]; // T
            pssm[i] = DnaProbVec(probs[0], probs[1], probs[2], probs[3]);
        }

        return pssm;
    }

private:
    DnaPSSM m_pssm;
    GenomeSeqFetch m_seqfetch;
    ScoringMode m_mode;
};

#endif // PWM_SCORER_H_