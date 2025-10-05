#ifndef PWM_SCORER_H_
#define PWM_SCORER_H_

#include <string>
#include <memory>
#include <vector>
#include "GenomeSeqScorer.h"
#include "DnaPSSM.h"

class PWMScorer : public GenomeSeqScorer
{
public:
    enum ScoringMode
    {
        TOTAL_LIKELIHOOD,  // For PWM function
        MAX_LIKELIHOOD,    // For PWM_MAX function
        MAX_LIKELIHOOD_POS // For PWM_MAX_POS function - returns position
    };

    PWMScorer(const DnaPSSM &pssm, const std::string &genome_root, bool extend = true,
              ScoringMode mode = TOTAL_LIKELIHOOD, char strand = 1,
              const std::vector<float>& spat_factor = std::vector<float>(),
              int spat_bin_size = 1);

    // Constructor with shared GenomeSeqFetch for caching
    PWMScorer(const DnaPSSM &pssm, GenomeSeqFetch* shared_seqfetch, bool extend = true,
              ScoringMode mode = TOTAL_LIKELIHOOD, char strand = 1,
              const std::vector<float>& spat_factor = std::vector<float>(),
              int spat_bin_size = 1);

    // Implement the virtual function from the base class
    float score_interval(const GInterval &interval, const GenomeChromKey &chromkey) override;

    // Static helper to create PWM from R matrix
    static DnaPSSM create_pssm_from_matrix(SEXP matrix);

private:
    DnaPSSM m_pssm;
    ScoringMode m_mode;

    // Optional spatial weighting
    bool m_use_spat = false;
    std::vector<float> m_spat_log_factors;
    int m_spat_bin_size = 1;
};

#endif // PWM_SCORER_H_