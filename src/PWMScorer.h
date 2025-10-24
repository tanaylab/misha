#ifndef PWM_SCORER_H_
#define PWM_SCORER_H_

#include <string>
#include <memory>
#include <vector>
#include <deque>
#include "GenomeSeqScorer.h"
#include "DnaPSSM.h"
#include "utils/RunningLogSumExp.h"
#include "utils/RunningMaxDeque.h"

class PWMScorer : public GenomeSeqScorer
{
public:
    enum ScoringMode
    {
        TOTAL_LIKELIHOOD,   // Integrated log-likelihood across all positions
        MAX_LIKELIHOOD,     // Maximum log-likelihood score
        MAX_LIKELIHOOD_POS, // Position of maximum log-likelihood
        MOTIF_COUNT         // Count of positions exceeding threshold
    };

    PWMScorer(const DnaPSSM &pssm, const std::string &genome_root, bool extend = true,
              ScoringMode mode = TOTAL_LIKELIHOOD, char strand = 1,
              const std::vector<float>& spat_factor = std::vector<float>(),
              int spat_bin_size = 1, float score_thresh = 0.0f);

    // Constructor with shared sequence fetcher for caching
    PWMScorer(const DnaPSSM &pssm, GenomeSeqFetch* shared_seqfetch, bool extend = true,
              ScoringMode mode = TOTAL_LIKELIHOOD, char strand = 1,
              const std::vector<float>& spat_factor = std::vector<float>(),
              int spat_bin_size = 1, float score_thresh = 0.0f);

    // Score a genomic interval using the PWM
    float score_interval(const GInterval &interval, const GenomeChromKey &chromkey) override;

    // Create PSSM from R matrix with columns A, C, G, T
    static DnaPSSM create_pssm_from_matrix(SEXP matrix);

    // Invalidate sliding window cache
    void invalidate_cache();

private:
    // Sliding window cache for contiguous intervals
    struct SlideCache {
        int chromid = -1;
        char strand_mode = 0;     // 1=plus, -1=minus, 0=both
        bool valid = false;
        int64_t last_interval_start = -1;
        int64_t last_interval_end = -1;
        size_t last_i_min = 0;
        size_t last_i_max = 0;
        size_t window_size = 0;
        size_t pos0 = 0;
        size_t stride = 0;

        // Mode-specific aggregators
        RunningLogSumExp rlse;
        RunningMaxDeque  rmax;
        std::deque<uint8_t> hits;
        int hit_count = 0;
    };

    // Sliding window methods
    float score_with_sliding_window(const std::string& target, const GInterval& original_interval,
                                     const GInterval& expanded_interval,
                                     size_t i_min, size_t i_max, size_t motif_len);
    float seed_sliding_window(const std::string& target, const GInterval& original_interval,
                              const GInterval& expanded_interval,
                              size_t i_min, size_t i_max, size_t motif_len);
    float try_slide_window(const std::string& target, const GInterval& original_interval,
                           const GInterval& expanded_interval,
                           size_t i_min, size_t i_max, size_t motif_len, size_t stride);
    
    // Scoring methods
    float score_without_spatial(const std::string& target, int64_t motif_length);
    float score_with_spatial(const std::string& target, int64_t motif_length);
    
    // Motif counting
    float count_motif_hits_no_spatial(const std::string& target, size_t motif_length);
    float count_motif_hits_with_spatial(const std::string& target, size_t motif_length);
    
    // Position finding
    float get_max_likelihood_pos_with_spatial(const std::string& target, size_t motif_length);
    float compute_position_result(size_t index, size_t target_length, size_t motif_length, int direction) const;
    
    // Utilities
    inline float get_spatial_log_factor(size_t pos_index) const;

    // Core members
    DnaPSSM m_pssm;
    ScoringMode m_mode;
    float m_score_thresh = 0.0f;

    // Spatial weighting
    bool m_use_spat = false;
    std::vector<float> m_spat_log_factors;
    int m_spat_bin_size = 1;

    // Cache
    SlideCache m_slide;
};

#endif // PWM_SCORER_H_
