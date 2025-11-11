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

    // Spatial sliding window cache with bin-aware delta updates
    struct SpatSlideCache {
        // Geometry / semantics
        int chromid = -1;
        char strand_mode = 0;
        int64_t last_interval_start = -1;
        int64_t last_interval_end = -1;
        size_t last_i_min = 0;
        size_t last_i_max = 0;
        size_t W = 0;          // window length (# evaluated positions)
        int B = 1;             // spatial bin size (bp)
        size_t bins = 0;       // number of spatial bins we maintain
        bool valid = false;

        // Ring buffer head: j=0 (new window) is ring index 'head'
        size_t head = 0;

        // Motif caches (contiguous, hot arrays)
        std::vector<float> motif_fwd;    // size W (log-like)
        std::vector<float> motif_rc;     // size W (log-like)
        std::vector<uint8_t> has_fwd;    // 0/1 per j
        std::vector<uint8_t> has_rc;     // 0/1

        // TOTAL (log-sum-exp) per-bin accumulators
        std::vector<double> bin_anchor;  // a_b
        std::vector<double> bin_sum_fwd; // s_b^f = sum exp(m_fwd - a_b)
        std::vector<double> bin_sum_rc;  // s_b^r = sum exp(m_rc  - a_b)
        std::vector<uint8_t> bin_dirty;  // 1 if this bin must be recomputed from scratch

        // MAX / MAX_POS per-bin max of motif-only
        struct BinMax {
            float val = -std::numeric_limits<float>::infinity();
            int idx = -1;     // ring index
            int dir = +1;     // +1 or -1
        };
        std::vector<BinMax> bin_max;

        // COUNT per-bin hit counts (position-level)
        std::vector<int> bin_hits_pos;
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

    // Spatial sliding window methods
    bool can_use_spatial_sliding(const GInterval& orig, const GInterval& expd,
                                 size_t i_min, size_t i_max, size_t motif_len,
                                 size_t& stride) const;

    void spat_seed(const std::string& target, const GInterval& expd,
                   size_t i_min, size_t i_max, size_t motif_len);

    void spat_slide_once(const std::string& target, const GInterval& expd,
                         size_t i_min, size_t i_max, size_t motif_len);

    // Answers (return the window's score for the current mode)
    float spat_answer_TOTAL();
    float spat_answer_MAX();
    float spat_answer_MAXPOS(const std::string& target, const GInterval& expd, size_t motif_len, size_t i_min);
    float spat_answer_COUNT() const;

    // Ring buffer and bin helpers
    inline size_t ring_idx_from_j(size_t j) const;
    inline size_t j_from_ring_idx(size_t ridx) const;
    inline size_t bin_of_j(size_t j) const;

    void compute_motif_at(const std::string& target, size_t i_in_target,
                          float& fwd, float& rc, uint8_t& has_fwd, uint8_t& has_rc);

    // TOTAL per-bin maintenance
    void bin_total_recompute(size_t b);
    void bin_total_add(size_t b, float m, bool is_rc);
    void bin_total_remove(size_t b, float m, bool is_rc);

    // MAX / MAX_POS per-bin maintenance
    void bin_max_maybe_recompute(size_t b);
    void bin_max_consider(size_t b, float m, int ridx, int dir);

    // COUNT helper
    inline bool is_hit(float m, size_t b) const;
    inline float combined_motif_score(size_t ridx) const;
    
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
    SpatSlideCache m_spat_slide;
};

#endif // PWM_SCORER_H_
