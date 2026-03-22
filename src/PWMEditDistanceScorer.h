#ifndef PWM_EDIT_DISTANCE_SCORER_H_
#define PWM_EDIT_DISTANCE_SCORER_H_

#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <limits>
#include "GenomeSeqScorer.h"
#include "DnaPSSM.h"

/**
 * PWMEditDistanceScorer: Computes edit-distance-related metrics for PWM windows.
 *
 * Supported return modes:
 *  - MIN_EDITS:        minimum edits across all positions in the interval
 *  - MIN_EDITS_POSITION: position (1-based, signed like pwm.max.pos) of the min-edits window
 *  - PWM_MAX_EDITS:    edits required at the pwm.max / pwm.max.pos location
 *
 * Each window evaluation uses an O(L) algorithm that leverages precomputed gain
 * bins to avoid per-window sorting when running in exact mode
 * (max_edits < 0). When max_edits >= 1, a fast heuristic is used that only
 * considers the top-k gains.
 *
 * score.min filtering: when m_score_min is not NaN, windows whose PWM score
 * is below m_score_min are skipped (edit distance returns NA for that window).
 * For MIN_EDITS/MIN_EDITS_POSITION modes, only windows passing score.min are
 * considered. For PWM_MAX_EDITS mode, the best-PWM window is found first
 * (regardless of score.min), then its edit distance is only computed if
 * its score >= score.min.
 */
class PWMEditDistanceScorer : public GenomeSeqScorer
{
public:
    enum class Mode {
        MIN_EDITS,
        MIN_EDITS_POSITION,
        PWM_MAX_EDITS
    };

    /**
     * Constructor
     * @param pssm The Position-Specific Scoring Matrix
     * @param shared_seqfetch Shared sequence fetcher for caching
     * @param threshold PWM score threshold to reach
     * @param max_edits Maximum edits to consider (-1 = exact, >=1 = heuristic)
     * @param extend Allow motif to extend beyond interval boundaries
     * @param strand Strand mode (0=both, 1=forward, -1=reverse)
     * @param mode Return mode
     * @param score_min Minimum PWM score filter (NaN = no filter)
     * @param max_indels Maximum number of insertions+deletions allowed (0 = substitutions only)
     */
    PWMEditDistanceScorer(const DnaPSSM& pssm,
                          GenomeSeqFetch* shared_seqfetch,
                          float threshold,
                          int max_edits = -1,
                          bool extend = true,
                          char strand = 0,
                          Mode mode = Mode::MIN_EDITS,
                          float score_min = std::numeric_limits<float>::quiet_NaN(),
                          int max_indels = 0);

    /**
     * Score a genomic interval - returns minimum edits needed
     * @return Number of edits (0+), or NaN if unreachable/exceeds max_edits
     */
    float score_interval(const GInterval& interval, const GenomeChromKey& chromkey) override;

    float get_last_min_edits() const { return m_last_metrics.min_edits; }
    float get_last_min_edits_pos() const { return m_last_metrics.min_edits_position; }
    float get_last_pwm_max_logp() const { return m_last_metrics.best_pwm_logp; }
    float get_last_pwm_max_edits() const { return m_last_metrics.best_pwm_edits; }
    float get_last_pwm_max_pos() const { return m_last_metrics.best_pwm_position; }
    int get_max_indels() const { return m_max_indels; }

private:
    struct ScanMetrics {
        float min_edits = std::numeric_limits<float>::quiet_NaN();
        float min_edits_position = std::numeric_limits<float>::quiet_NaN();
        size_t min_index = 0;
        int min_direction = 1;

        float best_pwm_logp = -std::numeric_limits<float>::infinity();
        float best_pwm_edits = std::numeric_limits<float>::quiet_NaN();
        float best_pwm_position = std::numeric_limits<float>::quiet_NaN();
        size_t best_pwm_index = 0;
        int best_pwm_direction = 1;
    };

    DnaPSSM m_pssm;
    float m_threshold;
    int m_max_edits;  // -1 = exact computation, >=1 = fast heuristic
    int m_max_indels; // 0 = substitutions only, >=1 = allow insertions/deletions via banded NW DP
    Mode m_mode;
    float m_score_min; // NaN = no filter; otherwise skip windows with PWM score < this
    ScanMetrics m_last_metrics;

    // Precomputed tables for exact mode
    std::vector<float> m_col_max_scores;     // s_max[i] - max score per column
    std::vector<float> m_gain_values;        // V[1..M] - sorted descending
    std::vector<std::vector<uint8_t>> m_bin_index;  // bin[i][b] - lookup table
    float m_S_max;                           // sum of column maxima

    // Reusable count vector for compute_exact (PERF-1: touched-list cleanup)
    std::vector<int> m_exact_count;
    std::vector<size_t> m_exact_touched;

    /**
     * Precompute gain value bins and lookup tables (called once in constructor)
     */
    void precompute_tables();

    ScanMetrics evaluate_windows(const std::string& seq, size_t interval_length);
    float encode_position(size_t index, size_t target_length, size_t motif_length, int direction) const;
    inline bool should_scan_forward() const;
    inline bool should_scan_reverse() const;
    inline char complement_base(char base) const;
    inline float compute_window_edits(const char* window_start, int seq_avail, bool reverse);

    /**
     * Compute PWM log-likelihood for a window (used for score.min filtering)
     */
    float compute_window_pwm_score(const char* window_start, bool reverse);

    /**
     * Compute exact minimum edits using histogram method
     * @param seq_ptr Sequence to evaluate (already extracted and possibly reverse-complemented)
     * @param reverse Whether to reverse complement
     * @return Number of edits needed, or NaN if unreachable
     */
    float compute_exact(const char* seq_ptr, bool reverse);

    /**
     * Compute minimum edits using fast heuristic (partial sort)
     * @param seq_ptr Sequence to evaluate
     * @param reverse Whether to reverse complement
     * @param max_k Maximum edits to consider
     * @return Number of edits needed (<=max_k), or NaN if exceeds max_k or unreachable
     */
    float compute_heuristic(const char* seq_ptr, bool reverse, int max_k);

    /**
     * Compute minimum edits using banded Needleman-Wunsch DP with indel support.
     * Tries all sequence window lengths from L-max_indels to L+max_indels,
     * aligning each against the motif with a DP that allows up to max_indels
     * total insertions + deletions.
     * @param seq_ptr Pointer to start of sequence window (must have at least L+max_indels bases available)
     * @param seq_len Total number of sequence bases available from seq_ptr
     * @param reverse Whether to reverse complement the sequence
     * @return Minimum edits needed (substitutions + indels) to reach threshold, or NaN if unreachable
     */
    float compute_with_indels(const char* seq_ptr, int seq_len, bool reverse);

    /**
     * Convert base character to index (A=0, C=1, G=2, T=3)
     * @return Base index, or 4 for unknown bases (N, etc.)
     */
    inline int base_to_index(char base) const {
        switch (base) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 4;  // Unknown base
        }
    }
};

#endif // PWM_EDIT_DISTANCE_SCORER_H_
