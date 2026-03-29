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
 * score.min/score.max filtering: when m_score_min (m_score_max) is not NaN,
 * windows whose PWM score is below m_score_min (above m_score_max) are
 * skipped (edit distance returns NA for that window).
 * For MIN_EDITS/MIN_EDITS_POSITION modes, only windows passing both filters
 * are considered. For PWM_MAX_EDITS mode, the best-PWM window is found first,
 * then its edit distance is only computed if its score passes both filters.
 */
class PWMEditDistanceScorer : public GenomeSeqScorer
{
public:
    enum class Mode {
        MIN_EDITS,
        MIN_EDITS_POSITION,
        PWM_MAX_EDITS
    };

    enum class Direction {
        ABOVE,  // minimum edits to bring score >= threshold (default)
        BELOW   // minimum edits to bring score <= threshold
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
     * @param score_max Maximum PWM score filter (NaN = no filter)
     * @param direction ABOVE = min edits to reach score >= threshold;
     *                  BELOW = min edits to bring score <= threshold
     */
    PWMEditDistanceScorer(const DnaPSSM& pssm,
                          GenomeSeqFetch* shared_seqfetch,
                          float threshold,
                          int max_edits = -1,
                          bool extend = true,
                          char strand = 0,
                          Mode mode = Mode::MIN_EDITS,
                          float score_min = std::numeric_limits<float>::quiet_NaN(),
                          int max_indels = 0,
                          float score_max = std::numeric_limits<float>::quiet_NaN(),
                          Direction direction = Direction::ABOVE);

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
    Direction m_direction;  // ABOVE = edits to raise score >= threshold; BELOW = edits to lower score <= threshold
    float m_score_min; // NaN = no filter; otherwise skip windows with PWM score < this
    float m_score_max; // NaN = no filter; otherwise skip windows with PWM score > this
    ScanMetrics m_last_metrics;

    // Precomputed tables for exact mode
    std::vector<float> m_col_max_scores;     // s_max[i] - max score per column
    std::vector<float> m_col_min_scores;     // s_min[i] - min score per column (for BELOW direction)
    std::vector<float> m_gain_values;        // V[1..M] - sorted descending (gains for ABOVE, losses for BELOW)
    std::vector<std::vector<uint8_t>> m_bin_index;  // bin[i][b] - lookup table
    float m_S_max;                           // sum of column maxima
    float m_S_min;                           // sum of column minima (for BELOW direction)
    std::vector<float> m_max_suffix_score;   // m_max_suffix_score[i] = sum of col maxima from i to L-1
    std::vector<float> m_max_gain_budget;    // m_max_gain_budget[k] = max total gain/loss from k substitutions (per-PSSM)

    // Flat precomputed PSSM lookup tables for cache-friendly access
    static constexpr int MAX_MOTIF_LEN_OPT = 64;
    float m_score_table[MAX_MOTIF_LEN_OPT][5];    // [motif_pos][base_index 0-3, 4=N]
    float m_gain_table[MAX_MOTIF_LEN_OPT][5];     // col_max - score
    bool m_mandatory_table[MAX_MOTIF_LEN_OPT][5];  // true if score is log-zero or non-finite

    // Reusable count vector for compute_exact (PERF-1: touched-list cleanup)
    std::vector<int> m_exact_count;
    std::vector<size_t> m_exact_touched;

    // Pigeonhole pre-filter: divide motif into (K+1) blocks; if a window matches
    // with at most K total edits, at least one block must match exactly (zero edits)
    // at some shift in {-D, ..., +D}.
    struct PrefilterBlock {
        int start;         // start column in motif
        int len;           // block length (number of columns)
        int num_entries;   // 4^len (size of viable bitset)
        std::vector<uint8_t> viable;  // viable[hash] = true if B-mer can match block with 0 edits
    };
    std::vector<PrefilterBlock> m_prefilter_blocks;
    bool m_use_prefilter;

    /**
     * Precompute gain value bins and lookup tables (called once in constructor)
     */
    void precompute_tables();

    ScanMetrics evaluate_windows(const std::string& seq, size_t interval_length);
    float encode_position(size_t index, size_t target_length, size_t motif_length, int direction) const;
    inline bool should_scan_forward() const;
    inline bool should_scan_reverse() const;
    inline char complement_base(char base) const;
    inline float compute_window_edits(const int* bidx, int seq_avail, bool reverse);

    /**
     * Get the effective PSSM score and gain for a motif position aligned with a sequence base.
     * Handles reverse complement, unknown bases (index==4), and log-zero PSSM entries.
     *
     * When the base has a log-zero (or -Inf) probability, this position is a mandatory
     * substitution: it must be edited regardless of budget. In that case:
     *   out_score = col_max  (score assuming the mandatory substitution is applied)
     *   out_gain  = 0.0f     (no additional gain available from this position)
     *   return value = true  (this position requires a mandatory edit)
     *
     * Normal case (base has finite, non-logzero probability):
     *   out_score = PSSM score for this base
     *   out_gain  = col_max - out_score
     *   return value = false
     *
     * Compatibility wrapper — kept for code outside the hot path.
     *
     * @param seq_ptr Pointer to sequence data
     * @param reverse Whether to reverse complement
     * @param motif_pos Motif column index
     * @param raw_seq_idx Index into seq_ptr for the aligned base
     * @param out_score Output: PSSM score at this alignment (col_max for mandatory edits)
     * @param out_gain Output: col_max - out_score (0 for mandatory edits)
     * @return true if this position requires a mandatory substitution
     */
    inline bool get_aligned_base_score(const char* seq_ptr, bool reverse,
                                       int motif_pos, int raw_seq_idx,
                                       float& out_score, float& out_gain) const;

    /**
     * Compute minimum total edits (indels + substitutions) to reach threshold.
     * Given an aligned score and a vector of per-column substitution gains, greedily
     * picks top gains to cover the deficit.
     * NOTE: Mutates gains via partial sort — caller must not reuse gains after this call.
     * @param aligned_score Sum of PSSM scores for the aligned columns
     * @param gains Per-column substitution gains (will be partially sorted in-place)
     * @param indels Number of indels used in this alignment family
     * @param best_edits_so_far Current best total edits (for pruning); NaN if none found yet
     * @return Total edits (indels + subs), or NaN if unreachable or pruned
     */
    float compute_min_edits_from_gains(double aligned_score, std::vector<float>& gains,
                                       int indels, float best_edits_so_far) const;

    /**
     * Quick O(1) deficit check: can this alignment family possibly beat best_edits_so_far?
     * Uses precomputed per-PSSM max gain budget to avoid expensive gains collection + sorting.
     * Call BEFORE collecting gains to skip hopeless candidates early.
     * @return true if the candidate is worth evaluating (might produce a result)
     */
    bool quick_deficit_check(double aligned_score, int indels, float best_edits_so_far) const;

    /**
     * Compute PWM log-likelihood for a window (used for score.min filtering)
     */
    float compute_window_pwm_score(const char* window_start, bool reverse);

    /**
     * Compute exact minimum edits using histogram method
     * @param bidx Precomputed base index array (forward or reverse-complemented)
     * @param reverse Whether to reverse index order (right-to-left)
     * @return Number of edits needed, or NaN if unreachable
     */
    float compute_exact(const int* bidx, bool reverse);

    /**
     * Compute minimum edits using fast heuristic (partial sort)
     * @param bidx Precomputed base index array (forward or reverse-complemented)
     * @param reverse Whether to reverse index order (right-to-left)
     * @param max_k Maximum edits to consider
     * @return Number of edits needed (<=max_k), or NaN if exceeds max_k or unreachable
     */
    float compute_heuristic(const int* bidx, bool reverse, int max_k);

    /**
     * Compute minimum edits using banded Needleman-Wunsch DP with indel support.
     * Tries all sequence window lengths from L-max_indels to L+max_indels,
     * aligning each against the motif with a DP that allows up to max_indels
     * total insertions + deletions.
     * @param bidx Precomputed base index array (forward or reverse-complemented)
     * @param seq_len Total number of sequence bases available from bidx
     * @param reverse Whether to reverse index order (right-to-left)
     * @return Minimum edits needed (substitutions + indels) to reach threshold, or NaN if unreachable
     */
    float compute_with_indels(const int* bidx, int seq_len, bool reverse);

    /**
     * Specialized exact solver for max_indels == 1.
     * Enumerates three alignment families (no-indel, one deletion, one insertion)
     * instead of the generic 3D banded DP, for better performance.
     * Produces identical results to compute_with_indels() when max_indels == 1.
     * @param bidx Precomputed base index array (forward or reverse-complemented)
     * @param seq_len Total number of sequence bases available from bidx
     * @param reverse Whether to reverse index order (right-to-left)
     * @return Minimum edits needed to reach threshold, or NaN if unreachable
     */
    float compute_with_one_indel(const int* bidx, int seq_len, bool reverse);

    /**
     * Specialized exact solver for max_indels == 2.
     * Enumerates six alignment families (no-indel, 1 del, 1 ins, 2 dels, 2 ins, 1 del + 1 ins)
     * instead of the generic 3D banded DP, for better performance.
     * Produces identical results to compute_with_indels() when max_indels == 2.
     * @param bidx Precomputed base index array (forward or reverse-complemented)
     * @param seq_len Total number of sequence bases available from bidx
     * @param reverse Whether to reverse index order (right-to-left)
     * @return Minimum edits needed to reach threshold, or NaN if unreachable
     */
    float compute_with_two_indels(const int* bidx, int seq_len, bool reverse);

    /**
     * Banded DP early-abandon filter for indel-enabled windows.
     * Runs a small stack-allocated banded Needleman-Wunsch DP row-by-row,
     * processing all window widths W in [L-D, L+D] simultaneously.
     * After each row, checks whether the best achievable score (current row max
     * plus max suffix score for remaining rows) can possibly reach the threshold.
     * If not, ALL alignment families are provably unreachable — skip the window.
     *
     * This is a PURE OPTIMIZATION — it never changes results. It returns true
     * (skip) ONLY when no alignment can possibly reach the threshold.
     *
     * @param seq_ptr Pointer to start of sequence window
     * @param seq_len Total number of sequence bases available from seq_ptr
     * @param reverse Whether to reverse complement the sequence
     * @return true if the window should be SKIPPED (provably unreachable)
     */
    bool early_abandon_banded_dp(const char* seq_ptr, int seq_len, bool reverse) const;

    /**
     * Pigeonhole pre-filter: checks whether ANY block of the motif can match
     * the sequence window exactly (zero edits) at some shift in {-D, ..., +D}.
     * If no block matches, the window provably cannot match within K total edits.
     * @param bidx Precomputed base index array (forward or reverse-complemented)
     * @param seq_avail Number of sequence bases available from bidx
     * @param reverse Whether reverse strand indexing is used
     * @return true if the window passes the filter (might match); false if safely skippable
     */
    bool passes_prefilter(const int* bidx, int seq_avail, bool reverse) const;

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
