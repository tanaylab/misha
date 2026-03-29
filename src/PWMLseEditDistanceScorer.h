#ifndef PWM_LSE_EDIT_DISTANCE_SCORER_H_
#define PWM_LSE_EDIT_DISTANCE_SCORER_H_

#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <queue>
#include "GenomeSeqScorer.h"
#include "DnaPSSM.h"
#include "PWMEditDistanceScorer.h"

/**
 * PWMLseEditDistanceScorer: Computes LSE-mode edit distance for PWM.
 *
 * Given an interval with sequence, the LSE score is F = log(sum(exp(S_p)))
 * where S_p is the PWM score at start position p. This scorer finds the
 * minimum number of base edits (substitutions) to bring F above or below
 * score.thresh, depending on the direction parameter.
 *
 * Algorithm:
 *  - Exhaustive search for k <= 2 (O(RL^2), provably optimal)
 *  - Greedy heuristic with partition function tracking for k >= 3
 *
 * Supported return modes:
 *  - LSE_EDIT_DISTANCE:     minimum edits to reach the threshold
 *  - LSE_EDIT_DISTANCE_POS: position of the most impactful single edit
 */
class PWMLseEditDistanceScorer : public GenomeSeqScorer
{
public:
    enum class Mode {
        LSE_EDIT_DISTANCE,
        LSE_EDIT_DISTANCE_POS
    };

    using Direction = PWMEditDistanceScorer::Direction;

    /**
     * Constructor
     * @param pssm The Position-Specific Scoring Matrix
     * @param shared_seqfetch Shared sequence fetcher for caching
     * @param threshold LSE score threshold to reach
     * @param max_edits Maximum edits to consider (-1 = unlimited)
     * @param extend Allow motif to extend beyond interval boundaries
     * @param strand Strand mode (0=both, 1=forward, -1=reverse)
     * @param mode Return mode
     * @param score_min Minimum LSE score filter (NaN = no filter)
     * @param score_max Maximum LSE score filter (NaN = no filter)
     * @param direction ABOVE = min edits to bring F >= threshold;
     *                  BELOW = min edits to bring F <= threshold
     */
    PWMLseEditDistanceScorer(const DnaPSSM& pssm,
                             GenomeSeqFetch* shared_seqfetch,
                             float threshold,
                             int max_edits = -1,
                             bool extend = true,
                             char strand = 0,
                             Mode mode = Mode::LSE_EDIT_DISTANCE,
                             float score_min = std::numeric_limits<float>::quiet_NaN(),
                             float score_max = std::numeric_limits<float>::quiet_NaN(),
                             Direction direction = Direction::ABOVE);

    /**
     * Score a genomic interval - returns minimum edits or position
     */
    float score_interval(const GInterval& interval, const GenomeChromKey& chromkey) override;

    float get_last_min_edits() const { return m_last_min_edits; }

private:
    DnaPSSM m_pssm;
    float m_threshold;
    int m_max_edits;
    Mode m_mode;
    Direction m_direction;
    float m_score_min;
    float m_score_max;
    float m_last_min_edits;

    // Precomputed tables
    std::vector<float> m_col_max_scores;  // max score per column
    std::vector<float> m_col_min_scores;  // min score per column (for BELOW direction)
    float m_S_max;                        // sum of column maxima
    float m_S_min;                        // sum of column minima (for BELOW direction)

    /**
     * Precompute column maxima (called once in constructor)
     */
    void precompute_tables();

    /**
     * Convert base character to index (A=0, C=1, G=2, T=3, unknown=4)
     */
    inline int base_to_index(char base) const {
        switch (base) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 4;
        }
    }

    inline char complement_base(char base) const {
        switch (base) {
            case 'A': case 'a': return 'T';
            case 'C': case 'c': return 'G';
            case 'G': case 'g': return 'C';
            case 'T': case 't': return 'A';
            default: return 'N';
        }
    }

    inline bool should_scan_forward() const {
        return m_pssm.is_bidirect() || m_strand >= 0;
    }

    inline bool should_scan_reverse() const {
        return m_pssm.is_bidirect() || m_strand <= 0;
    }

    /**
     * Encode position (1-based, signed by strand when bidirect)
     */
    float encode_position(size_t index, size_t target_length, size_t motif_length, int direction) const;

    /**
     * Compute per-window PWM score for a given sequence and direction
     * @param seq The full sequence
     * @param reverse Whether to reverse-complement
     * @return PWM score for this window
     */
    float compute_window_score(const char* window_start, bool reverse) const;

    struct LseResult {
        float min_edits = std::numeric_limits<float>::quiet_NaN();
        float best_pos = std::numeric_limits<float>::quiet_NaN();
        int best_edit_index = -1;    // genomic index of best single edit
        int best_edit_direction = 1; // direction for the best single edit position
    };

    /**
     * Core LSE edit distance algorithm for one strand direction.
     * @param seq The full sequence string
     * @param interval_length Length of the original interval
     * @param reverse Whether to reverse-complement
     * @return LseResult with min_edits and best position
     */
    LseResult compute_lse_edit_distance(const std::string& seq,
                                         size_t interval_length,
                                         bool reverse);

    /**
     * Initialize per-start scores and partition function
     * @param seq_bases Resolved sequence bases (already reverse-complemented if needed)
     * @param N Number of starts
     * @param[out] S_p Per-start scores
     * @param[out] A_p Shifted weights exp(S_p - m)
     * @param[out] Z Partition function sum
     * @param[out] m Shift (max S_p)
     * @param[out] F LSE score log(Z) + m
     */
    void init_partition(const std::vector<int>& seq_bases,
                        int N, int L,
                        std::vector<double>& S_p,
                        std::vector<double>& A_p,
                        double& Z, double& m, double& F) const;

    /**
     * Compute DeltaZ for a single edit at position q -> base y
     * @param q Genomic position (0-based in resolved sequence)
     * @param y_idx Target base index (0-3)
     * @param seq_bases Current resolved sequence bases
     * @param A_p Current shifted weights
     * @param N Number of starts
     * @param L Motif length
     * @return Change in Z from this edit
     */
    double compute_delta_z(int q, int y_idx,
                           const std::vector<int>& seq_bases,
                           const std::vector<double>& A_p,
                           int N, int L) const;

    /**
     * Apply an edit: update seq_bases, A_p, Z, S_p
     */
    void apply_edit(int q, int y_idx,
                    std::vector<int>& seq_bases,
                    std::vector<double>& S_p,
                    std::vector<double>& A_p,
                    double& Z, double& m,
                    int N, int L);

    /**
     * Try exhaustive k=1: check all single edits
     * @return number of edits (1) if found, or NaN if no single edit suffices
     */
    LseResult try_k1_exhaustive(const std::vector<int>& seq_bases,
                                 const std::vector<double>& A_p,
                                 double Z, double m, double F,
                                 int N, int L,
                                 size_t target_length, size_t motif_length,
                                 int direction) const;

    /**
     * Try exhaustive k=2: check all pairs of edits
     */
    LseResult try_k2_exhaustive(const std::vector<int>& seq_bases,
                                 const std::vector<double>& S_p,
                                 const std::vector<double>& A_p,
                                 double Z, double m, double F,
                                 int N, int L, int R,
                                 size_t target_length, size_t motif_length,
                                 int direction) const;

    /**
     * Greedy heuristic for k >= 3
     */
    LseResult greedy_search(std::vector<int>& seq_bases,
                            std::vector<double>& S_p,
                            std::vector<double>& A_p,
                            double& Z, double& m, double& F,
                            int N, int L, int R,
                            int first_edit_pos, int first_edit_base,
                            int second_edit_pos, int second_edit_base,
                            int edits_so_far,
                            size_t target_length, size_t motif_length,
                            int direction);
};

#endif // PWM_LSE_EDIT_DISTANCE_SCORER_H_
