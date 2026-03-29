/**
 * C_gseq_pwm_edits: Returns detailed edit information for PWM edit distance.
 *
 * For each input sequence, finds the optimal window and the specific base
 * changes (edits) needed to reach the threshold. Returns a long-format
 * data frame with one row per edit.
 *
 * Re-implements the O(L) histogram algorithm from PWMEditDistanceScorer
 * as a standalone version that additionally tracks which motif positions
 * and replacement bases are selected.
 */

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif

#include <cstdint>
#include "port.h"

#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>
#include <cctype>
#include <cfloat>

#include <R.h>
#include <Rinternals.h>

#include "DnaPSSM.h"
#include "PWMScorer.h"
#include "PwmCoreParams.h"
#include <set>
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>

namespace {

constexpr float kLogZeroThreshold = std::numeric_limits<float>::lowest() * 0.5f;

struct EditInfo {
    int motif_col;        // 1-based motif column (0 for deletions - will be NA in R)
    char ref_base;        // current/deleted base ('\0' for insertions)
    char alt_base;        // replacement/inserted base ('\0' for deletions)
    float gain;           // score improvement from this edit
    std::string edit_type; // "sub", "ins", "del"
};

struct WindowResult {
    int seq_idx;           // 0-based index into input sequences
    int strand;            // +1 or -1
    int window_start;      // 1-based position within sequence
    float score_before;
    float score_after;
    int n_edits;           // total edits needed (0 if already above threshold, -1 if unreachable)
    std::vector<EditInfo> edits;
    std::string window_seq;   // motif-length sequence at the window (as seen by PSSM)
    std::string mutated_seq;  // window_seq with edits applied
};

inline int base_to_index(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

inline char index_to_base(int idx) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    return (idx >= 0 && idx < 4) ? bases[idx] : 'N';
}

inline char complement_base(char base) {
    switch (base) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default: return 'N';
    }
}

/**
 * Compute the PWM score and detailed edits for a single window.
 *
 * This is the core algorithm: O(L) per window, same as compute_exact
 * in PWMEditDistanceScorer, but additionally returns the specific edits.
 */
WindowResult compute_window_edits_detailed(
    const char* seq_ptr, int L,
    const DnaPSSM& pssm,
    const std::vector<float>& col_max_scores,
    const std::vector<float>& col_min_scores,
    float S_max,
    float S_min,
    float threshold,
    int max_edits,    // -1 = no cap
    bool reverse,
    bool below,       // true = direction "below"
    // Precomputed tables for O(L) algorithm
    const std::vector<float>& gain_values,
    const std::vector<std::vector<uint8_t>>& bin_index)
{
    WindowResult result;
    result.n_edits = -1;  // unreachable by default
    result.score_before = 0.0f;
    result.score_after = 0.0f;

    // Per-position info for edit tracking
    struct PosInfo {
        int motif_col;        // 0-based motif column
        char ref_base;        // current base (after complement if reverse)
        int ref_idx;          // base index (0-3, or 4 for unknown)
        float base_score;     // current score contribution
        float gain;           // gain if switched to target base (positive = toward goal)
        float output_gain;    // gain value for output (positive for above, negative for below)
        int target_base_idx;  // index of target base (best for above, worst for below)
        bool mandatory;       // true if base is unknown/zero-prob (above only)
    };

    std::vector<PosInfo> positions(L);
    int mandatory_edits = 0;
    double true_score = 0.0;      // actual PWM score of the original sequence
    double adjusted_score = 0.0;  // score assuming mandatory edits already applied

    // Pass 1: compute scores and gains for each position
    for (int i = 0; i < L; i++) {
        PosInfo& p = positions[i];
        p.motif_col = i;

        int seq_idx = reverse ? (L - 1 - i) : i;
        char base = seq_ptr[seq_idx];
        if (reverse) {
            base = complement_base(base);
        }
        p.ref_base = base;
        p.ref_idx = base_to_index(base);

        if (below) {
            // "below" direction: target is worst base (min score)
            float worst_score = std::numeric_limits<float>::infinity();
            p.target_base_idx = 0;
            for (int b = 0; b < 4; b++) {
                float s = pssm[i].get_log_prob_from_code(b);
                if (s < worst_score) {
                    worst_score = s;
                    p.target_base_idx = b;
                }
            }

            if (p.ref_idx == 4) {
                // Unknown base: use mean log-probability
                float mean_logp = 0.0f;
                for (int b = 0; b < 4; b++) mean_logp += pssm[i].get_log_prob_from_code(b);
                mean_logp /= 4.0f;
                p.mandatory = false;  // no mandatory edits for "below"
                p.base_score = mean_logp;
                // loss = current - worst (positive means switching helps)
                p.gain = mean_logp - col_min_scores[i];
                p.output_gain = -(p.gain);  // negative in output
                true_score += static_cast<double>(mean_logp);
                adjusted_score += static_cast<double>(mean_logp);
            } else {
                float score = pssm[i].get_log_prob_from_code(p.ref_idx);
                p.mandatory = false;
                p.base_score = score;
                p.gain = score - col_min_scores[i];  // loss from switching to worst
                p.output_gain = -(p.gain);  // negative in output
                true_score += static_cast<double>(score);
                adjusted_score += static_cast<double>(score);
            }
        } else {
            // "above" direction: target is best base (max score)
            float best_score = -std::numeric_limits<float>::infinity();
            p.target_base_idx = 0;
            for (int b = 0; b < 4; b++) {
                float s = pssm[i].get_log_prob_from_code(b);
                if (s > best_score) {
                    best_score = s;
                    p.target_base_idx = b;
                }
            }

            if (p.ref_idx == 4) {
                // Unknown base: mandatory edit
                float mean_logp = 0.0f;
                for (int b = 0; b < 4; b++) mean_logp += pssm[i].get_log_prob_from_code(b);
                mean_logp /= 4.0f;

                p.mandatory = true;
                p.base_score = mean_logp;
                p.gain = col_max_scores[i] - mean_logp;
                p.output_gain = p.gain;
                mandatory_edits++;
                true_score += static_cast<double>(mean_logp);
                adjusted_score += static_cast<double>(col_max_scores[i]);
            } else {
                float score = pssm[i].get_log_prob_from_code(p.ref_idx);
                if (score <= kLogZeroThreshold || !std::isfinite(score)) {
                    p.mandatory = true;
                    p.base_score = score;
                    p.gain = col_max_scores[i] - score;
                    p.output_gain = p.gain;
                    mandatory_edits++;
                    true_score += static_cast<double>(score);
                    adjusted_score += static_cast<double>(col_max_scores[i]);
                } else {
                    p.mandatory = false;
                    p.base_score = score;
                    p.gain = col_max_scores[i] - score;
                    p.output_gain = p.gain;
                    true_score += static_cast<double>(score);
                    adjusted_score += static_cast<double>(score);
                }
            }
        }
    }

    result.score_before = static_cast<float>(true_score);

    // Build window_seq (sequence as seen by PSSM, i.e., reverse-complemented if reverse)
    result.window_seq.resize(L);
    for (int i = 0; i < L; i++) {
        result.window_seq[i] = positions[i].ref_base;
    }

    // For "below": surplus = adjusted_score - threshold (need to lose this much)
    // For "above": deficit = threshold - adjusted_score (need to gain this much)
    double gap = below
        ? (adjusted_score - static_cast<double>(threshold))
        : (static_cast<double>(threshold) - adjusted_score);

    // Already past threshold (after accounting for mandatory edits)?
    if (gap <= 0.0) {
        // Still apply max_edits cap: mandatory-only result must respect the budget.
        if (max_edits >= 0 && mandatory_edits > max_edits) {
            result.n_edits = -1;
            return result;
        }
        result.n_edits = mandatory_edits;
        // score_after = adjusted_score (which has mandatory edits applied)
        result.score_after = static_cast<float>(adjusted_score);
        result.mutated_seq = result.window_seq;
        // Add mandatory edits if any (above direction only)
        for (int i = 0; i < L; i++) {
            if (positions[i].mandatory) {
                EditInfo edit;
                edit.motif_col = i + 1;  // 1-based
                edit.ref_base = positions[i].ref_base;
                edit.alt_base = index_to_base(positions[i].target_base_idx);
                edit.gain = positions[i].output_gain;
                edit.edit_type = "sub";
                result.edits.push_back(edit);
                result.mutated_seq[i] = edit.alt_base;
            }
        }
        return result;
    }

    // Check reachability
    // For "above": max possible improvement = S_max - adjusted_score
    // For "below": max possible loss = adjusted_score - S_min
    double max_possible_delta = below
        ? (adjusted_score - static_cast<double>(S_min))
        : (static_cast<double>(S_max) - adjusted_score);
    if (max_possible_delta < gap) {
        result.n_edits = -1;  // unreachable
        return result;
    }

    // Sort non-mandatory positions by gain (descending) to find optimal edits
    std::vector<int> sorted_positions;
    sorted_positions.reserve(L);
    for (int i = 0; i < L; i++) {
        if (!positions[i].mandatory && positions[i].gain > 0.0f) {
            sorted_positions.push_back(i);
        }
    }
    std::sort(sorted_positions.begin(), sorted_positions.end(),
              [&positions](int a, int b) {
                  return positions[a].gain > positions[b].gain;
              });

    // Greedy: take largest gains until gap is covered
    double acc = 0.0;
    int edits = mandatory_edits;
    std::vector<EditInfo> edit_list;

    // First add mandatory edits (above direction only)
    for (int i = 0; i < L; i++) {
        if (positions[i].mandatory) {
            EditInfo edit;
            edit.motif_col = i + 1;
            edit.ref_base = positions[i].ref_base;
            edit.alt_base = index_to_base(positions[i].target_base_idx);
            edit.gain = positions[i].output_gain;
            edit.edit_type = "sub";
            edit_list.push_back(edit);
        }
    }

    // Then add greedy edits
    for (int idx : sorted_positions) {
        if (max_edits >= 0 && edits >= max_edits) {
            break;  // hit the cap
        }

        acc += static_cast<double>(positions[idx].gain);
        edits++;

        EditInfo edit;
        edit.motif_col = idx + 1;  // 1-based
        edit.ref_base = positions[idx].ref_base;
        edit.alt_base = index_to_base(positions[idx].target_base_idx);
        edit.gain = positions[idx].output_gain;
        edit.edit_type = "sub";
        edit_list.push_back(edit);

        if (acc >= gap) {
            // We've covered the gap
            double score_change = below ? -acc : acc;
            result.n_edits = edits;
            result.score_after = static_cast<float>(adjusted_score + score_change);
            result.edits = edit_list;
            // Build mutated_seq
            result.mutated_seq = result.window_seq;
            for (const auto& e : result.edits) {
                result.mutated_seq[e.motif_col - 1] = e.alt_base;
            }
            return result;
        }
    }

    // Couldn't reach threshold within max_edits
    result.n_edits = -1;
    return result;
}

/**
 * Compute detailed edits for a single window with indel support.
 *
 * Uses banded Needleman-Wunsch DP (same algorithm as
 * PWMEditDistanceScorer::compute_with_indels) but additionally tracks
 * per-edit info (type, position, base, gain).
 *
 * For each window length W in [L-D, L+D], aligns motif[0..L-1] against
 * seq[0..W-1] using a 3D DP table dp[i][j][k], then traces back the
 * alignment to identify specific edits.
 *
 * @param seq_ptr  pointer to the start of the candidate window (at least
 *                 L + max_indels bases available)
 * @param seq_avail number of bases available from seq_ptr
 * @param L        motif length
 * @param pssm     the PSSM
 * @param col_max_scores  per-column max scores
 * @param col_min_scores  per-column min scores
 * @param S_max    sum of col_max_scores
 * @param S_min    sum of col_min_scores
 * @param threshold target score
 * @param max_edits overall edit cap (-1 = no cap)
 * @param max_indels max indels allowed (D)
 * @param reverse  true if scanning reverse strand
 * @param below    true if direction is "below"
 */
WindowResult compute_window_edits_detailed_with_indels(
    const char* seq_ptr, int seq_avail, int L,
    const DnaPSSM& pssm,
    const std::vector<float>& col_max_scores,
    const std::vector<float>& col_min_scores,
    float S_max,
    float S_min,
    float threshold,
    int max_edits,
    int max_indels,
    bool reverse,
    bool below)
{
    const int D = max_indels;

    WindowResult best_result;
    best_result.n_edits = -1;
    best_result.score_before = 0.0f;
    best_result.score_after = 0.0f;

    // Try each window length W in [L-D, L+D]
    for (int W = std::max(1, L - D); W <= L + D; ++W) {
        if (W > seq_avail) break;

        const int rows = L + 1;
        const int cols = W + 1;
        const int indel_levels = D + 1;

        // Flattened 3D DP table
        // For "above" (maximize): init to -inf, pick max
        // For "below" (minimize): init to +inf, pick min
        const double dp_sentinel = below
            ? std::numeric_limits<double>::infinity()
            : -std::numeric_limits<double>::infinity();
        std::vector<double> dp(rows * cols * indel_levels, dp_sentinel);

        auto idx3 = [cols, indel_levels](int i, int j, int k) -> int {
            return i * cols * indel_levels + j * indel_levels + k;
        };

        auto is_valid = [below](double v) -> bool {
            if (below) return v < std::numeric_limits<double>::infinity() * 0.5;
            else return v > -std::numeric_limits<double>::infinity() * 0.5;
        };

        auto is_better = [below](double candidate, double current) -> bool {
            if (below) return candidate < current;
            else return candidate > current;
        };

        // Base case
        dp[idx3(0, 0, 0)] = 0.0;

        // First column: skip motif positions (insertions - each costs 1 indel)
        for (int i = 1; i <= std::min(L, D); ++i) {
            dp[idx3(i, 0, i)] = 0.0;
        }

        // First row: skip sequence positions (deletions - each costs 1 indel)
        for (int j = 1; j <= std::min(W, D); ++j) {
            dp[idx3(0, j, j)] = 0.0;
        }

        // Fill DP
        for (int i = 1; i <= L; ++i) {
            int j_min = std::max(1, i - D);
            int j_max = std::min(W, i + D);

            for (int j = j_min; j <= j_max; ++j) {
                // Get sequence base at position j-1
                int seq_idx = reverse ? (W - 1 - (j - 1)) : (j - 1);
                char base = seq_ptr[seq_idx];
                if (reverse) base = complement_base(base);
                int bidx = base_to_index(base);

                float base_score;
                if (bidx == 4) {
                    float min_s = std::numeric_limits<float>::infinity();
                    for (int b = 0; b < 4; b++) {
                        float s = pssm[i - 1].get_log_prob_from_code(b);
                        if (s < min_s) min_s = s;
                    }
                    base_score = min_s;
                } else {
                    base_score = pssm[i - 1].get_log_prob_from_code(bidx);
                }

                for (int k = 0; k <= D; ++k) {
                    // 1. Match/Substitution (diagonal)
                    if (std::abs((i - 1) - (j - 1)) <= D) {
                        double prev = dp[idx3(i - 1, j - 1, k)];
                        if (is_valid(prev)) {
                            double new_score = prev + static_cast<double>(base_score);
                            if (is_better(new_score, dp[idx3(i, j, k)])) {
                                dp[idx3(i, j, k)] = new_score;
                            }
                        }
                    }

                    if (k < D) {
                        // 2. Insertion: skip motif[i-1], advance motif not sequence
                        if (std::abs((i - 1) - j) <= D) {
                            double prev = dp[idx3(i - 1, j, k)];
                            if (is_valid(prev)) {
                                if (is_better(prev, dp[idx3(i, j, k + 1)])) {
                                    dp[idx3(i, j, k + 1)] = prev;
                                }
                            }
                        }

                        // 3. Deletion: skip seq[j-1], advance sequence not motif
                        if (std::abs(i - (j - 1)) <= D) {
                            double prev = dp[idx3(i, j - 1, k)];
                            if (is_valid(prev)) {
                                if (is_better(prev, dp[idx3(i, j, k + 1)])) {
                                    dp[idx3(i, j, k + 1)] = prev;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Extract results for each indel count k
        for (int k = 0; k <= D; ++k) {
            double score = dp[idx3(L, W, k)];
            if (!is_valid(score)) {
                continue;
            }

            // Traceback to find alignment and compute edits
            // We always traceback to get the full alignment, then determine
            // how many substitutions are needed on top of the k indels.

            // Track alignment operations during traceback
            // Operations: 'M' = match/sub (diagonal), 'I' = insertion (skip motif),
            //             'D' = deletion (skip seq)
            struct AlignOp {
                char op;        // 'M', 'I', 'D'
                int motif_pos;  // 0-based motif position (for M and I)
                int seq_pos;    // 0-based seq position (for M and D)
                float base_score;   // score at this aligned position (M only)
                int base_idx;       // base index (M only)
                char base_char;     // base character (M and D)
            };

            std::vector<AlignOp> alignment;
            alignment.reserve(L + D);

            int ti = L, tj = W, tk = k;
            bool traceback_ok = true;

            while (ti > 0 || tj > 0) {
                if (ti == 0 && tj > 0 && tk > 0) {
                    // Must be deletion
                    AlignOp op;
                    op.op = 'D';
                    op.motif_pos = -1;
                    op.seq_pos = tj - 1;
                    int s_idx = reverse ? (W - 1 - (tj - 1)) : (tj - 1);
                    char b = seq_ptr[s_idx];
                    if (reverse) b = complement_base(b);
                    op.base_char = b;
                    op.base_score = 0.0f;
                    op.base_idx = -1;
                    alignment.push_back(op);
                    tj--; tk--;
                    continue;
                }
                if (tj == 0 && ti > 0 && tk > 0) {
                    // Must be insertion
                    AlignOp op;
                    op.op = 'I';
                    op.motif_pos = ti - 1;
                    op.seq_pos = -1;
                    op.base_char = '\0';
                    op.base_score = 0.0f;
                    op.base_idx = -1;
                    alignment.push_back(op);
                    ti--; tk--;
                    continue;
                }
                if (ti == 0 || tj == 0) break;

                double cur = dp[idx3(ti, tj, tk)];

                // Try diagonal (match/substitution) first
                bool found = false;
                if (std::abs((ti - 1) - (tj - 1)) <= D) {
                    int s_idx = reverse ? (W - 1 - (tj - 1)) : (tj - 1);
                    char b = seq_ptr[s_idx];
                    if (reverse) b = complement_base(b);
                    int bi = base_to_index(b);

                    float bs;
                    if (bi == 4) {
                        float min_s = std::numeric_limits<float>::infinity();
                        for (int bb = 0; bb < 4; bb++) {
                            float s = pssm[ti - 1].get_log_prob_from_code(bb);
                            if (s < min_s) min_s = s;
                        }
                        bs = min_s;
                    } else {
                        bs = pssm[ti - 1].get_log_prob_from_code(bi);
                    }

                    double prev = dp[idx3(ti - 1, tj - 1, tk)];
                    if (is_valid(prev) &&
                        std::fabs((prev + static_cast<double>(bs)) - cur) < 1e-9 * std::max(1.0, std::fabs(cur))) {
                        AlignOp op;
                        op.op = 'M';
                        op.motif_pos = ti - 1;
                        op.seq_pos = tj - 1;
                        op.base_char = b;
                        op.base_score = bs;
                        op.base_idx = bi;
                        alignment.push_back(op);
                        ti--; tj--;
                        found = true;
                    }
                }

                if (!found && ti > 0 && tk > 0 && std::abs((ti - 1) - tj) <= D) {
                    // Try insertion (skip motif position)
                    double prev = dp[idx3(ti - 1, tj, tk - 1)];
                    if (is_valid(prev) &&
                        std::fabs(prev - cur) < 1e-9 * std::max(1.0, std::fabs(cur))) {
                        AlignOp op;
                        op.op = 'I';
                        op.motif_pos = ti - 1;
                        op.seq_pos = -1;
                        op.base_char = '\0';
                        op.base_score = 0.0f;
                        op.base_idx = -1;
                        alignment.push_back(op);
                        ti--; tk--;
                        found = true;
                    }
                }

                if (!found && tj > 0 && tk > 0 && std::abs(ti - (tj - 1)) <= D) {
                    // Try deletion (skip seq position)
                    double prev = dp[idx3(ti, tj - 1, tk - 1)];
                    if (is_valid(prev) &&
                        std::fabs(prev - cur) < 1e-9 * std::max(1.0, std::fabs(cur))) {
                        AlignOp op;
                        op.op = 'D';
                        op.motif_pos = -1;
                        op.seq_pos = tj - 1;
                        int s_idx = reverse ? (W - 1 - (tj - 1)) : (tj - 1);
                        char b = seq_ptr[s_idx];
                        if (reverse) b = complement_base(b);
                        op.base_char = b;
                        op.base_score = 0.0f;
                        op.base_idx = -1;
                        alignment.push_back(op);
                        tj--; tk--;
                        found = true;
                    }
                }

                if (!found) {
                    traceback_ok = false;
                    break;
                }
            }

            if (!traceback_ok) continue;

            // Reverse alignment (it was built backwards)
            std::reverse(alignment.begin(), alignment.end());

            // Collect gains from aligned (diagonal/M) positions for greedy sub selection
            struct AlignedPos {
                int align_idx;        // index into alignment vector
                float gain;           // improvement toward goal (always positive for sorting)
                float output_gain;    // gain for output (positive for above, negative for below)
                int target_base_idx;  // index of target base (best for above, worst for below)
                bool mandatory;       // unknown base or zero-prob (above only)
            };

            std::vector<AlignedPos> aligned_positions;
            aligned_positions.reserve(L);

            for (size_t a = 0; a < alignment.size(); ++a) {
                if (alignment[a].op != 'M') continue;

                int mc = alignment[a].motif_pos;

                if (below) {
                    // "below": target is worst base, gain = loss from switching
                    float worst_s = std::numeric_limits<float>::infinity();
                    int worst_b = 0;
                    for (int b = 0; b < 4; b++) {
                        float s = pssm[mc].get_log_prob_from_code(b);
                        if (s < worst_s) { worst_s = s; worst_b = b; }
                    }
                    float loss = alignment[a].base_score - col_min_scores[mc];

                    AlignedPos ap;
                    ap.align_idx = static_cast<int>(a);
                    ap.gain = loss;  // positive for sorting
                    ap.output_gain = -loss;  // negative in output
                    ap.target_base_idx = worst_b;
                    ap.mandatory = false;
                    aligned_positions.push_back(ap);
                } else {
                    // "above": target is best base
                    float best_s = -std::numeric_limits<float>::infinity();
                    int best_b = 0;
                    for (int b = 0; b < 4; b++) {
                        float s = pssm[mc].get_log_prob_from_code(b);
                        if (s > best_s) { best_s = s; best_b = b; }
                    }
                    float gain = col_max_scores[mc] - alignment[a].base_score;

                    bool mandatory_flag = (alignment[a].base_idx == 4) ||
                        (alignment[a].base_score <= kLogZeroThreshold ||
                         !std::isfinite(alignment[a].base_score));

                    AlignedPos ap;
                    ap.align_idx = static_cast<int>(a);
                    ap.gain = gain;
                    ap.output_gain = gain;
                    ap.target_base_idx = best_b;
                    ap.mandatory = mandatory_flag;
                    aligned_positions.push_back(ap);
                }
            }

            // Determine substitutions needed
            // For "above": account for mandatory subs (zero-prob / unknown)
            // For "below": no mandatory subs
            int mandatory_subs = 0;
            double adjusted = score;
            if (!below) {
                for (auto& ap : aligned_positions) {
                    if (ap.mandatory) {
                        adjusted += static_cast<double>(ap.gain);
                        mandatory_subs++;
                    }
                }
            }

            // Check if already past threshold after mandatory edits
            // For "above": adjusted >= threshold
            // For "below": adjusted <= threshold
            int subs_needed;
            bool already_past = below
                ? (adjusted <= static_cast<double>(threshold))
                : (adjusted >= static_cast<double>(threshold));
            if (already_past) {
                subs_needed = mandatory_subs;
            } else {
                // Need additional greedy subs
                // For "above": remaining gap = threshold - adjusted
                // For "below": remaining gap = adjusted - threshold
                double remaining_gap = below
                    ? (adjusted - static_cast<double>(threshold))
                    : (static_cast<double>(threshold) - adjusted);

                // Sort non-mandatory by gain descending
                std::vector<const AlignedPos*> optional;
                optional.reserve(aligned_positions.size());
                for (auto& ap : aligned_positions) {
                    if (!ap.mandatory && ap.gain > 1e-12f) {
                        optional.push_back(&ap);
                    }
                }
                std::sort(optional.begin(), optional.end(),
                          [](const AlignedPos* a, const AlignedPos* b) {
                              return a->gain > b->gain;
                          });

                double acc = 0.0;
                subs_needed = mandatory_subs;
                bool reachable = false;
                for (auto* ap : optional) {
                    if (max_edits >= 0 && (k + subs_needed) >= max_edits) break;
                    acc += static_cast<double>(ap->gain);
                    subs_needed++;
                    if (acc >= remaining_gap) {
                        reachable = true;
                        break;
                    }
                }
                if (!reachable) continue;
            }

            int total_edits = k + subs_needed;
            if (max_edits >= 0 && total_edits > max_edits) continue;

            // Is this better than our current best?
            if (best_result.n_edits >= 0 && total_edits >= best_result.n_edits) {
                // Not better (we want minimum edits)
                continue;
            }

            // Build detailed edit list and result
            // Mark which aligned positions need substitution
            std::set<int> sub_align_indices;
            // First, mandatory (above only)
            for (auto& ap : aligned_positions) {
                if (ap.mandatory) sub_align_indices.insert(ap.align_idx);
            }
            // Then greedy (re-sort and pick)
            if (subs_needed > mandatory_subs) {
                std::vector<const AlignedPos*> optional;
                for (auto& ap : aligned_positions) {
                    if (!ap.mandatory && ap.gain > 1e-12f) {
                        optional.push_back(&ap);
                    }
                }
                std::sort(optional.begin(), optional.end(),
                          [](const AlignedPos* a, const AlignedPos* b) {
                              return a->gain > b->gain;
                          });
                double acc = 0.0;
                double remaining_gap = below
                    ? (adjusted - static_cast<double>(threshold))
                    : (static_cast<double>(threshold) - adjusted);
                int extra = 0;
                for (auto* ap : optional) {
                    if (extra >= (subs_needed - mandatory_subs)) break;
                    acc += static_cast<double>(ap->gain);
                    sub_align_indices.insert(ap->align_idx);
                    extra++;
                    if (acc >= remaining_gap) break;
                }
            }

            // Build a lookup from align_idx to AlignedPos for target base info
            std::map<int, const AlignedPos*> ap_lookup;
            for (auto& ap : aligned_positions) {
                ap_lookup[ap.align_idx] = &ap;
            }

            // Build window_seq and mutated_seq as an alignment view:
            // Both strings have the same length, with '-' marking gaps.
            // - Deletion (seq base skipped): window_seq has the base, mutated_seq has '-'
            // - Insertion (motif pos has no seq base): window_seq has '-', mutated_seq has the inserted base
            // - Match/substitution: both have a base (same or different)
            std::string window_seq;
            std::string mutated_seq;
            std::vector<EditInfo> edit_list;
            double score_after = score;

            for (size_t a = 0; a < alignment.size(); ++a) {
                const auto& aop = alignment[a];

                if (aop.op == 'M') {
                    bool do_sub = (sub_align_indices.count(static_cast<int>(a)) > 0);

                    // Get target base from AlignedPos lookup
                    auto it = ap_lookup.find(static_cast<int>(a));
                    int target_b = 0;
                    float output_gain_val = 0.0f;
                    if (it != ap_lookup.end()) {
                        target_b = it->second->target_base_idx;
                        output_gain_val = it->second->output_gain;
                    }

                    window_seq += aop.base_char;
                    if (do_sub) {
                        EditInfo ei;
                        ei.motif_col = aop.motif_pos + 1; // 1-based
                        ei.ref_base = aop.base_char;
                        ei.alt_base = index_to_base(target_b);
                        ei.gain = output_gain_val;
                        ei.edit_type = "sub";
                        edit_list.push_back(ei);
                        mutated_seq += index_to_base(target_b);
                        // For "above": score goes up by gain
                        // For "below": score goes down by |gain| (output_gain is negative)
                        score_after += static_cast<double>(output_gain_val);
                    } else {
                        mutated_seq += aop.base_char;
                    }
                } else if (aop.op == 'I') {
                    // Insertion: motif position has no aligned seq base
                    // For "above": insert best base; for "below": insert worst base
                    int target_b;
                    float ins_gain;
                    if (below) {
                        float worst_s = std::numeric_limits<float>::infinity();
                        target_b = 0;
                        for (int b = 0; b < 4; b++) {
                            float s = pssm[aop.motif_pos].get_log_prob_from_code(b);
                            if (s < worst_s) { worst_s = s; target_b = b; }
                        }
                        ins_gain = col_min_scores[aop.motif_pos];  // negative contribution
                    } else {
                        float best_s = -std::numeric_limits<float>::infinity();
                        target_b = 0;
                        for (int b = 0; b < 4; b++) {
                            float s = pssm[aop.motif_pos].get_log_prob_from_code(b);
                            if (s > best_s) { best_s = s; target_b = b; }
                        }
                        ins_gain = col_max_scores[aop.motif_pos];
                    }

                    EditInfo ei;
                    ei.motif_col = aop.motif_pos + 1; // 1-based
                    ei.ref_base = '\0';
                    ei.alt_base = index_to_base(target_b);
                    ei.gain = ins_gain;
                    ei.edit_type = "ins";
                    edit_list.push_back(ei);
                    window_seq += '-';  // gap in original sequence
                    mutated_seq += index_to_base(target_b);
                    score_after += static_cast<double>(ins_gain);
                } else if (aop.op == 'D') {
                    // Deletion: seq base skipped (not aligned to any motif pos)
                    EditInfo ei;
                    ei.motif_col = 0;  // NA in R
                    ei.ref_base = aop.base_char;
                    ei.alt_base = '\0';
                    ei.gain = 0.0f;
                    ei.edit_type = "del";
                    edit_list.push_back(ei);
                    window_seq += aop.base_char;  // base present in original
                    mutated_seq += '-';            // gap: base removed
                }
            }

            // Compute score_before: the PWM score of the L-length window
            // at the same starting position (standard scoring, no alignment).
            // This matches the non-indel path where score_before is the true
            // PWM score of the original sequence.
            double sb = 0.0;
            for (int i = 0; i < L; i++) {
                int si = reverse ? (L - 1 - i) : i;
                char b = (si < seq_avail) ? seq_ptr[si] : 'N';
                if (reverse) b = complement_base(b);
                int bi = base_to_index(b);
                if (bi == 4) {
                    // Unknown: use mean log-probability
                    double mean_lp = 0.0;
                    for (int bb = 0; bb < 4; bb++)
                        mean_lp += static_cast<double>(pssm[i].get_log_prob_from_code(bb));
                    sb += mean_lp / 4.0;
                } else {
                    sb += static_cast<double>(pssm[i].get_log_prob_from_code(bi));
                }
            }
            float score_before_val = static_cast<float>(sb);

            best_result.n_edits = total_edits;
            best_result.score_before = score_before_val;
            best_result.score_after = static_cast<float>(score_after);
            best_result.edits = edit_list;
            best_result.window_seq = window_seq;
            best_result.mutated_seq = mutated_seq;
        }
    }

    return best_result;
}

/**
 * For a given sequence, find the best window and return its edit details.
 * Scans both strands if bidirect, returns the window with minimum edits.
 */
WindowResult find_best_window_edits(
    const std::string& seq,
    int roi_start_0, int roi_end_0,  // 0-based range of allowed window starts
    const DnaPSSM& pssm,
    const std::vector<float>& col_max_scores,
    const std::vector<float>& col_min_scores,
    float S_max,
    float S_min,
    float threshold,
    int max_edits,
    int max_indels,
    bool scan_forward, bool scan_reverse,
    float score_min, float score_max,
    bool below,
    const std::vector<float>& gain_values,
    const std::vector<std::vector<uint8_t>>& bin_index)
{
    const int L = pssm.length();
    const bool has_score_min = !std::isnan(score_min);
    const bool has_score_max = !std::isnan(score_max);
    // Disable score filter when indels are enabled: with indels, the DP alignment
    // can use L±D bases, so the L-length pre-score is not a valid pre-filter.
    // This matches PWMEditDistanceScorer::compute_interval behavior.
    const bool has_score_filter = (has_score_min || has_score_max) && (max_indels == 0);
    const int seqlen = static_cast<int>(seq.length());

    WindowResult best;
    best.n_edits = -1;
    best.strand = 0;
    best.window_start = 0;

    for (int s0 = roi_start_0; s0 <= roi_end_0; ++s0) {
        const char* window_start = seq.data() + s0;

        auto try_window = [&](bool reverse, int direction) {
            int seq_avail = seqlen - s0;
            // Score filter: compute L-length window score for filtering
            if (has_score_filter && seq_avail >= L) {
                float logp = 0.0f;
                for (int i = 0; i < L; i++) {
                    int si = reverse ? (L - 1 - i) : i;
                    char base = window_start[si];
                    if (reverse) base = complement_base(base);
                    int bidx = base_to_index(base);
                    if (bidx == 4) {
                        float sum = 0.0f;
                        for (int b = 0; b < 4; b++) sum += pssm[i].get_log_prob_from_code(b);
                        logp += sum / 4.0f;
                    } else {
                        logp += pssm[i].get_log_prob_from_code(bidx);
                    }
                }
                if (has_score_min && logp < score_min) return;
                if (has_score_max && logp > score_max) return;
            }

            WindowResult wr;
            if (max_indels > 0) {
                wr = compute_window_edits_detailed_with_indels(
                    window_start, seq_avail, L, pssm, col_max_scores, col_min_scores,
                    S_max, S_min, threshold, max_edits, max_indels, reverse, below);
            } else {
                wr = compute_window_edits_detailed(
                    window_start, L, pssm, col_max_scores, col_min_scores,
                    S_max, S_min, threshold, max_edits, reverse, below,
                    gain_values, bin_index);
            }

            if (wr.n_edits < 0) return;  // unreachable

            if (best.n_edits < 0 || wr.n_edits < best.n_edits ||
                (wr.n_edits == best.n_edits && wr.score_before > best.score_before)) {
                best = wr;
                best.strand = direction;
                best.window_start = s0 + 1;  // 1-based
            }
        };

        if (scan_forward) try_window(false, +1);
        if (scan_reverse) try_window(true, -1);
    }

    return best;
}

} // anonymous namespace


extern "C" {

/**
 * C_gseq_pwm_edits: R-callable function returning edit details.
 *
 * Parameters:
 *   r_seqs        - character vector of sequences
 *   r_pssm        - numeric matrix (columns A,C,G,T)
 *   r_score_thresh - target score threshold
 *   r_max_edits   - max edits (NULL or positive integer)
 *   r_max_indels  - max indels (NULL or non-negative integer; 0 = subs only)
 *   r_bidirect    - logical, scan both strands?
 *   r_strand_mode - integer (-1, 0, 1)
 *   r_prior       - numeric prior (pseudocount)
 *   r_roi_start   - integer vector of ROI starts (1-based), or NULL
 *   r_roi_end     - integer vector of ROI ends (1-based), or NULL
 *   r_extend      - logical or integer
 *   r_score_min   - numeric or NULL
 *   r_score_max   - numeric or NULL
 *
 * Returns: data.frame with columns:
 *   seq_idx, strand, window_start, score_before, score_after, n_edits,
 *   edit_num, motif_col, ref, alt, gain, edit_type, window_seq, mutated_seq
 */
SEXP C_gseq_pwm_edits(SEXP r_seqs, SEXP r_pssm, SEXP r_score_thresh,
                       SEXP r_max_edits, SEXP r_max_indels,
                       SEXP r_bidirect, SEXP r_strand_mode,
                       SEXP r_prior, SEXP r_roi_start, SEXP r_roi_end,
                       SEXP r_extend, SEXP r_score_min, SEXP r_score_max,
                       SEXP r_direction)
{
    try {
        if (!Rf_isString(r_seqs)) Rf_error("seqs must be a character vector");
        if (!Rf_isMatrix(r_pssm) || !Rf_isReal(r_pssm)) Rf_error("pssm must be a numeric matrix");

        int n_seqs = Rf_length(r_seqs);
        if (n_seqs == 0) {
            // Return empty data frame
            const int ncols = 14;
            SEXP result = PROTECT(Rf_allocVector(VECSXP, ncols));
            SEXP names = PROTECT(Rf_allocVector(STRSXP, ncols));
            const char* colnames[] = {"seq_idx", "strand", "window_start",
                                       "score_before", "score_after", "n_edits",
                                       "edit_num", "motif_col", "ref", "alt", "gain",
                                       "edit_type", "window_seq", "mutated_seq"};
            for (int i = 0; i < ncols; i++) {
                SET_STRING_ELT(names, i, Rf_mkChar(colnames[i]));
                SEXPTYPE tp = (i < 3 || i == 5 || i == 6 || i == 7) ? INTSXP :
                              (i == 8 || i == 9 || i == 11 || i == 12 || i == 13) ? STRSXP : REALSXP;
                SET_VECTOR_ELT(result, i, Rf_allocVector(tp, 0));
            }
            Rf_setAttrib(result, R_NamesSymbol, names);

            SEXP row_names = PROTECT(Rf_allocVector(INTSXP, 2));
            INTEGER(row_names)[0] = NA_INTEGER;
            INTEGER(row_names)[1] = 0;
            Rf_setAttrib(result, R_RowNamesSymbol, row_names);

            SEXP cls = PROTECT(Rf_mkString("data.frame"));
            Rf_setAttrib(result, R_ClassSymbol, cls);
            UNPROTECT(4);
            return result;
        }

        // Build PSSM
        PwmCoreParams core;
        core.pssm = PWMScorer::create_pssm_from_matrix(r_pssm);
        int w = core.pssm.size();
        core.prior = Rf_asReal(r_prior);
        core.apply_prior();

        float threshold = static_cast<float>(Rf_asReal(r_score_thresh));
        bool bidirect = Rf_asLogical(r_bidirect);
        int strand_mode = Rf_asInteger(r_strand_mode);
        if (bidirect) strand_mode = 0;

        bool scan_forward = bidirect || strand_mode >= 0;
        bool scan_reverse = bidirect || strand_mode <= 0;

        int max_edits = -1;
        if (!Rf_isNull(r_max_edits)) {
            max_edits = Rf_asInteger(r_max_edits);
        }

        int max_indels = 0;
        if (!Rf_isNull(r_max_indels)) {
            max_indels = Rf_asInteger(r_max_indels);
            if (max_indels < 0) max_indels = 0;
        }

        int extend_val;
        if (Rf_isLogical(r_extend)) {
            extend_val = Rf_asLogical(r_extend) ? (w - 1) : 0;
        } else {
            extend_val = Rf_asInteger(r_extend);
            if (extend_val < 0) extend_val = 0;
        }

        float score_min = std::numeric_limits<float>::quiet_NaN();
        if (!Rf_isNull(r_score_min)) {
            score_min = static_cast<float>(Rf_asReal(r_score_min));
        }
        float score_max_val = std::numeric_limits<float>::quiet_NaN();
        if (!Rf_isNull(r_score_max)) {
            score_max_val = static_cast<float>(Rf_asReal(r_score_max));
        }

        // Parse direction
        bool below = false;
        if (!Rf_isNull(r_direction) && Rf_isString(r_direction) && Rf_length(r_direction) > 0) {
            std::string dir_str = CHAR(STRING_ELT(r_direction, 0));
            below = (dir_str == "below");
        }

        // Precompute tables (same as PWMEditDistanceScorer::precompute_tables)
        std::vector<float> col_max_scores(w);
        std::vector<float> col_min_scores(w);
        float S_max = 0.0f;
        float S_min = 0.0f;
        for (int i = 0; i < w; i++) {
            float max_score = -std::numeric_limits<float>::infinity();
            float min_score = std::numeric_limits<float>::infinity();
            for (int b = 0; b < 4; b++) {
                float s = core.pssm[i].get_log_prob_from_code(b);
                if (s > max_score) max_score = s;
                if (s < min_score && std::isfinite(s)) min_score = s;
            }
            col_max_scores[i] = max_score;
            col_min_scores[i] = min_score;
            S_max += max_score;
            S_min += min_score;
        }

        // Gain/loss values sorted descending (direction-dependent)
        // ABOVE: delta[i][b] = col_max[i] - score[i][b]  (gain from switching to best)
        // BELOW: delta[i][b] = score[i][b] - col_min[i]   (loss from switching to worst)
        std::set<float, std::greater<float>> unique_gains;
        std::vector<std::vector<uint8_t>> bin_idx(w);
        for (int i = 0; i < w; i++) {
            bin_idx[i].resize(5);
            for (int b = 0; b < 4; b++) {
                float s = core.pssm[i].get_log_prob_from_code(b);
                float g = below ? (s - col_min_scores[i]) : (col_max_scores[i] - s);
                if (!std::isfinite(g) || g < 0.0f) g = 0.0f;
                unique_gains.insert(g);
            }
            // Unknown base (index 4): use max delta
            float max_delta = below ? (col_max_scores[i] - col_min_scores[i])
                                    : (col_max_scores[i] - col_min_scores[i]);
            if (!std::isfinite(max_delta) || max_delta < 0.0f) max_delta = 0.0f;
            unique_gains.insert(max_delta);
        }
        std::vector<float> gain_values(unique_gains.begin(), unique_gains.end());
        for (int i = 0; i < w; i++) {
            for (int b = 0; b < 5; b++) {
                float score;
                if (b < 4) {
                    score = core.pssm[i].get_log_prob_from_code(b);
                } else {
                    // Unknown base: use the worst-case score for this direction
                    if (below) {
                        score = -std::numeric_limits<float>::infinity();
                        for (int bb = 0; bb < 4; bb++) {
                            float s = core.pssm[i].get_log_prob_from_code(bb);
                            if (s > score) score = s;
                        }
                    } else {
                        score = std::numeric_limits<float>::infinity();
                        for (int bb = 0; bb < 4; bb++) {
                            float s = core.pssm[i].get_log_prob_from_code(bb);
                            if (s < score) score = s;
                        }
                    }
                }
                float g = below ? (score - col_min_scores[i]) : (col_max_scores[i] - score);
                if (!std::isfinite(g) || g < 0.0f) g = 0.0f;
                auto it = std::lower_bound(gain_values.begin(), gain_values.end(),
                                           g, std::greater<float>());
                bin_idx[i][b] = static_cast<uint8_t>(std::distance(gain_values.begin(), it));
            }
        }

        // Process each sequence and collect all edit rows
        struct EditRow {
            int seq_idx;
            int strand;
            int window_start;
            float score_before;
            float score_after;
            int n_edits;
            int edit_num;
            int motif_col;
            char ref;
            char alt;
            float gain;
            std::string edit_type;
            std::string window_seq;
            std::string mutated_seq;
        };

        std::vector<EditRow> all_rows;

        for (int si = 0; si < n_seqs; si++) {
            SEXP r_str = STRING_ELT(r_seqs, si);
            if (r_str == NA_STRING) continue;

            std::string seq = CHAR(r_str);
            // Convert to uppercase
            for (char& c : seq) c = toupper(c);

            int seqlen = static_cast<int>(seq.length());
            // Minimum sequence length: with indels, shortest window is max(1, w-D)
            int min_seq_len = (max_indels > 0) ? std::max(1, w - max_indels) : w;
            if (seqlen < min_seq_len) continue;

            // Compute ROI bounds
            int roi_start_1 = Rf_isNull(r_roi_start) ? 1 :
                (Rf_length(r_roi_start) == 1 ? INTEGER(r_roi_start)[0] : INTEGER(r_roi_start)[si]);
            int roi_end_1 = Rf_isNull(r_roi_end) ? seqlen :
                (Rf_length(r_roi_end) == 1 ? INTEGER(r_roi_end)[0] : INTEGER(r_roi_end)[si]);

            // Compute allowed window starts (0-based)
            // With indels, the DP tries window lengths W in [max(1,L-D), L+D].
            // The minimum bases needed from start is max(1, L-D) for the
            // shortest window. The DP handles longer windows up to seq_avail
            // via its own bounds check.
            int min_window = std::max(1, w - max_indels);
            int start_min0 = std::max(0, roi_start_1 - 1 - extend_val);
            int start_max0 = std::min(seqlen - min_window, roi_end_1 - w + extend_val);
            // For the non-indel path, still need at least L bases
            if (max_indels == 0) {
                start_max0 = std::min(start_max0, seqlen - w);
            }

            if (start_min0 > start_max0) continue;

            WindowResult wr = find_best_window_edits(
                seq, start_min0, start_max0,
                core.pssm, col_max_scores, col_min_scores, S_max, S_min,
                threshold, max_edits, max_indels,
                scan_forward, scan_reverse,
                score_min, score_max_val, below,
                gain_values, bin_idx);

            if (wr.n_edits == 0) {
                // Already above threshold — single row with edit_num=0
                EditRow row;
                row.seq_idx = si + 1;  // 1-based
                row.strand = wr.strand;
                row.window_start = wr.window_start;
                row.score_before = wr.score_before;
                row.score_after = wr.score_after;
                row.n_edits = 0;
                row.edit_num = 0;
                row.motif_col = NA_INTEGER;
                row.ref = 0;
                row.alt = 0;
                row.gain = 0.0f;
                row.edit_type = "";  // will be NA_STRING in R
                row.window_seq = wr.window_seq;
                row.mutated_seq = wr.mutated_seq;
                all_rows.push_back(row);
            } else if (wr.n_edits > 0 && !wr.edits.empty()) {
                for (size_t e = 0; e < wr.edits.size(); e++) {
                    EditRow row;
                    row.seq_idx = si + 1;
                    row.strand = wr.strand;
                    row.window_start = wr.window_start;
                    row.score_before = wr.score_before;
                    row.score_after = wr.score_after;
                    row.n_edits = wr.n_edits;
                    row.edit_num = static_cast<int>(e) + 1;
                    row.motif_col = wr.edits[e].motif_col;
                    row.ref = wr.edits[e].ref_base;
                    row.alt = wr.edits[e].alt_base;
                    row.gain = wr.edits[e].gain;
                    row.edit_type = wr.edits[e].edit_type;
                    row.window_seq = wr.window_seq;
                    row.mutated_seq = wr.mutated_seq;
                    all_rows.push_back(row);
                }
            }
            // If n_edits < 0 (unreachable), skip this sequence — no rows emitted
        }

        // Build R data frame
        int n_rows = static_cast<int>(all_rows.size());

        SEXP r_seq_idx = PROTECT(Rf_allocVector(INTSXP, n_rows));
        SEXP r_strand_out = PROTECT(Rf_allocVector(INTSXP, n_rows));
        SEXP r_wstart = PROTECT(Rf_allocVector(INTSXP, n_rows));
        SEXP r_sbefore = PROTECT(Rf_allocVector(REALSXP, n_rows));
        SEXP r_safter = PROTECT(Rf_allocVector(REALSXP, n_rows));
        SEXP r_nedits = PROTECT(Rf_allocVector(INTSXP, n_rows));
        SEXP r_editnum = PROTECT(Rf_allocVector(INTSXP, n_rows));
        SEXP r_mcol = PROTECT(Rf_allocVector(INTSXP, n_rows));
        SEXP r_ref = PROTECT(Rf_allocVector(STRSXP, n_rows));
        SEXP r_alt = PROTECT(Rf_allocVector(STRSXP, n_rows));
        SEXP r_gain = PROTECT(Rf_allocVector(REALSXP, n_rows));
        SEXP r_etype = PROTECT(Rf_allocVector(STRSXP, n_rows));
        SEXP r_wseq = PROTECT(Rf_allocVector(STRSXP, n_rows));
        SEXP r_mseq = PROTECT(Rf_allocVector(STRSXP, n_rows));

        for (int i = 0; i < n_rows; i++) {
            const EditRow& row = all_rows[i];
            INTEGER(r_seq_idx)[i] = row.seq_idx;
            INTEGER(r_strand_out)[i] = row.strand;
            INTEGER(r_wstart)[i] = row.window_start;
            REAL(r_sbefore)[i] = row.score_before;
            REAL(r_safter)[i] = row.score_after;
            INTEGER(r_nedits)[i] = row.n_edits;
            INTEGER(r_editnum)[i] = row.edit_num;
            // motif_col: 0 means NA (deletions have no motif position)
            INTEGER(r_mcol)[i] = (row.motif_col == 0) ? NA_INTEGER : row.motif_col;

            // ref and alt: handle '\0' as NA (insertions/deletions/no-edit rows)
            if (row.ref == 0) {
                SET_STRING_ELT(r_ref, i, NA_STRING);
            } else {
                char buf[2] = {row.ref, '\0'};
                SET_STRING_ELT(r_ref, i, Rf_mkChar(buf));
            }
            if (row.alt == 0) {
                SET_STRING_ELT(r_alt, i, NA_STRING);
            } else {
                char buf[2] = {row.alt, '\0'};
                SET_STRING_ELT(r_alt, i, Rf_mkChar(buf));
            }

            REAL(r_gain)[i] = row.gain;

            // edit_type: empty string means NA (n_edits==0 rows)
            if (row.edit_type.empty()) {
                SET_STRING_ELT(r_etype, i, NA_STRING);
            } else {
                SET_STRING_ELT(r_etype, i, Rf_mkChar(row.edit_type.c_str()));
            }

            SET_STRING_ELT(r_wseq, i, Rf_mkChar(row.window_seq.c_str()));
            SET_STRING_ELT(r_mseq, i, Rf_mkChar(row.mutated_seq.c_str()));
        }

        // Assemble data frame
        const int ncols = 14;
        SEXP result_df = PROTECT(Rf_allocVector(VECSXP, ncols));
        SET_VECTOR_ELT(result_df, 0, r_seq_idx);
        SET_VECTOR_ELT(result_df, 1, r_strand_out);
        SET_VECTOR_ELT(result_df, 2, r_wstart);
        SET_VECTOR_ELT(result_df, 3, r_sbefore);
        SET_VECTOR_ELT(result_df, 4, r_safter);
        SET_VECTOR_ELT(result_df, 5, r_nedits);
        SET_VECTOR_ELT(result_df, 6, r_editnum);
        SET_VECTOR_ELT(result_df, 7, r_mcol);
        SET_VECTOR_ELT(result_df, 8, r_ref);
        SET_VECTOR_ELT(result_df, 9, r_alt);
        SET_VECTOR_ELT(result_df, 10, r_gain);
        SET_VECTOR_ELT(result_df, 11, r_etype);
        SET_VECTOR_ELT(result_df, 12, r_wseq);
        SET_VECTOR_ELT(result_df, 13, r_mseq);

        SEXP names = PROTECT(Rf_allocVector(STRSXP, ncols));
        const char* colnames[] = {"seq_idx", "strand", "window_start",
                                   "score_before", "score_after", "n_edits",
                                   "edit_num", "motif_col", "ref", "alt", "gain",
                                   "edit_type", "window_seq", "mutated_seq"};
        for (int i = 0; i < ncols; i++) {
            SET_STRING_ELT(names, i, Rf_mkChar(colnames[i]));
        }
        Rf_setAttrib(result_df, R_NamesSymbol, names);

        SEXP row_names = PROTECT(Rf_allocVector(INTSXP, 2));
        INTEGER(row_names)[0] = NA_INTEGER;
        INTEGER(row_names)[1] = n_rows;
        Rf_setAttrib(result_df, R_RowNamesSymbol, row_names);

        SEXP cls = PROTECT(Rf_mkString("data.frame"));
        Rf_setAttrib(result_df, R_ClassSymbol, cls);

        UNPROTECT(18);
        return result_df;

    } catch (std::exception& e) {
        Rf_error("C_gseq_pwm_edits: %s", e.what());
    } catch (...) {
        Rf_error("C_gseq_pwm_edits: unknown error");
    }

    return R_NilValue;
}

} // extern "C"
