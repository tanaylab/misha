#include "PWMEditDistanceScorer.h"
#include "GenomeSeqFetch.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace {
constexpr float kLogZeroThreshold = std::numeric_limits<float>::lowest() * 0.5f;
}

PWMEditDistanceScorer::PWMEditDistanceScorer(const DnaPSSM& pssm,
                                             GenomeSeqFetch* shared_seqfetch,
                                             float threshold,
                                             int max_edits,
                                             bool extend,
                                             char strand,
                                             Mode mode,
                                             float score_min,
                                             int max_indels,
                                             float score_max)
    : GenomeSeqScorer(shared_seqfetch, extend, strand),
      m_pssm(pssm),
      m_threshold(threshold),
      m_max_edits(max_edits),
      m_max_indels(max_indels),
      m_mode(mode),
      m_score_min(score_min),
      m_score_max(score_max),
      m_S_max(0.0f)
{
    precompute_tables();
}

void PWMEditDistanceScorer::precompute_tables()
{
    int L = m_pssm.length();
    m_col_max_scores.resize(L);

    // Compute column maxima and S_max
    m_S_max = 0.0f;
    for (int i = 0; i < L; i++) {
        float max_score = -std::numeric_limits<float>::infinity();

        // Find max score across all bases for this column
        for (int b = 0; b < 4; b++) {
            float score = m_pssm[i].get_log_prob_from_code(b);
            if (score > max_score) {
                max_score = score;
            }
        }

        m_col_max_scores[i] = max_score;
        m_S_max += max_score;
    }

    // Collect all gain values
    std::set<float, std::greater<float>> unique_gains;  // Sorted descending
    m_bin_index.resize(L);

    for (int i = 0; i < L; i++) {
        m_bin_index[i].resize(5);  // 5 to handle unknown bases (index 4)

        for (int b = 0; b < 4; b++) {
            float score = m_pssm[i].get_log_prob_from_code(b);
            float gain = m_col_max_scores[i] - score;
            unique_gains.insert(gain);
        }

        // BUG-1b FIX: Handle unknown bases (N, etc.) - find minimum score
        // (not maximum) to compute worst-case (largest) gain
        float min_score = std::numeric_limits<float>::infinity();
        for (int b = 0; b < 4; b++) {
            float score = m_pssm[i].get_log_prob_from_code(b);
            if (score < min_score) {
                min_score = score;
            }
        }
        float gain_unknown = m_col_max_scores[i] - min_score;
        unique_gains.insert(gain_unknown);
    }

    // Copy to vector (already sorted descending)
    m_gain_values.assign(unique_gains.begin(), unique_gains.end());

    // PERF-1: Initialize reusable count vector
    m_exact_count.resize(m_gain_values.size(), 0);
    m_exact_touched.reserve(m_gain_values.size());

    // Build bin lookup table
    for (int i = 0; i < L; i++) {
        for (int b = 0; b < 5; b++) {  // Include unknown base index 4
            float score;
            if (b < 4) {
                score = m_pssm[i].get_log_prob_from_code(b);
            } else {
                // BUG-1b FIX: Unknown base - use minimum score (worst case)
                score = std::numeric_limits<float>::infinity();
                for (int bb = 0; bb < 4; bb++) {
                    float s = m_pssm[i].get_log_prob_from_code(bb);
                    if (s < score) {
                        score = s;
                    }
                }
            }

            float gain = m_col_max_scores[i] - score;

            // Find bin index (m_gain_values is sorted descending)
            auto it = std::lower_bound(m_gain_values.begin(), m_gain_values.end(),
                                      gain, std::greater<float>());
            m_bin_index[i][b] = static_cast<uint8_t>(std::distance(m_gain_values.begin(), it));
        }
    }
}

float PWMEditDistanceScorer::score_interval(const GInterval& interval,
                                            const GenomeChromKey& chromkey)
{
    m_last_metrics = ScanMetrics();

    const int motif_len = m_pssm.length();
    if (motif_len <= 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    // When indels are enabled, deletions from the motif consume extra sequence bases,
    // so we need to fetch max_indels additional bases beyond the normal expansion.
    GInterval exp_interval = calculate_expanded_interval(interval, chromkey, motif_len + m_max_indels);

    std::vector<char> seq_vec;
    try {
        m_seqfetch_ptr->read_interval(exp_interval, chromkey, seq_vec);
    } catch (...) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    std::string seq(seq_vec.begin(), seq_vec.end());
    if (seq.length() < static_cast<size_t>(motif_len)) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    int64_t interval_len = interval.end - interval.start;
    if (interval_len <= 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    ScanMetrics metrics = evaluate_windows(seq, static_cast<size_t>(interval_len));
    m_last_metrics = metrics;

    switch (m_mode) {
        case Mode::MIN_EDITS:
            return metrics.min_edits;
        case Mode::MIN_EDITS_POSITION:
            return metrics.min_edits_position;
        case Mode::PWM_MAX_EDITS:
            return metrics.best_pwm_edits;
        default:
            return std::numeric_limits<float>::quiet_NaN();
    }
}

float PWMEditDistanceScorer::compute_window_pwm_score(const char* window_start, bool reverse)
{
    const int L = m_pssm.length();
    float logp = 0.0f;

    for (int i = 0; i < L; i++) {
        int seq_idx = reverse ? (L - 1 - i) : i;
        char base = window_start[seq_idx];
        if (reverse) {
            base = complement_base(base);
        }

        int bidx = base_to_index(base);
        // BUG-1 FIX: check for unknown bases before calling get_log_prob_from_code
        if (bidx == 4) {
            // Unknown base: use mean log-probability (same as PWMScorer convention)
            float sum = 0.0f;
            for (int b = 0; b < 4; b++) {
                sum += m_pssm[i].get_log_prob_from_code(b);
            }
            logp += sum / 4.0f;
        } else {
            logp += m_pssm[i].get_log_prob_from_code(bidx);
        }
    }

    return logp;
}

float PWMEditDistanceScorer::compute_exact(const char* seq_ptr, bool reverse)
{
    const int L = m_pssm.length();
    if (L == 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    int mandatory_edits = 0;
    double adjusted_score = 0.0;

    // PERF-1: Use class-member count vector with touched-list cleanup
    // (m_exact_count is already zeroed from previous cleanup)
    m_exact_touched.clear();

    for (int i = 0; i < L; i++) {
        int seq_idx = reverse ? (L - 1 - i) : i;
        char base = seq_ptr[seq_idx];
        if (reverse) {
            base = complement_base(base);
        }

        int base_idx = base_to_index(base);

        // BUG-1 FIX: check for unknown bases before calling get_log_prob
        if (base_idx == 4) {
            // Unknown base: treat as mandatory edit with max score assumption
            mandatory_edits++;
            adjusted_score += static_cast<double>(m_col_max_scores[i]);
            continue;
        }

        float base_score = m_pssm[i].get_log_prob_from_code(base_idx);

        if (base_score <= kLogZeroThreshold || !std::isfinite(base_score)) {
            mandatory_edits++;
            adjusted_score += static_cast<double>(m_col_max_scores[i]);
            continue;
        }

        adjusted_score += static_cast<double>(base_score);

        int bin = m_bin_index[i][base_idx];
        if (m_exact_count[bin] == 0) {
            m_exact_touched.push_back(bin);
        }
        m_exact_count[bin]++;
    }

    double deficit = static_cast<double>(m_threshold) - adjusted_score;

    // Clean up count vector using touched list before any early return
    auto cleanup = [this]() {
        for (size_t idx : m_exact_touched) {
            m_exact_count[idx] = 0;
        }
    };

    if (deficit <= 0.0) {
        cleanup();
        return static_cast<float>(mandatory_edits);
    }

    double max_possible_gain = static_cast<double>(m_S_max) - adjusted_score;
    if (max_possible_gain < deficit) {
        cleanup();
        return std::numeric_limits<float>::quiet_NaN();
    }

    double acc = 0.0;
    int edits = mandatory_edits;

    for (size_t j = 0; j < m_gain_values.size(); j++) {
        int bin_count = m_exact_count[j];
        if (!bin_count) {
            continue;
        }

        double gain = static_cast<double>(m_gain_values[j]);
        if (gain <= 0.0) {
            continue;
        }

        double cap = gain * static_cast<double>(bin_count);

        if (acc + cap < deficit) {
            acc += cap;
            edits += bin_count;
        } else {
            double remaining = deficit - acc;
            if (remaining <= 0.0) {
                cleanup();
                return static_cast<float>(edits);
            }

            double need_raw = remaining / gain;
            if (!std::isfinite(need_raw)) {
                need_raw = static_cast<double>(bin_count);
            }

            int need = static_cast<int>(std::ceil(need_raw));
            need = std::max(1, std::min(need, bin_count));

            edits += need;
            cleanup();
            return static_cast<float>(edits);
        }
    }

    cleanup();
    return std::numeric_limits<float>::quiet_NaN();
}

float PWMEditDistanceScorer::compute_heuristic(const char* seq_ptr, bool reverse, int max_k)
{
    const int L = m_pssm.length();
    if (L == 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    int mandatory_edits = 0;
    double adjusted_score = 0.0;
    std::vector<float> gains;
    gains.reserve(L);

    for (int i = 0; i < L; i++) {
        int seq_idx = reverse ? (L - 1 - i) : i;
        char base = seq_ptr[seq_idx];
        if (reverse) {
            base = complement_base(base);
        }

        int base_idx = base_to_index(base);

        // BUG-1 FIX: check for unknown bases before calling get_log_prob
        if (base_idx == 4) {
            // Unknown base: treat as mandatory edit
            mandatory_edits++;
            adjusted_score += static_cast<double>(m_col_max_scores[i]);
            continue;
        }

        float base_score = m_pssm[i].get_log_prob_from_code(base_idx);

        if (base_score <= kLogZeroThreshold || !std::isfinite(base_score)) {
            mandatory_edits++;
            adjusted_score += static_cast<double>(m_col_max_scores[i]);
        } else {
            adjusted_score += static_cast<double>(base_score);
            gains.push_back(m_col_max_scores[i] - base_score);
        }
    }

    if (mandatory_edits > max_k) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    double deficit = static_cast<double>(m_threshold) - adjusted_score;
    if (deficit <= 0.0) {
        return static_cast<float>(mandatory_edits);
    }

    double max_possible_gain = static_cast<double>(m_S_max) - adjusted_score;
    if (max_possible_gain < deficit) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    int remaining_budget = max_k - mandatory_edits;
    if (remaining_budget <= 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    if (remaining_budget > static_cast<int>(gains.size())) {
        remaining_budget = static_cast<int>(gains.size());
    }

    if (remaining_budget <= 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    if (remaining_budget < static_cast<int>(gains.size())) {
        std::nth_element(gains.begin(), gains.begin() + remaining_budget, gains.end(),
                         std::greater<float>());
    }
    std::sort(gains.begin(), gains.begin() + remaining_budget, std::greater<float>());

    double acc = 0.0;
    for (int i = 0; i < remaining_budget; i++) {
        acc += static_cast<double>(gains[i]);
        if (acc >= deficit) {
            return static_cast<float>(mandatory_edits + i + 1);
        }
    }

    return std::numeric_limits<float>::quiet_NaN();
}

inline bool PWMEditDistanceScorer::should_scan_forward() const
{
    return m_pssm.is_bidirect() || m_strand >= 0;
}

inline bool PWMEditDistanceScorer::should_scan_reverse() const
{
    return m_pssm.is_bidirect() || m_strand <= 0;
}

inline char PWMEditDistanceScorer::complement_base(char base) const
{
    switch (base) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default: return 'N';
    }
}

float PWMEditDistanceScorer::encode_position(size_t index,
                                             size_t target_length,
                                             size_t motif_length,
                                             int direction) const
{
    float pos = static_cast<float>(index) + 1.0f;
    if (m_strand == -1) {
        pos = static_cast<float>(target_length) - pos - static_cast<float>(motif_length) + 1.0f;
    }
    if (m_pssm.is_bidirect()) {
        pos *= static_cast<float>(direction);
    }
    return pos;
}

float PWMEditDistanceScorer::compute_window_edits(const char* window_start, int seq_avail, bool reverse)
{
    // When indels are enabled, use banded NW DP
    if (m_max_indels > 0) {
        return compute_with_indels(window_start, seq_avail, reverse);
    }

    // Fast path: substitutions only
    if (m_max_edits < 0) {
        return compute_exact(window_start, reverse);
    }
    return compute_heuristic(window_start, reverse, m_max_edits);
}

PWMEditDistanceScorer::ScanMetrics PWMEditDistanceScorer::evaluate_windows(const std::string& seq,
                                                                           size_t interval_length)
{
    ScanMetrics metrics;

    const size_t motif_length = static_cast<size_t>(m_pssm.length());
    const size_t target_length = seq.length();
    if (motif_length == 0 || target_length < motif_length || interval_length == 0) {
        return metrics;
    }

    const bool scan_forward = should_scan_forward();
    const bool scan_reverse = should_scan_reverse();
    const bool need_min = (m_mode != Mode::PWM_MAX_EDITS);
    const bool need_min_pos = (m_mode == Mode::MIN_EDITS_POSITION);
    const bool need_pwm = (m_mode == Mode::PWM_MAX_EDITS);
    const bool has_score_min = !std::isnan(m_score_min);
    const bool has_score_max = !std::isnan(m_score_max);
    const bool has_score_filter = has_score_min || has_score_max;

    const size_t max_start = target_length - motif_length + 1;
    if (max_start == 0 || (!scan_forward && !scan_reverse)) {
        return metrics;
    }

    const char* seq_data = seq.data();
    bool pwm_found = false;

    auto maybe_update_min = [&](float edits, size_t idx, int direction) {
        if (!need_min || std::isnan(edits)) {
            return;
        }
        if (std::isnan(metrics.min_edits) ||
            edits < metrics.min_edits - 1e-6f ||
            (std::fabs(edits - metrics.min_edits) <= 1e-6f &&
             (idx < metrics.min_index ||
              (idx == metrics.min_index && direction > metrics.min_direction)))) {
            metrics.min_edits = edits;
            metrics.min_index = idx;
            metrics.min_direction = direction;
            if (need_min_pos) {
                metrics.min_edits_position = encode_position(idx, target_length, motif_length, direction);
            }
        }
    };

    auto maybe_update_pwm = [&](float logp, size_t idx, int direction) {
        if (!need_pwm) {
            return;
        }
        if (!pwm_found || logp > metrics.best_pwm_logp) {
            metrics.best_pwm_logp = logp;
            metrics.best_pwm_index = idx;
            metrics.best_pwm_direction = direction;
            pwm_found = true;
        }
    };

    for (size_t offset = 0; offset < max_start; ++offset) {
        const char* window_start = seq_data + offset;
        int seq_avail = static_cast<int>(target_length - offset);

        if (scan_forward) {
            if (need_min) {
                // score.min/score.max filtering: skip edit distance if PWM score is out of range
                bool pass_score_filter = true;
                if (has_score_filter) {
                    float logp = compute_window_pwm_score(window_start, /*reverse=*/false);
                    if (has_score_min && logp < m_score_min) pass_score_filter = false;
                    if (has_score_max && logp > m_score_max) pass_score_filter = false;
                }

                if (pass_score_filter) {
                    float edits = compute_window_edits(window_start, seq_avail, /*reverse=*/false);
                    maybe_update_min(edits, offset, +1);
                }
            }
            if (need_pwm) {
                float logp = 0.0f;
                std::string::const_iterator it = seq.begin() + offset;
                m_pssm.calc_like(it, logp);
                maybe_update_pwm(logp, offset, +1);
            }
        }

        if (scan_reverse) {
            if (need_min) {
                // score.min/score.max filtering: skip edit distance if PWM score is out of range
                bool pass_score_filter = true;
                if (has_score_filter) {
                    float logp = compute_window_pwm_score(window_start, /*reverse=*/true);
                    if (has_score_min && logp < m_score_min) pass_score_filter = false;
                    if (has_score_max && logp > m_score_max) pass_score_filter = false;
                }

                if (pass_score_filter) {
                    float edits = compute_window_edits(window_start, seq_avail, /*reverse=*/true);
                    maybe_update_min(edits, offset, -1);
                }
            }
            if (need_pwm) {
                float logp_rc = 0.0f;
                std::string::const_iterator it2 = seq.begin() + offset;
                m_pssm.calc_like_rc(it2, logp_rc);
                maybe_update_pwm(logp_rc, offset, -1);
            }
        }
    }

    if (need_pwm && pwm_found) {
        // For PWM_MAX_EDITS mode: compute edit distance at the best PWM window
        // Apply score.min/score.max filter
        if ((has_score_min && metrics.best_pwm_logp < m_score_min) ||
            (has_score_max && metrics.best_pwm_logp > m_score_max)) {
            metrics.best_pwm_edits = std::numeric_limits<float>::quiet_NaN();
        } else {
            const char* window_start = seq_data + metrics.best_pwm_index;
            int seq_avail = static_cast<int>(target_length - metrics.best_pwm_index);
            metrics.best_pwm_edits = compute_window_edits(window_start, seq_avail, metrics.best_pwm_direction < 0);
        }
        metrics.best_pwm_position = encode_position(metrics.best_pwm_index,
                                                    target_length,
                                                    motif_length,
                                                    metrics.best_pwm_direction);
    } else if (need_pwm) {
        metrics.best_pwm_edits = std::numeric_limits<float>::quiet_NaN();
        metrics.best_pwm_position = std::numeric_limits<float>::quiet_NaN();
    }

    return metrics;
}

float PWMEditDistanceScorer::compute_with_indels(const char* seq_ptr, int seq_len, bool reverse)
{
    const int L = m_pssm.length();
    if (L == 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    const int D = m_max_indels;
    float best_edits = std::numeric_limits<float>::quiet_NaN();

    // For each sequence window length W in [L-D, L+D], align the motif (length L)
    // against W sequence bases using a 3D DP:
    //   dp[i][j][k] = max PWM score aligning motif[0..i-1] with seq[0..j-1]
    //                 using exactly k indels (insertions + deletions)
    //
    // Band constraint: |i - j| <= D (prevents needing > D indels)
    //
    // For each final state dp[L][W][k], the score tells us the PWM log-likelihood
    // of the aligned positions with their actual sequence bases. To reach the
    // threshold, we may need additional substitutions: each substitution replaces
    // a mismatched base with the column-optimal base, gaining (col_max - current).
    // We traceback to find aligned positions and greedily pick the largest gains.
    //
    // Total edits = k (indels) + subs_needed.
    // We minimize this across all W and k.

    for (int W = std::max(1, L - D); W <= L + D; ++W) {
        if (W > seq_len) {
            break;
        }

        const int rows = L + 1;
        const int cols = W + 1;
        const int indel_levels = D + 1;  // k = 0, 1, ..., D

        // Flattened 3D DP table
        std::vector<double> dp(rows * cols * indel_levels, -std::numeric_limits<double>::infinity());

        auto idx3 = [cols, indel_levels](int i, int j, int k) -> int {
            return i * cols * indel_levels + j * indel_levels + k;
        };

        // Base case: empty alignment, 0 indels
        dp[idx3(0, 0, 0)] = 0.0;

        // First column: skip motif positions (each costs 1 indel)
        for (int i = 1; i <= std::min(L, D); ++i) {
            dp[idx3(i, 0, i)] = 0.0;
        }

        // First row: skip sequence positions (each costs 1 indel)
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
                if (reverse) {
                    base = complement_base(base);
                }
                int bidx = base_to_index(base);

                float base_score;
                if (bidx == 4) {
                    // Unknown base: use minimum score (worst case)
                    float min_s = std::numeric_limits<float>::infinity();
                    for (int b = 0; b < 4; b++) {
                        float s = m_pssm[i - 1].get_log_prob_from_code(b);
                        if (s < min_s) min_s = s;
                    }
                    base_score = min_s;
                } else {
                    base_score = m_pssm[i - 1].get_log_prob_from_code(bidx);
                }

                for (int k = 0; k <= D; ++k) {
                    // 1. Match/Substitution (diagonal): align motif[i-1] with seq[j-1]
                    if (std::abs((i - 1) - (j - 1)) <= D) {
                        double prev = dp[idx3(i - 1, j - 1, k)];
                        if (prev > -std::numeric_limits<double>::infinity() * 0.5) {
                            double new_score = prev + static_cast<double>(base_score);
                            if (new_score > dp[idx3(i, j, k)]) {
                                dp[idx3(i, j, k)] = new_score;
                            }
                        }
                    }

                    if (k < D) {
                        // 2. Insertion: skip motif[i-1], advance motif not sequence
                        if (std::abs((i - 1) - j) <= D) {
                            double prev = dp[idx3(i - 1, j, k)];
                            if (prev > -std::numeric_limits<double>::infinity() * 0.5) {
                                if (prev > dp[idx3(i, j, k + 1)]) {
                                    dp[idx3(i, j, k + 1)] = prev;
                                }
                            }
                        }

                        // 3. Deletion: skip seq[j-1], advance sequence not motif
                        if (std::abs(i - (j - 1)) <= D) {
                            double prev = dp[idx3(i, j - 1, k)];
                            if (prev > -std::numeric_limits<double>::infinity() * 0.5) {
                                if (prev > dp[idx3(i, j, k + 1)]) {
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
            if (score <= -std::numeric_limits<double>::infinity() * 0.5) {
                continue;
            }

            // If score already reaches threshold, no substitutions needed
            if (score >= static_cast<double>(m_threshold)) {
                float total_edits = static_cast<float>(k);
                if (std::isnan(best_edits) || total_edits < best_edits) {
                    best_edits = total_edits;
                }
                continue;
            }

            // Traceback to find aligned motif positions and their gains
            std::vector<float> aligned_gains;
            aligned_gains.reserve(L);
            int ti = L, tj = W, tk = k;

            while (ti > 0 || tj > 0) {
                if (ti == 0 && tj > 0 && tk > 0) {
                    // Must be deletion
                    tj--; tk--;
                    continue;
                }
                if (tj == 0 && ti > 0 && tk > 0) {
                    // Must be insertion
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
                            float s = m_pssm[ti - 1].get_log_prob_from_code(bb);
                            if (s < min_s) min_s = s;
                        }
                        bs = min_s;
                    } else {
                        bs = m_pssm[ti - 1].get_log_prob_from_code(bi);
                    }

                    double prev = dp[idx3(ti - 1, tj - 1, tk)];
                    if (prev > -std::numeric_limits<double>::infinity() * 0.5 &&
                        std::fabs((prev + static_cast<double>(bs)) - cur) < 1e-9) {
                        float gain = m_col_max_scores[ti - 1] - bs;
                        if (gain > 1e-12f) {
                            aligned_gains.push_back(gain);
                        }
                        ti--; tj--;
                        found = true;
                    }
                }

                if (!found && ti > 0 && tk > 0 && std::abs((ti - 1) - tj) <= D) {
                    // Try insertion
                    double prev = dp[idx3(ti - 1, tj, tk - 1)];
                    if (prev > -std::numeric_limits<double>::infinity() * 0.5 &&
                        std::fabs(prev - cur) < 1e-9) {
                        ti--; tk--;
                        found = true;
                    }
                }

                if (!found && tj > 0 && tk > 0 && std::abs(ti - (tj - 1)) <= D) {
                    // Try deletion
                    double prev = dp[idx3(ti, tj - 1, tk - 1)];
                    if (prev > -std::numeric_limits<double>::infinity() * 0.5 &&
                        std::fabs(prev - cur) < 1e-9) {
                        tj--; tk--;
                        found = true;
                    }
                }

                if (!found) break;  // DP inconsistency (shouldn't happen)
            }

            // Compute minimum substitutions needed
            double deficit = static_cast<double>(m_threshold) - score;
            if (deficit <= 0.0) {
                float total_edits = static_cast<float>(k);
                if (std::isnan(best_edits) || total_edits < best_edits) {
                    best_edits = total_edits;
                }
                continue;
            }

            // Check total available gain
            double total_gain = 0.0;
            for (float g : aligned_gains) {
                total_gain += static_cast<double>(g);
            }
            if (total_gain < deficit) {
                continue;  // Cannot reach threshold
            }

            // Sort gains descending, greedily pick largest
            std::sort(aligned_gains.begin(), aligned_gains.end(), std::greater<float>());

            double acc = 0.0;
            int subs = 0;
            for (float g : aligned_gains) {
                acc += static_cast<double>(g);
                subs++;
                if (acc >= deficit) break;
            }

            float total_edits = static_cast<float>(k + subs);
            if (std::isnan(best_edits) || total_edits < best_edits) {
                best_edits = total_edits;
            }
        }
    }

    // Apply max_edits cap if set
    if (!std::isnan(best_edits) && m_max_edits >= 1 && best_edits > static_cast<float>(m_max_edits)) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    return best_edits;
}
