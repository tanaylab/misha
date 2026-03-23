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

    // Precompute sorted column maxima (ascending) for insertion-family bounds
    m_sorted_col_max_asc = m_col_max_scores;
    std::sort(m_sorted_col_max_asc.begin(), m_sorted_col_max_asc.end());

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
    // When indels are enabled, use specialized or generic DP
    if (m_max_indels == 1) {
        return compute_with_one_indel(window_start, seq_avail, reverse);
    } else if (m_max_indels == 2) {
        return compute_with_two_indels(window_start, seq_avail, reverse);
    } else if (m_max_indels > 2) {
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
    const bool apply_score_filter = has_score_filter && (m_max_indels == 0);

    const size_t max_start = std::min(interval_length, target_length - motif_length + 1);
    if (max_start == 0 || (!scan_forward && !scan_reverse)) {
        return metrics;
    }

    const char* seq_data = seq.data();
    bool pwm_found = false;

    // Reachability pre-filter: when indels are enabled with an edit budget,
    // compute a cheap lower bound on minimum edits before calling the full solver.
    const bool use_indel_prefilter = (m_max_indels > 0 && m_max_edits >= 0);

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
                if (apply_score_filter) {
                    float logp = compute_window_pwm_score(window_start, /*reverse=*/false);
                    if (has_score_min && logp < m_score_min) pass_score_filter = false;
                    if (has_score_max && logp > m_score_max) pass_score_filter = false;
                }

                // Reachability pre-filter: skip if lower bound exceeds edit budget
                if (pass_score_filter && use_indel_prefilter) {
                    int lb = compute_indel_lower_bound(window_start, seq_avail, /*reverse=*/false);
                    if (lb > m_max_edits) {
                        pass_score_filter = false;
                    }
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
                if (apply_score_filter) {
                    float logp = compute_window_pwm_score(window_start, /*reverse=*/true);
                    if (has_score_min && logp < m_score_min) pass_score_filter = false;
                    if (has_score_max && logp > m_score_max) pass_score_filter = false;
                }

                // Reachability pre-filter: skip if lower bound exceeds edit budget
                if (pass_score_filter && use_indel_prefilter) {
                    int lb = compute_indel_lower_bound(window_start, seq_avail, /*reverse=*/true);
                    if (lb > m_max_edits) {
                        pass_score_filter = false;
                    }
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
        // Apply score.min/score.max filter (bypass when indels are enabled)
        if (apply_score_filter &&
            ((has_score_min && metrics.best_pwm_logp < m_score_min) ||
             (has_score_max && metrics.best_pwm_logp > m_score_max))) {
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
                        std::fabs((prev + static_cast<double>(bs)) - cur) < 1e-9 * std::max(1.0, std::fabs(cur))) {
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
                        std::fabs(prev - cur) < 1e-9 * std::max(1.0, std::fabs(cur))) {
                        ti--; tk--;
                        found = true;
                    }
                }

                if (!found && tj > 0 && tk > 0 && std::abs(ti - (tj - 1)) <= D) {
                    // Try deletion
                    double prev = dp[idx3(ti, tj - 1, tk - 1)];
                    if (prev > -std::numeric_limits<double>::infinity() * 0.5 &&
                        std::fabs(prev - cur) < 1e-9 * std::max(1.0, std::fabs(cur))) {
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

int PWMEditDistanceScorer::compute_indel_lower_bound(const char* window_start, int seq_avail, bool reverse)
{
    // Compute a lower bound on the minimum edits needed for ANY alignment family
    // (no-indel, 1-indel, 2-indel) to reach the threshold.
    //
    // Strategy:
    // 1. Compute the no-indel raw aligned score (O(L)).
    // 2. For the no-indel family (0 indels): compute min substitutions needed.
    //    This is a lower bound because we use partial sort to find the minimum.
    // 3. For D-indel families: the best possible score achievable is
    //    S_max (for deletions, all L columns get their best base) or
    //    S_max - sum_of_D_smallest_col_maxes (for insertions, we lose D columns).
    //    The minimum indels is D, so the lower bound for D-indel families is D
    //    (if the best achievable score >= threshold) or infinity (if not reachable).
    //
    // But the real per-window bound comes from the no-indel score:
    //   For the no-indel alignment, we know the exact score and gains.
    //   A D-indel alignment can shift at most D columns. The maximum score gain
    //   from shifting is bounded. But computing that tightly is expensive.
    //
    // PRACTICAL approach: use the no-indel raw score to compute a lower bound on
    // substitutions needed for the no-indel family, then check if indel families
    // could possibly do better enough to matter.
    //
    // The tightest cheap bound: for any alignment family with D indels,
    // total_edits >= D + subs_needed. The subs_needed for the best possible
    // D-indel alignment is 0 if S_max_achievable >= threshold. So the true
    // lower bound across all families is min(subs_no_indel, 1, 2, ...).
    // That's just 0 — not useful.
    //
    // Better: compute subs_no_indel (exact for 0-indel family). Then the overall
    // minimum is min(subs_no_indel, 1 + subs_1indel, 2 + subs_2indel, ...).
    // If subs_no_indel <= m_max_edits, we can't prune. So the useful case is
    // when the no-indel raw score is far below threshold.
    //
    // Key insight: the no-indel family's min_subs is a lower bound for the
    // 0-indel case. For D-indel families, the best achievable score is at most
    // S_max (deletions) or S_max - sum_of_D_smallest_col_max (insertions).
    // The min edits for D-indel family is D if reachable, otherwise infinity.
    //
    // So: lower_bound = min over D in {0..max_indels} of:
    //   D + max(0, ceil of subs needed if all aligned columns are optimal)
    //
    // For D=0: subs = subs_no_indel (computed from raw score + gains)
    // For D>=1 (deletions): if S_max >= threshold, lower_bound_D = D
    // For D>=1 (insertions): if S_max - sum_smallest_D >= threshold, lower_bound_D = D
    //
    // So the overall lower bound = min(subs_no_indel, min_D_reachable)
    // where min_D_reachable = smallest D such that some D-indel family is reachable.

    const int L = m_pssm.length();
    if (L == 0) return std::numeric_limits<int>::max();

    // Step 1: Compute no-indel raw score and per-column gains
    double raw_score = 0.0;
    std::vector<float> gains;
    gains.reserve(L);

    for (int i = 0; i < L; i++) {
        int seq_idx = reverse ? (L - 1 - i) : i;
        char base = window_start[seq_idx];
        if (reverse) {
            base = complement_base(base);
        }
        int bidx = base_to_index(base);
        float base_score;
        if (bidx == 4) {
            // Unknown base: use minimum score
            float min_s = std::numeric_limits<float>::infinity();
            for (int b = 0; b < 4; b++) {
                float s = m_pssm[i].get_log_prob_from_code(b);
                if (s < min_s) min_s = s;
            }
            base_score = min_s;
        } else {
            base_score = m_pssm[i].get_log_prob_from_code(bidx);
        }
        raw_score += static_cast<double>(base_score);
        float gain = m_col_max_scores[i] - base_score;
        if (gain > 1e-12f) {
            gains.push_back(gain);
        }
    }

    // Step 2: Compute min substitutions for no-indel family
    double deficit = static_cast<double>(m_threshold) - raw_score;
    int subs_no_indel;
    if (deficit <= 0.0) {
        subs_no_indel = 0;
    } else {
        // Check if reachable at all
        double total_gain = 0.0;
        for (float g : gains) total_gain += static_cast<double>(g);
        if (total_gain < deficit) {
            subs_no_indel = std::numeric_limits<int>::max(); // unreachable
        } else {
            // Use partial sort to find min subs needed
            // We need the smallest k such that sum of top-k gains >= deficit
            std::sort(gains.begin(), gains.end(), std::greater<float>());
            double acc = 0.0;
            subs_no_indel = 0;
            for (size_t i = 0; i < gains.size(); i++) {
                acc += static_cast<double>(gains[i]);
                subs_no_indel++;
                if (acc >= deficit) break;
            }
            if (acc < deficit) {
                subs_no_indel = std::numeric_limits<int>::max();
            }
        }
    }

    int best_lower_bound = subs_no_indel;

    // Step 3: Check D-indel families (D = 1 to m_max_indels)
    for (int D = 1; D <= m_max_indels; D++) {
        // Deletion family (W = L + D): all L motif columns are aligned.
        // Best achievable = S_max (all columns get their best base).
        if (static_cast<double>(m_S_max) >= static_cast<double>(m_threshold)) {
            // Reachable with D indels + 0 subs = D edits
            if (D < best_lower_bound) {
                best_lower_bound = D;
            }
        }

        // Insertion family (W = L - D): L - D motif columns are aligned.
        // Best achievable = S_max - sum of D smallest col_max values.
        if (D <= L - 1 && D <= static_cast<int>(m_sorted_col_max_asc.size())) {
            double ins_best = static_cast<double>(m_S_max);
            for (int j = 0; j < D; j++) {
                ins_best -= static_cast<double>(m_sorted_col_max_asc[j]);
            }
            if (ins_best >= static_cast<double>(m_threshold)) {
                if (D < best_lower_bound) {
                    best_lower_bound = D;
                }
            }
        }
    }

    return best_lower_bound;
}

inline void PWMEditDistanceScorer::get_aligned_base_score(const char* seq_ptr, bool reverse,
                                                          int motif_pos, int raw_seq_idx,
                                                          float& out_score, float& out_gain) const
{
    char base = seq_ptr[raw_seq_idx];
    if (reverse) {
        base = complement_base(base);
    }
    int bidx = base_to_index(base);
    if (bidx == 4) {
        // Unknown base: use minimum score (worst case), matching compute_with_indels
        float min_s = std::numeric_limits<float>::infinity();
        for (int b = 0; b < 4; b++) {
            float s = m_pssm[motif_pos].get_log_prob_from_code(b);
            if (s < min_s) min_s = s;
        }
        out_score = min_s;
    } else {
        out_score = m_pssm[motif_pos].get_log_prob_from_code(bidx);
    }
    out_gain = m_col_max_scores[motif_pos] - out_score;
}

float PWMEditDistanceScorer::compute_min_edits_from_gains(double aligned_score,
                                                          std::vector<float>& gains,
                                                          int indels,
                                                          float best_edits_so_far) const
{
    // Prune: if indels alone already can't beat best_edits, skip
    if (!std::isnan(best_edits_so_far) && static_cast<float>(indels) >= best_edits_so_far) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    // Prune: if indels alone exceed budget, skip
    if (m_max_edits >= 0 && indels > m_max_edits) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    double deficit = static_cast<double>(m_threshold) - aligned_score;
    if (deficit <= 0.0) {
        return static_cast<float>(indels);
    }

    // Compute max subs budget based on best_edits and m_max_edits
    int max_subs = static_cast<int>(gains.size());
    if (!std::isnan(best_edits_so_far)) {
        int ceiling = static_cast<int>(best_edits_so_far) - 1 - indels;
        if (ceiling < max_subs) max_subs = ceiling;
    }
    if (m_max_edits >= 0) {
        int budget_subs = m_max_edits - indels;
        if (budget_subs < max_subs) max_subs = budget_subs;
    }
    if (max_subs <= 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    // Check total available gain
    double total_gain = 0.0;
    for (float g : gains) {
        total_gain += static_cast<double>(g);
    }
    if (total_gain < deficit) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    // Use partial sort when max_subs < gains.size() for efficiency
    int n = static_cast<int>(gains.size());
    if (max_subs < n) {
        std::nth_element(gains.begin(), gains.begin() + max_subs, gains.end(),
                         std::greater<float>());
        std::sort(gains.begin(), gains.begin() + max_subs, std::greater<float>());
    } else {
        std::sort(gains.begin(), gains.end(), std::greater<float>());
    }

    double acc = 0.0;
    int subs = 0;
    int limit = std::min(max_subs, n);
    for (int i = 0; i < limit; i++) {
        acc += static_cast<double>(gains[i]);
        subs++;
        if (acc >= deficit) {
            return static_cast<float>(indels + subs);
        }
    }

    return std::numeric_limits<float>::quiet_NaN();
}

float PWMEditDistanceScorer::compute_with_one_indel(const char* seq_ptr, int seq_len, bool reverse)
{
    const int L = m_pssm.length();
    if (L == 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    float best_edits = std::numeric_limits<float>::quiet_NaN();

    // Shorthand for the shared helpers, capturing seq_ptr and reverse
    auto get_base_score = [&](int motif_pos, int raw_seq_idx, float& out_score, float& out_gain) {
        get_aligned_base_score(seq_ptr, reverse, motif_pos, raw_seq_idx, out_score, out_gain);
    };

    auto compute_total_edits = [&](double aligned_score, std::vector<float>& gains, int indels) -> float {
        return compute_min_edits_from_gains(aligned_score, gains, indels, best_edits);
    };

    // ===== Case A: No indel (W = L, k = 0) =====
    // Motif position i aligns with sequence position i (or L-1-i if reverse).
    // This is equivalent to what the DP computes for W=L, k=0.
    {
        double aligned_score = 0.0;
        std::vector<float> gains;
        gains.reserve(L);

        for (int i = 0; i < L; i++) {
            int raw_seq_idx = reverse ? (L - 1 - i) : i;
            float bs, gain;
            get_base_score(i, raw_seq_idx, bs, gain);
            aligned_score += static_cast<double>(bs);
            if (gain > 1e-12f) {
                gains.push_back(gain);
            }
        }

        float total = compute_total_edits(aligned_score, gains, 0);
        if (!std::isnan(total) && (std::isnan(best_edits) || total < best_edits)) {
            best_edits = total;
        }

        // Early return: 0 edits is optimal
        if (!std::isnan(best_edits) && best_edits <= 0.0f) {
            return best_edits;
        }
    }

    // ===== Case B: One deletion in sequence (W = L + 1, k = 1) =====
    // We have L+1 sequence bases and L motif positions. One sequence base is skipped.
    // The deletion at sequence position d means:
    //   motif[0..d-1]  aligns with seq[0..d-1]
    //   seq[d] is deleted (skipped)
    //   motif[d..L-1]  aligns with seq[d+1..L]
    //
    // Use prefix/suffix precomputation for O(1) per deletion position.
    if (L + 1 <= seq_len) {
        const int W = L + 1;

        // Precompute base scores and gains for all W sequence positions aligned with
        // motif columns in two configurations:
        //   prefix: motif[i] aligned with seq_raw[i] for i = 0..L-1
        //   suffix: motif[i] aligned with seq_raw[i+1] for i = 0..L-1
        // where seq_raw[j] = reverse ? (W-1-j) : j

        std::vector<float> prefix_scores(L);   // motif[i] <-> seq[i]
        std::vector<float> prefix_gains(L);
        std::vector<float> suffix_scores(L);   // motif[i] <-> seq[i+1]
        std::vector<float> suffix_gains(L);

        for (int i = 0; i < L; i++) {
            int raw_prefix = reverse ? (W - 1 - i) : i;
            get_base_score(i, raw_prefix, prefix_scores[i], prefix_gains[i]);

            int raw_suffix = reverse ? (W - 1 - (i + 1)) : (i + 1);
            get_base_score(i, raw_suffix, suffix_scores[i], suffix_gains[i]);
        }

        // Prefix sums for prefix alignment (motif[0..d-1] with seq[0..d-1])
        std::vector<double> prefix_cum_score(L + 1, 0.0);
        for (int i = 0; i < L; i++) {
            prefix_cum_score[i + 1] = prefix_cum_score[i] + static_cast<double>(prefix_scores[i]);
        }

        // Suffix sums for shifted alignment (motif[d..L-1] with seq[d+1..L])
        std::vector<double> suffix_cum_score(L + 1, 0.0);
        for (int i = L - 1; i >= 0; i--) {
            suffix_cum_score[i] = suffix_cum_score[i + 1] + static_cast<double>(suffix_scores[i]);
        }

        // Try each deletion position d in [0, L]
        // d=0: skip seq[0], motif[0..L-1] aligns with seq[1..L]
        // d=L: skip seq[L], motif[0..L-1] aligns with seq[0..L-1]
        for (int d = 0; d <= L; d++) {
            double aligned_score = prefix_cum_score[d] + suffix_cum_score[d];

            // Collect gains from aligned positions
            std::vector<float> gains;
            gains.reserve(L);
            for (int i = 0; i < d; i++) {
                if (prefix_gains[i] > 1e-12f) {
                    gains.push_back(prefix_gains[i]);
                }
            }
            for (int i = d; i < L; i++) {
                if (suffix_gains[i] > 1e-12f) {
                    gains.push_back(suffix_gains[i]);
                }
            }

            float total = compute_total_edits(aligned_score, gains, 1);
            if (!std::isnan(total) && (std::isnan(best_edits) || total < best_edits)) {
                best_edits = total;
                // Early exit: best_edits == 1 means just the indel, can't do better in this case
                if (best_edits <= 1.0f) break;
            }
        }
    }

    // Early return: if best_edits <= 1, no insertion case can beat it
    // (insertion also costs 1 indel)
    if (!std::isnan(best_edits) && best_edits <= 1.0f) {
        // Apply max_edits cap if set
        if (m_max_edits >= 1 && best_edits > static_cast<float>(m_max_edits)) {
            return std::numeric_limits<float>::quiet_NaN();
        }
        return best_edits;
    }

    // ===== Case C: One insertion in motif (W = L - 1, k = 1) =====
    // We have L-1 sequence bases and L motif positions. One motif position is skipped.
    // The insertion at motif position m means:
    //   motif[0..m-1]    aligns with seq[0..m-1]
    //   motif[m]          is skipped (insertion)
    //   motif[m+1..L-1]  aligns with seq[m..L-2]
    //
    // Need L-1 >= 1 (i.e. L >= 2) for this to apply.
    if (L >= 2) {
        const int W = L - 1;

        // Precompute base scores and gains for aligned columns:
        //   prefix: motif[i] aligned with seq_raw[i] for i = 0..L-2
        //   suffix: motif[i] aligned with seq_raw[i-1] for i = 1..L-1
        // where seq_raw[j] = reverse ? (W-1-j) : j

        std::vector<float> prefix_scores(L - 1);   // motif[i] <-> seq[i], i=0..L-2
        std::vector<float> prefix_gains(L - 1);
        std::vector<float> suffix_scores(L);        // motif[i] <-> seq[i-1], i=1..L-1
        std::vector<float> suffix_gains(L);

        for (int i = 0; i < L - 1; i++) {
            int raw_idx = reverse ? (W - 1 - i) : i;
            get_base_score(i, raw_idx, prefix_scores[i], prefix_gains[i]);
        }

        // suffix: motif[i] aligns with seq[i-1] for i=1..L-1
        for (int i = 1; i < L; i++) {
            int raw_idx = reverse ? (W - 1 - (i - 1)) : (i - 1);
            get_base_score(i, raw_idx, suffix_scores[i], suffix_gains[i]);
        }

        // Prefix cumulative sums (motif[0..m-1] with seq[0..m-1])
        std::vector<double> prefix_cum_score(L, 0.0);  // prefix_cum_score[m] = sum for i=0..m-1
        for (int i = 0; i < L - 1; i++) {
            prefix_cum_score[i + 1] = prefix_cum_score[i] + static_cast<double>(prefix_scores[i]);
        }

        // Suffix cumulative sums (motif[m+1..L-1] with seq[m..L-2])
        std::vector<double> suffix_cum_score(L, 0.0);  // suffix_cum_score[m] = sum for i=m+1..L-1
        for (int i = L - 1; i >= 1; i--) {
            suffix_cum_score[i - 1] = suffix_cum_score[i] + static_cast<double>(suffix_scores[i]);
        }

        // Try each insertion position m in [0, L-1]
        // m=0: skip motif[0], motif[1..L-1] aligns with seq[0..L-2]
        // m=L-1: skip motif[L-1], motif[0..L-2] aligns with seq[0..L-2]
        for (int m = 0; m < L; m++) {
            double aligned_score = prefix_cum_score[m] + suffix_cum_score[m];

            // Collect gains from aligned positions
            std::vector<float> gains;
            gains.reserve(L - 1);
            for (int i = 0; i < m; i++) {
                if (prefix_gains[i] > 1e-12f) {
                    gains.push_back(prefix_gains[i]);
                }
            }
            for (int i = m + 1; i < L; i++) {
                if (suffix_gains[i] > 1e-12f) {
                    gains.push_back(suffix_gains[i]);
                }
            }

            float total = compute_total_edits(aligned_score, gains, 1);
            if (!std::isnan(total) && (std::isnan(best_edits) || total < best_edits)) {
                best_edits = total;
                // Early exit: best_edits == 1 means just the indel, can't do better here
                if (best_edits <= 1.0f) break;
            }
        }
    }

    // Apply max_edits cap if set
    if (!std::isnan(best_edits) && m_max_edits >= 1 && best_edits > static_cast<float>(m_max_edits)) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    return best_edits;
}

float PWMEditDistanceScorer::compute_with_two_indels(const char* seq_ptr, int seq_len, bool reverse)
{
    const int L = m_pssm.length();
    if (L == 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    float best_edits = std::numeric_limits<float>::quiet_NaN();

    // Shorthand for the shared helpers, capturing seq_ptr and reverse
    auto get_base_score = [&](int motif_pos, int raw_seq_idx, float& out_score, float& out_gain) {
        get_aligned_base_score(seq_ptr, reverse, motif_pos, raw_seq_idx, out_score, out_gain);
    };

    auto compute_total_edits = [&](double aligned_score, std::vector<float>& gains, int indels) -> float {
        return compute_min_edits_from_gains(aligned_score, gains, indels, best_edits);
    };

    auto update_best = [&](float total) {
        if (!std::isnan(total) && (std::isnan(best_edits) || total < best_edits)) {
            best_edits = total;
        }
    };

    // ===== Case 1: No indel (W = L, k = 0) =====
    {
        double aligned_score = 0.0;
        std::vector<float> gains;
        gains.reserve(L);

        for (int i = 0; i < L; i++) {
            int raw_seq_idx = reverse ? (L - 1 - i) : i;
            float bs, gain;
            get_base_score(i, raw_seq_idx, bs, gain);
            aligned_score += static_cast<double>(bs);
            if (gain > 1e-12f) {
                gains.push_back(gain);
            }
        }

        update_best(compute_total_edits(aligned_score, gains, 0));

        // Early return: 0 edits is optimal
        if (!std::isnan(best_edits) && best_edits <= 0.0f) {
            return best_edits;
        }
    }

    // ===== Case 2: One deletion (W = L + 1, k = 1) =====
    if (L + 1 <= seq_len) {
        const int W = L + 1;

        std::vector<float> prefix_scores(L), prefix_gains(L);
        std::vector<float> suffix_scores(L), suffix_gains(L);

        for (int i = 0; i < L; i++) {
            int raw_prefix = reverse ? (W - 1 - i) : i;
            get_base_score(i, raw_prefix, prefix_scores[i], prefix_gains[i]);

            int raw_suffix = reverse ? (W - 1 - (i + 1)) : (i + 1);
            get_base_score(i, raw_suffix, suffix_scores[i], suffix_gains[i]);
        }

        std::vector<double> prefix_cum(L + 1, 0.0);
        for (int i = 0; i < L; i++) {
            prefix_cum[i + 1] = prefix_cum[i] + static_cast<double>(prefix_scores[i]);
        }

        std::vector<double> suffix_cum(L + 1, 0.0);
        for (int i = L - 1; i >= 0; i--) {
            suffix_cum[i] = suffix_cum[i + 1] + static_cast<double>(suffix_scores[i]);
        }

        for (int d = 0; d <= L; d++) {
            double aligned_score = prefix_cum[d] + suffix_cum[d];

            std::vector<float> gains;
            gains.reserve(L);
            for (int i = 0; i < d; i++) {
                if (prefix_gains[i] > 1e-12f) gains.push_back(prefix_gains[i]);
            }
            for (int i = d; i < L; i++) {
                if (suffix_gains[i] > 1e-12f) gains.push_back(suffix_gains[i]);
            }

            update_best(compute_total_edits(aligned_score, gains, 1));
            if (!std::isnan(best_edits) && best_edits <= 1.0f) break;
        }
    }

    // ===== Case 3: One insertion (W = L - 1, k = 1) =====
    // Skip entirely if best_edits <= 1 (can't beat a 1-indel result with another 1-indel)
    if (L >= 2 && (std::isnan(best_edits) || best_edits > 1.0f)) {
        const int W = L - 1;

        std::vector<float> prefix_scores(L - 1), prefix_gains(L - 1);
        std::vector<float> suffix_scores(L), suffix_gains(L);

        for (int i = 0; i < L - 1; i++) {
            int raw_idx = reverse ? (W - 1 - i) : i;
            get_base_score(i, raw_idx, prefix_scores[i], prefix_gains[i]);
        }

        for (int i = 1; i < L; i++) {
            int raw_idx = reverse ? (W - 1 - (i - 1)) : (i - 1);
            get_base_score(i, raw_idx, suffix_scores[i], suffix_gains[i]);
        }

        std::vector<double> prefix_cum(L, 0.0);
        for (int i = 0; i < L - 1; i++) {
            prefix_cum[i + 1] = prefix_cum[i] + static_cast<double>(prefix_scores[i]);
        }

        std::vector<double> suffix_cum(L, 0.0);
        for (int i = L - 1; i >= 1; i--) {
            suffix_cum[i - 1] = suffix_cum[i] + static_cast<double>(suffix_scores[i]);
        }

        for (int m = 0; m < L; m++) {
            double aligned_score = prefix_cum[m] + suffix_cum[m];

            std::vector<float> gains;
            gains.reserve(L - 1);
            for (int i = 0; i < m; i++) {
                if (prefix_gains[i] > 1e-12f) gains.push_back(prefix_gains[i]);
            }
            for (int i = m + 1; i < L; i++) {
                if (suffix_gains[i] > 1e-12f) gains.push_back(suffix_gains[i]);
            }

            update_best(compute_total_edits(aligned_score, gains, 1));
            if (!std::isnan(best_edits) && best_edits <= 1.0f) break;
        }
    }

    // Early return: if best_edits <= 1, no 2-indel case can beat it
    if (!std::isnan(best_edits) && best_edits <= 1.0f) {
        if (m_max_edits >= 1 && best_edits > static_cast<float>(m_max_edits)) {
            return std::numeric_limits<float>::quiet_NaN();
        }
        return best_edits;
    }

    // ===== Case 4: Two deletions (W = L + 2, k = 2) =====
    // Sequence has L+2 bases. Skip 2 bases at positions d1 < d2 in [0, L+1].
    // Segment layout:
    //   motif[0..d1-1]     <-> seq[0..d1-1]       (config0: no shift)
    //   motif[d1..d2-2]    <-> seq[d1+1..d2-1]    (config1: shift +1)
    //   motif[d2-1..L-1]   <-> seq[d2+1..L+1]     (config2: shift +2)
    // Skip if best_edits <= 2 (no 2-indel candidate can beat it)
    if (L + 2 <= seq_len && (std::isnan(best_edits) || best_edits > 2.0f)) {
        const int W = L + 2;

        std::vector<float> scores0(L), gains0(L);
        std::vector<float> scores1(L), gains1(L);
        std::vector<float> scores2(L), gains2(L);

        for (int i = 0; i < L; i++) {
            int raw0 = reverse ? (W - 1 - i) : i;
            get_base_score(i, raw0, scores0[i], gains0[i]);

            int raw1 = reverse ? (W - 1 - (i + 1)) : (i + 1);
            get_base_score(i, raw1, scores1[i], gains1[i]);

            int raw2 = reverse ? (W - 1 - (i + 2)) : (i + 2);
            get_base_score(i, raw2, scores2[i], gains2[i]);
        }

        // Prefix cumulative for config0
        std::vector<double> cum0(L + 1, 0.0);
        for (int i = 0; i < L; i++) {
            cum0[i + 1] = cum0[i] + static_cast<double>(scores0[i]);
        }

        // Suffix cumulative for config2
        std::vector<double> suf2(L + 1, 0.0);
        for (int i = L - 1; i >= 0; i--) {
            suf2[i] = suf2[i + 1] + static_cast<double>(scores2[i]);
        }

        // Prefix cumulative for config1 (middle segment)
        std::vector<double> cum1(L + 1, 0.0);
        for (int i = 0; i < L; i++) {
            cum1[i + 1] = cum1[i] + static_cast<double>(scores1[i]);
        }

        std::vector<float> gains_4;
        gains_4.reserve(L);
        for (int d1 = 0; d1 <= L + 1; d1++) {
            for (int d2 = d1 + 1; d2 <= L + 1; d2++) {
                // Seg 1: motif[0..d1-1] config0
                double seg1 = cum0[d1];

                // Seg 2: motif[d1..d2-2] config1
                double seg2 = 0.0;
                if (d2 - 1 > d1) {
                    seg2 = cum1[d2 - 1] - cum1[d1];
                }

                // Seg 3: motif[d2-1..L-1] config2
                double seg3 = 0.0;
                if (d2 - 1 < L) {
                    seg3 = suf2[d2 - 1];
                }

                double aligned_score = seg1 + seg2 + seg3;

                gains_4.clear();
                for (int i = 0; i < d1; i++) {
                    if (gains0[i] > 1e-12f) gains_4.push_back(gains0[i]);
                }
                for (int i = d1; i < d2 - 1; i++) {
                    if (gains1[i] > 1e-12f) gains_4.push_back(gains1[i]);
                }
                for (int i = d2 - 1; i < L; i++) {
                    if (gains2[i] > 1e-12f) gains_4.push_back(gains2[i]);
                }

                update_best(compute_total_edits(aligned_score, gains_4, 2));
                if (!std::isnan(best_edits) && best_edits <= 2.0f) break;
            }
            if (!std::isnan(best_edits) && best_edits <= 2.0f) break;
        }
    }

    // ===== Case 5: Two insertions (W = L - 2, k = 2) =====
    // Skip 2 motif columns m1 < m2 in [0, L-1].
    // Segment layout:
    //   motif[0..m1-1]     <-> seq[0..m1-1]       (config0: no shift)
    //   motif[m1+1..m2-1]  <-> seq[m1..m2-2]      (config1: motif shifted, seq[i-1])
    //   motif[m2+1..L-1]   <-> seq[m2-1..L-3]     (config2: seq[i-2])
    // Skip if best_edits <= 2
    if (L >= 3 && (std::isnan(best_edits) || best_edits > 2.0f)) {
        const int W = L - 2;

        std::vector<float> scores0(L), gains0_v(L);
        std::vector<float> scores1(L), gains1_v(L);
        std::vector<float> scores2(L), gains2_v(L);

        for (int i = 0; i < L; i++) {
            if (i < W) {
                int raw0 = reverse ? (W - 1 - i) : i;
                get_base_score(i, raw0, scores0[i], gains0_v[i]);
            }
            if (i >= 1 && i - 1 < W) {
                int raw1 = reverse ? (W - 1 - (i - 1)) : (i - 1);
                get_base_score(i, raw1, scores1[i], gains1_v[i]);
            }
            if (i >= 2 && i - 2 < W) {
                int raw2 = reverse ? (W - 1 - (i - 2)) : (i - 2);
                get_base_score(i, raw2, scores2[i], gains2_v[i]);
            }
        }

        // Prefix cumulative for config0
        std::vector<double> cum0(L + 1, 0.0);
        for (int i = 0; i < std::min(L, W); i++) {
            cum0[i + 1] = cum0[i] + static_cast<double>(scores0[i]);
        }

        // Suffix cumulative for config2
        std::vector<double> suf2(L + 1, 0.0);
        for (int i = L - 1; i >= 2; i--) {
            if (i - 2 < W) {
                suf2[i] = suf2[i + 1] + static_cast<double>(scores2[i]);
            } else {
                suf2[i] = suf2[i + 1];
            }
        }

        // Prefix cumulative for config1
        std::vector<double> cum1(L + 1, 0.0);
        for (int i = 1; i < L; i++) {
            if (i - 1 < W) {
                cum1[i + 1] = cum1[i] + static_cast<double>(scores1[i]);
            } else {
                cum1[i + 1] = cum1[i];
            }
        }

        std::vector<float> gains_5;
        gains_5.reserve(L - 2);
        for (int m1 = 0; m1 < L; m1++) {
            for (int m2 = m1 + 1; m2 < L; m2++) {
                double seg1 = cum0[m1];

                double seg2 = 0.0;
                if (m2 > m1 + 1) {
                    seg2 = cum1[m2] - cum1[m1 + 1];
                }

                double seg3 = 0.0;
                if (m2 + 1 < L) {
                    seg3 = suf2[m2 + 1];
                }

                double aligned_score = seg1 + seg2 + seg3;

                gains_5.clear();
                for (int i = 0; i < m1; i++) {
                    if (gains0_v[i] > 1e-12f) gains_5.push_back(gains0_v[i]);
                }
                for (int i = m1 + 1; i < m2; i++) {
                    if (gains1_v[i] > 1e-12f) gains_5.push_back(gains1_v[i]);
                }
                for (int i = m2 + 1; i < L; i++) {
                    if (gains2_v[i] > 1e-12f) gains_5.push_back(gains2_v[i]);
                }

                update_best(compute_total_edits(aligned_score, gains_5, 2));
                if (!std::isnan(best_edits) && best_edits <= 2.0f) break;
            }
            if (!std::isnan(best_edits) && best_edits <= 2.0f) break;
        }
    }

    // ===== Case 6: One deletion + one insertion (W = L, k = 2) =====
    // Delete seq[d], skip motif[m]. Remaining L-1 pairs align 1:1.
    //
    // d < m: motif[0..d-1] <-> seq[0..d-1],
    //        motif[d..m-1] <-> seq[d+1..m],
    //        motif[m+1..L-1] <-> seq[m+1..L-1]
    //
    // d > m: motif[0..m-1] <-> seq[0..m-1],
    //        motif[m+1..d] <-> seq[m..d-1],
    //        motif[d+1..L-1] <-> seq[d+1..L-1]
    //
    // d == m: motif[0..d-1] <-> seq[0..d-1],
    //         motif[d+1..L-1] <-> seq[d+1..L-1]
    // Skip if best_edits <= 2
    if (std::isnan(best_edits) || best_edits > 2.0f) {
        const int W = L;

        std::vector<float> scores_ns(L), gains_ns(L);   // no shift
        std::vector<float> scores_sp1(L), gains_sp1(L); // seq shifted +1
        std::vector<float> scores_sm1(L), gains_sm1(L); // seq shifted -1

        for (int i = 0; i < L; i++) {
            int raw_ns = reverse ? (W - 1 - i) : i;
            get_base_score(i, raw_ns, scores_ns[i], gains_ns[i]);

            if (i + 1 < W) {
                int raw_sp1 = reverse ? (W - 1 - (i + 1)) : (i + 1);
                get_base_score(i, raw_sp1, scores_sp1[i], gains_sp1[i]);
            }

            if (i >= 1) {
                int raw_sm1 = reverse ? (W - 1 - (i - 1)) : (i - 1);
                get_base_score(i, raw_sm1, scores_sm1[i], gains_sm1[i]);
            }
        }

        // Prefix sums for no-shift
        std::vector<double> cum_ns(L + 1, 0.0);
        for (int i = 0; i < L; i++) {
            cum_ns[i + 1] = cum_ns[i] + static_cast<double>(scores_ns[i]);
        }

        // Prefix sums for seq+1 shift
        std::vector<double> cum_sp1(L + 1, 0.0);
        for (int i = 0; i < L - 1; i++) {
            cum_sp1[i + 1] = cum_sp1[i] + static_cast<double>(scores_sp1[i]);
        }

        // Prefix sums for seq-1 shift (valid for i >= 1)
        std::vector<double> cum_sm1(L + 1, 0.0);
        for (int i = 1; i < L; i++) {
            cum_sm1[i + 1] = cum_sm1[i] + static_cast<double>(scores_sm1[i]);
        }

        // Suffix sums for no-shift
        std::vector<double> suf_ns(L + 1, 0.0);
        for (int i = L - 1; i >= 0; i--) {
            suf_ns[i] = suf_ns[i + 1] + static_cast<double>(scores_ns[i]);
        }

        std::vector<float> gains_6;
        gains_6.reserve(L - 1);
        for (int d = 0; d < L; d++) {
            for (int m = 0; m < L; m++) {
                // d == m is always dominated by the no-indel case (Case 1):
                // it skips both seq[d] and motif[d], aligning L-2 pairs with 2 indels,
                // while Case 1 aligns all L pairs with 0 indels.
                if (d == m) continue;

                double aligned_score = 0.0;
                gains_6.clear();

                if (d < m) {
                    aligned_score += cum_ns[d];
                    for (int i = 0; i < d; i++) {
                        if (gains_ns[i] > 1e-12f) gains_6.push_back(gains_ns[i]);
                    }

                    aligned_score += cum_sp1[m] - cum_sp1[d];
                    for (int i = d; i < m; i++) {
                        if (gains_sp1[i] > 1e-12f) gains_6.push_back(gains_sp1[i]);
                    }

                    if (m + 1 < L) {
                        aligned_score += suf_ns[m + 1];
                        for (int i = m + 1; i < L; i++) {
                            if (gains_ns[i] > 1e-12f) gains_6.push_back(gains_ns[i]);
                        }
                    }
                } else {
                    // d > m (d == m is skipped above)
                    aligned_score += cum_ns[m];
                    for (int i = 0; i < m; i++) {
                        if (gains_ns[i] > 1e-12f) gains_6.push_back(gains_ns[i]);
                    }

                    aligned_score += cum_sm1[d + 1] - cum_sm1[m + 1];
                    for (int i = m + 1; i <= d; i++) {
                        if (gains_sm1[i] > 1e-12f) gains_6.push_back(gains_sm1[i]);
                    }

                    if (d + 1 < L) {
                        aligned_score += suf_ns[d + 1];
                        for (int i = d + 1; i < L; i++) {
                            if (gains_ns[i] > 1e-12f) gains_6.push_back(gains_ns[i]);
                        }
                    }
                }

                update_best(compute_total_edits(aligned_score, gains_6, 2));
                if (!std::isnan(best_edits) && best_edits <= 2.0f) break;
            }
            if (!std::isnan(best_edits) && best_edits <= 2.0f) break;
        }
    }

    // Apply max_edits cap if set
    if (!std::isnan(best_edits) && m_max_edits >= 1 && best_edits > static_cast<float>(m_max_edits)) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    return best_edits;
}
