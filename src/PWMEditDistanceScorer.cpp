#include "PWMEditDistanceScorer.h"
#include "GenomeSeqFetch.h"
#include <algorithm>
#include <cmath>
#include <cstring>
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
                                             float score_max,
                                             Direction direction)
    : GenomeSeqScorer(shared_seqfetch, extend, strand),
      m_pssm(pssm),
      m_threshold(threshold),
      m_max_edits(max_edits),
      m_max_indels(max_indels),
      m_mode(mode),
      m_direction(direction),
      m_score_min(score_min),
      m_score_max(score_max),
      m_S_max(0.0f),
      m_S_min(0.0f)
{
    precompute_tables();
}

float PWMEditDistanceScorer::compute_column_ic(int col) const {
    float entropy = 0.0f;
    for (int b = 0; b < 4; b++) {
        float log_prob = m_pssm[col].get_log_prob_from_code(b);
        if (log_prob > -1e10f && std::isfinite(log_prob)) {
            float p = std::exp(log_prob);
            if (p > 0.0f) {
                entropy -= p * std::log2(p);
            }
        }
    }
    return 2.0f - entropy;  // IC = log2(4) - H
}

void PWMEditDistanceScorer::precompute_tables()
{
    int L = m_pssm.length();
    const bool below = (m_direction == Direction::BELOW);

    m_col_max_scores.resize(L);
    m_col_min_scores.resize(L);

    // Compute column maxima, minima, S_max, and S_min
    m_S_max = 0.0f;
    m_S_min = 0.0f;
    for (int i = 0; i < L; i++) {
        float max_score = -std::numeric_limits<float>::infinity();
        float min_score = std::numeric_limits<float>::infinity();

        for (int b = 0; b < 4; b++) {
            float score = m_pssm[i].get_log_prob_from_code(b);
            if (score > max_score) max_score = score;
            if (score < min_score) min_score = score;
        }

        m_col_max_scores[i] = max_score;
        m_col_min_scores[i] = min_score;
        m_S_max += max_score;
        m_S_min += min_score;
    }

    // Precompute suffix sums for early-abandon pruning.
    // Uses target_score() which returns col_max for ABOVE, col_min for BELOW.
    m_max_suffix_score.resize(L + 1, 0.0f);
    for (int i = L - 1; i >= 0; i--) {
        m_max_suffix_score[i] = m_max_suffix_score[i + 1] + target_score(i);
    }

    // IC-ordered column processing: sort columns by information content descending
    // so early-abandon kicks in faster in compute_heuristic (subs-only mode).
    m_use_ic_order = (m_max_indels == 0 && L <= MAX_MOTIF_LEN_OPT);
    if (m_use_ic_order) {
        std::vector<std::pair<float, int>> col_ics(L);
        for (int i = 0; i < L; i++) {
            col_ics[i] = {compute_column_ic(i), i};
        }
        std::sort(col_ics.begin(), col_ics.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });
        for (int i = 0; i < L; i++) {
            m_ic_col_order[i] = col_ics[i].second;
        }
        // Suffix target scores in IC-sorted order
        m_ic_suffix_target[L] = 0.0f;
        for (int i = L - 1; i >= 0; i--) {
            m_ic_suffix_target[i] = m_ic_suffix_target[i + 1] + target_score(m_ic_col_order[i]);
        }
    }

    // Precompute maximum possible delta (gain for ABOVE, loss for BELOW) from k substitutions.
    // ABOVE: max_gain_budget[k] = sum of top k (col_max - col_min)
    // BELOW: max_gain_budget[k] = sum of top k (col_max - col_min)
    // Both directions use the same per-column spread; the interpretation changes but the
    // values are identical: col_max - col_min is the maximum delta achievable at column i.
    {
        std::vector<float> col_max_deltas(L);
        for (int i = 0; i < L; i++) {
            col_max_deltas[i] = m_col_max_scores[i] - m_col_min_scores[i];
        }
        std::sort(col_max_deltas.begin(), col_max_deltas.end(), std::greater<float>());
        m_max_gain_budget.resize(L + 1, 0.0f);
        for (int k = 0; k < L; k++) {
            m_max_gain_budget[k + 1] = m_max_gain_budget[k] + col_max_deltas[k];
        }
    }

    // Collect all delta values (gain for ABOVE, loss for BELOW).
    // ABOVE: delta[i][b] = col_max[i] - PSSM[i][b]  (score improvement from editing position i)
    // BELOW: delta[i][b] = PSSM[i][b] - col_min[i]  (score reduction from editing position i)
    std::set<float, std::greater<float>> unique_deltas;  // Sorted descending
    m_bin_index.resize(L);

    for (int i = 0; i < L; i++) {
        m_bin_index[i].resize(5);  // 5 to handle unknown bases (index 4)

        for (int b = 0; b < 4; b++) {
            float score = m_pssm[i].get_log_prob_from_code(b);
            unique_deltas.insert(compute_position_delta(score, i));
        }

        // Unknown bases (N): worst case for the direction → max possible delta
        unique_deltas.insert(m_col_max_scores[i] - m_col_min_scores[i]);
    }

    // Copy to vector (already sorted descending)
    m_gain_values.assign(unique_deltas.begin(), unique_deltas.end());

    // PERF-1: Initialize reusable count vector
    m_exact_count.resize(m_gain_values.size(), 0);
    m_exact_touched.reserve(m_gain_values.size());

    // Build bin lookup table
    for (int i = 0; i < L; i++) {
        for (int b = 0; b < 5; b++) {  // Include unknown base index 4
            float score = (b < 4)
                ? m_pssm[i].get_log_prob_from_code(b)
                : (below ? m_col_max_scores[i] : m_col_min_scores[i]);  // worst case for direction

            float delta = compute_position_delta(score, i);
            auto it = std::lower_bound(m_gain_values.begin(), m_gain_values.end(),
                                      delta, std::greater<float>());
            m_bin_index[i][b] = static_cast<uint8_t>(std::distance(m_gain_values.begin(), it));
        }
    }

    // Populate flat PSSM lookup tables for cache-friendly access.
    // For ABOVE: mandatory = log-zero → assume col_max score, gain = 0
    // For BELOW: log-zero makes score -Inf → already below any threshold → not mandatory
    for (int i = 0; i < L && i < MAX_MOTIF_LEN_OPT; i++) {
        for (int b = 0; b < 4; b++) {
            float raw = m_pssm[i].get_log_prob_from_code(b);
            bool is_logzero = (raw <= kLogZeroThreshold || !std::isfinite(raw));
            if (below) {
                // BELOW: log-zero → raw -Inf stored, adjusted_score becomes -Inf,
                // deficit ≤ 0, returns 0 edits. Loss = 0 (can't decrease below -Inf).
                m_mandatory_table[i][b] = false;
                m_score_table[i][b] = raw;
                m_gain_table[i][b] = is_logzero ? 0.0f : (raw - m_col_min_scores[i]);
            } else {
                m_mandatory_table[i][b] = is_logzero;
                m_score_table[i][b] = is_logzero ? m_col_max_scores[i] : raw;
                m_gain_table[i][b] = is_logzero ? 0.0f : (m_col_max_scores[i] - raw);
            }
        }
        // N-base (index 4): ABOVE = mandatory (worst case), BELOW = not mandatory (worst case opposite)
        m_mandatory_table[i][4] = !below;
        m_score_table[i][4] = m_col_max_scores[i];
        m_gain_table[i][4] = below ? (m_col_max_scores[i] - m_col_min_scores[i]) : 0.0f;
    }

    // Build pigeonhole pre-filter blocks (must come after m_mandatory_table is populated).
    // Pigeonhole only applies to ABOVE direction; for BELOW, the invariant
    // ("at least one block must match exactly") doesn't hold in the same way
    // since we're trying to disrupt matches, not find them.
    m_use_prefilter = false;
    m_prefilter_blocks.clear();
    if (!below) {
        int K = m_max_edits;
        // For exact mode (K < 0), K is unbounded — can't form finite block count.
        // For K == 0, we already check exact match only, no pre-filter needed.
        if (K >= 1 && L >= 3 * (K + 1) && L <= MAX_MOTIF_LEN_OPT) {
            int num_blocks = K + 1;
            m_prefilter_blocks.resize(num_blocks);

            if (m_max_indels == 0) {
                // === SUBS-ONLY: IC-sorted non-contiguous column groups ===
                std::vector<std::pair<float, int>> col_ics(L);
                for (int i = 0; i < L; i++) {
                    col_ics[i] = {compute_column_ic(i), i};
                }
                std::sort(col_ics.begin(), col_ics.end(),
                          [](const auto& a, const auto& b) { return a.first > b.first; });
                for (int b = 0; b < num_blocks; b++) {
                    m_prefilter_blocks[b].columns.clear();
                }
                int block_size = L / num_blocks;
                for (int rank = 0; rank < L; rank++) {
                    int block_idx = rank / block_size;
                    if (block_idx >= num_blocks) block_idx = num_blocks - 1;
                    m_prefilter_blocks[block_idx].columns.push_back(col_ics[rank].second);
                }
            } else {
                // === INDEL MODE: contiguous blocks (as before), stored as columns ===
                for (int b = 0; b < num_blocks; b++) {
                    int block_start = b * L / num_blocks;
                    int block_end = (b + 1) * L / num_blocks;
                    m_prefilter_blocks[b].columns.clear();
                    for (int i = block_start; i < block_end; i++) {
                        m_prefilter_blocks[b].columns.push_back(i);
                    }
                }
            }

            // Build viable tables and compute avg_ic for each block.
            // For subs-only ABOVE mode, use score-aware viability: a hash is viable
            // only if block_score + outside_col_max >= threshold. This is correct
            // because with 0 edits in the block, its score is fixed, and K edits
            // on outside columns can at best bring each to col_max.
            bool score_aware = (m_max_indels == 0 && !below);

            for (int b = 0; b < num_blocks; b++) {
                PrefilterBlock& blk = m_prefilter_blocks[b];
                int block_len = (int)blk.columns.size();
                blk.num_entries = 1 << (2 * block_len);
                blk.viable.assign(blk.num_entries, false);

                float total_ic = 0.0f;
                for (int col : blk.columns) {
                    total_ic += compute_column_ic(col);
                }
                blk.avg_ic = total_ic / block_len;

                // Compute sum of col_max for columns outside this block
                float outside_col_max = 0.0f;
                if (score_aware) {
                    std::vector<bool> in_block(L, false);
                    for (int col : blk.columns) in_block[col] = true;
                    for (int i = 0; i < L; i++) {
                        if (!in_block[i]) outside_col_max += m_col_max_scores[i];
                    }
                }

                for (int h = 0; h < blk.num_entries; h++) {
                    bool ok = true;
                    float block_score = 0.0f;
                    for (int j = 0; j < block_len && ok; j++) {
                        int base = (h >> (2 * j)) & 3;
                        if (m_mandatory_table[blk.columns[j]][base]) {
                            ok = false;
                        } else if (score_aware) {
                            block_score += m_score_table[blk.columns[j]][base];
                        }
                    }
                    if (ok && score_aware) {
                        ok = (block_score + outside_col_max >= m_threshold);
                    }
                    blk.viable[h] = ok;
                }
            }

            // Sort blocks by avg_ic descending (highest IC checked first)
            std::sort(m_prefilter_blocks.begin(), m_prefilter_blocks.end(),
                      [](const PrefilterBlock& a, const PrefilterBlock& b) {
                          return a.avg_ic > b.avg_ic;
                      });

            m_use_prefilter = true;
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

float PWMEditDistanceScorer::compute_exact(const int* bidx, bool reverse)
{
    const int L = m_pssm.length();
    if (L == 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    if (!is_globally_reachable()) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    int mandatory_edits = 0;
    double adjusted_score = 0.0;

    // PERF-1: Use class-member count vector with touched-list cleanup
    m_exact_touched.clear();

    for (int i = 0; i < L; i++) {
        int seq_idx = reverse ? (L - 1 - i) : i;
        int base_idx = bidx[seq_idx];

        if (m_mandatory_table[i][base_idx]) {
            mandatory_edits++;
            adjusted_score += static_cast<double>(m_score_table[i][base_idx]);
            continue;
        }

        adjusted_score += static_cast<double>(m_score_table[i][base_idx]);

        int bin = m_bin_index[i][base_idx];
        if (m_exact_count[bin] == 0) {
            m_exact_touched.push_back(bin);
        }
        m_exact_count[bin]++;
    }

    double deficit = compute_deficit(adjusted_score);

    auto cleanup = [this]() {
        for (size_t idx : m_exact_touched) {
            m_exact_count[idx] = 0;
        }
    };

    if (deficit <= 0.0) {
        cleanup();
        return static_cast<float>(mandatory_edits);
    }

    if (max_possible_delta(adjusted_score) < deficit) {
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

        double delta_val = static_cast<double>(m_gain_values[j]);
        if (delta_val <= 0.0) {
            continue;
        }

        double cap = delta_val * static_cast<double>(bin_count);

        if (acc + cap < deficit) {
            acc += cap;
            edits += bin_count;
        } else {
            double remaining = deficit - acc;
            if (remaining <= 0.0) {
                cleanup();
                return static_cast<float>(edits);
            }

            double need_raw = remaining / delta_val;
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

float PWMEditDistanceScorer::compute_heuristic(const int* bidx, bool reverse, int max_k)
{
    const int L = m_pssm.length();
    if (L == 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    if (!is_globally_reachable()) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    const bool below = is_below();

    int mandatory_edits = 0;
    double adjusted_score = 0.0;
    m_heur_deltas.clear();

    // Flat top-K delta tracker: for K <= 3, a sorted array of K floats is faster
    // than std::push_heap/pop_heap. For larger K, falls back to insertion sort
    // which is still efficient for small K.
    double top_deltas_sum = 0.0;
    float top_k[4] = {};  // up to K=4 supported; stack-allocated, no heap
    int n_top_k = 0;

    // Choose column processing order and matching suffix table:
    // - IC-sorted order (high-IC first) for subs-only: early-abandon kicks in faster
    // - Positional order for indel mode (alignment-dependent)
    const int* col_order = m_use_ic_order ? m_ic_col_order : nullptr;
    const float* suffix_target = m_use_ic_order ? m_ic_suffix_target : m_max_suffix_score.data();
    const int capped_k = std::min(max_k, 4);

    for (int i = 0; i < L; i++) {
        int col = col_order ? col_order[i] : i;
        int seq_idx = reverse ? (L - 1 - col) : col;
        int base_idx = bidx[seq_idx];

        // Use precomputed lookup tables (complement already applied in bidx)
        if (m_mandatory_table[col][base_idx]) {
            mandatory_edits++;
            adjusted_score += static_cast<double>(m_score_table[col][base_idx]);
        } else {
            float base_score = m_score_table[col][base_idx];
            adjusted_score += static_cast<double>(base_score);
            float delta = m_gain_table[col][base_idx];
            m_heur_deltas.push_back(delta);

            // Update flat top-K tracker (insertion sort, K elements max)
            if (delta > 0.0f) {
                if (n_top_k < capped_k) {
                    // Still filling: insert in sorted position (descending)
                    int pos = n_top_k;
                    while (pos > 0 && top_k[pos - 1] < delta) {
                        top_k[pos] = top_k[pos - 1];
                        pos--;
                    }
                    top_k[pos] = delta;
                    n_top_k++;
                    top_deltas_sum += static_cast<double>(delta);
                } else if (delta > top_k[n_top_k - 1]) {
                    // Replace smallest in top-K
                    top_deltas_sum -= static_cast<double>(top_k[n_top_k - 1]);
                    int pos = n_top_k - 1;
                    while (pos > 0 && top_k[pos - 1] < delta) {
                        top_k[pos] = top_k[pos - 1];
                        pos--;
                    }
                    top_k[pos] = delta;
                    top_deltas_sum += static_cast<double>(delta);
                }
            }
        }

        // Column-by-column suffix-bound early-abandon:
        // ABOVE: adjusted_score + top_k_gains + suffix_max < threshold → unreachable
        // BELOW: adjusted_score - top_k_losses + suffix_min > threshold → unreachable
        if (below) {
            if (adjusted_score - top_deltas_sum + static_cast<double>(suffix_target[i + 1])
                > static_cast<double>(m_threshold)) {
                return std::numeric_limits<float>::quiet_NaN();
            }
        } else {
            if (adjusted_score + top_deltas_sum + static_cast<double>(suffix_target[i + 1])
                < static_cast<double>(m_threshold)) {
                return std::numeric_limits<float>::quiet_NaN();
            }
        }
    }

    if (mandatory_edits > max_k) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    double deficit = compute_deficit(adjusted_score);
    if (deficit <= 0.0) {
        return static_cast<float>(mandatory_edits);
    }

    if (max_possible_delta(adjusted_score) < deficit) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    int remaining_budget = max_k - mandatory_edits;
    if (remaining_budget <= 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    if (remaining_budget > static_cast<int>(m_heur_deltas.size())) {
        remaining_budget = static_cast<int>(m_heur_deltas.size());
    }

    if (remaining_budget <= 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    if (remaining_budget < static_cast<int>(m_heur_deltas.size())) {
        std::nth_element(m_heur_deltas.begin(), m_heur_deltas.begin() + remaining_budget,
                         m_heur_deltas.end(), std::greater<float>());
    }
    std::sort(m_heur_deltas.begin(), m_heur_deltas.begin() + remaining_budget,
              std::greater<float>());

    double acc = 0.0;
    for (int i = 0; i < remaining_budget; i++) {
        acc += static_cast<double>(m_heur_deltas[i]);
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

bool PWMEditDistanceScorer::passes_prefilter(const int* bidx, int seq_avail, bool reverse) const
{
    const int L = m_pssm.length();

    for (const auto& block : m_prefilter_blocks) {
        int block_len = (int)block.columns.size();

        for (int shift = -m_max_indels; shift <= m_max_indels; shift++) {
            int hash = 0;
            bool valid = true;

            for (int j = 0; j < block_len && valid; j++) {
                int motif_col = block.columns[j];
                int seq_idx;
                if (!reverse) {
                    seq_idx = motif_col + shift;
                } else {
                    seq_idx = L - 1 - motif_col + shift;
                }

                if (seq_idx < 0 || seq_idx >= seq_avail) {
                    valid = false;
                    break;
                }

                int b = bidx[seq_idx];
                if (b >= 4) {
                    valid = false;
                    break;
                }

                hash += b << (2 * j);
            }

            if (valid && block.viable[hash]) {
                return true;
            }
        }
    }

    return false;
}

float PWMEditDistanceScorer::compute_window_edits(const int* bidx, int seq_avail, bool reverse)
{
    // When indels are enabled, use specialized solvers.
    // Apply pigeonhole pre-filter first to skip windows that provably cannot
    // match within the edit budget.
    if (m_max_indels == 1) {
        if (m_use_prefilter && !passes_prefilter(bidx, seq_avail, reverse)) {
            return std::numeric_limits<float>::quiet_NaN();
        }
        return compute_with_one_indel(bidx, seq_avail, reverse);
    } else if (m_max_indels == 2) {
        if (m_use_prefilter && !passes_prefilter(bidx, seq_avail, reverse)) {
            return std::numeric_limits<float>::quiet_NaN();
        }
        return compute_with_two_indels(bidx, seq_avail, reverse);
    } else if (m_max_indels > 2) {
        if (m_use_prefilter && !passes_prefilter(bidx, seq_avail, reverse)) {
            return std::numeric_limits<float>::quiet_NaN();
        }
        return compute_with_indels(bidx, seq_avail, reverse);
    }

    // Fast path: substitutions only
    if (m_max_edits < 0) {
        return compute_exact(bidx, reverse);
    }
    // Pigeonhole pre-filter also helps pure substitution mode (shifts = {0})
    if (m_use_prefilter && !passes_prefilter(bidx, seq_avail, reverse)) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    return compute_heuristic(bidx, reverse, m_max_edits);
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

    // Precompute base indices for both strands to eliminate per-call
    // complement_base() and base_to_index() overhead from the inner loop.
    // fwd_bidx[j] = base_to_index(seq[j]) — forward strand
    // rev_bidx[j] = base_to_index(complement_base(seq[j])) — reverse strand
    std::vector<int> fwd_bidx(target_length);
    std::vector<int> rev_bidx(target_length);
    for (size_t j = 0; j < target_length; j++) {
        char base = seq_data[j];
        fwd_bidx[j] = base_to_index(base);
        rev_bidx[j] = base_to_index(complement_base(base));
    }

    // Sliding-window N-count for fast skip of N-heavy regions.
    // N bases force mandatory edits; if a window has more than max_edits Ns,
    // it's unreachable. We maintain a running count, O(1) per window step.
    bool use_n_skip = (m_max_edits >= 0 && max_start > 0);
    int n_count = 0;
    if (use_n_skip) {
        for (size_t j = 0; j < motif_length && j < target_length; j++) {
            n_count += (fwd_bidx[j] >= 4) ? 1 : 0;
        }
    }

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
        // Sliding N-count maintenance (for offset > 0, slide the window by 1)
        if (use_n_skip && offset > 0) {
            // Remove the base that just left the window
            n_count -= (fwd_bidx[offset - 1] >= 4) ? 1 : 0;
            // Add the base that just entered the window
            if (offset + motif_length - 1 < target_length) {
                n_count += (fwd_bidx[offset + motif_length - 1] >= 4) ? 1 : 0;
            }
        }

        // Skip windows with too many N-bases (each N forces a mandatory edit)
        if (use_n_skip && n_count > m_max_edits) {
            continue;
        }

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

                if (pass_score_filter) {
                    const int* bidx = fwd_bidx.data() + offset;
                    float edits = compute_window_edits(bidx, seq_avail, /*reverse=*/false);
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

                if (pass_score_filter) {
                    const int* bidx = rev_bidx.data() + offset;
                    float edits = compute_window_edits(bidx, seq_avail, /*reverse=*/true);
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
            int seq_avail = static_cast<int>(target_length - metrics.best_pwm_index);
            bool is_reverse = metrics.best_pwm_direction < 0;
            const int* bidx = is_reverse
                ? rev_bidx.data() + metrics.best_pwm_index
                : fwd_bidx.data() + metrics.best_pwm_index;
            metrics.best_pwm_edits = compute_window_edits(bidx, seq_avail, is_reverse);
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

float PWMEditDistanceScorer::compute_with_indels(const int* bidx_arr, int seq_len, bool reverse)
{
    const int L = m_pssm.length();
    if (L == 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    const int D = m_max_indels;
    const bool below = (m_direction == Direction::BELOW);
    float best_edits = std::numeric_limits<float>::quiet_NaN();

    // For each sequence window length W in [L-D, L+D], align the motif (length L)
    // against W sequence bases using a 3D DP:
    //   ABOVE: dp[i][j][k] = max PWM score aligning motif[0..i-1] with seq[0..j-1]
    //   BELOW: dp[i][j][k] = min PWM score aligning motif[0..i-1] with seq[0..j-1]
    //                 using exactly k indels (insertions + deletions)
    //
    // Band constraint: |i - j| <= D (prevents needing > D indels)
    //
    // For each final state dp[L][W][k], the score tells us the PWM log-likelihood
    // of the aligned positions with their actual sequence bases. To reach the
    // threshold, we may need additional substitutions: each substitution replaces
    // a base with the column-optimal base for the given direction.
    // ABOVE: gain = col_max - current (raise score toward threshold)
    // BELOW: gain = current - col_min (lower score toward threshold)
    // We traceback to find aligned positions and greedily pick the largest gains.
    //
    // Total edits = k (indels) + subs_needed.
    // We minimize this across all W and k.

    // DP initialization sentinel: -inf for ABOVE (maximize), +inf for BELOW (minimize)
    const double dp_sentinel = below
        ? std::numeric_limits<double>::infinity()
        : -std::numeric_limits<double>::infinity();

    for (int W = std::max(1, L - D); W <= L + D; ++W) {
        if (W > seq_len) {
            break;
        }

        const int rows = L + 1;
        const int cols = W + 1;
        const int indel_levels = D + 1;  // k = 0, 1, ..., D

        // Flattened 3D DP table
        std::vector<double> dp(rows * cols * indel_levels, dp_sentinel);

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
                // Get precomputed base index at position j-1
                int seq_idx = reverse ? (W - 1 - (j - 1)) : (j - 1);
                int bi = bidx_arr[seq_idx];

                // Use precomputed score from PSSM (complement already applied in bidx_arr).
                // For the generic DP we need the raw PSSM score (not the mandatory-adjusted one),
                // because the DP traceback handles gains separately.
                float base_score;
                if (bi == 4) {
                    // Unknown base: worst case for the direction
                    // ABOVE: minimum score (hardest to reach threshold from above)
                    // BELOW: maximum score (hardest to get below threshold)
                    if (below) {
                        float max_s = -std::numeric_limits<float>::infinity();
                        for (int b = 0; b < 4; b++) {
                            float s = m_pssm[i - 1].get_log_prob_from_code(b);
                            if (s > max_s) max_s = s;
                        }
                        base_score = max_s;
                    } else {
                        float min_s = std::numeric_limits<float>::infinity();
                        for (int b = 0; b < 4; b++) {
                            float s = m_pssm[i - 1].get_log_prob_from_code(b);
                            if (s < min_s) min_s = s;
                        }
                        base_score = min_s;
                    }
                } else {
                    base_score = m_pssm[i - 1].get_log_prob_from_code(bi);
                }

                for (int k = 0; k <= D; ++k) {
                    // 1. Match/Substitution (diagonal): align motif[i-1] with seq[j-1]
                    if (std::abs((i - 1) - (j - 1)) <= D) {
                        double prev = dp[idx3(i - 1, j - 1, k)];
                        if (below) {
                            // BELOW: minimize score
                            if (prev < std::numeric_limits<double>::infinity() * 0.5) {
                                double new_score = prev + static_cast<double>(base_score);
                                if (new_score < dp[idx3(i, j, k)]) {
                                    dp[idx3(i, j, k)] = new_score;
                                }
                            }
                        } else {
                            // ABOVE: maximize score
                            if (prev > -std::numeric_limits<double>::infinity() * 0.5) {
                                double new_score = prev + static_cast<double>(base_score);
                                if (new_score > dp[idx3(i, j, k)]) {
                                    dp[idx3(i, j, k)] = new_score;
                                }
                            }
                        }
                    }

                    if (k < D) {
                        // 2. Insertion: skip motif[i-1], advance motif not sequence
                        if (std::abs((i - 1) - j) <= D) {
                            double prev = dp[idx3(i - 1, j, k)];
                            if (below) {
                                if (prev < std::numeric_limits<double>::infinity() * 0.5) {
                                    if (prev < dp[idx3(i, j, k + 1)]) {
                                        dp[idx3(i, j, k + 1)] = prev;
                                    }
                                }
                            } else {
                                if (prev > -std::numeric_limits<double>::infinity() * 0.5) {
                                    if (prev > dp[idx3(i, j, k + 1)]) {
                                        dp[idx3(i, j, k + 1)] = prev;
                                    }
                                }
                            }
                        }

                        // 3. Deletion: skip seq[j-1], advance sequence not motif
                        if (std::abs(i - (j - 1)) <= D) {
                            double prev = dp[idx3(i, j - 1, k)];
                            if (below) {
                                if (prev < std::numeric_limits<double>::infinity() * 0.5) {
                                    if (prev < dp[idx3(i, j, k + 1)]) {
                                        dp[idx3(i, j, k + 1)] = prev;
                                    }
                                }
                            } else {
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
        }

        // Extract results for each indel count k
        for (int k = 0; k <= D; ++k) {
            double score = dp[idx3(L, W, k)];
            // Check for sentinel (unreachable state)
            bool is_sentinel = below
                ? (score >= std::numeric_limits<double>::infinity() * 0.5)
                : (score <= -std::numeric_limits<double>::infinity() * 0.5);
            if (is_sentinel) {
                continue;
            }

            // If score already satisfies threshold, no substitutions needed
            if (threshold_satisfied(score)) {
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

                // Helper: check if a DP value is a valid (non-sentinel) state
                auto is_valid = [&](double val) -> bool {
                    return below
                        ? (val < std::numeric_limits<double>::infinity() * 0.5)
                        : (val > -std::numeric_limits<double>::infinity() * 0.5);
                };

                // Try diagonal (match/substitution) first
                bool found = false;
                if (std::abs((ti - 1) - (tj - 1)) <= D) {
                    int s_idx = reverse ? (W - 1 - (tj - 1)) : (tj - 1);
                    int bi = bidx_arr[s_idx];

                    float bs;
                    if (bi == 4) {
                        // Unknown base: use same logic as DP fill
                        if (below) {
                            float max_s = -std::numeric_limits<float>::infinity();
                            for (int bb = 0; bb < 4; bb++) {
                                float s = m_pssm[ti - 1].get_log_prob_from_code(bb);
                                if (s > max_s) max_s = s;
                            }
                            bs = max_s;
                        } else {
                            float min_s = std::numeric_limits<float>::infinity();
                            for (int bb = 0; bb < 4; bb++) {
                                float s = m_pssm[ti - 1].get_log_prob_from_code(bb);
                                if (s < min_s) min_s = s;
                            }
                            bs = min_s;
                        }
                    } else {
                        bs = m_pssm[ti - 1].get_log_prob_from_code(bi);
                    }

                    double prev = dp[idx3(ti - 1, tj - 1, tk)];
                    if (is_valid(prev) &&
                        std::fabs((prev + static_cast<double>(bs)) - cur) < 1e-9 * std::max(1.0, std::fabs(cur))) {
                        // ABOVE: gain = col_max - bs (how much we can raise the score)
                        // BELOW: gain = bs - col_min (how much we can lower the score)
                        float gain = below
                            ? (bs - m_col_min_scores[ti - 1])
                            : (m_col_max_scores[ti - 1] - bs);
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
                    if (is_valid(prev) &&
                        std::fabs(prev - cur) < 1e-9 * std::max(1.0, std::fabs(cur))) {
                        ti--; tk--;
                        found = true;
                    }
                }

                if (!found && tj > 0 && tk > 0 && std::abs(ti - (tj - 1)) <= D) {
                    // Try deletion
                    double prev = dp[idx3(ti, tj - 1, tk - 1)];
                    if (is_valid(prev) &&
                        std::fabs(prev - cur) < 1e-9 * std::max(1.0, std::fabs(cur))) {
                        tj--; tk--;
                        found = true;
                    }
                }

                if (!found) break;  // DP inconsistency (shouldn't happen)
            }

            double deficit = compute_deficit(score);
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

bool PWMEditDistanceScorer::early_abandon_banded_dp(const char* seq_ptr, int seq_len, bool reverse) const
{
    // Banded DP early-abandon filter for indel-enabled windows.
    //
    // Key insight: the edit distance metric allows substitutions at any matched
    // position. Each substitution brings a column's score to its optimal value
    // for the given direction.
    //
    // ABOVE: each matched column achieves col_max (upper-bound on score).
    //   dp computes max achievable score. Abandon if row_max + suffix_max < threshold.
    // BELOW: each matched column achieves col_min (lower-bound on score).
    //   dp computes min achievable score. Abandon if row_min + suffix_min > threshold.
    //
    // Band constraint: |i - j| <= D.
    //
    // For insertion-heavy families (many motif columns skipped), the achievable
    // score changes because fewer columns contribute. This is where the filter
    // provides value.

    const int L = m_pssm.length();
    const int D = m_max_indels;
    const bool below = (m_direction == Direction::BELOW);

    if (L == 0 || D < 0) return false;

    // Stack allocation limits
    static constexpr int MAX_MOTIF_LEN = 40;
    static constexpr int MAX_INDELS = 3;
    if (L > MAX_MOTIF_LEN || D > MAX_INDELS) return false;

    const int max_j = std::min(seq_len, L + D);
    if (max_j <= 0) return false;

    static constexpr int MAX_J = MAX_MOTIF_LEN + MAX_INDELS + 1;
    static constexpr int KD = MAX_INDELS + 1;

    double dp_prev[MAX_J][KD];
    double dp_cur[MAX_J][KD];

    // ABOVE: sentinel = -inf (maximizing), BELOW: sentinel = +inf (minimizing)
    const double SENTINEL = below ? 1e300 : -1e300;
    const double threshold_d = static_cast<double>(m_threshold);

    // Initialize dp_prev (row 0) to sentinel
    for (int j = 0; j <= max_j; j++) {
        for (int k = 0; k <= D; k++) {
            dp_prev[j][k] = SENTINEL;
        }
    }

    // Row 0: dp[0][0][0] = 0 (empty alignment)
    dp_prev[0][0] = 0.0;
    // dp[0][j][j] = 0 for j=1..D (skip j sequence bases = j deletions)
    for (int j = 1; j <= std::min(D, max_j); j++) {
        dp_prev[j][j] = 0.0;
    }

    // Early abandon check for row 0
    // ABOVE: 0 + suffix_max < threshold => skip
    // BELOW: 0 + suffix_min > threshold => skip
    if (below) {
        if (0.0 + static_cast<double>(m_max_suffix_score[0]) > threshold_d) {
            return true;
        }
    } else {
        if (0.0 + static_cast<double>(m_max_suffix_score[0]) < threshold_d) {
            return true;
        }
    }

    // Row-by-row fill
    for (int i = 1; i <= L; i++) {
        const int j_lo = std::max(0, i - D);
        const int j_hi = std::min(max_j, i + D);

        // Initialize dp_cur for this row's range
        for (int j = j_lo; j <= j_hi; j++) {
            for (int k = 0; k <= D; k++) {
                dp_cur[j][k] = SENTINEL;
            }
        }

        // j=0: only insertion (skip motif[i-1]) is possible
        if (j_lo == 0) {
            for (int k = 1; k <= D; k++) {
                double prev = dp_prev[0][k - 1];
                bool is_valid = below ? (prev < SENTINEL) : (prev > SENTINEL);
                bool is_better = below ? (prev < dp_cur[0][k]) : (prev > dp_cur[0][k]);
                if (is_valid && is_better) {
                    dp_cur[0][k] = prev;
                }
            }
        }

        // ABOVE: use col_max (substitutions can raise score to maximum)
        // BELOW: use col_min (substitutions can lower score to minimum)
        const double col_score = below
            ? static_cast<double>(m_col_min_scores[i - 1])
            : static_cast<double>(m_col_max_scores[i - 1]);

        for (int j = std::max(1, j_lo); j <= j_hi; j++) {
            for (int k = 0; k <= D; k++) {
                double val = SENTINEL;

                // 1. Match (diagonal): dp_prev[j-1][k] + col_score
                //    Substitution is free in edit-distance terms (counted separately).
                if (std::abs((i - 1) - (j - 1)) <= D) {
                    double prev = dp_prev[j - 1][k];
                    bool is_valid = below ? (prev < SENTINEL) : (prev > SENTINEL);
                    if (is_valid) {
                        double cand = prev + col_score;
                        bool is_better = below ? (cand < val) : (cand > val);
                        if (is_better) val = cand;
                    }
                }

                // 2. Insertion (skip motif[i-1]): dp_prev[j][k-1]
                if (k >= 1 && std::abs((i - 1) - j) <= D) {
                    double prev = dp_prev[j][k - 1];
                    bool is_valid = below ? (prev < SENTINEL) : (prev > SENTINEL);
                    bool is_better = below ? (prev < val) : (prev > val);
                    if (is_valid && is_better) {
                        val = prev;
                    }
                }

                // 3. Deletion (skip seq[j-1]): dp_cur[j-1][k-1]
                if (k >= 1) {
                    double prev = dp_cur[j - 1][k - 1];
                    bool is_valid = below ? (prev < SENTINEL) : (prev > SENTINEL);
                    bool is_better = below ? (prev < val) : (prev > val);
                    if (is_valid && is_better) {
                        val = prev;
                    }
                }

                bool val_better = below ? (val < dp_cur[j][k]) : (val > dp_cur[j][k]);
                if (val_better) {
                    dp_cur[j][k] = val;
                }
            }
        }

        // Early abandon: find row optimum across all valid (j, k)
        // ABOVE: row_max + suffix_max < threshold => skip
        // BELOW: row_min + suffix_min > threshold => skip
        double row_opt = SENTINEL;
        for (int j = j_lo; j <= j_hi; j++) {
            for (int k = 0; k <= D; k++) {
                bool is_better = below
                    ? (dp_cur[j][k] < row_opt)
                    : (dp_cur[j][k] > row_opt);
                if (is_better) {
                    row_opt = dp_cur[j][k];
                }
            }
        }

        if (below) {
            if (row_opt + static_cast<double>(m_max_suffix_score[i]) > threshold_d) {
                return true;
            }
        } else {
            if (row_opt + static_cast<double>(m_max_suffix_score[i]) < threshold_d) {
                return true;
            }
        }

        // Swap rows
        for (int j = 0; j <= max_j; j++) {
            for (int k = 0; k <= D; k++) {
                dp_prev[j][k] = SENTINEL;
            }
        }
        for (int j = j_lo; j <= j_hi; j++) {
            for (int k = 0; k <= D; k++) {
                dp_prev[j][k] = dp_cur[j][k];
            }
        }
    }

    return false;
}

inline bool PWMEditDistanceScorer::get_aligned_base_score(const char* seq_ptr, bool reverse,
                                                          int motif_pos, int raw_seq_idx,
                                                          float& out_score, float& out_gain) const
{
    char base = seq_ptr[raw_seq_idx];
    if (reverse) {
        base = complement_base(base);
    }
    int bidx = base_to_index(base);
    out_score = m_score_table[motif_pos][bidx];
    out_gain = m_gain_table[motif_pos][bidx];
    return m_mandatory_table[motif_pos][bidx];
}

bool PWMEditDistanceScorer::quick_deficit_check(double aligned_score, int indels, float best_edits_so_far) const
{
    // O(1) check: can this alignment family possibly produce a result better than best_edits_so_far?
    // Uses precomputed per-PSSM max gain budget to avoid expensive gains collection + sorting.

    // Prune: if indels alone already can't beat best_edits, skip
    if (!std::isnan(best_edits_so_far) && static_cast<float>(indels) >= best_edits_so_far) {
        return false;
    }

    // Prune: if indels alone exceed budget, skip
    if (m_max_edits >= 0 && indels > m_max_edits) {
        return false;
    }

    double deficit = compute_deficit(aligned_score);
    if (deficit <= 0.0) {
        return true;  // Already at or past threshold
    }

    // Compute max substitution budget
    int max_subs = m_pssm.length();  // upper bound
    if (!std::isnan(best_edits_so_far)) {
        int ceiling = static_cast<int>(best_edits_so_far) - 1 - indels;
        if (ceiling < max_subs) max_subs = ceiling;
    }
    if (m_max_edits >= 0) {
        int budget_subs = m_max_edits - indels;
        if (budget_subs < max_subs) max_subs = budget_subs;
    }
    if (max_subs <= 0) {
        return false;
    }

    // Quick check using precomputed per-PSSM max gain budget
    int budget_idx = std::min(max_subs, static_cast<int>(m_max_gain_budget.size()) - 1);
    if (budget_idx <= 0) {
        return false;
    }
    return static_cast<double>(m_max_gain_budget[budget_idx]) >= deficit;
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

    double deficit = compute_deficit(aligned_score);
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

float PWMEditDistanceScorer::compute_with_one_indel(const int* bidx_arr, int seq_len, bool reverse)
{
    const int L = m_pssm.length();
    if (L == 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    // Fall back to generic DP solver for motifs exceeding stack array limit
    if (L > MAX_MOTIF_LEN_OPT) {
        return compute_with_indels(bidx_arr, seq_len, reverse);
    }

    float best_edits = std::numeric_limits<float>::quiet_NaN();

    // get_base_score: direct table lookup using precomputed base indices.
    // Returns true if the position requires a mandatory substitution.
    auto get_base_score = [&](int motif_pos, int raw_seq_idx,
                               float& out_score, float& out_gain) -> bool {
        int bi = bidx_arr[raw_seq_idx];
        out_score = m_score_table[motif_pos][bi];
        out_gain = m_gain_table[motif_pos][bi];
        return m_mandatory_table[motif_pos][bi];
    };

    // compute_total_edits: pass indels + mandatory_subs as the base fixed-edit count.
    auto compute_total_edits = [&](double aligned_score, std::vector<float>& gains,
                                   int indels, int mandatory_subs) -> float {
        return compute_min_edits_from_gains(aligned_score, gains,
                                            indels + mandatory_subs, best_edits);
    };

    // ===== Case A: No indel (W = L, k = 0) =====
    // Motif position i aligns with sequence position i (or L-1-i if reverse).
    // This is equivalent to what the DP computes for W=L, k=0.
    {
        double aligned_score = 0.0;
        float gains_arr[MAX_MOTIF_LEN_OPT];
        int n_gains = 0;
        int mandatory_subs = 0;

        for (int i = 0; i < L; i++) {
            int raw_seq_idx = reverse ? (L - 1 - i) : i;
            float bs, gain;
            bool mandatory = get_base_score(i, raw_seq_idx, bs, gain);
            aligned_score += static_cast<double>(bs);
            if (mandatory) {
                mandatory_subs++;
            } else if (gain > 1e-12f) {
                gains_arr[n_gains++] = gain;
            }
        }

        std::vector<float> gains(gains_arr, gains_arr + n_gains);
        float total = compute_total_edits(aligned_score, gains, 0, mandatory_subs);
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

        // Precompute base scores, gains, and mandatory flags for all W sequence
        // positions aligned with motif columns in two configurations:
        //   prefix: motif[i] aligned with seq_raw[i] for i = 0..L-1
        //   suffix: motif[i] aligned with seq_raw[i+1] for i = 0..L-1
        // where seq_raw[j] = reverse ? (W-1-j) : j

        float prefix_scores[MAX_MOTIF_LEN_OPT], prefix_gains[MAX_MOTIF_LEN_OPT];
        float suffix_scores[MAX_MOTIF_LEN_OPT], suffix_gains[MAX_MOTIF_LEN_OPT];
        bool prefix_mandatory[MAX_MOTIF_LEN_OPT], suffix_mandatory[MAX_MOTIF_LEN_OPT];
        memset(prefix_mandatory, 0, sizeof(bool) * L);
        memset(suffix_mandatory, 0, sizeof(bool) * L);

        for (int i = 0; i < L; i++) {
            int raw_prefix = reverse ? (W - 1 - i) : i;
            prefix_mandatory[i] = get_base_score(i, raw_prefix,
                                                  prefix_scores[i], prefix_gains[i]);

            int raw_suffix = reverse ? (W - 1 - (i + 1)) : (i + 1);
            suffix_mandatory[i] = get_base_score(i, raw_suffix,
                                                  suffix_scores[i], suffix_gains[i]);
        }

        // Prefix sums for prefix alignment (motif[0..d-1] with seq[0..d-1])
        double prefix_cum_score[MAX_MOTIF_LEN_OPT + 1];
        int prefix_cum_mandatory[MAX_MOTIF_LEN_OPT + 1];
        prefix_cum_score[0] = 0.0;
        prefix_cum_mandatory[0] = 0;
        for (int i = 0; i < L; i++) {
            prefix_cum_score[i + 1] = prefix_cum_score[i] + static_cast<double>(prefix_scores[i]);
            prefix_cum_mandatory[i + 1] = prefix_cum_mandatory[i] + (prefix_mandatory[i] ? 1 : 0);
        }

        // Suffix sums for shifted alignment (motif[d..L-1] with seq[d+1..L])
        double suffix_cum_score[MAX_MOTIF_LEN_OPT + 1];
        int suffix_cum_mandatory[MAX_MOTIF_LEN_OPT + 1];
        suffix_cum_score[L] = 0.0;
        suffix_cum_mandatory[L] = 0;
        for (int i = L - 1; i >= 0; i--) {
            suffix_cum_score[i] = suffix_cum_score[i + 1] + static_cast<double>(suffix_scores[i]);
            suffix_cum_mandatory[i] = suffix_cum_mandatory[i + 1] + (suffix_mandatory[i] ? 1 : 0);
        }

        // Try each deletion position d in [0, L]
        // d=0: skip seq[0], motif[0..L-1] aligns with seq[1..L]
        // d=L: skip seq[L], motif[0..L-1] aligns with seq[0..L-1]
        float gains_arr[MAX_MOTIF_LEN_OPT];
        for (int d = 0; d <= L; d++) {
            double aligned_score = prefix_cum_score[d] + suffix_cum_score[d];

            // Quick O(1) deficit check: skip expensive gains collection if hopeless
            if (!quick_deficit_check(aligned_score, 1, best_edits)) continue;

            int mandatory_subs = prefix_cum_mandatory[d] + suffix_cum_mandatory[d];

            // Collect gains from aligned positions (skip mandatory positions)
            int n_gains = 0;
            for (int i = 0; i < d; i++) {
                if (!prefix_mandatory[i] && prefix_gains[i] > 1e-12f) {
                    gains_arr[n_gains++] = prefix_gains[i];
                }
            }
            for (int i = d; i < L; i++) {
                if (!suffix_mandatory[i] && suffix_gains[i] > 1e-12f) {
                    gains_arr[n_gains++] = suffix_gains[i];
                }
            }

            std::vector<float> gains(gains_arr, gains_arr + n_gains);
            float total = compute_total_edits(aligned_score, gains, 1, mandatory_subs);
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

        // Precompute base scores, gains, and mandatory flags for aligned columns:
        //   prefix: motif[i] aligned with seq_raw[i] for i = 0..L-2
        //   suffix: motif[i] aligned with seq_raw[i-1] for i = 1..L-1
        // where seq_raw[j] = reverse ? (W-1-j) : j

        float prefix_scores[MAX_MOTIF_LEN_OPT], prefix_gains[MAX_MOTIF_LEN_OPT];
        float suffix_scores[MAX_MOTIF_LEN_OPT], suffix_gains[MAX_MOTIF_LEN_OPT];
        bool prefix_mandatory[MAX_MOTIF_LEN_OPT], suffix_mandatory[MAX_MOTIF_LEN_OPT];
        memset(prefix_mandatory, 0, sizeof(bool) * (L - 1));
        memset(suffix_mandatory, 0, sizeof(bool) * L);

        for (int i = 0; i < L - 1; i++) {
            int raw_idx = reverse ? (W - 1 - i) : i;
            prefix_mandatory[i] = get_base_score(i, raw_idx,
                                                  prefix_scores[i], prefix_gains[i]);
        }

        // suffix: motif[i] aligns with seq[i-1] for i=1..L-1
        for (int i = 1; i < L; i++) {
            int raw_idx = reverse ? (W - 1 - (i - 1)) : (i - 1);
            suffix_mandatory[i] = get_base_score(i, raw_idx,
                                                  suffix_scores[i], suffix_gains[i]);
        }

        // Prefix cumulative sums (motif[0..m-1] with seq[0..m-1])
        double prefix_cum_score[MAX_MOTIF_LEN_OPT];
        int prefix_cum_mandatory[MAX_MOTIF_LEN_OPT];
        prefix_cum_score[0] = 0.0;
        prefix_cum_mandatory[0] = 0;
        for (int i = 0; i < L - 1; i++) {
            prefix_cum_score[i + 1] = prefix_cum_score[i] + static_cast<double>(prefix_scores[i]);
            prefix_cum_mandatory[i + 1] = prefix_cum_mandatory[i] + (prefix_mandatory[i] ? 1 : 0);
        }

        // Suffix cumulative sums (motif[m+1..L-1] with seq[m..L-2])
        double suffix_cum_score[MAX_MOTIF_LEN_OPT];
        int suffix_cum_mandatory[MAX_MOTIF_LEN_OPT];
        suffix_cum_score[L - 1] = 0.0;
        suffix_cum_mandatory[L - 1] = 0;
        for (int i = L - 1; i >= 1; i--) {
            suffix_cum_score[i - 1] = suffix_cum_score[i] + static_cast<double>(suffix_scores[i]);
            suffix_cum_mandatory[i - 1] = suffix_cum_mandatory[i] + (suffix_mandatory[i] ? 1 : 0);
        }

        // Try each insertion position m in [0, L-1]
        // m=0: skip motif[0], motif[1..L-1] aligns with seq[0..L-2]
        // m=L-1: skip motif[L-1], motif[0..L-2] aligns with seq[0..L-2]
        float gains_arr[MAX_MOTIF_LEN_OPT];
        for (int m = 0; m < L; m++) {
            double aligned_score = prefix_cum_score[m] + suffix_cum_score[m];

            // Quick O(1) deficit check: skip expensive gains collection if hopeless
            if (!quick_deficit_check(aligned_score, 1, best_edits)) continue;

            int mandatory_subs = prefix_cum_mandatory[m] + suffix_cum_mandatory[m];

            // Collect gains from aligned positions (skip mandatory positions)
            int n_gains = 0;
            for (int i = 0; i < m; i++) {
                if (!prefix_mandatory[i] && prefix_gains[i] > 1e-12f) {
                    gains_arr[n_gains++] = prefix_gains[i];
                }
            }
            for (int i = m + 1; i < L; i++) {
                if (!suffix_mandatory[i] && suffix_gains[i] > 1e-12f) {
                    gains_arr[n_gains++] = suffix_gains[i];
                }
            }

            std::vector<float> gains(gains_arr, gains_arr + n_gains);
            float total = compute_total_edits(aligned_score, gains, 1, mandatory_subs);
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

float PWMEditDistanceScorer::compute_with_two_indels(const int* bidx_arr, int seq_len, bool reverse)
{
    const int L = m_pssm.length();
    if (L == 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    // Fall back to generic DP solver for motifs exceeding stack array limit
    if (L > MAX_MOTIF_LEN_OPT) {
        return compute_with_indels(bidx_arr, seq_len, reverse);
    }

    float best_edits = std::numeric_limits<float>::quiet_NaN();

    // get_base_score: direct table lookup using precomputed base indices.
    // Returns true if the position requires a mandatory substitution.
    auto get_base_score = [&](int motif_pos, int raw_seq_idx,
                               float& out_score, float& out_gain) -> bool {
        int bi = bidx_arr[raw_seq_idx];
        out_score = m_score_table[motif_pos][bi];
        out_gain = m_gain_table[motif_pos][bi];
        return m_mandatory_table[motif_pos][bi];
    };

    // compute_total_edits: pass indels + mandatory_subs as the base fixed-edit count.
    auto compute_total_edits = [&](double aligned_score, std::vector<float>& gains,
                                   int indels, int mandatory_subs) -> float {
        return compute_min_edits_from_gains(aligned_score, gains,
                                            indels + mandatory_subs, best_edits);
    };

    auto update_best = [&](float total) {
        if (!std::isnan(total) && (std::isnan(best_edits) || total < best_edits)) {
            best_edits = total;
        }
    };

    // ===== Case 1: No indel (W = L, k = 0) =====
    {
        double aligned_score = 0.0;
        float gains_arr[MAX_MOTIF_LEN_OPT];
        int n_gains = 0;
        int mandatory_subs = 0;

        for (int i = 0; i < L; i++) {
            int raw_seq_idx = reverse ? (L - 1 - i) : i;
            float bs, gain;
            bool mandatory = get_base_score(i, raw_seq_idx, bs, gain);
            aligned_score += static_cast<double>(bs);
            if (mandatory) {
                mandatory_subs++;
            } else if (gain > 1e-12f) {
                gains_arr[n_gains++] = gain;
            }
        }

        std::vector<float> gains(gains_arr, gains_arr + n_gains);
        update_best(compute_total_edits(aligned_score, gains, 0, mandatory_subs));

        // Early return: 0 edits is optimal
        if (!std::isnan(best_edits) && best_edits <= 0.0f) {
            return best_edits;
        }
    }

    // Helper: collect gains and mandatory count for a set of (motif_col, score, gain, mandatory) entries.
    // Fills gains vector (non-mandatory positions with gain > 1e-12) and returns mandatory count.
    // Used by the multi-segment cases below.

    // ===== Case 2: One deletion (W = L + 1, k = 1) =====
    if (L + 1 <= seq_len) {
        const int W = L + 1;

        float prefix_scores[MAX_MOTIF_LEN_OPT], prefix_gains_arr[MAX_MOTIF_LEN_OPT];
        float suffix_scores[MAX_MOTIF_LEN_OPT], suffix_gains_arr[MAX_MOTIF_LEN_OPT];
        bool prefix_mandatory[MAX_MOTIF_LEN_OPT], suffix_mandatory[MAX_MOTIF_LEN_OPT];
        memset(prefix_mandatory, 0, sizeof(bool) * L);
        memset(suffix_mandatory, 0, sizeof(bool) * L);

        for (int i = 0; i < L; i++) {
            int raw_prefix = reverse ? (W - 1 - i) : i;
            prefix_mandatory[i] = get_base_score(i, raw_prefix,
                                                  prefix_scores[i], prefix_gains_arr[i]);

            int raw_suffix = reverse ? (W - 1 - (i + 1)) : (i + 1);
            suffix_mandatory[i] = get_base_score(i, raw_suffix,
                                                  suffix_scores[i], suffix_gains_arr[i]);
        }

        double prefix_cum[MAX_MOTIF_LEN_OPT + 1];
        int prefix_cum_mand[MAX_MOTIF_LEN_OPT + 1];
        prefix_cum[0] = 0.0;
        prefix_cum_mand[0] = 0;
        for (int i = 0; i < L; i++) {
            prefix_cum[i + 1] = prefix_cum[i] + static_cast<double>(prefix_scores[i]);
            prefix_cum_mand[i + 1] = prefix_cum_mand[i] + (prefix_mandatory[i] ? 1 : 0);
        }

        double suffix_cum[MAX_MOTIF_LEN_OPT + 1];
        int suffix_cum_mand[MAX_MOTIF_LEN_OPT + 1];
        suffix_cum[L] = 0.0;
        suffix_cum_mand[L] = 0;
        for (int i = L - 1; i >= 0; i--) {
            suffix_cum[i] = suffix_cum[i + 1] + static_cast<double>(suffix_scores[i]);
            suffix_cum_mand[i] = suffix_cum_mand[i + 1] + (suffix_mandatory[i] ? 1 : 0);
        }

        float gains_arr[MAX_MOTIF_LEN_OPT];
        for (int d = 0; d <= L; d++) {
            double aligned_score = prefix_cum[d] + suffix_cum[d];

            if (!quick_deficit_check(aligned_score, 1, best_edits)) continue;

            int mandatory_subs = prefix_cum_mand[d] + suffix_cum_mand[d];
            int n_gains = 0;
            for (int i = 0; i < d; i++) {
                if (!prefix_mandatory[i] && prefix_gains_arr[i] > 1e-12f)
                    gains_arr[n_gains++] = prefix_gains_arr[i];
            }
            for (int i = d; i < L; i++) {
                if (!suffix_mandatory[i] && suffix_gains_arr[i] > 1e-12f)
                    gains_arr[n_gains++] = suffix_gains_arr[i];
            }

            std::vector<float> gains(gains_arr, gains_arr + n_gains);
            update_best(compute_total_edits(aligned_score, gains, 1, mandatory_subs));
            if (!std::isnan(best_edits) && best_edits <= 1.0f) break;
        }
    }

    // ===== Case 3: One insertion (W = L - 1, k = 1) =====
    // Skip entirely if best_edits <= 1 (can't beat a 1-indel result with another 1-indel)
    if (L >= 2 && (std::isnan(best_edits) || best_edits > 1.0f)) {
        const int W = L - 1;

        float prefix_scores[MAX_MOTIF_LEN_OPT], prefix_gains_arr[MAX_MOTIF_LEN_OPT];
        float suffix_scores[MAX_MOTIF_LEN_OPT], suffix_gains_arr[MAX_MOTIF_LEN_OPT];
        bool prefix_mandatory[MAX_MOTIF_LEN_OPT], suffix_mandatory[MAX_MOTIF_LEN_OPT];
        memset(prefix_mandatory, 0, sizeof(bool) * (L - 1));
        memset(suffix_mandatory, 0, sizeof(bool) * L);

        for (int i = 0; i < L - 1; i++) {
            int raw_idx = reverse ? (W - 1 - i) : i;
            prefix_mandatory[i] = get_base_score(i, raw_idx,
                                                  prefix_scores[i], prefix_gains_arr[i]);
        }

        for (int i = 1; i < L; i++) {
            int raw_idx = reverse ? (W - 1 - (i - 1)) : (i - 1);
            suffix_mandatory[i] = get_base_score(i, raw_idx,
                                                  suffix_scores[i], suffix_gains_arr[i]);
        }

        double prefix_cum[MAX_MOTIF_LEN_OPT];
        int prefix_cum_mand[MAX_MOTIF_LEN_OPT];
        prefix_cum[0] = 0.0;
        prefix_cum_mand[0] = 0;
        for (int i = 0; i < L - 1; i++) {
            prefix_cum[i + 1] = prefix_cum[i] + static_cast<double>(prefix_scores[i]);
            prefix_cum_mand[i + 1] = prefix_cum_mand[i] + (prefix_mandatory[i] ? 1 : 0);
        }

        double suffix_cum[MAX_MOTIF_LEN_OPT];
        int suffix_cum_mand[MAX_MOTIF_LEN_OPT];
        suffix_cum[L - 1] = 0.0;
        suffix_cum_mand[L - 1] = 0;
        for (int i = L - 1; i >= 1; i--) {
            suffix_cum[i - 1] = suffix_cum[i] + static_cast<double>(suffix_scores[i]);
            suffix_cum_mand[i - 1] = suffix_cum_mand[i] + (suffix_mandatory[i] ? 1 : 0);
        }

        float gains_arr[MAX_MOTIF_LEN_OPT];
        for (int m = 0; m < L; m++) {
            double aligned_score = prefix_cum[m] + suffix_cum[m];

            if (!quick_deficit_check(aligned_score, 1, best_edits)) continue;

            int mandatory_subs = prefix_cum_mand[m] + suffix_cum_mand[m];
            int n_gains = 0;
            for (int i = 0; i < m; i++) {
                if (!prefix_mandatory[i] && prefix_gains_arr[i] > 1e-12f)
                    gains_arr[n_gains++] = prefix_gains_arr[i];
            }
            for (int i = m + 1; i < L; i++) {
                if (!suffix_mandatory[i] && suffix_gains_arr[i] > 1e-12f)
                    gains_arr[n_gains++] = suffix_gains_arr[i];
            }

            std::vector<float> gains(gains_arr, gains_arr + n_gains);
            update_best(compute_total_edits(aligned_score, gains, 1, mandatory_subs));
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

        float scores0[MAX_MOTIF_LEN_OPT], gains0[MAX_MOTIF_LEN_OPT];
        float scores1[MAX_MOTIF_LEN_OPT], gains1[MAX_MOTIF_LEN_OPT];
        float scores2[MAX_MOTIF_LEN_OPT], gains2[MAX_MOTIF_LEN_OPT];
        bool mand0[MAX_MOTIF_LEN_OPT], mand1[MAX_MOTIF_LEN_OPT], mand2[MAX_MOTIF_LEN_OPT];
        memset(mand0, 0, sizeof(bool) * L);
        memset(mand1, 0, sizeof(bool) * L);
        memset(mand2, 0, sizeof(bool) * L);

        for (int i = 0; i < L; i++) {
            int raw0 = reverse ? (W - 1 - i) : i;
            mand0[i] = get_base_score(i, raw0, scores0[i], gains0[i]);

            int raw1 = reverse ? (W - 1 - (i + 1)) : (i + 1);
            mand1[i] = get_base_score(i, raw1, scores1[i], gains1[i]);

            int raw2 = reverse ? (W - 1 - (i + 2)) : (i + 2);
            mand2[i] = get_base_score(i, raw2, scores2[i], gains2[i]);
        }

        // Prefix cumulative for config0
        double cum0[MAX_MOTIF_LEN_OPT + 1];
        int cum0_mand[MAX_MOTIF_LEN_OPT + 1];
        cum0[0] = 0.0;
        cum0_mand[0] = 0;
        for (int i = 0; i < L; i++) {
            cum0[i + 1] = cum0[i] + static_cast<double>(scores0[i]);
            cum0_mand[i + 1] = cum0_mand[i] + (mand0[i] ? 1 : 0);
        }

        // Suffix cumulative for config2
        double suf2[MAX_MOTIF_LEN_OPT + 1];
        int suf2_mand[MAX_MOTIF_LEN_OPT + 1];
        suf2[L] = 0.0;
        suf2_mand[L] = 0;
        for (int i = L - 1; i >= 0; i--) {
            suf2[i] = suf2[i + 1] + static_cast<double>(scores2[i]);
            suf2_mand[i] = suf2_mand[i + 1] + (mand2[i] ? 1 : 0);
        }

        // Prefix cumulative for config1 (middle segment)
        double cum1[MAX_MOTIF_LEN_OPT + 1];
        int cum1_mand[MAX_MOTIF_LEN_OPT + 1];
        cum1[0] = 0.0;
        cum1_mand[0] = 0;
        for (int i = 0; i < L; i++) {
            cum1[i + 1] = cum1[i] + static_cast<double>(scores1[i]);
            cum1_mand[i + 1] = cum1_mand[i] + (mand1[i] ? 1 : 0);
        }

        float gains_arr[MAX_MOTIF_LEN_OPT];
        for (int d1 = 0; d1 <= L + 1; d1++) {
            for (int d2 = d1 + 1; d2 <= L + 1; d2++) {
                // Seg 1: motif[0..d1-1] config0
                double seg1 = cum0[d1];
                int mand1_count = cum0_mand[d1];

                // Seg 2: motif[d1..d2-2] config1
                double seg2 = 0.0;
                int mand2_count = 0;
                if (d2 - 1 > d1) {
                    seg2 = cum1[d2 - 1] - cum1[d1];
                    mand2_count = cum1_mand[d2 - 1] - cum1_mand[d1];
                }

                // Seg 3: motif[d2-1..L-1] config2
                double seg3 = 0.0;
                int mand3_count = 0;
                if (d2 - 1 < L) {
                    seg3 = suf2[d2 - 1];
                    mand3_count = suf2_mand[d2 - 1];
                }

                double aligned_score = seg1 + seg2 + seg3;

                if (!quick_deficit_check(aligned_score, 2, best_edits)) continue;

                int mandatory_subs = mand1_count + mand2_count + mand3_count;
                int n_gains = 0;
                for (int i = 0; i < d1; i++) {
                    if (!mand0[i] && gains0[i] > 1e-12f) gains_arr[n_gains++] = gains0[i];
                }
                for (int i = d1; i < d2 - 1; i++) {
                    if (!mand1[i] && gains1[i] > 1e-12f) gains_arr[n_gains++] = gains1[i];
                }
                for (int i = d2 - 1; i < L; i++) {
                    if (!mand2[i] && gains2[i] > 1e-12f) gains_arr[n_gains++] = gains2[i];
                }

                std::vector<float> gains(gains_arr, gains_arr + n_gains);
                update_best(compute_total_edits(aligned_score, gains, 2, mandatory_subs));
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

        float scores0[MAX_MOTIF_LEN_OPT], gains0_v[MAX_MOTIF_LEN_OPT];
        float scores1[MAX_MOTIF_LEN_OPT], gains1_v[MAX_MOTIF_LEN_OPT];
        float scores2[MAX_MOTIF_LEN_OPT], gains2_v[MAX_MOTIF_LEN_OPT];
        bool mand0[MAX_MOTIF_LEN_OPT], mand1[MAX_MOTIF_LEN_OPT], mand2[MAX_MOTIF_LEN_OPT];
        memset(mand0, 0, sizeof(bool) * L);
        memset(mand1, 0, sizeof(bool) * L);
        memset(mand2, 0, sizeof(bool) * L);

        for (int i = 0; i < L; i++) {
            if (i < W) {
                int raw0 = reverse ? (W - 1 - i) : i;
                mand0[i] = get_base_score(i, raw0, scores0[i], gains0_v[i]);
            }
            if (i >= 1 && i - 1 < W) {
                int raw1 = reverse ? (W - 1 - (i - 1)) : (i - 1);
                mand1[i] = get_base_score(i, raw1, scores1[i], gains1_v[i]);
            }
            if (i >= 2 && i - 2 < W) {
                int raw2 = reverse ? (W - 1 - (i - 2)) : (i - 2);
                mand2[i] = get_base_score(i, raw2, scores2[i], gains2_v[i]);
            }
        }

        // Prefix cumulative for config0
        double cum0[MAX_MOTIF_LEN_OPT + 1];
        int cum0_mand[MAX_MOTIF_LEN_OPT + 1];
        cum0[0] = 0.0;
        cum0_mand[0] = 0;
        for (int i = 0; i < std::min(L, W); i++) {
            cum0[i + 1] = cum0[i] + static_cast<double>(scores0[i]);
            cum0_mand[i + 1] = cum0_mand[i] + (mand0[i] ? 1 : 0);
        }

        // Suffix cumulative for config2
        double suf2[MAX_MOTIF_LEN_OPT + 1];
        int suf2_mand[MAX_MOTIF_LEN_OPT + 1];
        suf2[L] = 0.0;
        suf2_mand[L] = 0;
        for (int i = L - 1; i >= 2; i--) {
            if (i - 2 < W) {
                suf2[i] = suf2[i + 1] + static_cast<double>(scores2[i]);
                suf2_mand[i] = suf2_mand[i + 1] + (mand2[i] ? 1 : 0);
            } else {
                suf2[i] = suf2[i + 1];
                suf2_mand[i] = suf2_mand[i + 1];
            }
        }

        // Prefix cumulative for config1
        double cum1[MAX_MOTIF_LEN_OPT + 1];
        int cum1_mand[MAX_MOTIF_LEN_OPT + 1];
        cum1[0] = 0.0;
        cum1_mand[0] = 0;
        // cum1 starts accumulating from index 1
        cum1[1] = 0.0;
        cum1_mand[1] = 0;
        for (int i = 1; i < L; i++) {
            if (i - 1 < W) {
                cum1[i + 1] = cum1[i] + static_cast<double>(scores1[i]);
                cum1_mand[i + 1] = cum1_mand[i] + (mand1[i] ? 1 : 0);
            } else {
                cum1[i + 1] = cum1[i];
                cum1_mand[i + 1] = cum1_mand[i];
            }
        }

        float gains_arr[MAX_MOTIF_LEN_OPT];
        for (int m1 = 0; m1 < L; m1++) {
            for (int m2 = m1 + 1; m2 < L; m2++) {
                double seg1 = cum0[m1];
                int mand1_count = cum0_mand[m1];

                double seg2 = 0.0;
                int mand2_count = 0;
                if (m2 > m1 + 1) {
                    seg2 = cum1[m2] - cum1[m1 + 1];
                    mand2_count = cum1_mand[m2] - cum1_mand[m1 + 1];
                }

                double seg3 = 0.0;
                int mand3_count = 0;
                if (m2 + 1 < L) {
                    seg3 = suf2[m2 + 1];
                    mand3_count = suf2_mand[m2 + 1];
                }

                double aligned_score = seg1 + seg2 + seg3;

                if (!quick_deficit_check(aligned_score, 2, best_edits)) continue;

                int mandatory_subs = mand1_count + mand2_count + mand3_count;
                int n_gains = 0;
                for (int i = 0; i < m1; i++) {
                    if (!mand0[i] && gains0_v[i] > 1e-12f) gains_arr[n_gains++] = gains0_v[i];
                }
                for (int i = m1 + 1; i < m2; i++) {
                    if (!mand1[i] && gains1_v[i] > 1e-12f) gains_arr[n_gains++] = gains1_v[i];
                }
                for (int i = m2 + 1; i < L; i++) {
                    if (!mand2[i] && gains2_v[i] > 1e-12f) gains_arr[n_gains++] = gains2_v[i];
                }

                std::vector<float> gains(gains_arr, gains_arr + n_gains);
                update_best(compute_total_edits(aligned_score, gains, 2, mandatory_subs));
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

        float scores_ns[MAX_MOTIF_LEN_OPT], gains_ns[MAX_MOTIF_LEN_OPT];   // no shift
        float scores_sp1[MAX_MOTIF_LEN_OPT], gains_sp1[MAX_MOTIF_LEN_OPT]; // seq shifted +1
        float scores_sm1[MAX_MOTIF_LEN_OPT], gains_sm1[MAX_MOTIF_LEN_OPT]; // seq shifted -1
        bool mand_ns[MAX_MOTIF_LEN_OPT], mand_sp1[MAX_MOTIF_LEN_OPT], mand_sm1[MAX_MOTIF_LEN_OPT];
        memset(mand_ns, 0, sizeof(bool) * L);
        memset(mand_sp1, 0, sizeof(bool) * L);
        memset(mand_sm1, 0, sizeof(bool) * L);

        for (int i = 0; i < L; i++) {
            int raw_ns = reverse ? (W - 1 - i) : i;
            mand_ns[i] = get_base_score(i, raw_ns, scores_ns[i], gains_ns[i]);

            if (i + 1 < W) {
                int raw_sp1 = reverse ? (W - 1 - (i + 1)) : (i + 1);
                mand_sp1[i] = get_base_score(i, raw_sp1, scores_sp1[i], gains_sp1[i]);
            }

            if (i >= 1) {
                int raw_sm1 = reverse ? (W - 1 - (i - 1)) : (i - 1);
                mand_sm1[i] = get_base_score(i, raw_sm1, scores_sm1[i], gains_sm1[i]);
            }
        }

        // Prefix sums for no-shift
        double cum_ns[MAX_MOTIF_LEN_OPT + 1];
        int cum_ns_mand[MAX_MOTIF_LEN_OPT + 1];
        cum_ns[0] = 0.0;
        cum_ns_mand[0] = 0;
        for (int i = 0; i < L; i++) {
            cum_ns[i + 1] = cum_ns[i] + static_cast<double>(scores_ns[i]);
            cum_ns_mand[i + 1] = cum_ns_mand[i] + (mand_ns[i] ? 1 : 0);
        }

        // Prefix sums for seq+1 shift
        double cum_sp1[MAX_MOTIF_LEN_OPT + 1];
        int cum_sp1_mand[MAX_MOTIF_LEN_OPT + 1];
        cum_sp1[0] = 0.0;
        cum_sp1_mand[0] = 0;
        for (int i = 0; i < L - 1; i++) {
            cum_sp1[i + 1] = cum_sp1[i] + static_cast<double>(scores_sp1[i]);
            cum_sp1_mand[i + 1] = cum_sp1_mand[i] + (mand_sp1[i] ? 1 : 0);
        }

        // Prefix sums for seq-1 shift (valid for i >= 1)
        double cum_sm1[MAX_MOTIF_LEN_OPT + 1];
        int cum_sm1_mand[MAX_MOTIF_LEN_OPT + 1];
        cum_sm1[0] = 0.0;
        cum_sm1_mand[0] = 0;
        cum_sm1[1] = 0.0;
        cum_sm1_mand[1] = 0;
        for (int i = 1; i < L; i++) {
            cum_sm1[i + 1] = cum_sm1[i] + static_cast<double>(scores_sm1[i]);
            cum_sm1_mand[i + 1] = cum_sm1_mand[i] + (mand_sm1[i] ? 1 : 0);
        }

        // Suffix sums for no-shift
        double suf_ns[MAX_MOTIF_LEN_OPT + 1];
        int suf_ns_mand[MAX_MOTIF_LEN_OPT + 1];
        suf_ns[L] = 0.0;
        suf_ns_mand[L] = 0;
        for (int i = L - 1; i >= 0; i--) {
            suf_ns[i] = suf_ns[i + 1] + static_cast<double>(scores_ns[i]);
            suf_ns_mand[i] = suf_ns_mand[i + 1] + (mand_ns[i] ? 1 : 0);
        }

        float gains_arr[MAX_MOTIF_LEN_OPT];
        for (int d = 0; d < L; d++) {
            for (int m = 0; m < L; m++) {
                // d == m is always dominated by the no-indel case (Case 1):
                // it skips both seq[d] and motif[d], aligning L-2 pairs with 2 indels,
                // while Case 1 aligns all L pairs with 0 indels.
                if (d == m) continue;

                // Compute aligned_score and mandatory count first for quick deficit check
                double aligned_score = 0.0;
                int mandatory_subs = 0;
                if (d < m) {
                    aligned_score = cum_ns[d] + (cum_sp1[m] - cum_sp1[d]);
                    mandatory_subs = cum_ns_mand[d] + (cum_sp1_mand[m] - cum_sp1_mand[d]);
                    if (m + 1 < L) {
                        aligned_score += suf_ns[m + 1];
                        mandatory_subs += suf_ns_mand[m + 1];
                    }
                } else {
                    aligned_score = cum_ns[m] + (cum_sm1[d + 1] - cum_sm1[m + 1]);
                    mandatory_subs = cum_ns_mand[m] + (cum_sm1_mand[d + 1] - cum_sm1_mand[m + 1]);
                    if (d + 1 < L) {
                        aligned_score += suf_ns[d + 1];
                        mandatory_subs += suf_ns_mand[d + 1];
                    }
                }

                if (!quick_deficit_check(aligned_score, 2, best_edits)) continue;

                // Collect gains only for promising candidates (skip mandatory positions)
                int n_gains = 0;
                if (d < m) {
                    for (int i = 0; i < d; i++) {
                        if (!mand_ns[i] && gains_ns[i] > 1e-12f) gains_arr[n_gains++] = gains_ns[i];
                    }
                    for (int i = d; i < m; i++) {
                        if (!mand_sp1[i] && gains_sp1[i] > 1e-12f) gains_arr[n_gains++] = gains_sp1[i];
                    }
                    for (int i = m + 1; i < L; i++) {
                        if (!mand_ns[i] && gains_ns[i] > 1e-12f) gains_arr[n_gains++] = gains_ns[i];
                    }
                } else {
                    for (int i = 0; i < m; i++) {
                        if (!mand_ns[i] && gains_ns[i] > 1e-12f) gains_arr[n_gains++] = gains_ns[i];
                    }
                    for (int i = m + 1; i <= d; i++) {
                        if (!mand_sm1[i] && gains_sm1[i] > 1e-12f) gains_arr[n_gains++] = gains_sm1[i];
                    }
                    for (int i = d + 1; i < L; i++) {
                        if (!mand_ns[i] && gains_ns[i] > 1e-12f) gains_arr[n_gains++] = gains_ns[i];
                    }
                }

                std::vector<float> gains(gains_arr, gains_arr + n_gains);
                update_best(compute_total_edits(aligned_score, gains, 2, mandatory_subs));
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
