#include "PWMLseEditDistanceScorer.h"
#include "GenomeSeqFetch.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <queue>

PWMLseEditDistanceScorer::PWMLseEditDistanceScorer(const DnaPSSM& pssm,
                                                     GenomeSeqFetch* shared_seqfetch,
                                                     float threshold,
                                                     int max_edits,
                                                     bool extend,
                                                     char strand,
                                                     Mode mode,
                                                     float score_min)
    : GenomeSeqScorer(shared_seqfetch, extend, strand),
      m_pssm(pssm),
      m_threshold(threshold),
      m_max_edits(max_edits),
      m_mode(mode),
      m_score_min(score_min),
      m_last_min_edits(std::numeric_limits<float>::quiet_NaN()),
      m_S_max(0.0f)
{
    precompute_tables();
}

void PWMLseEditDistanceScorer::precompute_tables()
{
    int L = m_pssm.length();
    m_col_max_scores.resize(L);
    m_S_max = 0.0f;

    for (int i = 0; i < L; i++) {
        float max_score = -std::numeric_limits<float>::infinity();
        for (int b = 0; b < 4; b++) {
            float score = m_pssm[i].get_log_prob_from_code(b);
            if (score > max_score) {
                max_score = score;
            }
        }
        m_col_max_scores[i] = max_score;
        m_S_max += max_score;
    }
}

float PWMLseEditDistanceScorer::score_interval(const GInterval& interval,
                                                const GenomeChromKey& chromkey)
{
    m_last_min_edits = std::numeric_limits<float>::quiet_NaN();

    const int motif_len = m_pssm.length();
    if (motif_len <= 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    GInterval exp_interval = calculate_expanded_interval(interval, chromkey, motif_len);

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

    const bool scan_forward = should_scan_forward();
    const bool scan_reverse = should_scan_reverse();

    LseResult best_result;

    if (scan_forward) {
        LseResult fwd = compute_lse_edit_distance(seq, static_cast<size_t>(interval_len), false);
        if (!std::isnan(fwd.min_edits)) {
            if (std::isnan(best_result.min_edits) || fwd.min_edits < best_result.min_edits) {
                best_result = fwd;
            }
        }
    }

    if (scan_reverse) {
        LseResult rev = compute_lse_edit_distance(seq, static_cast<size_t>(interval_len), true);
        if (!std::isnan(rev.min_edits)) {
            if (std::isnan(best_result.min_edits) || rev.min_edits < best_result.min_edits) {
                best_result = rev;
            }
        }
    }

    m_last_min_edits = best_result.min_edits;

    switch (m_mode) {
        case Mode::LSE_EDIT_DISTANCE:
            return best_result.min_edits;
        case Mode::LSE_EDIT_DISTANCE_POS:
            return best_result.best_pos;
        default:
            return std::numeric_limits<float>::quiet_NaN();
    }
}

float PWMLseEditDistanceScorer::encode_position(size_t index,
                                                 size_t target_length,
                                                 size_t motif_length,
                                                 int direction) const
{
    // Convert 0-based index to 1-based position
    float pos = static_cast<float>(index) + 1.0f;
    if (m_strand == -1) {
        pos = static_cast<float>(target_length) - pos + 1.0f;
    }
    if (m_pssm.is_bidirect()) {
        pos *= static_cast<float>(direction);
    }
    return pos;
}

float PWMLseEditDistanceScorer::compute_window_score(const char* window_start, bool reverse) const
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
        if (bidx == 4) {
            // Unknown base: use mean log-probability
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

void PWMLseEditDistanceScorer::init_partition(const std::vector<int>& seq_bases,
                                               int N, int L,
                                               std::vector<double>& S_p,
                                               std::vector<double>& A_p,
                                               double& Z, double& m, double& F) const
{
    S_p.resize(N);
    A_p.resize(N);

    // Compute per-start scores
    for (int p = 0; p < N; p++) {
        double score = 0.0;
        for (int i = 0; i < L; i++) {
            int b = seq_bases[p + i];
            if (b == 4) {
                // Unknown base: use mean log-probability
                float sum = 0.0f;
                for (int bb = 0; bb < 4; bb++) {
                    sum += m_pssm[i].get_log_prob_from_code(bb);
                }
                score += static_cast<double>(sum / 4.0f);
            } else {
                score += static_cast<double>(m_pssm[i].get_log_prob_from_code(b));
            }
        }
        S_p[p] = score;
    }

    // Find max for numerical stability
    m = -std::numeric_limits<double>::infinity();
    for (int p = 0; p < N; p++) {
        if (S_p[p] > m) m = S_p[p];
    }

    // Compute shifted weights and partition function
    Z = 0.0;
    for (int p = 0; p < N; p++) {
        A_p[p] = std::exp(S_p[p] - m);
        Z += A_p[p];
    }

    F = std::log(Z) + m;
}

double PWMLseEditDistanceScorer::compute_delta_z(int q, int y_idx,
                                                  const std::vector<int>& seq_bases,
                                                  const std::vector<double>& A_p,
                                                  int N, int L) const
{
    int cur_base = seq_bases[q];
    if (cur_base == y_idx) return 0.0;  // No change

    double delta_z = 0.0;

    // For each start p that covers position q: p in [q-L+1, q], clipped to [0, N-1]
    int p_min = std::max(0, q - L + 1);
    int p_max = std::min(N - 1, q);

    for (int p = p_min; p <= p_max; p++) {
        int i = q - p;  // column index in the motif

        // Score change at column i: old base -> new base
        double old_score, new_score;
        if (cur_base == 4) {
            float sum = 0.0f;
            for (int b = 0; b < 4; b++) {
                sum += m_pssm[i].get_log_prob_from_code(b);
            }
            old_score = static_cast<double>(sum / 4.0f);
        } else {
            old_score = static_cast<double>(m_pssm[i].get_log_prob_from_code(cur_base));
        }
        new_score = static_cast<double>(m_pssm[i].get_log_prob_from_code(y_idx));

        double delta_s = new_score - old_score;
        // DeltaZ contribution = A_p * expm1(delta_s)
        // Using expm1 for numerical stability when delta_s is near zero
        delta_z += A_p[p] * std::expm1(delta_s);
    }

    return delta_z;
}

void PWMLseEditDistanceScorer::apply_edit(int q, int y_idx,
                                           std::vector<int>& seq_bases,
                                           std::vector<double>& S_p,
                                           std::vector<double>& A_p,
                                           double& Z, double& m,
                                           int N, int L)
{
    int old_base = seq_bases[q];
    seq_bases[q] = y_idx;

    int p_min = std::max(0, q - L + 1);
    int p_max = std::min(N - 1, q);

    for (int p = p_min; p <= p_max; p++) {
        int i = q - p;  // column index

        double old_score;
        if (old_base == 4) {
            float sum = 0.0f;
            for (int b = 0; b < 4; b++) {
                sum += m_pssm[i].get_log_prob_from_code(b);
            }
            old_score = static_cast<double>(sum / 4.0f);
        } else {
            old_score = static_cast<double>(m_pssm[i].get_log_prob_from_code(old_base));
        }
        double new_score = static_cast<double>(m_pssm[i].get_log_prob_from_code(y_idx));
        double delta_s = new_score - old_score;

        double A_old = A_p[p];
        S_p[p] += delta_s;
        A_p[p] = A_old * std::exp(delta_s);
        Z += A_old * std::expm1(delta_s);
    }

    // Renormalize if drift is too large (threshold for double: ~700)
    double max_sp = -std::numeric_limits<double>::infinity();
    for (int p = p_min; p <= p_max; p++) {
        if (S_p[p] > max_sp) max_sp = S_p[p];
    }

    if (max_sp - m > 500.0) {
        // Full renormalization
        double new_m = -std::numeric_limits<double>::infinity();
        for (int p = 0; p < N; p++) {
            if (S_p[p] > new_m) new_m = S_p[p];
        }
        double scale = std::exp(m - new_m);
        for (int p = 0; p < N; p++) {
            A_p[p] *= scale;
        }
        Z *= scale;
        m = new_m;
    }
}

PWMLseEditDistanceScorer::LseResult
PWMLseEditDistanceScorer::try_k1_exhaustive(const std::vector<int>& seq_bases,
                                             const std::vector<double>& A_p,
                                             double Z, double m, double F,
                                             int N, int L,
                                             size_t target_length, size_t motif_length,
                                             int direction) const
{
    LseResult result;
    double T = static_cast<double>(m_threshold);
    int R = static_cast<int>(seq_bases.size());

    double best_F_after = -std::numeric_limits<double>::infinity();
    int best_q = -1;

    for (int q = 0; q < R; q++) {
        int cur = seq_bases[q];
        for (int y = 0; y < 4; y++) {
            if (y == cur) continue;

            double dz = compute_delta_z(q, y, seq_bases, A_p, N, L);
            double new_Z = Z + dz;
            if (new_Z <= 0.0) continue;

            double new_F = std::log(new_Z) + m;
            if (new_F >= T) {
                // Found a valid 1-edit solution
                if (new_F > best_F_after) {
                    best_F_after = new_F;
                    best_q = q;
                }
            }
        }
    }

    if (best_q >= 0) {
        result.min_edits = 1.0f;
        result.best_edit_index = best_q;
        result.best_edit_direction = direction;
        result.best_pos = encode_position(static_cast<size_t>(best_q),
                                          target_length, motif_length, direction);
    }

    return result;
}

PWMLseEditDistanceScorer::LseResult
PWMLseEditDistanceScorer::try_k2_exhaustive(const std::vector<int>& seq_bases,
                                             const std::vector<double>& S_p,
                                             const std::vector<double>& A_p,
                                             double Z, double m, double F,
                                             int N, int L, int R,
                                             size_t target_length, size_t motif_length,
                                             int direction) const
{
    LseResult result;
    double T = static_cast<double>(m_threshold);

    double best_F_after = -std::numeric_limits<double>::infinity();
    int best_q1 = -1;

    // For each first edit (q1, y1)
    for (int q1 = 0; q1 < R; q1++) {
        int cur1 = seq_bases[q1];

        for (int y1 = 0; y1 < 4; y1++) {
            if (y1 == cur1) continue;

            // Compute the effect of first edit on partition function
            // We need to temporarily apply it, then check all second edits

            // Compute updated A_p values for starts affected by q1
            int p_min1 = std::max(0, q1 - L + 1);
            int p_max1 = std::min(N - 1, q1);

            // Store changes from first edit
            std::vector<double> A_p_temp(A_p);
            std::vector<double> S_p_temp(S_p);
            double Z_temp = Z;

            for (int p = p_min1; p <= p_max1; p++) {
                int i = q1 - p;
                double old_s, new_s;
                if (cur1 == 4) {
                    float sum = 0.0f;
                    for (int b = 0; b < 4; b++) sum += m_pssm[i].get_log_prob_from_code(b);
                    old_s = static_cast<double>(sum / 4.0f);
                } else {
                    old_s = static_cast<double>(m_pssm[i].get_log_prob_from_code(cur1));
                }
                new_s = static_cast<double>(m_pssm[i].get_log_prob_from_code(y1));
                double ds = new_s - old_s;

                double A_old = A_p_temp[p];
                S_p_temp[p] += ds;
                A_p_temp[p] = A_old * std::exp(ds);
                Z_temp += A_old * std::expm1(ds);
            }

            // Now try all second edits (q2, y2) with q2 > q1
            // (q2 < q1 was already covered when q1 was smaller)
            for (int q2 = q1 + 1; q2 < R; q2++) {
                // Use original base at q2 (first edit was at q1, not q2)
                int cur2 = seq_bases[q2];

                for (int y2 = 0; y2 < 4; y2++) {
                    if (y2 == cur2) continue;

                    // Compute DeltaZ for second edit using updated A_p
                    double dz2 = 0.0;
                    int p_min2 = std::max(0, q2 - L + 1);
                    int p_max2 = std::min(N - 1, q2);

                    for (int p = p_min2; p <= p_max2; p++) {
                        int i = q2 - p;
                        double old_s2;
                        if (cur2 == 4) {
                            float sum = 0.0f;
                            for (int b = 0; b < 4; b++) sum += m_pssm[i].get_log_prob_from_code(b);
                            old_s2 = static_cast<double>(sum / 4.0f);
                        } else {
                            old_s2 = static_cast<double>(m_pssm[i].get_log_prob_from_code(cur2));
                        }
                        double new_s2 = static_cast<double>(m_pssm[i].get_log_prob_from_code(y2));
                        double ds2 = new_s2 - old_s2;
                        dz2 += A_p_temp[p] * std::expm1(ds2);
                    }

                    double new_Z = Z_temp + dz2;
                    if (new_Z <= 0.0) continue;

                    double new_F = std::log(new_Z) + m;
                    if (new_F >= T) {
                        if (new_F > best_F_after) {
                            best_F_after = new_F;
                            best_q1 = q1;
                        }
                    }
                }
            }
        }
    }

    if (best_q1 >= 0) {
        result.min_edits = 2.0f;
        result.best_edit_index = best_q1;
        result.best_edit_direction = direction;
        result.best_pos = encode_position(static_cast<size_t>(best_q1),
                                          target_length, motif_length, direction);
    }

    return result;
}

PWMLseEditDistanceScorer::LseResult
PWMLseEditDistanceScorer::greedy_search(std::vector<int>& seq_bases,
                                         std::vector<double>& S_p,
                                         std::vector<double>& A_p,
                                         double& Z, double& m, double& F,
                                         int N, int L, int R,
                                         int first_edit_pos, int first_edit_base,
                                         int second_edit_pos, int second_edit_base,
                                         int edits_so_far,
                                         size_t target_length, size_t motif_length,
                                         int direction)
{
    LseResult result;
    double T = static_cast<double>(m_threshold);

    // Apply prior edits (from exhaustive k=2 that failed)
    if (first_edit_pos >= 0) {
        apply_edit(first_edit_pos, first_edit_base, seq_bases, S_p, A_p, Z, m, N, L);
        F = std::log(Z) + m;
    }
    if (second_edit_pos >= 0) {
        apply_edit(second_edit_pos, second_edit_base, seq_bases, S_p, A_p, Z, m, N, L);
        F = std::log(Z) + m;
    }

    if (F >= T) {
        result.min_edits = static_cast<float>(edits_so_far);
        if (first_edit_pos >= 0) {
            result.best_edit_index = first_edit_pos;
            result.best_edit_direction = direction;
            result.best_pos = encode_position(static_cast<size_t>(first_edit_pos),
                                              target_length, motif_length, direction);
        }
        return result;
    }

    int first_greedy_edit_pos = -1;
    int k = edits_so_far;
    int max_k = (m_max_edits > 0) ? m_max_edits : R;  // Practical limit

    // Greedy: use a max-heap with stale-key checking
    // Heap entry: (delta_z, position, best_base)
    struct HeapEntry {
        double delta_z;
        int pos;
        int best_base;
        bool operator<(const HeapEntry& o) const { return delta_z < o.delta_z; }
    };

    // Initialize heap with DeltaZ for all positions
    std::priority_queue<HeapEntry> heap;
    std::vector<double> stored_dz(R, 0.0);

    for (int q = 0; q < R; q++) {
        int cur = seq_bases[q];
        double best_dz = -std::numeric_limits<double>::infinity();
        int best_y = -1;

        for (int y = 0; y < 4; y++) {
            if (y == cur) continue;
            double dz = compute_delta_z(q, y, seq_bases, A_p, N, L);
            if (dz > best_dz) {
                best_dz = dz;
                best_y = y;
            }
        }

        if (best_y >= 0 && best_dz > 0.0) {
            stored_dz[q] = best_dz;
            heap.push({best_dz, q, best_y});
        }
    }

    while (F < T && !heap.empty() && k < max_k) {
        HeapEntry top = heap.top();
        heap.pop();

        int q = top.pos;

        // Stale-key check: recompute DeltaZ
        int cur = seq_bases[q];
        double best_dz = -std::numeric_limits<double>::infinity();
        int best_y = -1;

        for (int y = 0; y < 4; y++) {
            if (y == cur) continue;
            double dz = compute_delta_z(q, y, seq_bases, A_p, N, L);
            if (dz > best_dz) {
                best_dz = dz;
                best_y = y;
            }
        }

        if (best_dz <= 0.0) {
            continue;  // No improvement possible at this position
        }

        // Check if key is stale (recomputed value significantly different)
        if (best_dz < top.delta_z - 1e-10 && !heap.empty() && best_dz < heap.top().delta_z) {
            // Push back with updated key
            stored_dz[q] = best_dz;
            heap.push({best_dz, q, best_y});
            continue;
        }

        // Apply edit
        apply_edit(q, best_y, seq_bases, S_p, A_p, Z, m, N, L);
        F = std::log(Z) + m;
        k++;

        if (first_greedy_edit_pos < 0) {
            first_greedy_edit_pos = q;
        }

        // Recompute DeltaZ in the (2L-1)-wide diamond around q
        int q_lo = std::max(0, q - L + 1);
        int q_hi = std::min(R - 1, q + L - 1);

        for (int qp = q_lo; qp <= q_hi; qp++) {
            int curp = seq_bases[qp];
            double best_dzp = -std::numeric_limits<double>::infinity();
            int best_yp = -1;

            for (int y = 0; y < 4; y++) {
                if (y == curp) continue;
                double dz = compute_delta_z(qp, y, seq_bases, A_p, N, L);
                if (dz > best_dzp) {
                    best_dzp = dz;
                    best_yp = y;
                }
            }

            if (best_yp >= 0 && best_dzp > 0.0) {
                stored_dz[qp] = best_dzp;
                heap.push({best_dzp, qp, best_yp});
            }
        }
    }

    if (F >= T) {
        result.min_edits = static_cast<float>(k);
        // Position: report the first edit position
        int report_pos = (first_edit_pos >= 0) ? first_edit_pos : first_greedy_edit_pos;
        if (report_pos >= 0) {
            result.best_edit_index = report_pos;
            result.best_edit_direction = direction;
            result.best_pos = encode_position(static_cast<size_t>(report_pos),
                                              target_length, motif_length, direction);
        }
    }

    return result;
}

PWMLseEditDistanceScorer::LseResult
PWMLseEditDistanceScorer::compute_lse_edit_distance(const std::string& seq,
                                                     size_t interval_length,
                                                     bool reverse)
{
    LseResult result;

    const int L = m_pssm.length();
    const size_t target_length = seq.length();
    if (L <= 0 || target_length < static_cast<size_t>(L) || interval_length == 0) {
        return result;
    }

    const size_t max_start = std::min(interval_length, target_length - static_cast<size_t>(L) + 1);
    if (max_start == 0) return result;

    const int N = static_cast<int>(max_start);  // Number of starts
    const int R = static_cast<int>(target_length);  // Sequence length
    const int direction = reverse ? -1 : 1;

    // Resolve sequence to base indices (handling reverse complement)
    std::vector<int> seq_bases(R);
    for (int q = 0; q < R; q++) {
        if (reverse) {
            char base = complement_base(seq[R - 1 - q]);
            seq_bases[q] = base_to_index(base);
        } else {
            seq_bases[q] = base_to_index(seq[q]);
        }
    }

    // Initialize partition function
    std::vector<double> S_p, A_p;
    double Z, m, F;
    init_partition(seq_bases, N, L, S_p, A_p, Z, m, F);

    // score.min filter: if current LSE score is below score.min, skip
    if (!std::isnan(m_score_min) && F < static_cast<double>(m_score_min)) {
        return result;
    }

    double T = static_cast<double>(m_threshold);

    // Already above threshold? 0 edits needed
    if (F >= T) {
        result.min_edits = 0.0f;
        result.best_pos = std::numeric_limits<float>::quiet_NaN();  // No edit needed
        return result;
    }

    // Quick unreachability check: S_max + log(N) is the maximum possible F
    double max_possible_F = static_cast<double>(m_S_max) + std::log(static_cast<double>(N));
    if (T > max_possible_F + 1e-9) {
        return result;  // Unreachable even with all bases optimal
    }

    // Check max_edits cap
    if (m_max_edits == 0) {
        return result;  // No edits allowed but F < T
    }

    // Try k=1 exhaustive
    LseResult k1 = try_k1_exhaustive(seq_bases, A_p, Z, m, F,
                                      N, L, target_length, L, direction);
    if (!std::isnan(k1.min_edits)) {
        return k1;
    }

    if (m_max_edits == 1) {
        return result;  // Only 1 edit allowed
    }

    // Try k=2 exhaustive
    LseResult k2 = try_k2_exhaustive(seq_bases, S_p, A_p, Z, m, F,
                                      N, L, R, target_length, L, direction);
    if (!std::isnan(k2.min_edits)) {
        return k2;
    }

    if (m_max_edits == 2) {
        return result;  // Only 2 edits allowed
    }

    // For k >= 3, use greedy heuristic
    // We need to find the best pair of first two edits to seed the greedy.
    // The greedy starts from scratch (no prior exhaustive edits applied),
    // since the exhaustive search for k=2 already failed.

    // Make copies for greedy (it modifies state)
    std::vector<int> seq_bases_copy = seq_bases;
    std::vector<double> S_p_copy = S_p;
    std::vector<double> A_p_copy = A_p;
    double Z_copy = Z, m_copy = m, F_copy = F;

    LseResult greedy = greedy_search(seq_bases_copy, S_p_copy, A_p_copy,
                                      Z_copy, m_copy, F_copy,
                                      N, L, R,
                                      -1, -1, -1, -1, 0,
                                      target_length, L, direction);

    return greedy;
}
