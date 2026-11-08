#include "PWMLseEditDistanceScorer.h"
#include "GenomeSeqFetch.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <queue>
#include <cstdint>

PWMLseEditDistanceScorer::PWMLseEditDistanceScorer(const DnaPSSM& pssm,
                                                     GenomeSeqFetch* shared_seqfetch,
                                                     float threshold,
                                                     int max_edits,
                                                     bool extend,
                                                     char strand,
                                                     Mode mode,
                                                     float score_min,
                                                     float score_max,
                                                     Direction direction)
    : GenomeSeqScorer(shared_seqfetch, extend, strand),
      m_pssm(pssm),
      m_threshold(threshold),
      m_max_edits(max_edits),
      m_mode(mode),
      m_direction(direction),
      m_score_min(score_min),
      m_score_max(score_max),
      m_last_min_edits(std::numeric_limits<float>::quiet_NaN()),
      m_S_max(0.0f),
      m_S_min(0.0f)
{
    precompute_tables();
}

void PWMLseEditDistanceScorer::precompute_tables()
{
    int L = m_pssm.length();
    m_col_max_scores.resize(L);
    m_col_min_scores.resize(L);
    m_S_max = 0.0f;
    m_S_min = 0.0f;
    m_logp.assign(L * 5, 0.0);

    for (int i = 0; i < L; i++) {
        float sum = 0.0f;
        for (int b = 0; b < 4; b++) {
            m_logp[i * 5 + b] = static_cast<double>(m_pssm[i].get_log_prob_from_code(b));
            sum += m_pssm[i].get_log_prob_from_code(b);
        }
        m_logp[i * 5 + 4] = static_cast<double>(sum / 4.0f);

        float max_score = -std::numeric_limits<float>::infinity();
        float min_score = std::numeric_limits<float>::infinity();
        for (int b = 0; b < 4; b++) {
            float score = m_pssm[i].get_log_prob_from_code(b);
            if (score > max_score) {
                max_score = score;
            }
            if (score < min_score) {
                min_score = score;
            }
        }
        m_col_max_scores[i] = max_score;
        m_col_min_scores[i] = min_score;
        m_S_max += max_score;
        m_S_min += min_score;
    }
}

double PWMLseEditDistanceScorer::score_interval(const GInterval& interval,
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

    if (scan_forward && scan_reverse) {
        // Both strands: the target is the strand-combined LSE (the pwm vtrack's value), and one substitution
        // changes both strands, so the two are searched together rather than each on its own.
        best_result = compute_lse_edit_distance(seq, static_cast<size_t>(interval_len), false, true);
    } else if (scan_forward) {
        LseResult fwd = compute_lse_edit_distance(seq, static_cast<size_t>(interval_len), false);
        if (!std::isnan(fwd.min_edits)) {
            if (std::isnan(best_result.min_edits) || fwd.min_edits < best_result.min_edits) {
                best_result = fwd;
            }
        }
    } else if (scan_reverse) {
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

double PWMLseEditDistanceScorer::encode_position(size_t index,
                                                  size_t target_length,
                                                  size_t motif_length,
                                                  int direction) const
{
    // The reverse pass scores a reversed-and-complemented COPY of the target
    // (see compute_lse_edits: seq_bases[q] = complement(seq[R-1-q])), so `index`
    // is an RC-array index there and maps to forward 1-based R - index. Gate on
    // the pass direction, not on m_strand: the reverse pass also runs under
    // bidirect = TRUE, where m_strand is not -1. `index` is a single edit
    // position, not a window start, so motif_length does not enter.
    (void)motif_length;
    double pos = static_cast<double>(index) + 1.0;
    if (direction == -1) {
        pos = static_cast<double>(target_length) - static_cast<double>(index);
    }
    if (m_pssm.is_bidirect()) {
        pos *= static_cast<double>(direction);
    }
    return pos;
}

void PWMLseEditDistanceScorer::init_partition(const std::vector<int>& seq_bases,
                                               int N, int L,
                                               std::vector<double>& S_p,
                                               std::vector<double>& A_p,
                                               double& Z, double& m, double& F) const
{
    const int W = m_ndirs * N;
    S_p.resize(W);
    A_p.resize(W);

    // Compute per-window scores. A window holding a non-ACGT base is not scored, as in the pwm vtrack: its
    // weight is 0 and stays 0 under any edit, since every edit multiplies it.
    for (int w = 0; w < W; w++) {
        const int p = w % N;
        double score = 0.0;
        for (int i = 0; i < L; i++) {
            if (seq_bases[p + i] == 4) {
                score = -std::numeric_limits<double>::infinity();
                break;
            }
            score += win_score(w, p + i, seq_bases[p + i], N, L);
        }
        S_p[w] = score;
    }

    // Find max for numerical stability
    m = -std::numeric_limits<double>::infinity();
    for (int w = 0; w < W; w++) {
        if (S_p[w] > m) m = S_p[w];
    }

    // No scorable window: F = -inf, and the caller reports NaN
    if (!std::isfinite(m)) {
        std::fill(A_p.begin(), A_p.end(), 0.0);
        Z = 0.0;
        F = m;
        return;
    }

    // Compute shifted weights and partition function
    Z = 0.0;
    for (int w = 0; w < W; w++) {
        A_p[w] = std::exp(S_p[w] - m);
        Z += A_p[w];
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

    // Every window that covers q: starts p in [q-L+1, q], clipped to [0, N-1], on each scanned strand
    const int p_min = std::max(0, q - L + 1);
    const int p_max = std::min(N - 1, q);

    // Forward windows: column q - p. DeltaZ contribution = A_w * expm1(delta_s); expm1 is stable for
    // delta_s near zero.
    const double *lp = m_logp.data();
    for (int p = p_min; p <= p_max; p++) {
        const double *col = lp + (q - p) * 5;
        delta_z += A_p[p] * std::expm1(col[y_idx] - col[cur_base]);
    }
    // Reverse-strand windows of a joint search: column L - 1 - (q - p), complemented base
    if (m_ndirs == 2) {
        const int cy = 3 - y_idx;
        const int cc = cur_base == 4 ? 4 : 3 - cur_base;
        for (int p = p_min; p <= p_max; p++) {
            const double *col = lp + (L - 1 - (q - p)) * 5;
            delta_z += A_p[N + p] * std::expm1(col[cy] - col[cc]);
        }
    }

    return delta_z;
}

double PWMLseEditDistanceScorer::edit_position(int q, int y, const std::vector<int>& seq_bases,
                                               const std::vector<double>& A_p, int N, int L, int direction) const
{
    const size_t R = seq_bases.size();
    if (direction != 0)
        return encode_position(static_cast<size_t>(q), R, static_cast<size_t>(L), direction);

    double dz[2] = {0.0, 0.0};
    const int p_min = std::max(0, q - L + 1);
    const int p_max = std::min(N - 1, q);
    for (int d = 0; d < m_ndirs; d++) {
        for (int p = p_min; p <= p_max; p++) {
            const int w = d * N + p;
            dz[d] += A_p[w] * std::expm1(win_score(w, q, y, N, L) - win_score(w, q, seq_bases[q], N, L));
        }
    }
    const double pos = static_cast<double>(q) + 1.0;
    return std::fabs(dz[1]) > std::fabs(dz[0]) ? -pos : pos;
}

void PWMLseEditDistanceScorer::apply_edit(int q, int y_idx,
                                           std::vector<int>& seq_bases,
                                           std::vector<double>& S_p,
                                           std::vector<double>& A_p,
                                           double& Z, double& m,
                                           int N, int L)
{
    const int old_base = seq_bases[q];
    seq_bases[q] = y_idx;

    const int p_min = std::max(0, q - L + 1);
    const int p_max = std::min(N - 1, q);

    for (int d = 0; d < m_ndirs; d++) {
        for (int p = p_min; p <= p_max; p++) {
            const int w = d * N + p;
            const double delta_s = win_score(w, q, y_idx, N, L) - win_score(w, q, old_base, N, L);
            const double A_old = A_p[w];
            S_p[w] += delta_s;
            A_p[w] = A_old * std::exp(delta_s);
            Z += A_old * std::expm1(delta_s);
        }
    }

    // Renormalize if drift is too large (threshold for double: ~700)
    double max_sp = -std::numeric_limits<double>::infinity();
    for (int d = 0; d < m_ndirs; d++) {
        for (int p = p_min; p <= p_max; p++) {
            if (S_p[d * N + p] > max_sp) max_sp = S_p[d * N + p];
        }
    }

    if (max_sp - m > 500.0) {
        // Full renormalization
        const int W = m_ndirs * N;
        double new_m = -std::numeric_limits<double>::infinity();
        for (int w = 0; w < W; w++) {
            if (S_p[w] > new_m) new_m = S_p[w];
        }
        double scale = std::exp(m - new_m);
        for (int w = 0; w < W; w++) {
            A_p[w] *= scale;
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
    const bool is_below = (m_direction == Direction::BELOW);

    // For ABOVE: track best (highest) F after edit
    // For BELOW: track best (lowest) F after edit
    double best_F_after = is_below ? std::numeric_limits<double>::infinity()
                                   : -std::numeric_limits<double>::infinity();
    int best_q = -1;
    int best_y = -1;

    for (int q = 0; q < R; q++) {
        int cur = seq_bases[q];
        for (int y = 0; y < 4; y++) {
            if (y == cur) continue;

            double dz = compute_delta_z(q, y, seq_bases, A_p, N, L);
            double new_Z = Z + dz;
            if (new_Z <= 0.0) continue;

            double new_F = std::log(new_Z) + m;
            bool meets_threshold = is_below ? (new_F <= T) : (new_F >= T);
            if (meets_threshold) {
                // Found a valid 1-edit solution
                bool is_better = is_below ? (new_F < best_F_after) : (new_F > best_F_after);
                if (is_better) {
                    best_F_after = new_F;
                    best_q = q;
                    best_y = y;
                }
            }
        }
    }

    if (best_q >= 0) {
        result.min_edits = 1.0f;
        result.best_edit_index = best_q;
        result.best_edit_direction = direction;
        result.best_pos = edit_position(best_q, best_y, seq_bases, A_p, N, L, direction);
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
    const bool is_below = (m_direction == Direction::BELOW);

    double best_F_after = is_below ? std::numeric_limits<double>::infinity()
                                   : -std::numeric_limits<double>::infinity();
    int best_q1 = -1;
    int best_y1 = -1;

    std::vector<double> A_p_temp;
    std::vector<int> bases_temp(seq_bases);

    // For each first edit (q1, y1)
    for (int q1 = 0; q1 < R; q1++) {
        int cur1 = seq_bases[q1];

        for (int y1 = 0; y1 < 4; y1++) {
            if (y1 == cur1) continue;

            // Apply the first edit to a copy of the weights, then check all second edits
            const int p_min1 = std::max(0, q1 - L + 1);
            const int p_max1 = std::min(N - 1, q1);

            A_p_temp = A_p;
            double Z_temp = Z;

            for (int d = 0; d < m_ndirs; d++) {
                for (int p = p_min1; p <= p_max1; p++) {
                    const int w = d * N + p;
                    const double ds = win_score(w, q1, y1, N, L) - win_score(w, q1, cur1, N, L);
                    const double A_old = A_p_temp[w];
                    A_p_temp[w] = A_old * std::exp(ds);
                    Z_temp += A_old * std::expm1(ds);
                }
            }

            // Now try all second edits (q2, y2) with q2 > q1
            // (q2 < q1 was already covered when q1 was smaller); q2 != q1, so the
            // base at q2 is still the original one
            bases_temp[q1] = y1;
            for (int q2 = q1 + 1; q2 < R; q2++) {
                for (int y2 = 0; y2 < 4; y2++) {
                    if (y2 == seq_bases[q2]) continue;

                    double dz2 = compute_delta_z(q2, y2, bases_temp, A_p_temp, N, L);
                    double new_Z = Z_temp + dz2;
                    if (new_Z <= 0.0) continue;

                    double new_F = std::log(new_Z) + m;
                    bool meets_threshold = is_below ? (new_F <= T) : (new_F >= T);
                    if (meets_threshold) {
                        bool is_better = is_below ? (new_F < best_F_after) : (new_F > best_F_after);
                        if (is_better) {
                            best_F_after = new_F;
                            best_q1 = q1;
                            best_y1 = y1;
                        }
                    }
                }
            }
            bases_temp[q1] = cur1;
        }
    }

    if (best_q1 >= 0) {
        result.min_edits = 2.0f;
        result.best_edit_index = best_q1;
        result.best_edit_direction = direction;
        result.best_pos = edit_position(best_q1, best_y1, seq_bases, A_p, N, L, direction);
    }

    return result;
}

namespace {
// Hash of one substitution (position q to base y), XORed into the hash of an edit set
inline uint64_t edit_hash(int q, int y)
{
    uint64_t z = static_cast<uint64_t>(q) * 4u + static_cast<uint64_t>(y) + 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
}

PWMLseEditDistanceScorer::LseResult
PWMLseEditDistanceScorer::beam_search(const std::vector<int>& seq_bases,
                                       const std::vector<double>& S_p,
                                       const std::vector<double>& A_p,
                                       double Z, double m,
                                       int N, int L, int R,
                                       int direction)
{
    LseResult result;
    const double T = static_cast<double>(m_threshold);
    const bool is_below = (m_direction == Direction::BELOW);
    // F moves towards the threshold when sgn * F grows
    const double sgn = is_below ? -1.0 : 1.0;
    auto passes = [&](double F) { return is_below ? (F <= T) : (F >= T); };
    const int max_k = (m_max_edits > 0) ? m_max_edits : R;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // Edit sets of d substitutions at d distinct positions, the BEAM_WIDTH best by F at each d. Unlike a
    // single greedy path, the beam keeps sets that score worse than their parent, which is how it gets
    // past a local optimum (windows share bases, and under bidirect a substitution that helps one strand
    // can hurt the other).
    static constexpr int BEAM_WIDTH = 16;

    // dz[q * 4 + y]: change in Z of substituting y at q; NaN where no edit is allowed (the position was
    // already edited, holds a non-ACGT base whose windows are unscored and stay so, or y is its base).
    // An edit at q changes the weights of the windows covering q only, so a child inherits its parent's
    // table and recomputes the positions within L - 1 of the edit.
    struct State {
        std::vector<int> bases;
        std::vector<double> S, A, dz;
        double Z, m;
        uint64_t hash;
        double first_pos;  // reported position of its first edit
    };

    auto fill_dz = [&](State& st, int q) {
        for (int y = 0; y < 4; y++) {
            st.dz[q * 4 + y] = (st.bases[q] != seq_bases[q] || seq_bases[q] == 4 || y == seq_bases[q])
                                   ? nan
                                   : compute_delta_z(q, y, st.bases, st.A, N, L);
        }
    };

    std::vector<State> beam(1);
    {
        State& st = beam[0];
        st.bases = seq_bases;
        st.S = S_p;
        st.A = A_p;
        st.Z = Z;
        st.m = m;
        st.hash = 0;
        st.first_pos = nan;
        st.dz.assign(static_cast<size_t>(R) * 4, nan);
        for (int q = 0; q < R; q++) fill_dz(st, q);
    }

    struct Child {
        double key;  // sgn * F
        int parent;
        int q, y;
        uint64_t hash;
    };
    std::vector<Child> children;
    auto by_key = [](const Child& a, const Child& b) { return a.key > b.key; };

    double best_key = sgn * (std::log(Z) + m);
    int last_improved = 0;

    for (int d = 1; d <= max_k; d++) {
        children.clear();
        for (int si = 0; si < static_cast<int>(beam.size()); si++) {
            const State& st = beam[si];
            for (int e = 0; e < R * 4; e++) {
                const double dz = st.dz[e];
                if (std::isnan(dz)) continue;
                const double new_Z = st.Z + dz;
                if (new_Z <= 0.0) continue;
                const int q = e / 4, y = e % 4;
                children.push_back({sgn * (std::log(new_Z) + st.m), si, q, y, st.hash ^ edit_hash(q, y)});
            }
        }
        if (children.empty()) break;

        // The best few, enough to find BEAM_WIDTH distinct sets (the same set can come from two parents)
        size_t n_top = std::min(children.size(), static_cast<size_t>(4 * BEAM_WIDTH));
        std::partial_sort(children.begin(), children.begin() + n_top, children.end(), by_key);

        // The best child is the best set of d edits found: done if it reaches the threshold
        const Child& top = children.front();
        if (passes(sgn * top.key)) {
            const State& par = beam[top.parent];
            result.min_edits = static_cast<float>(d);
            result.best_edit_direction = direction;
            if (std::isnan(par.first_pos)) {
                result.best_edit_index = top.q;
                result.best_pos = edit_position(top.q, top.y, par.bases, par.A, N, L, direction);
            } else {
                result.best_pos = par.first_pos;
            }
            return result;
        }

        // Stop when the best set has not moved F towards the threshold for L more edits: no window holds
        // more than L bases, so a run that long without progress is past what the search can still cross
        if (top.key > best_key + 1e-12) {
            best_key = top.key;
            last_improved = d;
        } else if (d - last_improved >= L) {
            break;
        }

        // Next beam: the best distinct sets
        std::vector<State> next;
        next.reserve(BEAM_WIDTH);
        std::vector<uint64_t> seen;
        for (size_t ci = 0; ci < children.size() && static_cast<int>(next.size()) < BEAM_WIDTH; ci++) {
            if (ci == n_top) {
                // more duplicates than expected: order the rest too
                std::sort(children.begin() + n_top, children.end(), by_key);
                n_top = children.size();
            }
            const Child& c = children[ci];
            if (std::find(seen.begin(), seen.end(), c.hash) != seen.end()) continue;
            seen.push_back(c.hash);
            const State& par = beam[c.parent];
            State st = par;
            if (std::isnan(st.first_pos)) {
                st.first_pos = edit_position(c.q, c.y, par.bases, par.A, N, L, direction);
            }
            const double m_before = st.m;
            apply_edit(c.q, c.y, st.bases, st.S, st.A, st.Z, st.m, N, L);
            st.hash = c.hash;
            // apply_edit may rescale every weight (a new shift m); the inherited dz scale with them
            if (st.m != m_before) {
                const double scale = std::exp(m_before - st.m);
                for (double& v : st.dz) v *= scale;
            }
            const int q_lo = std::max(0, c.q - L + 1);
            const int q_hi = std::min(R - 1, c.q + L - 1);
            for (int q = q_lo; q <= q_hi; q++) fill_dz(st, q);
            next.push_back(std::move(st));
        }
        beam.swap(next);
    }

    return result;
}

PWMLseEditDistanceScorer::LseResult
PWMLseEditDistanceScorer::compute_lse_edit_distance(const std::string& seq,
                                                     size_t interval_length,
                                                     bool reverse,
                                                     bool joint)
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
    // A joint search reads both strands on the forward sequence (direction 0); a one-strand search reads
    // its strand as the forward strand of the sequence or of its reverse complement.
    m_ndirs = joint ? 2 : 1;
    const int direction = joint ? 0 : (reverse ? -1 : 1);
    if (joint) reverse = false;

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

    // No window without a non-ACGT base: nothing to edit towards or away from the threshold
    if (!std::isfinite(F)) {
        return result;
    }

    // score.min/score.max filter: if current LSE score is out of range, skip
    if (!std::isnan(m_score_min) && F < static_cast<double>(m_score_min)) {
        return result;
    }
    if (!std::isnan(m_score_max) && F > static_cast<double>(m_score_max)) {
        return result;
    }

    double T = static_cast<double>(m_threshold);
    const bool is_below = (m_direction == Direction::BELOW);

    // Already at target? 0 edits needed
    if (is_below ? (F <= T) : (F >= T)) {
        result.min_edits = 0.0f;
        result.best_pos = std::numeric_limits<double>::quiet_NaN();  // No edit needed
        return result;
    }

    // Quick unreachability check. Under any substitutions every scorable window scores between S_min and
    // S_max, and a window holding a non-ACGT base keeps weight 0, so with n scorable windows
    // F = log(sum(exp(S_p))) lies in [S_min + log(n), S_max + log(n)].
    int n_scorable = 0;
    for (double sp : S_p) {
        if (std::isfinite(sp)) n_scorable++;
    }
    const double log_n = std::log(static_cast<double>(n_scorable));
    if (is_below) {
        double min_possible_F = static_cast<double>(m_S_min) + log_n;
        if (T < min_possible_F - 1e-6) {
            return result;  // Unreachable: F can never go below S_min + log(n)
        }
    } else {
        double max_possible_F = static_cast<double>(m_S_max) + log_n;
        if (T > max_possible_F + 1e-6) {
            return result;  // Unreachable even with all bases optimal
        }
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

    // For k >= 3, a beam search over edit sets (heuristic: the count is an upper bound on the minimum,
    // and NaN means it found no solution)
    return beam_search(seq_bases, S_p, A_p, Z, m, N, L, R, direction);
}
