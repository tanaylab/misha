#ifndef POTTS_MODEL_H_
#define POTTS_MODEL_H_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

// DnaPSSM.h uses bare vector/string/ostream and assumes "using namespace std"
// is already active - normally supplied by port.h, which config.h (the header
// DnaPSSM.h pulls in itself) does NOT provide. Any translation unit that
// reaches DnaPSSM.h only through this header (as GseqPotts.cpp and the
// TrackExpressionVars potts branch do) needs it declared here, before the
// include, not after.
using namespace std;

#include "DnaPSSM.h" // DnaLookupTables::BASE_ENCODE

// A pairwise (Potts) energy model over a fixed-width window:
//
//   S(x) = intercept + sum_i e[i][x_i] + sum_k J_k[x_{p1k}][x_{p2k}]
//
// score_codes() must agree with the tests' pure-R oracle implementation of
// this scoring rule to float precision. Tables and accumulator are double:
// at W = 20 the J table is 24 KB and still L2-resident, and 211 float adds at
// magnitude ~20 would cost ~4e-4 of error for no measured speed gain.
//
// THE REVERSE STRAND IS A SECOND MODEL, not a reversed sequence. Scoring the
// reverse strand at an anchor means scoring the model on comp(x_{W-1-i}), and
// reindexing that gives a Potts of the same shape - see rc(). Building it once
// makes the reverse strand cost exactly what the forward strand costs, where
// reverse-complementing per anchor would dominate the loop.
//
// PottsModel has two kernels for the same score:
//
//   score_codes_naive()   - one table lookup per single-site term and one per
//                            pair term: 1 + W + npair lookups.
//   score_codes_blocked() - positions packed two at a time into a 0..15
//                            "block code" (once per window, for both
//                            strands, by potts_block_codes()); one table per
//                            block folds its own singles (and the pair
//                            within the block, if any)
//                            into one lookup, and one table per BLOCK PAIR
//                            that some coupling links (every block pair,
//                            from 80% linked) folds its up to 4 cross-block
//                            couplings into one lookup: ceil(W/2) + (block
//                            pairs with a table) lookups. At
//                            W = 20, full pairwise, that is 10 + 45 = 55
//                            instead of 211; for a sum of four side-by-side
//                            models at W = 41 (253 pairs), 21 + 73 = 94
//                            instead of 295.
//
// score_codes() picks one of the two per model (m_use_blocked, set once in
// build_blocked_tables()). The blocked kernel is 2-5x faster on a model whose
// couplings are local or dense, and faster on an order-1 model, but up to
// 1.5x slower on a wide model whose couplings are scattered so that each
// needs a block pair of its own - see the gate in build_blocked_tables().
// The naive kernel is also the fallback for W > 256 and the equivalence
// tests' reference.
//
// build_blocked_tables() generalizes to odd W (a trailing size-1 block) and
// to sparse or absent pairs (a missing coupling contributes 0 to its block
// pair's table, and a block pair that no coupling links has no table, or an
// all-zero one from 80% linked), so score_codes_blocked() agrees with
// score_codes_naive() for ANY
// PottsModel - not only the dense, even-W case it is designed to win on. For
// W > 256 (MAX_BLOCKS = 128), build_blocked_tables() leaves the blocked
// tables unbuilt (blocked_available() is false) and score_codes() falls back
// to the naive kernel.
class PottsModel {
public:
    PottsModel() = default;

    // e_rowmajor: W*4 doubles, row i = position i, columns A, C, G, T.
    // j_rowmajor: npair*16 doubles, row k = J_k with element [a][b] at b*4 + a
    //             (column-major per 4x4 block, i.e. what as.numeric() of the
    //             block gives in R).
    // p1, p2:     0-BASED position indices, p1[k] < p2[k]. A (p1, p2) pair may
    //             repeat - see build_blocked_tables() - and every kernel below
    //             treats a repeat as an ADDITIONAL coupling term, summed in,
    //             exactly as score_codes_naive()'s loop over k naturally does.
    PottsModel(int W, double intercept,
               const double *e_rowmajor, const double *j_rowmajor,
               const std::vector<int> &p1, const std::vector<int> &p2)
        : m_W(W), m_intercept(intercept), m_p1(p1), m_p2(p2)
    {
        m_e.assign(e_rowmajor, e_rowmajor + (std::size_t)W * 4);
        // An order-1 model (no pairs - a per-position table with no couplings, and
        // test 1's npair_mode = "none" exercises it) passes p1.empty() and,
        // from PottsParams.h, a null j_rowmajor. Forming a [ptr, ptr+0) range
        // out of a null pointer is formally UB even though libstdc++ treats it
        // as a benign no-op, so this is guarded rather than relied upon.
        if (!p1.empty())
            m_J.assign(j_rowmajor, j_rowmajor + p1.size() * 16);
        build_blocked_tables();
        build_pair_index();
    }

    int width() const { return m_W; }
    std::size_t npair() const { return m_p1.size(); }
    double intercept() const { return m_intercept; }

    // Whether build_blocked_tables() actually built usable tables (false only
    // for W > 256 - see the class comment). score_codes_blocked() requires
    // this; score_codes() checks it internally, any other caller must check
    // it itself before calling score_codes_blocked() directly.
    bool blocked_available() const { return m_nblocks > 0; }

    // Whether score_codes() runs the blocked kernel for this model - the
    // gate in build_blocked_tables(). C_potts_score_codes_cmp reports it so
    // that tests can check the gate rather than restate it.
    bool use_blocked() const { return m_use_blocked; }

    // The reference kernel: one table lookup per single-site term and one per
    // pair term, in whatever order pairs were given. 1 + W + npair lookups.
    inline double score_codes_naive(const int8_t *c) const
    {
        double s = m_intercept;
        for (int i = 0; i < m_W; ++i)
            s += m_e[(std::size_t)i * 4 + c[i]];
        const std::size_t np = m_p1.size();
        for (std::size_t k = 0; k < np; ++k)
            s += m_J[k * 16 + (std::size_t)c[m_p2[k]] * 4 + (std::size_t)c[m_p1[k]]];
        return s;
    }

    // The dinucleotide-blocked kernel - see the class comment. bc points at
    // the window's first block code from potts_block_codes(): block b's code is
    // bc[2 * b]. The caller must guarantee blocked_available() (score_codes()
    // checks it; a direct caller, such as the equivalence test's
    // C_potts_score_codes_cmp, must check it itself).
    inline double score_codes_blocked(const int8_t *bc) const
    {
        double s = m_intercept;
        for (int b = 0; b < m_nblocks; ++b)
            s += m_block_table[(std::size_t)b * 16 + (std::size_t)bc[2 * b]];

        // Only the block pairs some coupling links, in (u, v) order; pair k's
        // 16 x 16 table starts at k * 256. Any other block pair's table would
        // be all zeros, and adding 0 leaves the sum as it was. When every
        // block pair has a table (see build_blocked_tables()), u and v come
        // from the loop counters instead of m_pair_u / m_pair_v - the same
        // terms in the same order, and measured 6-30% faster on full
        // pairwise models (W = 20 to 64).
        if (m_pair_u.size() == (std::size_t)m_nblocks * (std::size_t)(m_nblocks - 1) / 2) {
            std::size_t k = 0;
            for (int u = 0; u < m_nblocks; ++u) {
                const std::size_t cu = (std::size_t)bc[2 * u] * 16;
                for (int v = u + 1; v < m_nblocks; ++v, ++k)
                    s += m_pair_table[k * 256 + cu + (std::size_t)bc[2 * v]];
            }
        } else {
            for (std::size_t k = 0; k < m_pair_u.size(); ++k)
                s += m_pair_table[k * 256 + (std::size_t)bc[m_pair_u[k]] * 16 + (std::size_t)bc[m_pair_v[k]]];
        }
        return s;
    }

    // Production entry point - GseqPotts.cpp and PottsParams' consumers call
    // this and never the two kernels above directly. Runs the kernel that
    // build_blocked_tables()'s gate picked for THIS model. c and bc are the
    // window's codes and block codes (potts_encode(), potts_block_codes()).
    inline double score_codes(const int8_t *c, const int8_t *bc) const
    {
        return m_use_blocked ? score_codes_blocked(bc) : score_codes_naive(c);
    }

    // score_codes() of the four anchors at c, c + 1, c + 2 and c + 3 (and bc
    // likewise), every window scorable. Each anchor adds the same terms in the
    // same order as score_codes(), so out[k] is bit-identical to
    // score_codes(c + k, bc + k); the four sums only run side by side. One
    // anchor's sum is a single chain of dependent additions, each waiting for
    // the one before it; four chains fill those waits.
    inline void score_codes4(const int8_t *c, const int8_t *bc, double out[4]) const
    {
        double s0 = m_intercept, s1 = m_intercept, s2 = m_intercept, s3 = m_intercept;
        if (m_use_blocked) {
            for (int b = 0; b < m_nblocks; ++b) {
                const double *t = &m_block_table[(std::size_t)b * 16];
                const int8_t *p = bc + 2 * b;
                s0 += t[p[0]];
                s1 += t[p[1]];
                s2 += t[p[2]];
                s3 += t[p[3]];
            }
            if (m_pair_u.size() == (std::size_t)m_nblocks * (std::size_t)(m_nblocks - 1) / 2) {
                std::size_t k = 0;
                for (int u = 0; u < m_nblocks; ++u) {
                    const int8_t *pu = bc + 2 * u;
                    const std::size_t c0 = (std::size_t)pu[0] * 16, c1 = (std::size_t)pu[1] * 16,
                                      c2 = (std::size_t)pu[2] * 16, c3 = (std::size_t)pu[3] * 16;
                    for (int v = u + 1; v < m_nblocks; ++v, ++k) {
                        const double *t = &m_pair_table[k * 256];
                        const int8_t *pv = bc + 2 * v;
                        s0 += t[c0 + pv[0]];
                        s1 += t[c1 + pv[1]];
                        s2 += t[c2 + pv[2]];
                        s3 += t[c3 + pv[3]];
                    }
                }
            } else {
                for (std::size_t k = 0; k < m_pair_u.size(); ++k) {
                    const double *t = &m_pair_table[k * 256];
                    const int8_t *pu = bc + m_pair_u[k];
                    const int8_t *pv = bc + m_pair_v[k];
                    s0 += t[(std::size_t)pu[0] * 16 + pv[0]];
                    s1 += t[(std::size_t)pu[1] * 16 + pv[1]];
                    s2 += t[(std::size_t)pu[2] * 16 + pv[2]];
                    s3 += t[(std::size_t)pu[3] * 16 + pv[3]];
                }
            }
        } else {
            for (int i = 0; i < m_W; ++i) {
                const double *e = &m_e[(std::size_t)i * 4];
                s0 += e[c[i]];
                s1 += e[c[i + 1]];
                s2 += e[c[i + 2]];
                s3 += e[c[i + 3]];
            }
            const std::size_t np = m_p1.size();
            for (std::size_t k = 0; k < np; ++k) {
                const double *J = &m_J[k * 16];
                const int8_t *a = c + m_p1[k];
                const int8_t *b = c + m_p2[k];
                s0 += J[(std::size_t)b[0] * 4 + a[0]];
                s1 += J[(std::size_t)b[1] * 4 + a[1]];
                s2 += J[(std::size_t)b[2] * 4 + a[2]];
                s3 += J[(std::size_t)b[3] * 4 + a[3]];
            }
        }
        out[0] = s0;
        out[1] = s1;
        out[2] = s2;
        out[3] = s3;
    }

    // Change in score_codes(c) when position j is changed from c[j] to base b
    // (0..3), every other position as in c. Touches only the field at j and
    // the pairs that contain j (m_pairs_at), so it costs 1 + deg(j) lookups
    // instead of a full rescore.
    inline double delta_codes(const int8_t *c, int j, int b) const
    {
        const int x = c[j];
        double d = m_e[(std::size_t)j * 4 + b] - m_e[(std::size_t)j * 4 + x];
        for (const int k : m_pairs_at[j]) {
            const double *Jk = &m_J[(std::size_t)k * 16];
            if (m_p1[k] == j) {
                const std::size_t row = (std::size_t)c[m_p2[k]] * 4;
                d += Jk[row + b] - Jk[row + x];
            } else {
                const int a = c[m_p1[k]];
                d += Jk[(std::size_t)b * 4 + a] - Jk[(std::size_t)x * 4 + a];
            }
        }
        return d;
    }

    // The largest delta_codes(c, j, b) over every position j and every base b
    // other than c[j]; -inf when W == 0. Each (j, b) is computed with the same
    // terms, in the same order and with the same arithmetic as delta_codes(),
    // only four bases at a time and without its per-pair branch, so the value
    // is bit-identical to the largest delta_codes(). Since rounding f + d is
    // monotone in d, f + best_delta(c) is also the largest f + delta_codes()
    // exactly - what potts.max with max_edits = 1 needs on one strand.
    //
    // Two positions are summed side by side, eight running sums in registers:
    // one position alone makes every addition wait for the one before it.
    // Each sum still adds its pairs in m_pairs_at[j]'s order.
    inline double best_delta(const int8_t *c) const
    {
        double best = -std::numeric_limits<double>::infinity();
        int j = 0;
        for (; j + 1 < m_W; j += 2) {
            const int x = c[j];
            const int y = c[j + 1];
            const double *ej = &m_e[(std::size_t)j * 4];
            const double *ek = &m_e[(std::size_t)(j + 1) * 4];
            const double ex = ej[x], ey = ek[y];
            double a0 = ej[0] - ex, a1 = ej[1] - ex, a2 = ej[2] - ex, a3 = ej[3] - ex;
            double b0 = ek[0] - ey, b1 = ek[1] - ey, b2 = ek[2] - ey, b3 = ek[3] - ey;
            std::size_t p = m_field_start[j];
            const std::size_t pe = m_field_start[(std::size_t)j + 1];
            std::size_t q = pe;
            const std::size_t qe = m_field_start[(std::size_t)j + 2];
            for (; p < pe && q < qe; ++p, ++q) {
                const double *r = &m_field_tab[p * 16 + (std::size_t)c[m_field_other[p]] * 4];
                const double *s = &m_field_tab[q * 16 + (std::size_t)c[m_field_other[q]] * 4];
                const double rx = r[x], sy = s[y];
                a0 += r[0] - rx;
                a1 += r[1] - rx;
                a2 += r[2] - rx;
                a3 += r[3] - rx;
                b0 += s[0] - sy;
                b1 += s[1] - sy;
                b2 += s[2] - sy;
                b3 += s[3] - sy;
            }
            for (; p < pe; ++p) {
                const double *r = &m_field_tab[p * 16 + (std::size_t)c[m_field_other[p]] * 4];
                const double rx = r[x];
                a0 += r[0] - rx;
                a1 += r[1] - rx;
                a2 += r[2] - rx;
                a3 += r[3] - rx;
            }
            for (; q < qe; ++q) {
                const double *s = &m_field_tab[q * 16 + (std::size_t)c[m_field_other[q]] * 4];
                const double sy = s[y];
                b0 += s[0] - sy;
                b1 += s[1] - sy;
                b2 += s[2] - sy;
                b3 += s[3] - sy;
            }
            best = max_other(best, x, a0, a1, a2, a3);
            best = max_other(best, y, b0, b1, b2, b3);
        }
        if (j < m_W) {
            const int x = c[j];
            const double *ej = &m_e[(std::size_t)j * 4];
            const double ex = ej[x];
            double a0 = ej[0] - ex, a1 = ej[1] - ex, a2 = ej[2] - ex, a3 = ej[3] - ex;
            const std::size_t pe = m_field_start[(std::size_t)j + 1];
            for (std::size_t p = m_field_start[j]; p < pe; ++p) {
                const double *r = &m_field_tab[p * 16 + (std::size_t)c[m_field_other[p]] * 4];
                const double rx = r[x];
                a0 += r[0] - rx;
                a1 += r[1] - rx;
                a2 += r[2] - rx;
                a3 += r[3] - rx;
            }
            best = max_other(best, x, a0, a1, a2, a3);
        }
        return best;
    }

    // The complemented twin: rc().score_codes(w) == score_codes(revcomp(w)).
    //
    //   e'[i][b]         = e[W-1-i][comp(b)]
    //   pairs'           = (W-1-p2, W-1-p1), which keeps u < v since p1 < p2
    //   J'_{(u,v)}[a][b] = J_k[comp(b)][comp(a)]
    //
    // Pairs are NOT re-sorted: score_codes() iterates the pair list in whatever
    // order it is given, so the order only has to be self-consistent with m_J.
    PottsModel rc() const
    {
        static const int comp[4] = {3, 2, 1, 0}; // A<->T, C<->G in A,C,G,T order
        PottsModel o;
        o.m_W = m_W;
        o.m_intercept = m_intercept;
        o.m_e.assign((std::size_t)m_W * 4, 0.0);
        for (int i = 0; i < m_W; ++i)
            for (int b = 0; b < 4; ++b)
                o.m_e[(std::size_t)i * 4 + b] = m_e[(std::size_t)(m_W - 1 - i) * 4 + comp[b]];

        const std::size_t np = m_p1.size();
        o.m_p1.resize(np);
        o.m_p2.resize(np);
        o.m_J.assign(np * 16, 0.0);
        for (std::size_t k = 0; k < np; ++k) {
            o.m_p1[k] = m_W - 1 - m_p2[k];
            o.m_p2[k] = m_W - 1 - m_p1[k];
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    o.m_J[k * 16 + (std::size_t)b * 4 + a] =
                        m_J[k * 16 + (std::size_t)comp[a] * 4 + comp[b]];
        }
        o.build_blocked_tables();
        o.build_pair_index();
        return o;
    }

private:
    // The largest of best and the d_b with b != x.
    static inline double max_other(double best, int x, double d0, double d1, double d2, double d3)
    {
        if (x != 0 && d0 > best)
            best = d0;
        if (x != 1 && d1 > best)
            best = d1;
        if (x != 2 && d2 > best)
            best = d2;
        if (x != 3 && d3 > best)
            best = d3;
        return best;
    }

    // Supports W up to 256. build_blocked_tables() disables the blocked
    // kernel beyond that rather than build tables that grow as W^2 (jsum
    // below is 8 MB at W = 256) - no PottsModel using it today gets remotely
    // close (every test width, and a realistic motif width, is well under
    // 100).
    static constexpr int MAX_BLOCKS = 128;

    // m_pairs_at[i]: indices k of the pairs that contain position i. Also the
    // tables best_delta() reads: for position j, one entry per pair of
    // m_pairs_at[j], in that order, holding the pair's other position and its
    // 4x4 table turned so that row (other base) * 4 + (base at j) is the
    // coupling, whichever end of the pair j is.
    void build_pair_index()
    {
        m_pairs_at.assign((std::size_t)m_W, std::vector<int>());
        for (std::size_t k = 0; k < m_p1.size(); ++k) {
            m_pairs_at[m_p1[k]].push_back((int)k);
            m_pairs_at[m_p2[k]].push_back((int)k);
        }
        m_field_start.assign((std::size_t)m_W + 1, 0);
        m_field_other.clear();
        m_field_tab.clear();
        for (int j = 0; j < m_W; ++j) {
            m_field_start[j] = m_field_other.size();
            for (const int k : m_pairs_at[j]) {
                const double *Jk = &m_J[(std::size_t)k * 16];
                const bool first = (m_p1[k] == j);
                m_field_other.push_back(first ? m_p2[k] : m_p1[k]);
                for (int o = 0; o < 4; ++o)
                    for (int b = 0; b < 4; ++b)
                        // delta_codes() reads Jk[o * 4 + b] when j is p1 (o the base at p2)
                        // and Jk[b * 4 + o] when j is p2 (o the base at p1)
                        m_field_tab.push_back(first ? Jk[o * 4 + b] : Jk[b * 4 + o]);
            }
        }
        m_field_start[m_W] = m_field_other.size();

    }

    void build_blocked_tables()
    {
        m_nblocks = (m_W + 1) / 2; // ceil(W / 2)
        if (m_nblocks == 0 || m_nblocks > MAX_BLOCKS) {
            m_nblocks = 0; // blocked_available() is false; score_codes() uses naive
            return;
        }

        // Sum every k's 4x4 block into jsum[(i*W+j)*16 + b*4+a], i < j. A
        // (p1, p2) pair repeated across multiple k (rejected by R's
        // .coerce_potts_model() today, but PottsParams::parse() is a .Call
        // trust boundary and should not lean on that) accumulates here
        // exactly as score_codes_naive()'s loop over k does - there is no
        // "does a pair exist" flag to get out of sync with that sum, only a
        // table that starts at 0 and has every k's contribution added once.
        vector<double> jsum((std::size_t)m_W * (std::size_t)m_W * 16, 0.0);
        for (std::size_t k = 0; k < m_p1.size(); ++k) {
            const std::size_t dst = ((std::size_t)m_p1[k] * m_W + (std::size_t)m_p2[k]) * 16;
            const std::size_t src = k * 16;
            for (int f = 0; f < 16; ++f)
                jsum[dst + f] += m_J[src + f];
        }
        auto pair_j = [&](int i, int j, int xi, int xj) -> double {
            return jsum[((std::size_t)i * m_W + j) * 16 + (std::size_t)xj * 4 + xi];
        };

        // Per-block tables, 16 entries each: the block's own single-site
        // terms, plus the within-block pair if the model has one. A block
        // code is 4 * (first base) + (second base). The trailing block of an
        // odd W has one position, and its table ignores the second base
        // (whatever follows the window, or 0 - see potts_block_codes()).
        m_block_table.assign((std::size_t)m_nblocks * 16, 0.0);
        for (int b = 0; b < m_nblocks; ++b) {
            const int p0 = 2 * b;
            const bool paired = p0 + 1 < m_W;
            for (int code = 0; code < 16; ++code) {
                const int x0 = code >> 2, x1 = code & 3;
                double v = m_e[(std::size_t)p0 * 4 + x0];
                if (paired) {
                    v += m_e[(std::size_t)(p0 + 1) * 4 + x1];
                    v += pair_j(p0, p0 + 1, x0, x1);
                }
                m_block_table[(std::size_t)b * 16 + (std::size_t)code] = v;
            }
        }

        // Per-block-pair tables: up to 4 cross terms between the (1 or 2)
        // positions of block u and the (1 or 2) positions of block v. Unless
        // most block pairs are linked (below), only a block pair that some
        // coupling links (a position of u with a position of v) gets a
        // table: any other block pair's table would be all zeros. A model
        // whose couplings are local links few of its block pairs - a sum of
        // four side-by-side models at W = 41 links 73 of its 210 - and the
        // zero tables were most of the kernel's lookups. Every table is
        // 16 x 16; when v is the trailing size-1 block, its column ignores
        // the second base, as the block table does.
        vector<char> linked((std::size_t)m_nblocks * (std::size_t)m_nblocks, 0);
        std::size_t nlinked = 0;
        for (std::size_t k = 0; k < m_p1.size(); ++k) {
            const int bu = m_p1[k] / 2, bv = m_p2[k] / 2;
            char &l = linked[(std::size_t)bu * m_nblocks + bv];
            if (bu != bv && !l) {
                l = 1;
                ++nlinked;
            }
        }
        // Once most block pairs are linked, every block pair gets a table
        // (all zeros for an unlinked one) so that the kernel takes u and v
        // from its loop counters. Reading them from m_pair_u / m_pair_v costs
        // more per lookup, and from about 80% linked down, the fewer lookups
        // no longer make up for it (measured at W = 20, 41 and 64).
        const std::size_t npairs_blocks = (std::size_t)m_nblocks * (std::size_t)(m_nblocks - 1) / 2;
        if (5 * nlinked >= 4 * npairs_blocks) {
            for (int u = 0; u < m_nblocks; ++u)
                for (int v = u + 1; v < m_nblocks; ++v)
                    linked[(std::size_t)u * m_nblocks + v] = 1;
            nlinked = npairs_blocks;
        }

        m_pair_u.clear();
        m_pair_v.clear();
        m_pair_table.assign(nlinked * 256, 0.0);
        for (int u = 0; u < m_nblocks; ++u) {
            const int pu0 = 2 * u;
            const int pu1 = (pu0 + 1 < m_W) ? pu0 + 1 : -1;
            for (int v = u + 1; v < m_nblocks; ++v) {
                if (!linked[(std::size_t)u * m_nblocks + v])
                    continue;
                const int pv0 = 2 * v;
                const int pv1 = (pv0 + 1 < m_W) ? pv0 + 1 : -1;
                const std::size_t off = m_pair_u.size() * 256;
                m_pair_u.push_back(2 * u);
                m_pair_v.push_back(2 * v);

                for (int cu = 0; cu < 16; ++cu) {
                    const int xu0 = cu >> 2, xu1 = cu & 3;
                    for (int cv = 0; cv < 16; ++cv) {
                        const int xv0 = cv >> 2, xv1 = cv & 3;

                        double s = pair_j(pu0, pv0, xu0, xv0);
                        if (pv1 >= 0)
                            s += pair_j(pu0, pv1, xu0, xv1);
                        if (pu1 >= 0)
                            s += pair_j(pu1, pv0, xu1, xv0);
                        if (pu1 >= 0 && pv1 >= 0)
                            s += pair_j(pu1, pv1, xu1, xv1);

                        m_pair_table[off + (std::size_t)cu * 16 + (std::size_t)cv] = s;
                    }
                }
            }
        }

        // The gate, in naive lookups: blocked when a * nlinked + 0.5 *
        // m_nblocks - 4 < W + npair, where nlinked counts the block pairs with
        // a table and a block pair's lookup costs a = 1 in the loop over every
        // block pair, 1.5 through the m_pair_u / m_pair_v list. A block pair's
        // table is 2 KB where the naive kernel reads 128 bytes per coupling,
        // so couplings scattered one per block pair cost the blocked kernel
        // more than the naive one, and wide models with them run naive.
        // Fitted to gseq.potts() times of both kernels on 245 models (W 6-256;
        // order 1, bands, random, side-by-side blocks, full pairwise and
        // scattered couplings, some at 81% of block pairs; one strand and
        // both): the kernel this picks was at most 12% slower than the other
        // one, 0.15% on average. The two kernels add the same terms in
        // different orders, so their scores can differ in the last bits; for
        // odd W, rc() links a different set of block pairs, so the two
        // strands of one model can run different kernels.
        const std::size_t a2 = (nlinked == npairs_blocks) ? 2 : 3; // 2a
        m_use_blocked = a2 * nlinked + (std::size_t)m_nblocks < 2 * ((std::size_t)m_W + m_p1.size()) + 8;
    }

    int m_W = 0;
    double m_intercept = 0.0;
    std::vector<double> m_e; // W*4,      [i*4 + base]
    std::vector<double> m_J; // npair*16, [k*16 + b*4 + a]
    std::vector<int> m_p1, m_p2;
    std::vector<std::vector<int>> m_pairs_at;
    // best_delta()'s tables, built by build_pair_index(): position j's entries
    // are [m_field_start[j], m_field_start[j + 1]).
    std::vector<std::size_t> m_field_start;
    std::vector<int> m_field_other;
    std::vector<double> m_field_tab;

    // Blocked-kernel tables, built by build_blocked_tables(). m_nblocks == 0
    // means "not built" (W too wide for MAX_BLOCKS) - blocked_available()
    // reports this. m_use_blocked is the gate's per-model choice; the tables
    // are built whenever W allows, so score_codes_blocked() stays callable for
    // the equivalence test whichever kernel production runs.
    int m_nblocks = 0;
    bool m_use_blocked = false;
    std::vector<double> m_block_table; // nblocks * 16, block b's at b * 16
    // The block pairs with a table - the linked ones, or every one from 80%
    // linked - in (u, v) order, and their 16 x 16 tables (pair k's at
    // k * 256). m_pair_u and m_pair_v hold each block's first position
    // (2u, 2v), which indexes the block codes directly: storing the block
    // index and doubling it per lookup took 14% longer through gseq.potts()
    // on a W = 41 band of 8 (1.23 vs 1.08 s).
    std::vector<int> m_pair_u, m_pair_v;
    std::vector<double> m_pair_table;
};

// Encode a target once per interval. codes[p] is 0..3 or -1 for any non-ACGT
// (BASE_ENCODE folds case). nbad[p] counts non-ACGT bases in [0, p), so anchor
// i of width W is scorable iff nbad[i + W] - nbad[i] == 0 - O(1) per anchor
// instead of re-walking W bases. nbad has size codes.size() + 1.
inline void potts_encode(const std::string &target,
                         std::vector<int8_t> &codes,
                         std::vector<int32_t> &nbad)
{
    const std::size_t n = target.size();
    codes.resize(n);
    nbad.resize(n + 1);
    nbad[0] = 0;
    for (std::size_t p = 0; p < n; ++p) {
        const int8_t c = DnaLookupTables::BASE_ENCODE[(unsigned char)target[p]];
        codes[p] = c;
        nbad[p + 1] = nbad[p] + (c < 0 ? 1 : 0);
    }
}

// block_codes[p] = 4 * codes[p] + codes[p + 1], the code of a block that
// starts at p, for p in [lo, hi) - what score_codes_blocked() reads. Built
// once per window for both strands rather than by the kernel per anchor, and
// only over the windows about to be scored: a vtrack that slides by one base
// encodes its whole target but scores one anchor, and a pass over the whole
// target made that 35-55% slower. block_codes is sized to codes. Past the
// end, or before a non-ACGT base, the second base counts as 0: a scorable
// window reads that only for the trailing size-1 block of an odd W, whose
// tables ignore the second base. A non-ACGT first base gives a negative code,
// never read, since a window with that base is not scored.
inline void potts_block_codes(const std::vector<int8_t> &codes,
                              std::vector<int8_t> &block_codes,
                              std::size_t lo, std::size_t hi)
{
    const std::size_t n = codes.size();
    block_codes.resize(n);
    if (hi > n)
        hi = n;
    // Every p but the last has a next base; the last gets 0.
    const std::size_t end = (hi == n && lo < hi) ? hi - 1 : hi;
    for (std::size_t p = lo; p < end; ++p) {
        const int8_t next = codes[p + 1];
        block_codes[p] = (int8_t)(4 * codes[p] + (next < 0 ? 0 : next));
    }
    if (end < hi)
        block_codes[end] = (int8_t)(4 * codes[end]);
}

#endif // POTTS_MODEL_H_
