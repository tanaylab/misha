#ifndef POTTS_MODEL_H_
#define POTTS_MODEL_H_

#include <cstddef>
#include <cstdint>
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
// score_codes() has two implementations, chosen once per model in
// build_blocked_tables() rather than per anchor:
//
//   score_codes_naive()   - one table lookup per single-site term and one per
//                            pair term: 1 + W + npair lookups.
//   score_codes_blocked() - positions packed two at a time into a 0..15
//                            "block code"; one table per block folds its own
//                            singles (and the pair within the block, if any)
//                            into one lookup, and one table per BLOCK PAIR
//                            folds the up to 4 cross-block couplings into one
//                            lookup: ceil(W/2) + C(ceil(W/2), 2) lookups. At
//                            W = 20, full pairwise, that is 10 + 45 = 55
//                            instead of 211.
//
// score_codes_blocked()'s lookup count does NOT depend on how many pairs the
// model actually has (it is fixed by W alone), while score_codes_naive()'s
// does. So the blocked kernel only wins when the model is dense enough that
// 1 + W + npair exceeds ceil(W/2) + C(ceil(W/2), 2) - measured through
// gseq.potts() on an idle host at W = 20: ~3.4x faster full pairwise (190
// pairs), a real but smaller win at 64 pairs, and 4-6x SLOWER at 0 pairs
// (an order-1 model, i.e. one with no couplings at all), where it would spend 55
// lookups computing what 21 could. score_codes() therefore picks whichever
// kernel has the smaller lookup count for THIS model (m_use_blocked, set once
// in build_blocked_tables()), not the blocked kernel unconditionally.
//
// build_blocked_tables() generalizes to odd W (a trailing size-1 block) and
// to sparse or absent pairs (a missing coupling contributes 0 to its slot),
// so score_codes_blocked() agrees with score_codes_naive() for ANY
// PottsModel - not only the dense, even-W case it is designed to win on. For
// a W wide enough that ceil(W/2) would overflow the fixed block-code scratch
// array (MAX_BLOCKS = 128, i.e. W > 256), build_blocked_tables() leaves the
// blocked tables unbuilt (blocked_available() is false) and score_codes()
// falls back to the naive kernel.
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
    }

    int width() const { return m_W; }
    std::size_t npair() const { return m_p1.size(); }
    double intercept() const { return m_intercept; }

    // Whether build_blocked_tables() actually built usable tables (false only
    // for W > 256 - see the class comment). score_codes_blocked() requires
    // this; score_codes() checks it internally, any other caller must check
    // it itself before calling score_codes_blocked() directly.
    bool blocked_available() const { return m_nblocks > 0; }

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

    // The dinucleotide-blocked kernel - see the class comment. The caller
    // must guarantee blocked_available() (score_codes() checks it; a direct
    // caller, such as the equivalence test's C_potts_score_codes_cmp, must
    // check it itself).
    inline double score_codes_blocked(const int8_t *c) const
    {
        int8_t code[MAX_BLOCKS];
        for (int b = 0; b < m_nblocks; ++b) {
            const int p0 = 2 * b;
            code[b] = (m_block_size[b] == 2) ? (int8_t)(4 * c[p0] + c[p0 + 1]) : c[p0];
        }

        double s = m_intercept;
        for (int b = 0; b < m_nblocks; ++b)
            s += m_block_table[m_block_offset[b] + (std::size_t)code[b]];

        std::size_t pp = 0;
        for (int u = 0; u < m_nblocks; ++u)
            for (int v = u + 1; v < m_nblocks; ++v, ++pp)
                s += m_pair_table[m_pair_offset[pp] +
                                   (std::size_t)code[u] * (std::size_t)m_card[v] + (std::size_t)code[v]];
        return s;
    }

    // Production entry point - GseqPotts.cpp and PottsParams' consumers call
    // this and never the two kernels above directly. Picks whichever kernel
    // has the smaller lookup count for THIS model - see the class comment.
    inline double score_codes(const int8_t *c) const
    {
        return m_use_blocked ? score_codes_blocked(c) : score_codes_naive(c);
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
        return o;
    }

private:
    // Supports W up to 256. build_blocked_tables() disables the blocked
    // kernel rather than overflow the score_codes_blocked() scratch array if
    // ceil(W/2) exceeds this - no PottsModel using it today gets remotely
    // close (test widths top out at 21; a realistic motif width tops out
    // well under 100).
    static constexpr int MAX_BLOCKS = 128;

    void build_blocked_tables()
    {
        m_nblocks = (m_W + 1) / 2; // ceil(W / 2)
        if (m_nblocks == 0 || m_nblocks > MAX_BLOCKS) {
            m_nblocks = 0; // blocked_available() is false; score_codes() uses naive
            return;
        }

        m_block_size.assign(m_nblocks, 2);
        m_card.assign(m_nblocks, 16);
        if (m_W % 2 == 1) {
            m_block_size.back() = 1;
            m_card.back() = 4;
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

        // Per-block tables: the block's own single-site terms, plus the
        // within-block pair if the model has one.
        m_block_offset.assign((std::size_t)m_nblocks + 1, 0);
        for (int b = 0; b < m_nblocks; ++b)
            m_block_offset[b + 1] = m_block_offset[b] + (std::size_t)m_card[b];
        m_block_table.assign(m_block_offset[m_nblocks], 0.0);

        for (int b = 0; b < m_nblocks; ++b) {
            const int p0 = 2 * b;
            const bool paired = m_block_size[b] == 2;
            const int p1_ = paired ? p0 + 1 : -1;
            for (int code = 0; code < m_card[b]; ++code) {
                const int x0 = paired ? (code >> 2) : code;
                const int x1 = paired ? (code & 3) : 0;
                double v = m_e[(std::size_t)p0 * 4 + x0];
                if (paired) {
                    v += m_e[(std::size_t)p1_ * 4 + x1];
                    v += pair_j(p0, p1_, x0, x1);
                }
                m_block_table[m_block_offset[b] + (std::size_t)code] = v;
            }
        }

        // Per-block-pair tables: up to 4 cross terms between the (1 or 2)
        // positions of block u and the (1 or 2) positions of block v.
        const int npairs_blocks = m_nblocks * (m_nblocks - 1) / 2;
        m_pair_offset.assign((std::size_t)npairs_blocks + 1, 0);
        {
            int pp = 0;
            for (int u = 0; u < m_nblocks; ++u)
                for (int v = u + 1; v < m_nblocks; ++v, ++pp)
                    m_pair_offset[pp + 1] =
                        m_pair_offset[pp] + (std::size_t)m_card[u] * (std::size_t)m_card[v];
        }
        m_pair_table.assign(m_pair_offset[npairs_blocks], 0.0);

        int pp = 0;
        for (int u = 0; u < m_nblocks; ++u) {
            const int pu0 = 2 * u;
            const int pu1 = (m_block_size[u] == 2) ? pu0 + 1 : -1;
            for (int v = u + 1; v < m_nblocks; ++v, ++pp) {
                const int pv0 = 2 * v;
                const int pv1 = (m_block_size[v] == 2) ? pv0 + 1 : -1;
                for (int cu = 0; cu < m_card[u]; ++cu) {
                    const int xu0 = (pu1 >= 0) ? (cu >> 2) : cu;
                    const int xu1 = (pu1 >= 0) ? (cu & 3) : 0;
                    for (int cv = 0; cv < m_card[v]; ++cv) {
                        const int xv0 = (pv1 >= 0) ? (cv >> 2) : cv;
                        const int xv1 = (pv1 >= 0) ? (cv & 3) : 0;

                        double s = pair_j(pu0, pv0, xu0, xv0);
                        if (pv1 >= 0)
                            s += pair_j(pu0, pv1, xu0, xv1);
                        if (pu1 >= 0)
                            s += pair_j(pu1, pv0, xu1, xv0);
                        if (pu1 >= 0 && pv1 >= 0)
                            s += pair_j(pu1, pv1, xu1, xv1);

                        m_pair_table[m_pair_offset[pp] +
                                      (std::size_t)cu * (std::size_t)m_card[v] + (std::size_t)cv] = s;
                    }
                }
            }
        }

        // The gate: use the blocked kernel only when it actually has fewer
        // lookups for THIS model. Its lookup count (m_nblocks + npairs_blocks)
        // is fixed by W alone; the naive kernel's (1 + W + npair) shrinks with
        // the model's pair count, so a sparse or order-1 model can make naive
        // the faster choice even though blocked wins the dense case this
        // kernel was designed for. See the class comment for measured numbers.
        const std::size_t blocked_lookups = (std::size_t)m_nblocks + (std::size_t)npairs_blocks;
        const std::size_t naive_lookups = 1 + (std::size_t)m_W + m_p1.size();
        m_use_blocked = blocked_lookups < naive_lookups;
    }

    int m_W = 0;
    double m_intercept = 0.0;
    std::vector<double> m_e; // W*4,      [i*4 + base]
    std::vector<double> m_J; // npair*16, [k*16 + b*4 + a]
    std::vector<int> m_p1, m_p2;

    // Blocked-kernel tables, built by build_blocked_tables(). m_nblocks == 0
    // means "not built" (W too wide for MAX_BLOCKS) - blocked_available()
    // reports this. m_use_blocked is the separate, per-model dispatch
    // decision score_codes() acts on (see build_blocked_tables()'s gate) -
    // the tables are always built when W allows, even when m_use_blocked is
    // false, so score_codes_blocked() stays callable (and correct) for the
    // equivalence test regardless of which kernel production picks.
    int m_nblocks = 0;
    bool m_use_blocked = false;
    std::vector<int> m_block_size;      // per block, 1 or 2 positions
    std::vector<int> m_card;            // per block, 4^size
    std::vector<double> m_block_table;
    std::vector<std::size_t> m_block_offset; // nblocks+1, prefix sums into m_block_table
    std::vector<double> m_pair_table;
    std::vector<std::size_t> m_pair_offset;  // npairs_blocks+1, prefix sums into m_pair_table
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

#endif // POTTS_MODEL_H_
