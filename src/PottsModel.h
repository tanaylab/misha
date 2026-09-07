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
// This is motifmodel's Potts, and score_codes() must agree with its
// potts_score() kernel to float precision. Tables and accumulator are double:
// at W = 20 the J table is 24 KB and still L2-resident, and 211 float adds at
// magnitude ~20 would cost ~4e-4 of error for no measured speed gain.
//
// THE REVERSE STRAND IS A SECOND MODEL, not a reversed sequence. Scoring the
// reverse strand at an anchor means scoring the model on comp(x_{W-1-i}), and
// reindexing that gives a Potts of the same shape - see rc(). Building it once
// makes the reverse strand cost exactly what the forward strand costs, where
// reverse-complementing per anchor would dominate the loop.
class PottsModel {
public:
    PottsModel() = default;

    // e_rowmajor: W*4 doubles, row i = position i, columns A, C, G, T.
    // j_rowmajor: npair*16 doubles, row k = J_k with element [a][b] at b*4 + a
    //             (motifmodel's column-major 4x4 flattening).
    // p1, p2:     0-BASED position indices, p1[k] < p2[k].
    PottsModel(int W, double intercept,
               const double *e_rowmajor, const double *j_rowmajor,
               const std::vector<int> &p1, const std::vector<int> &p2)
        : m_W(W), m_intercept(intercept), m_p1(p1), m_p2(p2)
    {
        m_e.assign(e_rowmajor, e_rowmajor + (std::size_t)W * 4);
        // An order-1 model (no pairs - fit_motif(order = 1) produces one, and
        // test 1's npair_mode = "none" exercises it) passes p1.empty() and,
        // from PottsParams.h, a null j_rowmajor. Forming a [ptr, ptr+0) range
        // out of a null pointer is formally UB even though libstdc++ treats it
        // as a benign no-op, so this is guarded rather than relied upon.
        if (!p1.empty())
            m_J.assign(j_rowmajor, j_rowmajor + p1.size() * 16);
    }

    int width() const { return m_W; }
    std::size_t npair() const { return m_p1.size(); }
    double intercept() const { return m_intercept; }

    // Score a window of W base codes in 0..3. The caller guarantees the length
    // and that no code is negative; potts_encode()'s nbad prefix sum is the gate.
    inline double score_codes(const int8_t *c) const
    {
        double s = m_intercept;
        for (int i = 0; i < m_W; ++i)
            s += m_e[(std::size_t)i * 4 + c[i]];
        const std::size_t np = m_p1.size();
        for (std::size_t k = 0; k < np; ++k)
            s += m_J[k * 16 + (std::size_t)c[m_p2[k]] * 4 + (std::size_t)c[m_p1[k]]];
        return s;
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
        return o;
    }

private:
    int m_W = 0;
    double m_intercept = 0.0;
    std::vector<double> m_e; // W*4,      [i*4 + base]
    std::vector<double> m_J; // npair*16, [k*16 + b*4 + a]
    std::vector<int> m_p1, m_p2;
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
