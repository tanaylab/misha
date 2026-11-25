// PwmFastScan.h - skip PWM windows that provably cannot reach a score cut.
//
// Shared byte-identical by R misha and PyMisha. Standard library only.
//
// A thresholded PWM scan (count, hits) or a maximum (max, max.pos) only needs
// the windows that can reach a cut: the threshold, or the best score so far.
// BoundScanner proves cheaply that most windows cannot and hands the rest to
// the caller, which scores them with its own unchanged code. A skipped window
// is one whose score, as the caller would have computed it in float or
// double, is provably below the cut, so every result stays bit-identical to
// scoring every window (MISHA_DISABLE_PWM_SKIP=1 in the callers).
//
// Technique: lookahead scoring with precomputed tables over groups of letters,
// as published in
//   [1] MOODS: fast search for position weight matrix matches in DNA
//       sequences. Bioinformatics 25(23):3181-3182 (2009).
//   [2] Finding significant matches of position weight matrices in linear
//       time. IEEE/ACM Trans Comput Biol Bioinform 8(1):69-79 (2011).
//
//   Codes.  The scanned letters are read as {A, C, G, T, other} (either case)
//           and every position p gets the base-5 code of the kQ letters
//           starting at p, first letter most significant.
//   Blocks. Each strand's weights, in window order, are cut into blocks of kQ
//           consecutive columns starting at column 0; the last block may be
//           shorter. A block has a table of its summed weight for every code
//           (a shorter block ignores the trailing letters of the code), and
//           its maximum entry.
//   Order.  A strand evaluates its blocks in decreasing order of range (the
//           largest minus the smallest finite entry), the block that tells
//           windows apart best first.
//   Bound.  After k blocks a window's partial sum plus the maxima of the
//           blocks not yet added bounds its strand score from above. The
//           window is dropped from that strand as soon as this falls below the
//           cut (less a rounding margin). It is handed to the caller if it
//           survives every block on a strand the caller reads.
//
// Soundness. Let a window's m weights on a strand be w_i (finite; a window with
// a -inf weight scores -inf in every caller and is below any cut the scanner
// tests, since cuts of -inf or NaN skip nothing), S = sum w_i and A = sum |w_i|.
//
//   1. Callers add the m float weights one by one, in float (calc_like) or in
//      double (gseq), in either column order. Recursive summation of m terms
//      with unit roundoff u returns at most S + gamma * A, with
//      gamma = m u / (1 - m u); float's u = 2^-24 covers double too. No
//      partial sum overflows: init() refuses weights whose column maxima of
//      |w| add up past 1e30.
//   2. Every table entry is built from the inflated weights w + 2 gamma |w|
//      in double and rounded up to float, so it is at least the block's sum of
//      w + gamma |w| (the double rounding is below 2^-50 |w| per term, far
//      inside the extra gamma |w| >= 2^-24 |w|). Summed over the blocks, the
//      entries a window picks are >= S + gamma A >= the caller's score. Block
//      maxima are entries too, so the bound after k blocks is >= that score.
//   3. The scan adds at most nb (blocks per strand) table entries in double,
//      and the cut is lowered by the maxima of the blocks still to come, also
//      in double. Each of those roundings is at most 2^-53 times a value no
//      larger than Mag + |cut|, where Mag is the larger over the strands of
//      sum over blocks of max |finite entry|. The window is dropped only when
//      the computed bound is below cut - (nb + 4) 2^-48 (Mag + |cut|), which
//      exceeds those (2 nb + 4) roundings sixteen times over. So a dropped
//      window's exact bound, and with it the caller's strand score, is below
//      the cut.
//   4. strand_cut() turns a threshold on the caller's final score into a cut
//      on each strand's score; see the derivation there.
//
// The scanner is never wrong in the other direction: a window it hands over is
// simply scored. Every margin errs that way.

#ifndef PWM_FAST_SCAN_H_
#define PWM_FAST_SCAN_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace pwm_fast_scan {

// Below this many windows callers score every window. Measured on random
// motifs of width 8-20: a running maximum breaks even at about 6 windows per
// scan and a threshold count at under 4; 16 leaves room for motifs that carry
// little information, whose windows are harder to rule out.
static const size_t kMinWindows = 16;

namespace detail {

static const int kQ = 4;                  // letters per code
static const int kAlpha = 5;              // A, C, G, T, other
static const int kOther = 4;
static const int kCodes = 625;            // kAlpha^kQ
static const unsigned kTopPow = 125;      // kAlpha^(kQ - 1)
static const size_t kChunk = 4096;        // windows encoded at a time
static const size_t kBatch = 64;          // windows per batched first-block test

// 0..3 for A, C, G, T in either case, kOther for every other byte
inline const unsigned char *letter_codes()
{
    struct Table {
        unsigned char v[256];
        Table()
        {
            for (int i = 0; i < 256; ++i)
                v[i] = (unsigned char)kOther;
            v['A'] = v['a'] = 0;
            v['C'] = v['c'] = 1;
            v['G'] = v['g'] = 2;
            v['T'] = v['t'] = 3;
        }
    };
    static const Table table;
    return table.v;
}

// The smallest float >= e, for e in [-inf, FLT_MAX]. When the conversion
// rounded down, step one float up: a positive float's next value up has the
// next bit pattern, a negative one's the previous.
inline float float_at_least(double e)
{
    float f = (float)e;
    if ((double)f < e) {
        uint32_t bits;
        std::memcpy(&bits, &f, sizeof bits);
        if (f == 0.0f)
            bits = 1;   // the smallest positive float
        else if (f > 0.0f)
            ++bits;
        else
            --bits;
        std::memcpy(&f, &bits, sizeof f);
    }
    return f;
}

// Writes the codes of positions p0 .. p0 + np - 1 of s, which holds nlet
// readable letters. Letters past the end read as A: only a shorter last block
// reaches past a window, and its table ignores those letters.
inline void encode(const char *s, size_t nlet, size_t p0, size_t np, uint16_t *out)
{
    const unsigned char *lc = letter_codes();
    unsigned c = 0;
    size_t a = p0;    // next letter to shift in
    for (int i = 0; i < kQ - 1; ++i, ++a)
        c = c * kAlpha + (a < nlet ? lc[(unsigned char)s[a]] : 0u);
    size_t p = 0;
    for (; p < np && a < nlet; ++p, ++a) {
        c = (c % kTopPow) * kAlpha + lc[(unsigned char)s[a]];
        out[p] = (uint16_t)c;
    }
    for (; p < np; ++p) {
        c = (c % kTopPow) * kAlpha;
        out[p] = (uint16_t)c;
    }
}

} // namespace detail

// A cut on each strand's raw score (the sum of its column weights) such that a
// window whose used strands all score below it has a final score below
// `thresh`. Also used for maxima, as strand_cut(best_so_far, ...).
//
//   n_lse      2 when the final score is the log-sum-exp of the two strands,
//              1 when it is one strand or the larger of the two.
//   extra_max  the largest value added to a strand score afterwards (the
//              largest spatial log factor, which may be negative), or 0.
//
// Derivation. The caller's final score is built from the strand scores (all
// below the cut c) by at most three roundings: the log-sum-exp
// hi + log1p(exp(lo - hi)), whose correction lies in [0, ln n_lse] up to
// library error (well under 1e-6), and the addition of a spatial term <= e.
// Rounding is monotone, so the final score is at most what these operations
// give for c + ln n_lse + e, which is that value plus three roundings of at most
// 2^-24 relative to numbers no larger than |thresh| + |e| + 1, plus the library
// error. The slack 1e-5 (1 + |thresh| + |e|) covers this many times over.
inline double strand_cut(double thresh, int n_lse, double extra_max)
{
    const double inf = std::numeric_limits<double>::infinity();
    if (std::isnan(thresh) || std::isnan(extra_max))
        return std::numeric_limits<double>::quiet_NaN();   // skip nothing
    if (thresh == -inf || extra_max == inf)
        return -inf;                                        // anything may pass
    if (thresh == inf || extra_max == -inf)
        return inf;                                         // no finite score passes
    const double lse = n_lse == 2 ? 0.69314718055994531 : n_lse > 2 ? std::log((double)n_lse) : 0.0;
    const double slack = 1e-5 * (1.0 + std::fabs(thresh) + std::fabs(extra_max));
    return thresh - extra_max - lse - slack;
}

class BoundScanner {
public:
    BoundScanner() : m_(0), nb_(0), usable_(false), mag_(0.0) {}

    // logp: m rows of 4 weights (A, C, G, T), the float values the caller adds
    // on the forward strand. The reverse strand of the window at the same start
    // reads column m-1-k, complemented, at its k-th letter. other: the weight of
    // any other letter, per column (nullptr: -inf, such windows never score);
    // the reverse strand charges column m-1-k's at its k-th letter.
    void init(const float *logp, int m, const float *other = nullptr);

    // False when m <= 0, a weight is NaN or +inf, or the weights are too large
    // for the rounding argument (see the top of the file). Callers then score
    // every window.
    bool usable() const { return usable_; }

    int width() const { return m_; }

    // Calls fn(x), in ascending order, for every window start x in [0, nwin)
    // that may score at least `cut` on a used strand (use_fwd: the forward
    // weights, use_rc: the reverse ones); a skipped window scores below it on
    // every used strand. s holds nwin + m - 1 readable letters. `cut` is read
    // again after every call to fn, so a running maximum may raise it. A cut
    // of -inf or NaN skips nothing; +inf skips everything. Without a usable
    // scanner, or with no strand used, fn is called for every x.
    //
    // const and allocation-local: one scanner may serve several threads.
    template <class Fn>
    void scan(const char *s, size_t nwin, bool use_fwd, bool use_rc, const double &cut, Fn &&fn) const;

private:
    struct Strand {
        std::vector<float> tab;    // nb tables of kCodes entries, in lookahead order
        std::vector<int> off;      // first column of each block, in that order
        std::vector<double> rest;  // rest[k]: sum of the maxima of the blocks after the k-th
        double total;              // sum of all block maxima
        double mag;                // sum over blocks of max |finite entry|
        Strand() : total(0.0), mag(0.0) {}
    };

    void build(Strand &st, const std::vector<double> &w);
    // Whether the window whose codes start at cx passes every block of st,
    // with lim[k] the partial sum it needs after the k-th block
    bool passes(const Strand &st, const uint16_t *cx, const double *lim) const;

    int m_;
    int nb_;
    bool usable_;
    double mag_;
    Strand st_[2];   // forward, reverse
};

inline void BoundScanner::init(const float *logp, int m, const float *other)
{
    using namespace detail;
    *this = BoundScanner();
    m_ = m;
    if (m <= 0 || !logp)
        return;
    const double u = std::ldexp(1.0, -24);
    if (m * u >= 0.01)
        return;
    const double inf = std::numeric_limits<double>::infinity();

    // weight of letter b (0..3, kOther) in column k
    auto weight = [&](int k, int b) -> double {
        if (b < 4)
            return (double)logp[(size_t)k * 4 + b];
        return other ? (double)other[k] : -inf;
    };
    double abs_sum = 0.0;
    for (int k = 0; k < m; ++k) {
        double col = 0.0;
        for (int b = 0; b < kAlpha; ++b) {
            const double w = weight(k, b);
            if (std::isnan(w) || w == inf)
                return;
            if (w != -inf)
                col = std::max(col, std::fabs(w));
        }
        abs_sum += col;
    }
    if (!(abs_sum <= 1e30))
        return;

    // point 2 at the top: inflate every finite weight by 2 gamma |w|
    const double gamma = m * u / (1.0 - m * u);
    auto inflate = [&](double w) -> double { return w == -inf ? w : w + 2.0 * gamma * std::fabs(w); };
    std::vector<double> fwd((size_t)m * kAlpha), rev((size_t)m * kAlpha);
    for (int k = 0; k < m; ++k) {
        for (int b = 0; b < kAlpha; ++b) {
            const int rb = b < 4 ? 3 - b : b;   // complement; other stays other
            fwd[(size_t)k * kAlpha + b] = inflate(weight(k, b));
            rev[(size_t)k * kAlpha + b] = inflate(weight(m - 1 - k, rb));
        }
    }
    nb_ = (m + kQ - 1) / kQ;
    build(st_[0], fwd);
    build(st_[1], rev);
    mag_ = std::max(st_[0].mag, st_[1].mag);
    usable_ = true;
}

inline void BoundScanner::build(Strand &st, const std::vector<double> &w)
{
    using namespace detail;
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<float> tab((size_t)nb_ * kCodes);
    std::vector<double> range(nb_), maxv(nb_), absmax(nb_), sums(kCodes);
    for (int j = 0; j < nb_; ++j) {
        const int off = j * kQ;
        const int len = std::min(kQ, m_ - off);
        // sums[code] = the block's weights summed over the code's letters,
        // first letter first; built one letter at a time, in place (the code
        // with prefix c and next letter d is c * kAlpha + d >= c)
        sums[0] = 0.0;
        size_t cnt = 1;
        for (int i = 0; i < kQ; ++i, cnt *= kAlpha) {
            for (size_t c = cnt; c-- > 0;) {
                const double base = sums[c];
                for (int d = kAlpha - 1; d >= 0; --d)
                    sums[c * kAlpha + d] = i < len ? base + w[(size_t)(off + i) * kAlpha + d] : base;
            }
        }
        float hi = -std::numeric_limits<float>::infinity();
        float lo = std::numeric_limits<float>::infinity();
        double amax = 0.0;
        for (int c = 0; c < kCodes; ++c) {
            const float t = float_at_least(sums[c]);
            tab[(size_t)j * kCodes + c] = t;
            hi = std::max(hi, t);
            if (t != -std::numeric_limits<float>::infinity()) {
                lo = std::min(lo, t);
                amax = std::max(amax, std::fabs((double)t));
            }
        }
        maxv[j] = hi;
        absmax[j] = amax;
        // a block with no finite entry rules out every window: test it first
        range[j] = lo <= hi ? (double)hi - (double)lo : inf;
    }

    std::vector<int> order(nb_);
    for (int j = 0; j < nb_; ++j)
        order[j] = j;
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return range[a] > range[b]; });

    st.tab.assign((size_t)nb_ * kCodes, 0.0f);
    st.off.assign(nb_, 0);
    st.rest.assign(nb_, 0.0);
    st.mag = 0.0;
    for (int k = 0; k < nb_; ++k) {
        const int j = order[k];
        std::copy(tab.begin() + (size_t)j * kCodes, tab.begin() + (size_t)(j + 1) * kCodes,
                  st.tab.begin() + (size_t)k * kCodes);
        st.off[k] = j * kQ;
        st.mag += absmax[j];
    }
    double rest = 0.0;
    for (int k = nb_ - 1; k >= 0; --k) {
        st.rest[k] = rest;
        rest += maxv[order[k]];
    }
    st.total = rest;
}

inline bool BoundScanner::passes(const Strand &st, const uint16_t *cx, const double *lim) const
{
    const float *tab = st.tab.data();
    const int *off = st.off.data();
    double p = 0.0;
    for (int k = 0; k < nb_; ++k, tab += detail::kCodes) {
        p += tab[cx[off[k]]];
        if (p < lim[k])
            return false;
    }
    return true;
}

template <class Fn>
void BoundScanner::scan(const char *s, size_t nwin, bool use_fwd, bool use_rc, const double &cut, Fn &&fn) const
{
    using namespace detail;
    if (nwin == 0)
        return;
    if (!usable_ || (!use_fwd && !use_rc)) {
        for (size_t x = 0; x < nwin; ++x)
            fn(x);
        return;
    }
    const double inf = std::numeric_limits<double>::infinity();
    const bool use[2] = {use_fwd, use_rc};
    const size_t nlet = nwin + (size_t)m_ - 1;
    const size_t span = (size_t)(nb_ - 1) * kQ;   // offset of the last block
    const size_t chunk_max = std::min(nwin, kChunk);
    std::vector<uint16_t> codes(chunk_max + span);
    std::vector<double> lim(2 * (size_t)nb_);

    // lim[t * nb + k]: the partial sum strand t needs after its k-th block,
    // (cut - margin) - rest[k] (point 3 at the top). on[t]: strand t is read
    // and can reach the cut at all. state: every window is scored (a cut of
    // -inf or NaN), windows are tested, or none can pass.
    enum { kAll, kTest, kNone };
    int state = kAll;
    bool built = false, on[2] = {false, false};
    double built_for = 0.0;
    auto stale = [&]() { return !built || !(cut == built_for); };
    auto rebuild = [&]() {
        const double c = cut;
        built = true;
        built_for = c;
        if (!(c > -inf)) {
            state = kAll;
            return;
        }
        if (c == inf) {
            state = kNone;
            return;
        }
        const double ce = c - (nb_ + 4) * std::ldexp(mag_ + std::fabs(c), -48);
        for (int t = 0; t < 2; ++t) {
            on[t] = use[t] && !(0.0 < ce - st_[t].total);
            for (int k = 0; k < nb_; ++k)
                lim[(size_t)t * nb_ + k] = ce - st_[t].rest[k];
        }
        state = (on[0] || on[1]) ? kTest : kNone;
    };

    const float *first_tab[2] = {st_[0].tab.data(), st_[1].tab.data()};
    const int first_off[2] = {st_[0].off[0], st_[1].off[0]};
    uint16_t keep[kBatch];

    for (size_t x0 = 0; x0 < nwin;) {
        const size_t n = std::min(chunk_max, nwin - x0);
        encode(s, nlet, x0, n + span, codes.data());
        size_t i = 0;
        while (i < n) {
            if (stale())
                rebuild();
            if (state == kNone)
                return;   // no window can pass, and only fn could change the cut
            if (state == kAll) {
                fn(x0 + i);
                ++i;
                continue;
            }
            // The first block of both strands over a batch of windows, without
            // branches (the first step of passes(): 0 + entry < lim[0]); only
            // the windows that survive it are tested in full.
            const size_t bn = std::min(kBatch, n - i);
            const double lim0 = on[0] ? lim[0] : inf, lim1 = on[1] ? lim[nb_] : inf;
            const uint16_t *cx = codes.data() + i;
            size_t nk = 0;
            for (size_t j = 0; j < bn; ++j) {
                const bool pass0 = !((double)first_tab[0][cx[j + first_off[0]]] < lim0);
                const bool pass1 = !((double)first_tab[1][cx[j + first_off[1]]] < lim1);
                keep[nk] = (uint16_t)j;
                nk += (size_t)(pass0 | pass1);
            }
            const double filtered_for = built_for;
            size_t next = i + bn;
            for (size_t k = 0; k < nk; ++k) {
                const size_t j = i + keep[k];
                const uint16_t *cj = codes.data() + j;
                if (!(on[0] && passes(st_[0], cj, lim.data())) &&
                    !(on[1] && passes(st_[1], cj, lim.data() + nb_)))
                    continue;
                fn(x0 + j);
                if (stale()) {
                    rebuild();
                    // A raised cut keeps the batch's survivors a superset of
                    // the windows that can pass; any other change restarts
                    // after this window.
                    if (state != kTest || !(built_for >= filtered_for)) {
                        next = j + 1;
                        break;
                    }
                }
            }
            i = next;
        }
        x0 += n;
    }
}

// One window that passes seq_hits()' threshold
struct SeqHit {
    int64_t start;   // 0-based start of the window
    int strand;      // 1 or -1: the strand read alone, or under bidirect the one
                     // that scores higher on its own (+1 on a tie)
    double score;    // the thresholded score (the strands' log-sum-exp under bidirect)
};

// Appends to `out`, in ascending order, the windows starting in [lo, hi]
// (0-based) of seq (len letters) that score at least thresh, scored exactly as
// misha's gseq.pwm(mode = "count") scores a sequence without gaps, spatial
// weights or neutral characters: per strand the double sum of the float
// weights in logp (forward: columns 0..m-1 over the window's letters; reverse:
// the k-th letter complemented against column m-1-k), and a window with any
// other letter does not score. Under bidirect the strands are combined by the
// same log-sum-exp; otherwise strand -1 reads the reverse strand and any other
// value the forward one. bound: nullptr, or a scanner initialised on logp.
inline void seq_hits(BoundScanner *bound, const float *logp, int m, const char *seq, size_t len,
                     int64_t lo, int64_t hi, bool bidirect, int strand, double thresh,
                     std::vector<SeqHit> &out)
{
    if (m <= 0 || lo < 0 || hi < lo || (uint64_t)hi + (uint64_t)m > (uint64_t)len)
        return;
    const double neg_inf = -std::numeric_limits<double>::infinity();
    const unsigned char *lc = detail::letter_codes();
    const bool use_fwd = bidirect || strand != -1;
    const bool use_rc = bidirect || strand == -1;

    auto score_at = [&](int64_t x) {
        const char *win = seq + x;
        double f = neg_inf, r = neg_inf;
        if (use_fwd) {
            double sum = 0.0;
            for (int k = 0; k < m; ++k) {
                const unsigned b = lc[(unsigned char)win[k]];
                if (b > 3) {
                    sum = neg_inf;
                    break;
                }
                sum += logp[(size_t)k * 4 + b];
            }
            f = sum;
        }
        if (use_rc) {
            double sum = 0.0;
            for (int k = 0; k < m; ++k) {
                const unsigned b = lc[(unsigned char)win[k]];
                if (b > 3) {
                    sum = neg_inf;
                    break;
                }
                sum += logp[(size_t)(m - 1 - k) * 4 + (3 - b)];
            }
            r = sum;
        }
        double v;
        int st;
        if (bidirect) {
            // gseq.pwm's log-sum-exp of the two strands, operation for operation
            if (f == neg_inf)
                v = r;
            else if (r == neg_inf)
                v = f;
            else if (f > r)
                v = f + std::log1p(std::exp(r - f));
            else
                v = r + std::log1p(std::exp(f - r));
            st = r > f ? -1 : 1;
        } else {
            v = strand == -1 ? r : f;
            st = strand == -1 ? -1 : 1;
        }
        if (v >= thresh) {
            SeqHit h;
            h.start = x;
            h.strand = st;
            h.score = v;
            out.push_back(h);
        }
    };

    const size_t n = (size_t)(hi - lo) + 1;
    if (bound && bound->usable() && bound->width() == m && n >= kMinWindows) {
        const double cut = strand_cut(thresh, bidirect ? 2 : 1, 0.0);
        bound->scan(seq + lo, n, use_fwd, use_rc, cut, [&](size_t i) { score_at(lo + (int64_t)i); });
    } else {
        for (int64_t x = lo; x <= hi; ++x)
            score_at(x);
    }
}

} // namespace pwm_fast_scan

#endif // PWM_FAST_SCAN_H_
