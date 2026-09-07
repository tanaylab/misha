#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#include "PottsScorer.h"
#include "util.h" // log_sum_log

PottsScorer::PottsScorer(const PottsModel &model, const std::string &genome_root,
                         bool extend, ScoringMode mode, bool bidirect, char strand,
                         double score_thresh)
    : GenomeSeqScorer(genome_root, extend, strand), m_model(model), m_rc(model.rc()),
      m_mode(mode), m_bidirect(bidirect), m_score_thresh(score_thresh)
{
}

PottsScorer::PottsScorer(const PottsModel &model, GenomeSeqFetch *shared_seqfetch,
                         bool extend, ScoringMode mode, bool bidirect, char strand,
                         double score_thresh)
    : GenomeSeqScorer(shared_seqfetch, extend, strand), m_model(model), m_rc(model.rc()),
      m_mode(mode), m_bidirect(bidirect), m_score_thresh(score_thresh)
{
}

void PottsScorer::invalidate_cache()
{
    // No cache yet - the sliding-window one lands separately, deliberately
    // after the uncached path is tested, so a cache bug cannot hide behind an
    // uncached-path bug. Declared and called from both
    // TrackExpressionVars::start_chrom() overloads from the start, rather than
    // having the call sites added later and one of the two forgotten.
}

double PottsScorer::anchor_value(const int8_t *c, bool union_max, int &dir) const
{
    // When m_strand == -1 the fetched target is already reverse-complemented,
    // so the ORIGINAL forward strand is what the twin reads off it, and vice
    // versa. Same inversion as PWMScorer's score_{forward,reverse}_original().
    const bool check_fwd = m_bidirect || m_strand != -1;
    const bool check_rev = m_bidirect || m_strand == -1;

    double f = -std::numeric_limits<double>::infinity();
    double r = -std::numeric_limits<double>::infinity();
    if (check_fwd)
        f = (m_strand == -1) ? m_rc.score_codes(c) : m_model.score_codes(c);
    if (check_rev)
        r = (m_strand == -1) ? m_model.score_codes(c) : m_rc.score_codes(c);

    dir = 1;
    if (check_fwd && check_rev) {
        if (union_max) {
            if (r > f) {
                dir = -1;
                return r;
            }
            return f;
        }
        double u = f;
        log_sum_log(u, r);
        return u;
    }
    if (check_fwd)
        return f;
    dir = -1;
    return r;
}

float PottsScorer::compute_position_result(size_t index, size_t target_length,
                                           size_t motif_length, int direction) const
{
    if (target_length < motif_length)
        return std::numeric_limits<float>::quiet_NaN();

    float pos_result = float(index) + 1.0f; // 1-based

    if (m_strand == -1) {
        // The target is reverse-complemented: a window at target index `index`
        // covers forward-strand 0-based [target_length - index - motif_length, ...),
        // so the 1-based forward offset is that plus one. Signed arithmetic
        // because an unsigned underflow here would read as a huge position.
        pos_result = float(std::ptrdiff_t(target_length) - std::ptrdiff_t(index) -
                           std::ptrdiff_t(motif_length)) + 1.0f;
    }

    if (m_bidirect)
        pos_result = pos_result * direction;

    return pos_result;
}

float PottsScorer::score_interval(const GInterval &interval, const GenomeChromKey &chromkey)
{
    m_last_max_score = -std::numeric_limits<double>::infinity();

    const int64_t motif_length = m_model.width();
    GInterval expanded = calculate_expanded_interval(interval, chromkey, motif_length);
    expanded.strand = m_strand;

    if (!m_extend && (expanded.end - expanded.start) < motif_length)
        return std::numeric_limits<float>::quiet_NaN();

    try {
        std::vector<char> seq;
        m_seqfetch_ptr->read_interval(expanded, chromkey, seq);
        const std::string target(seq.begin(), seq.end());

        const size_t tlen = target.size();
        const size_t motif_len = (size_t)motif_length;
        if (tlen < motif_len)
            return (m_mode == MOTIF_COUNT) ? 0.0f
                                           : std::numeric_limits<float>::quiet_NaN();

        // Clamp to anchors whose STARTS fall inside the iterator interval.
        // Transliterated from PWMScorer::score_interval(); extension is
        // END-only, so extra_left is 0 whenever extend is true.
        size_t i_min = 0;
        size_t i_max = tlen - motif_len;
        const int64_t interval_len = interval.end - interval.start;
        if (interval_len > 0) {
            const int64_t extra_left = std::max<int64_t>(0, interval.start - expanded.start);
            const int64_t max_valid = (int64_t)(tlen - motif_len);
            auto clamp_index = [&](int64_t idx) -> size_t {
                if (idx < 0)
                    return 0;
                if (idx > max_valid)
                    return (size_t)std::max<int64_t>(0, max_valid);
                return (size_t)idx;
            };
            const size_t lo = clamp_index(extra_left);
            const size_t hi = clamp_index(extra_left + interval_len - 1);
            if (lo <= hi) {
                i_min = std::max(i_min, lo);
                i_max = std::min(i_max, hi);
            } else {
                // Unreachable as written: clamp_index() is monotone, so
                // interval_len > 0 forces lo <= hi. Kept because the sibling
                // has it.
                i_min = i_max = lo;
            }
        }
        // Also unreachable: lo <= max_valid, so i_min <= i_max above. If it
        // ever went live it would reset the scan to index 0 and score anchors
        // OUTSIDE the interval, which is why it is worth naming rather than
        // trusting - the sibling's identical line is the only reason it is
        // still here.
        if (i_min > i_max)
            i_min = 0;

        potts_encode(target, m_codes, m_nbad);

        const bool union_max = (m_mode == MAX_LIKELIHOOD_POS);
        double acc = -std::numeric_limits<double>::infinity();
        bool have_acc = false;
        double best = -std::numeric_limits<double>::infinity();
        size_t best_i = 0;
        int best_dir = 1;
        int count = 0;
        bool any = false;

        for (size_t i = i_min; i <= i_max; ++i) {
            // A Potts has no prior, so an ambiguous base leaves the anchor
            // unscorable rather than charged some fallback energy. The prefix
            // sum makes that an O(1) test.
            if (m_nbad[i + motif_len] - m_nbad[i] != 0)
                continue;
            int dir = 1;
            const double u = anchor_value(&m_codes[i], union_max, dir);
            any = true;
            if (u > best) {
                best = u;
                best_i = i;
                best_dir = dir;
            }
            if (m_mode == TOTAL_LIKELIHOOD) {
                // Seeded from the first scorable anchor rather than from -inf:
                // util.h's double log_sum_log() (util.h:57) has no isinf()
                // guard, unlike the float overload at util.h:16, so
                // log_sum_log(-inf, -inf) computes exp(NaN) = NaN. Every u here
                // is finite by construction, but seeding this way means the
                // accumulator can never touch that path. Same guard as
                // C_gseq_potts.
                if (!have_acc) {
                    acc = u;
                    have_acc = true;
                } else {
                    log_sum_log(acc, u);
                }
            } else if (m_mode == MOTIF_COUNT && u >= m_score_thresh) {
                ++count;
            }
        }

        if (any)
            m_last_max_score = best;

        switch (m_mode) {
        case MOTIF_COUNT:
            return (float)count;
        case TOTAL_LIKELIHOOD:
            return any ? (float)acc : std::numeric_limits<float>::quiet_NaN();
        case MAX_LIKELIHOOD:
            return any ? (float)best : std::numeric_limits<float>::quiet_NaN();
        case MAX_LIKELIHOOD_POS:
            return any ? compute_position_result(best_i, tlen, motif_len, best_dir)
                       : std::numeric_limits<float>::quiet_NaN();
        }
        return std::numeric_limits<float>::quiet_NaN();
    } catch (TGLException &e) {
        // Exactly what PWMScorer::score_interval() catches: a sequence read
        // that fails on one interval reports NaN rather than aborting the scan,
        // but nothing wider than a TGLException is swallowed.
        return std::numeric_limits<float>::quiet_NaN();
    }
}
