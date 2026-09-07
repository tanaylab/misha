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
    // Called from both TrackExpressionVars::start_chrom() overloads. The
    // chromid guard in score_with_sliding_window() would catch a chromosome
    // change on its own; this also drops the deques, so a scan does not carry
    // one chromosome's window's worth of memory into the next.
    m_slide.valid = false;
    m_slide.populated = false;
    m_slide.stride = 0;
    m_slide.rlse.clear();
    m_slide.rmax.clear();
    m_slide.hits.clear();
    m_slide.hit_count = 0;

    // The seed scratch is one entry per anchor in the last seeded window, which
    // is one iterator interval's worth - a whole chromosome, if that is what the
    // iterator was. Released here rather than kept the way m_codes is, so a
    // scan that opens with one huge interval does not hold that buffer for the
    // life of the scorer. It costs one allocation per chromosome.
    std::vector<float>().swap(m_win_val);
    std::vector<int8_t>().swap(m_win_dir);
    std::vector<uint8_t>().swap(m_win_hit);
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

// Task 6's anchor loop, unchanged but for MAX_LIKELIHOOD_POS's tie-break, which
// is described where it happens below.
float PottsScorer::score_direct(size_t i_min, size_t i_max, size_t motif_len, size_t tlen,
                                bool fill_window)
{
    const bool union_max = (m_mode == MAX_LIKELIHOOD_POS);
    double acc = -std::numeric_limits<double>::infinity();
    bool have_acc = false;
    double best = -std::numeric_limits<double>::infinity();
    int count = 0;
    bool any = false;

    // MAX_LIKELIHOOD_POS's argmax is tracked separately from `best`, and
    // deliberately not the same way.
    //
    // `best` is a double, and its tie-break is "first anchor scanned wins",
    // which is the lowest TARGET index - the lowest genomic start going forward
    // and the HIGHEST in the reverse orientation, where the fetched target is
    // reverse-complemented. RunningMaxDeque cannot reproduce either: it stores
    // float, and its front is the earliest-pushed among equals, which is the
    // lowest genomic start in both orientations.
    //
    // Two rules for one question is how a cache returns a plausible wrong
    // answer, so this loop uses the deque's: maximise the FLOAT value, break a
    // tie on the lowest window slot. Slot 0 is the lowest genomic start, so
    // that is the deque's rule exactly, and the seeded, the slid and the
    // un-populated paths then all name the same anchor. Measured before this
    // was made single-valued: an all-tied model scanned with iterator = 10 and
    // no iterator shift, on the reverse orientation, reported position 10 from
    // the contiguous scan and 1 from the same intervals fetched one at a
    // time - 19 of 20 bins.
    //
    // It costs one float and one size_t comparison per anchor, against ~110
    // table lookups. The alternative - populating the window whatever the
    // stride, so the deque always has the answer - reintroduces the ~17%
    // overhead the un-populated path exists to avoid.
    float pos_best = -std::numeric_limits<float>::infinity();
    size_t pos_slot = 0;
    size_t pos_i = 0;
    int pos_dir = 1;

    if (fill_window) {
        const size_t W = i_max - i_min + 1;
        // An unscorable anchor keeps the identity these slots are filled with:
        // -INFINITY for a log-sum-exp and for a maximum, 0 for a hit count.
        // The running structures are fixed-size windows, so such an anchor
        // cannot simply be omitted without breaking the stride arithmetic.
        if (m_mode == MOTIF_COUNT) {
            m_win_hit.assign(W, 0);
        } else {
            m_win_val.assign(W, -std::numeric_limits<float>::infinity());
            if (m_mode == MAX_LIKELIHOOD_POS)
                m_win_dir.assign(W, 1);
        }
    }

    for (size_t i = i_min; i <= i_max; ++i) {
        // A Potts has no prior, so an ambiguous base leaves the anchor
        // unscorable rather than charged some fallback energy. The prefix
        // sum makes that an O(1) test.
        if (m_nbad[i + motif_len] - m_nbad[i] != 0)
            continue;
        int dir = 1;
        const double u = anchor_value(&m_codes[i], union_max, dir);
        const size_t s = (m_strand == -1) ? (i_max - i) : (i - i_min);
        if (u > best)
            best = u;
        if (m_mode == MAX_LIKELIHOOD_POS) {
            const float v = (float)u;
            if (!any || v > pos_best || (v == pos_best && s < pos_slot)) {
                pos_best = v;
                pos_slot = s;
                pos_i = i;
                pos_dir = dir;
            }
        }
        any = true;
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
        if (fill_window) {
            if (m_mode == MOTIF_COUNT) {
                m_win_hit[s] = (u >= m_score_thresh) ? 1 : 0;
            } else {
                m_win_val[s] = (float)u;
                if (m_mode == MAX_LIKELIHOOD_POS)
                    m_win_dir[s] = (int8_t)dir;
            }
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
        return any ? compute_position_result(pos_i, tlen, motif_len, pos_dir)
                   : std::numeric_limits<float>::quiet_NaN();
    }
    return std::numeric_limits<float>::quiet_NaN();
}

float PottsScorer::slid_answer(const GInterval &expanded_interval, size_t motif_len,
                               size_t tlen)
{
    // BOTH aggregators report -inf, not NaN, for a window whose every anchor
    // was pushed as -INFINITY: RunningLogSumExp::value() short-circuits on a
    // non-finite max, and RunningMaxDeque::value() returns -inf on an empty
    // deque. An interval with no scorable anchor is NaN - 0 for potts.count -
    // which is what the direct path returns, so the mapping has to happen
    // here. An -Inf leaking out instead would survive arithmetic rather than
    // poison it, which is the harder failure to notice; the model tables are
    // validated finite, so -inf here can only mean "nothing scorable".
    switch (m_mode) {
    case MOTIF_COUNT:
        // No running maximum in this mode, so m_last_max_score is left alone.
        return (float)m_slide.hit_count;
    case TOTAL_LIKELIHOOD: {
        const double v = m_slide.rlse.value();
        if (!m_slide.rlse.maxdq.empty())
            m_last_max_score = m_slide.rlse.maxdq.front();
        return std::isfinite(v) ? (float)v : std::numeric_limits<float>::quiet_NaN();
    }
    case MAX_LIKELIHOOD: {
        const double v = m_slide.rmax.value();
        m_last_max_score = v;
        return std::isfinite(v) ? (float)v : std::numeric_limits<float>::quiet_NaN();
    }
    case MAX_LIKELIHOOD_POS: {
        const double v = m_slide.rmax.value();
        m_last_max_score = v;
        if (!std::isfinite(v))
            return std::numeric_limits<float>::quiet_NaN();
        // The deque holds absolute genomic starts, and the reported position is
        // an index into the CURRENT fetched target, so it has to be mapped back
        // through this call's expanded.start - the whole reason
        // RunningMaxDeque stores coordinates rather than indices.
        const int64_t off =
            m_slide.rmax.argmax_genomic_position() - int64_t(expanded_interval.start);
        const size_t i = (m_strand == -1)
                             ? (tlen - motif_len - (size_t)off)
                             : (size_t)off;
        return compute_position_result(i, tlen, motif_len,
                                       m_slide.rmax.argmax_direction());
    }
    }
    return std::numeric_limits<float>::quiet_NaN();
}

float PottsScorer::seed_sliding_window(const GInterval &original_interval,
                                       const GInterval &expanded_interval,
                                       size_t i_min, size_t i_max, size_t motif_len,
                                       size_t tlen, bool populate)
{
    const size_t W = i_max - i_min + 1;
    const float direct = score_direct(i_min, i_max, motif_len, tlen, populate);

    m_slide.valid = true;
    m_slide.populated = populate;
    m_slide.chromid = original_interval.chromid;
    m_slide.strand_mode = m_strand;
    m_slide.last_interval_start = original_interval.start;
    m_slide.last_interval_end = original_interval.end;
    m_slide.last_expanded_start = expanded_interval.start;
    m_slide.last_tlen = tlen;
    m_slide.last_i_min = i_min;
    m_slide.last_i_max = i_max;
    m_slide.window_size = W;
    m_slide.stride = 0;
    // Slot s of the window sits at win_lo_key + s, in either orientation.
    m_slide.win_lo_key =
        anchor_key(expanded_interval.start, slot_index(0, i_min, i_max), motif_len, tlen);

    if (!populate)
        return direct;

    switch (m_mode) {
    case TOTAL_LIKELIHOOD:
        m_slide.rlse.init(m_win_val); // exactly W long, in window-slot order
        return direct;

    case MAX_LIKELIHOOD:
    case MAX_LIKELIHOOD_POS:
        m_slide.rmax.clear();
        for (size_t s = 0; s < W; ++s) {
            const int dir = (m_mode == MAX_LIKELIHOOD_POS) ? int(m_win_dir[s]) : 1;
            m_slide.rmax.push(m_win_val[s], dir, m_slide.win_lo_key + int64_t(s));
        }
        // MAX_LIKELIHOOD's cached answer is the maximum of the same floats the
        // direct reduction maximised as doubles, and float() is monotone, so
        // the two are the same number and the direct one is returned.
        //
        // Both modes return score_direct()'s answer, and for both it is the
        // same number the deque would give.
        //
        // MAX_LIKELIHOOD: the cached value is the maximum of the same floats
        // the direct reduction maximised as doubles, and float() is monotone,
        // so max-of-floats == float(max-of-doubles).
        //
        // MAX_LIKELIHOOD_POS: score_direct() breaks its tie on the lowest
        // window slot precisely so this holds - see the comment on `pos_best`
        // there. Relative to Task 6 that means a tie names a different, equally
        // maximal anchor: measured over 72,000 scanned positions in repeat-rich
        // sequence at W = 4, 6, 8, 4,087 positions name a different anchor and
        // 0 report a different SCORE. Task 6's own tests already treat the
        // argmax index as undefined under a tie, and "lowest genomic start
        // wins" is one rule in both orientations where "lowest target index
        // wins" silently flips direction with the fetch - the artefact
        // .potts_params() works around by clamping strand to 1 under bidirect.
        return direct;

    case MOTIF_COUNT:
        m_slide.hits.assign(m_win_hit.begin(), m_win_hit.end());
        m_slide.hit_count = 0;
        for (size_t s = 0; s < W; ++s)
            m_slide.hit_count += m_win_hit[s];
        return direct;
    }
    return std::numeric_limits<float>::quiet_NaN();
}

float PottsScorer::try_slide_window(const GInterval &original_interval,
                                    const GInterval &expanded_interval,
                                    size_t i_min, size_t i_max, size_t motif_len,
                                    size_t tlen, size_t stride, bool &slid)
{
    slid = false;
    const size_t W = m_slide.window_size;

    // Every incoming target index lies in [i_min, i_max], so one bound check
    // covers them all - and it happens before anything is mutated, so a give-up
    // cannot leave the aggregators half-advanced.
    if (i_max + motif_len > tlen || stride >= W)
        return std::numeric_limits<float>::quiet_NaN();

    // The `stride` anchors leaving are the ones with the lowest genomic starts,
    // the `stride` arriving are the highest: slots [W - stride, W) of the NEW
    // window, whose genomic starts are win_lo_key + W ... win_lo_key + W +
    // stride - 1. In the forward orientation those are the target indices just
    // below i_max; in the reverse, where the target is reverse-complemented,
    // the ones just above i_min.
    const bool union_max = (m_mode == MAX_LIKELIHOOD_POS);

    switch (m_mode) {
    case TOTAL_LIKELIHOOD: {
        for (size_t k = 0; k < stride; ++k)
            m_slide.rlse.pop_front();

        // RunningLogSumExp::push(-INFINITY) while its running maximum is
        // already -inf evaluates exp(-inf + inf) = NaN into sum_scaled, and the
        // NaN then survives into the first FINITE push - so the window reports
        // NaN for the rest of its life although it has scorable anchors. It is
        // not reachable from the pwm family, whose per-anchor value is -inf
        // only for a zero-prior PSSM, which is why the shared utility has never
        // needed a guard; it is reached here by an unscorable anchor, which is
        // ordinary genomic N.
        //
        // Guarded rather than fixed in RunningLogSumExp.h, which pwm shares:
        // the batch that pushes -Inf and then a real value while the window
        // holds nothing scorable gives up and re-seeds, which rebuilds the
        // accumulator from scratch. A window that stays entirely unscorable is
        // left alone - value() short-circuits on the non-finite maximum and
        // never reads sum_scaled - so a scan of an assembly gap still slides.
        const bool nothing_scorable = !std::isfinite(m_slide.rlse.M);
        bool incoming_scorable = false;
        for (size_t k = 0; k < stride; ++k) {
            const size_t i = slot_index(W - stride + k, i_min, i_max);
            int dir = 1;
            const bool bad = (m_nbad[i + motif_len] - m_nbad[i] != 0);
            if (!bad)
                incoming_scorable = true;
            m_slide.rlse.push(bad ? -std::numeric_limits<float>::infinity()
                                  : (float)anchor_value(&m_codes[i], union_max, dir));
        }
        if (nothing_scorable && incoming_scorable)
            return std::numeric_limits<float>::quiet_NaN(); // slid stays false
        break;
    }

    case MAX_LIKELIHOOD:
    case MAX_LIKELIHOOD_POS:
        m_slide.rmax.pop_front(m_slide.win_lo_key + int64_t(stride));
        for (size_t k = 0; k < stride; ++k) {
            const size_t i = slot_index(W - stride + k, i_min, i_max);
            int dir = 1;
            const bool bad = (m_nbad[i + motif_len] - m_nbad[i] != 0);
            const float v = bad ? -std::numeric_limits<float>::infinity()
                                : (float)anchor_value(&m_codes[i], union_max, dir);
            m_slide.rmax.push(v, dir, m_slide.win_lo_key + int64_t(W) + int64_t(k));
        }
        break;

    case MOTIF_COUNT:
        for (size_t k = 0; k < stride && !m_slide.hits.empty(); ++k) {
            m_slide.hit_count -= m_slide.hits.front();
            m_slide.hits.pop_front();
        }
        for (size_t k = 0; k < stride; ++k) {
            const size_t i = slot_index(W - stride + k, i_min, i_max);
            int dir = 1;
            const bool bad = (m_nbad[i + motif_len] - m_nbad[i] != 0);
            // The threshold comparison stays in double, as it is in
            // score_direct(), so a cached count and a re-seeded one cannot
            // differ over an anchor sitting on the threshold.
            const uint8_t hit =
                (!bad && anchor_value(&m_codes[i], union_max, dir) >= m_score_thresh) ? 1
                                                                                     : 0;
            m_slide.hits.push_back(hit);
            m_slide.hit_count += hit;
        }
        break;
    }

    m_slide.win_lo_key += int64_t(stride);
    m_slide.last_interval_start = original_interval.start;
    m_slide.last_interval_end = original_interval.end;
    m_slide.last_expanded_start = expanded_interval.start;
    m_slide.stride = stride;
    slid = true;
    return slid_answer(expanded_interval, motif_len, tlen);
}

float PottsScorer::score_with_sliding_window(const GInterval &original_interval,
                                             const GInterval &expanded_interval,
                                             size_t i_min, size_t i_max, size_t motif_len,
                                             size_t tlen)
{
    size_t stride = 0;
    if (m_slide.valid) {
        const int64_t step_start = original_interval.start - m_slide.last_interval_start;
        const int64_t step_end = original_interval.end - m_slide.last_interval_end;
        if (step_start > 0 && step_start == step_end)
            stride = (size_t)step_start;
    }

    const bool can_slide =
        m_slide.valid && m_slide.populated &&
        m_slide.chromid == original_interval.chromid &&
        m_slide.strand_mode == m_strand && stride > 0 &&
        // STRICTLY less, where PWMScorer allows stride == window_size. At
        // stride == window_size consecutive windows share no anchor, so a
        // "slide" evaluates the whole window again and pays the aggregator
        // churn on top. That is not a corner case: it is every scan with no
        // gvtrack.iterator shift, where the iterator intervals tile the genome
        // instead of overlapping. So the rule this guard now encodes is simply
        // "slide only where the windows overlap".
        stride < m_slide.window_size &&
        (m_slide.stride == 0 || stride == m_slide.stride) &&
        i_min == m_slide.last_i_min && i_max == m_slide.last_i_max &&
        // Neither of these is in PWMScorer's guard, and the genomic-key
        // arithmetic rests on both: the target has to be the same length, since
        // the reverse orientation keys an anchor on its distance from the
        // target's END, and it has to have moved by exactly `stride`, since
        // that is the amount every cached key is taken to have shifted by.
        // Both hold for a uniform iterator and neither is free to assume near a
        // chromosome edge.
        tlen == m_slide.last_tlen &&
        int64_t(expanded_interval.start) - m_slide.last_expanded_start == int64_t(stride);

    if (can_slide) {
        bool slid = false;
        const float result = try_slide_window(original_interval, expanded_interval, i_min,
                                              i_max, motif_len, tlen, stride, slid);
        if (slid)
            return result;
    }

    // Populate the aggregators only where a later slide could reuse an anchor.
    // A stride at or beyond the window size shares nothing with its
    // predecessor, so keeping the window would be pure overhead - measured at
    // ~17% on a no-shift iterator = 500 scan of hg38 chr1 before this test was
    // added. stride == 0 means "no usable history yet" (the first call of a
    // scan, or a non-monotone step), where there is nothing to predict from and
    // the window has to be built to have a chance of paying.
    const size_t W = i_max - i_min + 1;
    const bool populate = (stride == 0) || (stride < W);
    return seed_sliding_window(original_interval, expanded_interval, i_min, i_max,
                               motif_len, tlen, populate);
}

float PottsScorer::score_interval(const GInterval &interval, const GenomeChromKey &chromkey)
{
    m_last_max_score = -std::numeric_limits<double>::infinity();

    const int64_t motif_length = m_model.width();
    GInterval expanded = calculate_expanded_interval(interval, chromkey, motif_length);
    expanded.strand = m_strand;

    if (!m_extend && (expanded.end - expanded.start) < motif_length) {
        // Every early return drops the cache. The guard below cannot tell a
        // skipped interval from a contiguous step - expanded.start tracks
        // interval.start exactly, because extension is END-only - so a call
        // that answers without touching the window must not leave one behind
        // for its successor to slide off.
        invalidate_cache();
        return std::numeric_limits<float>::quiet_NaN();
    }

    try {
        std::vector<char> seq;
        m_seqfetch_ptr->read_interval(expanded, chromkey, seq);
        const std::string target(seq.begin(), seq.end());

        const size_t tlen = target.size();
        const size_t motif_len = (size_t)motif_length;
        if (tlen < motif_len) {
            invalidate_cache();
            return (m_mode == MOTIF_COUNT) ? 0.0f
                                           : std::numeric_limits<float>::quiet_NaN();
        }

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

        return score_with_sliding_window(interval, expanded, i_min, i_max, motif_len, tlen);
    } catch (TGLException &e) {
        // Exactly what PWMScorer::score_interval() catches: a sequence read
        // that fails on one interval reports NaN rather than aborting the scan,
        // but nothing wider than a TGLException is swallowed.
        invalidate_cache();
        return std::numeric_limits<float>::quiet_NaN();
    }
}
