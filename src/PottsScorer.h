#ifndef POTTS_SCORER_H_
#define POTTS_SCORER_H_

#include <cstdint>
#include <deque>
#include <limits>
#include <string>
#include <vector>

#include "GenomeSeqScorer.h"
#include "PottsModel.h"
#include "utils/RunningLogSumExp.h"
#include "utils/RunningMaxDeque.h"

// Scores a genomic interval under a pairwise (Potts) energy model, reducing
// over every anchor whose start falls inside the interval.
//
// A SIBLING OF PWMScorer, NOT A SUBCLASS. The reductions and the anchor
// clamping match it exactly, but PWMScorer is 1495 lines with two open bug-fix
// branches in its spatial half, and templating a per-anchor kernel out of it
// would put this change on top of that. The duplication is the cheaper risk.
//
// No spatial weighting: the potts family does not accept spat_* at all.
class PottsScorer : public GenomeSeqScorer {
public:
    enum ScoringMode {
        TOTAL_LIKELIHOOD,   // log-sum-exp over anchors
        MAX_LIKELIHOOD,     // best anchor
        MAX_LIKELIHOOD_POS, // 1-based position of the best anchor
        MOTIF_COUNT         // anchors at or above score_thresh
    };

    // score_thresh is double, unlike PWMScorer's float: the tables, the
    // accumulator and the comparison in score_direct() are all double, and
    // this value comes straight from the user, so there is no reason to round
    // it to ~7 digits on the way in.
    PottsScorer(const PottsModel &model, const std::string &genome_root,
                bool extend, ScoringMode mode, bool bidirect, char strand,
                double score_thresh);

    PottsScorer(const PottsModel &model, GenomeSeqFetch *shared_seqfetch,
                bool extend, ScoringMode mode, bool bidirect, char strand,
                double score_thresh);

    float score_interval(const GInterval &interval, const GenomeChromKey &chromkey) override;

    void invalidate_cache();

    // The best per-anchor score from the last score_interval() call, or -inf if
    // it scored nothing. SequenceVarProcessor::score_potts_var() (in
    // SequenceVarProcessor.cpp) calls this to aggregate MAX_LIKELIHOOD_POS
    // across a filter's unmasked parts - the position alone does not say which
    // part's position to keep. Mirrors PWMEditDistanceScorer::get_last_min_edits().
    //
    // On a call answered by a SLID window it comes out of the aggregator, so it
    // carries the float precision the aggregators store rather than the double
    // the direct reduction computes; the only current caller reads it under
    // MAX_LIKELIHOOD_POS, so this can affect which of two nearly-tied parts its
    // filter aggregation prefers, but not correctness beyond that. For
    // MOTIF_COUNT - the one mode that keeps no running maximum - it is left at
    // -inf; not read there today, since a filtered potts.count sums scores, not
    // this.
    double get_last_max_score() const { return m_last_max_score; }

private:
    // The per-anchor value in the ORIGINAL genome's orientation. `union_max`
    // picks the strand union: false is log-sum-exp (TOTAL_LIKELIHOOD,
    // MAX_LIKELIHOOD, MOTIF_COUNT), true is the maximum (MAX_LIKELIHOOD_POS,
    // which has to name a strand). Same asymmetry as the pwm family, and as
    // C_gseq_potts, on purpose.
    //
    // Always finite for a scorable anchor: .coerce_potts_model() rejects a
    // model whose worst-case window magnitude - W * max|e| + npair * max|J| +
    // |intercept| - does not fit in a float, so no anchor can overflow the
    // float this returns through and -inf out of the aggregators below can
    // only mean "no scorable anchor".
    double anchor_value(const int8_t *c, bool union_max, int &dir) const;

    // 1-based position of `index`, in forward-strand orientation, signed by
    // `direction` when bidirect. Transliterated from
    // PWMScorer::compute_position_result().
    float compute_position_result(size_t index, size_t target_length,
                                  size_t motif_length, int direction) const;

    // Task 6's anchor loop, extracted unchanged. With `fill_window` it also
    // stashes each anchor into the seed scratch below, indexed by window slot,
    // so seeding evaluates every anchor exactly once instead of once for the
    // answer and once for the aggregator.
    float score_direct(size_t i_min, size_t i_max, size_t motif_len, size_t tlen,
                       bool fill_window);

    // Sliding-window cache for contiguous intervals. Mirrors PWMScorer's
    // non-spatial SlideCache and reuses both of its aggregators as they are;
    // nothing spatial, because the potts family has no spat_* parameters.
    //
    // It pays more here than it does there. When the iterator advances by a
    // stride smaller than the scan window - gvtrack.iterator(sshift = -250,
    // eshift = 250) with iterator = 1 is a 501-anchor window moving 1 bp - the
    // cache evaluates `stride` anchors instead of `window_size`, and each Potts
    // anchor costs ~10x a PWM anchor.
    //
    // The window is held in ASCENDING GENOMIC ORDER in all three structures,
    // not in target-index order, which is where this departs from PWMScorer.
    // In the reverse orientation the fetched target is reverse-complemented, so
    // the anchor with the LOWEST genomic start sits at the HIGHEST target
    // index. Keying on the genome rather than on the target is what lets one
    // set of pop-front/push-back calls serve both orientations - no
    // pop_back/push_front branch, and so no O(window) maxdq rebuild on every
    // reverse-orientation step - and RunningMaxDeque's position-keyed eviction
    // needs a monotonically increasing key anyway.
    struct SlideCache {
        int chromid = -1;
        char strand_mode = 0;
        // `valid` says the geometry below describes the immediately preceding
        // call, which is what the stride is measured against. `populated` says
        // the aggregators actually hold that window. They come apart because
        // the geometry is tracked even when the window is not worth keeping -
        // see the `populate` decision in score_with_sliding_window().
        bool valid = false;
        bool populated = false;
        int64_t last_interval_start = -1;
        int64_t last_interval_end = -1;
        int64_t last_expanded_start = -1;
        size_t last_tlen = 0;
        size_t last_i_min = 0;
        size_t last_i_max = 0;
        size_t window_size = 0;
        size_t stride = 0;
        // Genomic start of the window's lowest anchor, tracked here rather than
        // read back from RunningMaxDeque::base_genomic_pos: push() resets that
        // member whenever the deque happens to be empty, and the eviction
        // arithmetic must not depend on whether it was.
        int64_t win_lo_key = 0;

        RunningLogSumExp rlse;
        RunningMaxDeque rmax;
        std::deque<uint8_t> hits;
        int hit_count = 0;
    };

    float score_with_sliding_window(const GInterval &original_interval,
                                    const GInterval &expanded_interval,
                                    size_t i_min, size_t i_max, size_t motif_len,
                                    size_t tlen);
    // With `populate` false it records the geometry and answers from
    // score_direct() without touching the aggregators, so the next call cannot
    // slide but also paid nothing for the option.
    float seed_sliding_window(const GInterval &original_interval,
                              const GInterval &expanded_interval,
                              size_t i_min, size_t i_max, size_t motif_len,
                              size_t tlen, bool populate);
    // `slid` comes back false when the slide had to give up, in which case the
    // caller re-seeds. It is an out-parameter and not a NaN return, as it is in
    // PWMScorer: NaN is a legitimate potts answer - the window whose every
    // anchor is unscorable - so overloading it would re-seed the whole window
    // at every step of an assembly gap.
    float try_slide_window(const GInterval &original_interval,
                           const GInterval &expanded_interval,
                           size_t i_min, size_t i_max, size_t motif_len,
                           size_t tlen, size_t stride, bool &slid);
    // The current window's answer, read out of the aggregators.
    float slid_answer(const GInterval &expanded_interval, size_t motif_len, size_t tlen);

    // Target index of window slot `s` (slot 0 is the lowest genomic start), and
    // the genomic start of the anchor at target index `i`. Both fold in the
    // reverse orientation's index reversal, and they are inverses.
    inline size_t slot_index(size_t s, size_t i_min, size_t i_max) const
    {
        return (m_strand == -1) ? (i_max - s) : (i_min + s);
    }
    inline int64_t anchor_key(int64_t expanded_start, size_t i, size_t motif_len,
                              size_t tlen) const
    {
        return (m_strand == -1)
                   ? expanded_start + int64_t(tlen) - int64_t(motif_len) - int64_t(i)
                   : expanded_start + int64_t(i);
    }

    PottsModel m_model;
    PottsModel m_rc;
    ScoringMode m_mode;
    bool m_bidirect;
    double m_score_thresh;
    double m_last_max_score = -std::numeric_limits<double>::infinity();

    // Scratch, reused across calls so a per-bp iterator does not reallocate.
    std::vector<int8_t> m_codes;
    std::vector<int32_t> m_nbad;

    SlideCache m_slide;
    // Seed scratch in window-slot order: the value for every mode but
    // MOTIF_COUNT, the winning strand for MAX_LIKELIHOOD_POS, the hit bit for
    // MOTIF_COUNT. Members for the same reason m_codes is.
    std::vector<float> m_win_val;
    std::vector<int8_t> m_win_dir;
    std::vector<uint8_t> m_win_hit;
};

#endif // POTTS_SCORER_H_
