#ifndef POTTS_SCORER_H_
#define POTTS_SCORER_H_

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "GenomeSeqScorer.h"
#include "PottsModel.h"

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
    // accumulator and the comparison in score_interval() are all double, and
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
    // it scored nothing. SequenceVarProcessor needs it to aggregate
    // MAX_LIKELIHOOD_POS across a filter's unmasked parts - the position alone
    // does not say which part's position to keep. Mirrors
    // PWMEditDistanceScorer::get_last_min_edits().
    double get_last_max_score() const { return m_last_max_score; }

private:
    // The per-anchor value in the ORIGINAL genome's orientation. `union_max`
    // picks the strand union: false is log-sum-exp (TOTAL_LIKELIHOOD,
    // MAX_LIKELIHOOD, MOTIF_COUNT), true is the maximum (MAX_LIKELIHOOD_POS,
    // which has to name a strand). Same asymmetry as the pwm family, and as
    // C_gseq_potts, on purpose.
    double anchor_value(const int8_t *c, bool union_max, int &dir) const;

    // 1-based position of `index`, in forward-strand orientation, signed by
    // `direction` when bidirect. Transliterated from
    // PWMScorer::compute_position_result().
    float compute_position_result(size_t index, size_t target_length,
                                  size_t motif_length, int direction) const;

    PottsModel m_model;
    PottsModel m_rc;
    ScoringMode m_mode;
    bool m_bidirect;
    double m_score_thresh;
    double m_last_max_score = -std::numeric_limits<double>::infinity();

    // Scratch, reused across calls so a per-bp iterator does not reallocate.
    std::vector<int8_t> m_codes;
    std::vector<int32_t> m_nbad;
};

#endif // POTTS_SCORER_H_
