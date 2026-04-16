#ifndef GLMVARPROCESSOR_H_
#define GLMVARPROCESSOR_H_

#include "TrackExpressionVars.h"
#include "GenomeTrack1D.h"
#include "GenomeTrackFixedBin.h"
#include "GenomeTrackSparse.h"
#include <vector>
#include <cstdint>
#include <cmath>

class GlmVarProcessor {
public:
    GlmVarProcessor(rdb::IntervUtils &iu) : m_iu(iu) {}
    ~GlmVarProcessor() = default;

    // Pre-build filtered GLM var list (call once at start_chrom)
    void prepare_batch(TrackExpressionVars::Track_vars &track_vars);

    void process_glm_vars(
        TrackExpressionVars::Track_vars &track_vars,
        const GInterval &interval,
        unsigned idx
    );

private:
    rdb::IntervUtils &m_iu;
    std::vector<TrackExpressionVars::Track_var*> m_glm_vars;

    // Scratch buffers
    std::vector<float> m_raw_bins;
    std::vector<bool> m_inter_referenced;  // per-entry: true if used by any interaction
    const void *m_inter_referenced_owner{nullptr};  // tracks which var's interactions built the bitmap

    void process_single_glm_var(
        TrackExpressionVars::Track_var &var,
        const GInterval &interval,
        unsigned idx
    );

    // Aggregate a track window [start, end) using sum or lse
    double aggregate_window(
        const TrackExpressionVars::Track_var::GlmEntry &entry,
        int64_t start, int64_t end,
        bool is_lse
    );

    // Read a single bin value at a genome position
    double read_single_bin(
        const TrackExpressionVars::Track_var::GlmEntry &entry,
        int64_t pos
    );

    // Inline helpers
    static inline double apply_scaling(
        double raw,
        const TrackExpressionVars::Track_var::GlmScaling &s,
        double scale_factor)
    {
        if (!s.enabled) return raw;
        if (s.simple_cap) {
            // Simple cap-and-divide: min(raw, max_cap) / dis_from_cap
            double capped = std::min(raw, s.max_cap);
            double scaled = capped / s.dis_from_cap;
            return std::isfinite(scaled) ? scaled : 0.0;
        }
        double ceiled = std::min(raw - s.max_cap, 0.0);
        double floored = std::max(ceiled, -s.dis_from_cap);
        double scaled = scale_factor * (floored + s.dis_from_cap) / s.dis_from_cap;
        return std::isfinite(scaled) ? scaled : 0.0;
    }

    // Fast exp approximation for bounded inputs (GLM logistic range ~[-20, 20]).
    // Splits x = n*ln2 + r (r in [0, ln2)), uses ldexp(2^n) and a degree-5
    // Taylor series for exp(r). Max relative error < 0.01% in [-20, 20].
    static inline double fast_exp(double x)
    {
        if (x < -20.0) return 0.0;
        if (x > 20.0) return std::exp(x);

        const double LN2 = 0.6931471805599453;
        const double LOG2E = 1.4426950408889634;
        double y = x * LOG2E;
        double n = std::floor(y);
        double r = x - n * LN2;  // r in [0, ln2)

        // Taylor exp(r) for r in [0, 0.693]: 1 + r + r²/2 + r³/6 + r⁴/24 + r⁵/120
        // Horner form: 1 + r*(1 + r*(1/2 + r*(1/6 + r*(1/24 + r/120))))
        double p = 1.0 + r * (1.0 + r * (0.5 + r * (1.0/6.0 + r * (1.0/24.0 + r * (1.0/120.0)))));
        return std::ldexp(p, (int)n);
    }

    static inline double apply_transform(
        double x,
        const TrackExpressionVars::Track_var::GlmTransform &t)
    {
        double input = x + t.pre_shift;
        return t.L / (1.0 + fast_exp(-t.k * (input - t.x_0))) + t.post_shift;
    }
};

#endif /* GLMVARPROCESSOR_H_ */
