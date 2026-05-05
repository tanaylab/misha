#ifndef GLMFEATUREEXTRACTOR_H_
#define GLMFEATUREEXTRACTOR_H_

#include "rdbutils.h"
#include "GenomeTrackFixedBin.h"
#include "GenomeTrackSparse.h"
#include "GenomeTrack.h"

#include <vector>
#include <string>
#include <memory>
#include <cstdint>
#include <cmath>
#include <algorithm>

struct LmTransformConfig {
    double L;
    double k;
    double x_0;
    double pre_shift;
    double post_shift;
};

struct LmScalingConfig {
    double max_cap;
    double dis_from_cap;
    double scale_factor;
};

class GlmFeatureExtractor {
public:
    GlmFeatureExtractor(rdb::IntervUtils &iu) : m_iu(iu) {}

    void extract(
        const std::vector<std::string> &track_names,
        const int *peak_chromids,
        const int64_t *peak_starts,
        const int64_t *peak_ends,
        int n_peaks,
        int tile_size,
        int flank_size,
        const std::vector<LmScalingConfig> &scaling,
        const std::vector<LmTransformConfig> &transforms,
        const std::string &gc_track_name,
        double gc_scale_factor,
        double *output,
        int n_cols,
        int n_threads = 1
    );

private:
    rdb::IntervUtils &m_iu;

    struct TrackHandle {
        std::string track_dir;
        GenomeTrack::Type type;
        std::shared_ptr<GenomeTrack> track;
        GenomeTrackFixedBin *fixedbin = nullptr;
        GenomeTrackSparse *sparse = nullptr;
        unsigned bin_size = 0;
    };

    void open_track(TrackHandle &handle, int chromid);
    // Worker-callable version: takes the chromkey directly so multiple
    // threads don't race on m_iu.get_chromkey() (which is itself harmless
    // but keeps the worker free of any shared mutable state).
    static void open_track_static(TrackHandle &handle, int chromid,
                                  const GenomeChromKey &chromkey);

    static double aggregate_lse(const TrackHandle &handle, int64_t start, int64_t end);
    static double aggregate_sum(const TrackHandle &handle, int64_t start, int64_t end);

    // Binary search: find first interval index where intervals[i].end > pos
    static size_t sparse_lower_bound(const GIntervals &intervals, int64_t pos);

    static inline double apply_scaling(double raw, const LmScalingConfig &s) {
        double ceiled = std::min(raw - s.max_cap, 0.0);
        double floored = std::max(ceiled, -s.dis_from_cap);
        double scaled = s.scale_factor * (floored + s.dis_from_cap) / s.dis_from_cap;
        return std::isfinite(scaled) ? scaled : 0.0;
    }

    static inline double apply_transform(double x, const LmTransformConfig &t) {
        double input = x + t.pre_shift;
        return t.L / (1.0 + std::exp(-t.k * (input - t.x_0))) + t.post_shift;
    }
};

#endif /* GLMFEATUREEXTRACTOR_H_ */
