#ifndef TRACKVARPROCESSOR_H_
#define TRACKVARPROCESSOR_H_

#include "TrackExpressionVars.h"
#include "GenomeTrack1D.h"
#include "GenomeTrack2D.h"
#include "GInterval.h"
#include "GInterval2D.h"
#include "DiagonalBand.h"
#include "rdbutils.h"
#include <vector>

class TrackVarProcessor {
public:
	TrackVarProcessor(rdb::IntervUtils &iu) : m_iu(iu) {}

	void process_track_vars(
		TrackExpressionVars::Track_vars &track_vars,
		const GInterval &interval,
		const GInterval2D &interval2d,
		const DiagonalBand &band,
		unsigned idx
	);

	// Weighted sum over the filter's unmasked parts (shared with the value-based vtracks)
	static double aggregate_weighted_sum_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);

	// Quantile over the filter's unmasked parts (shared with the value-based vtracks)
	static double aggregate_quantile_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts, double percentile,
	                                             rdb::IntervUtils &iu);

	// A value (pos = false) or its position (pos = true) drawn uniformly from the values of the filter's unmasked
	// parts, each value once; the track must compute SIZE (shared with the value-based vtracks)
	static double aggregate_sample_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts, bool pos);

	// "nearest" when the filter's unmasked parts hold no value: the value of the sparse, array or value-based track
	// whose unmasked stretch is nearest to the parts, the left one of two at the same distance; NaN if the
	// chromosome has no such value and for a dense track (shared with the value-based vtracks)
	static double nearest_unmasked_value(GenomeTrack1D &track, const Genome1DFilter &filter, const std::vector<GInterval> &parts);

private:
	rdb::IntervUtils &m_iu;

	void process_single_track_var_1d(
		TrackExpressionVars::Track_var &var,
		const GInterval &interval,
		unsigned idx
	);

	void process_single_track_var_2d(
		TrackExpressionVars::Track_var &var,
		const GInterval2D &interval,
		const DiagonalBand &band,
		unsigned idx
	);

	// Helpers to deduplicate switch logic in process_single_track_var_1d
	double extract_single_value(const TrackExpressionVars::Track_var &var, GenomeTrack1D &track, int64_t base_start);
	double aggregate_multi_value(const TrackExpressionVars::Track_var &var, GenomeTrack1D &track,
	                             const std::vector<GInterval> &parts, int64_t base_start);

	// Filter aggregation helpers for 1D tracks
	double aggregate_avg_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_sum_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_min_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_max_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_max_pos_abs_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_max_pos_rel_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts, int64_t base_start);
	double aggregate_min_pos_abs_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_min_pos_rel_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts, int64_t base_start);
	bool find_best_max_pos_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts, double &best_pos);
	bool find_best_min_pos_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts, double &best_pos);
	double aggregate_stddev_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_exists_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_size_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_first_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts, bool pos);
	double aggregate_last_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts, bool pos);
	double aggregate_lse_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
};

#endif /* TRACKVARPROCESSOR_H_ */

