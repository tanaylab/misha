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
	double aggregate_quantile_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts, double percentile);
	double aggregate_exists_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_size_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_sample_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_sample_pos_abs_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts);
	double aggregate_sample_pos_rel_with_filter(GenomeTrack1D &track, const std::vector<GInterval> &parts, int64_t base_start);
};

#endif /* TRACKVARPROCESSOR_H_ */

