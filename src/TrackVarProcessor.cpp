#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

#include "TrackVarProcessor.h"
#include "TrackExpressionVars.h"
#include "GenomeTrack.h"
#include "GenomeTrack1D.h"
#include "GenomeTrack2D.h"
#include "StreamPercentiler.h"
#include "rdbutils.h"

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

using namespace rdb;
using namespace std;

void TrackVarProcessor::process_track_vars(
	TrackExpressionVars::Track_vars &track_vars,
	const GInterval &interval,
	const GInterval2D &interval2d,
	const DiagonalBand &band,
	unsigned idx)
{
	for (TrackExpressionVars::Track_vars::iterator ivar = track_vars.begin(); ivar != track_vars.end(); ++ivar) {
		// Skip sequence-based vtracks (already processed separately)
		if (TrackExpressionVars::is_sequence_based_function(ivar->val_func)) {
			continue;
		}

		if (GenomeTrack::is_1d(ivar->track_n_imdf->type)) {
			process_single_track_var_1d(*ivar, interval, idx);
		} else {
			process_single_track_var_2d(*ivar, interval2d, band, idx);
		}
	}
}

void TrackVarProcessor::process_single_track_var_1d(
	TrackExpressionVars::Track_var &var,
	const GInterval &interval,
	unsigned idx)
{
	GenomeTrack1D &track = *(GenomeTrack1D *)var.track_n_imdf->track;
	const GInterval &base_interval = var.track_n_imdf->imdf1d ?
		var.track_n_imdf->imdf1d->interval : interval;
	const int64_t base_start = base_interval.start;

	if (var.track_n_imdf->imdf1d && var.track_n_imdf->imdf1d->out_of_range) {
		var.var[idx] = numeric_limits<double>::quiet_NaN();
		return;
	}

	// Check if filter applies and get unmasked parts
	std::vector<GInterval> unmasked_parts;
	bool has_filter = false;

	if (var.filter) {
		const GInterval &eval_interval = var.track_n_imdf->imdf1d ?
			var.track_n_imdf->imdf1d->interval : interval;

		var.filter->subtract(eval_interval, unmasked_parts);
		has_filter = true;

		if (unmasked_parts.empty()) {
			// Completely masked - return NaN
			var.var[idx] = numeric_limits<double>::quiet_NaN();
			return;
		}
	}

	// If filter exists and resulted in multiple unmasked parts, aggregate across them
	// If only one part (or filter but no filtering needed), use normal path below
	if (has_filter && unmasked_parts.size() > 1) {
		// Aggregate over unmasked parts using helper methods
		double result = std::numeric_limits<double>::quiet_NaN();

		switch (var.val_func) {
		case TrackExpressionVars::Track_var::REG:
		case TrackExpressionVars::Track_var::PV:
			result = aggregate_avg_with_filter(track, unmasked_parts);
			break;
		case TrackExpressionVars::Track_var::STDDEV:
			result = aggregate_stddev_with_filter(track, unmasked_parts);
			break;
		case TrackExpressionVars::Track_var::SUM:
			result = aggregate_sum_with_filter(track, unmasked_parts);
			break;
		case TrackExpressionVars::Track_var::REG_MIN:
		case TrackExpressionVars::Track_var::PV_MIN:
			result = aggregate_min_with_filter(track, unmasked_parts);
			break;
		case TrackExpressionVars::Track_var::REG_MAX:
		case TrackExpressionVars::Track_var::PV_MAX:
			result = aggregate_max_with_filter(track, unmasked_parts);
			break;
		case TrackExpressionVars::Track_var::MAX_POS_ABS:
			result = aggregate_max_pos_abs_with_filter(track, unmasked_parts);
			break;
		case TrackExpressionVars::Track_var::MAX_POS_REL:
			result = aggregate_max_pos_rel_with_filter(track, unmasked_parts, base_start);
			break;
		case TrackExpressionVars::Track_var::MIN_POS_ABS:
			result = aggregate_min_pos_abs_with_filter(track, unmasked_parts);
			break;
		case TrackExpressionVars::Track_var::MIN_POS_REL:
			result = aggregate_min_pos_rel_with_filter(track, unmasked_parts, base_start);
			break;
		case TrackExpressionVars::Track_var::REG_NEAREST:
			// For nearest, use the first unmasked part (closest to start)
			track.read_interval(unmasked_parts[0]);
			result = track.last_nearest();
			break;
		case TrackExpressionVars::Track_var::QUANTILE:
			result = aggregate_quantile_with_filter(track, unmasked_parts, var.percentile);
			break;
		case TrackExpressionVars::Track_var::EXISTS:
			result = aggregate_exists_with_filter(track, unmasked_parts);
			break;
		case TrackExpressionVars::Track_var::SIZE:
			result = aggregate_size_with_filter(track, unmasked_parts);
			break;
		case TrackExpressionVars::Track_var::SAMPLE:
			result = aggregate_sample_with_filter(track, unmasked_parts);
			break;
		case TrackExpressionVars::Track_var::SAMPLE_POS_ABS:
			result = aggregate_sample_pos_abs_with_filter(track, unmasked_parts);
			break;
		case TrackExpressionVars::Track_var::SAMPLE_POS_REL:
			result = aggregate_sample_pos_rel_with_filter(track, unmasked_parts, base_start);
			break;
		case TrackExpressionVars::Track_var::FIRST:
			// For first, use the first unmasked part
			track.read_interval(unmasked_parts[0]);
			result = track.last_first();
			break;
		case TrackExpressionVars::Track_var::FIRST_POS_ABS:
			// For first position, use the first unmasked part
			track.read_interval(unmasked_parts[0]);
			result = track.last_first_pos();
			break;
		case TrackExpressionVars::Track_var::FIRST_POS_REL:
			// For first position relative, use the first unmasked part
			track.read_interval(unmasked_parts[0]);
			result = track.last_first_pos() - base_start;
			break;
		case TrackExpressionVars::Track_var::LAST:
			// For last, use the last unmasked part
			track.read_interval(unmasked_parts.back());
			result = track.last_last();
			break;
		case TrackExpressionVars::Track_var::LAST_POS_ABS:
			// For last position, use the last unmasked part
			track.read_interval(unmasked_parts.back());
			result = track.last_last_pos();
			break;
		case TrackExpressionVars::Track_var::LAST_POS_REL:
			// For last position relative, use the last unmasked part
			track.read_interval(unmasked_parts.back());
			result = track.last_last_pos() - base_start;
			break;
		// Sequence-based functions are already handled above
		default:
			if (!TrackExpressionVars::is_sequence_based_function(var.val_func))
				verror("Internal error: unsupported function %d", var.val_func);
			break;
		}

		var.var[idx] = result;

		// Restore track state after filter processing
		// The aggregate functions call track.read_interval() which modifies the shared track state.
		// Other vtracks sharing this Track_n_imdf need the original interval restored.
		track.read_interval(base_interval);
	} else if (has_filter && unmasked_parts.size() == 1) {
		// Single unmasked part - read just that part
		track.read_interval(unmasked_parts[0]);

		switch (var.val_func) {
		case TrackExpressionVars::Track_var::REG:
		case TrackExpressionVars::Track_var::PV:
			var.var[idx] = track.last_avg();
			break;
		case TrackExpressionVars::Track_var::REG_MIN:
		case TrackExpressionVars::Track_var::PV_MIN:
			var.var[idx] = track.last_min();
			break;
		case TrackExpressionVars::Track_var::REG_MAX:
		case TrackExpressionVars::Track_var::PV_MAX:
			var.var[idx] = track.last_max();
			break;
		case TrackExpressionVars::Track_var::MAX_POS_ABS:
			var.var[idx] = track.last_max_pos();
			break;
		case TrackExpressionVars::Track_var::MAX_POS_REL:
			var.var[idx] = track.last_max_pos() - base_start;
			break;
		case TrackExpressionVars::Track_var::MIN_POS_ABS:
			var.var[idx] = track.last_min_pos();
			break;
		case TrackExpressionVars::Track_var::MIN_POS_REL:
			var.var[idx] = track.last_min_pos() - base_start;
			break;
		case TrackExpressionVars::Track_var::REG_NEAREST:
			var.var[idx] = track.last_nearest();
			break;
		case TrackExpressionVars::Track_var::STDDEV:
			var.var[idx] = track.last_stddev();
			break;
		case TrackExpressionVars::Track_var::SUM:
			var.var[idx] = track.last_sum();
			break;
		case TrackExpressionVars::Track_var::QUANTILE:
			var.var[idx] = track.last_quantile(var.percentile);
			break;
		case TrackExpressionVars::Track_var::EXISTS:
			var.var[idx] = track.last_exists();
			break;
		case TrackExpressionVars::Track_var::SIZE:
			var.var[idx] = track.last_size();
			break;
		case TrackExpressionVars::Track_var::SAMPLE:
			var.var[idx] = track.last_sample();
			break;
		case TrackExpressionVars::Track_var::SAMPLE_POS_ABS:
			var.var[idx] = track.last_sample_pos();
			break;
		case TrackExpressionVars::Track_var::SAMPLE_POS_REL:
			var.var[idx] = track.last_sample_pos() - base_start;
			break;
		case TrackExpressionVars::Track_var::FIRST:
			var.var[idx] = track.last_first();
			break;
		case TrackExpressionVars::Track_var::FIRST_POS_ABS:
			var.var[idx] = track.last_first_pos();
			break;
		case TrackExpressionVars::Track_var::FIRST_POS_REL:
			var.var[idx] = track.last_first_pos() - base_start;
			break;
		case TrackExpressionVars::Track_var::LAST:
			var.var[idx] = track.last_last();
			break;
		case TrackExpressionVars::Track_var::LAST_POS_ABS:
			var.var[idx] = track.last_last_pos();
			break;
		case TrackExpressionVars::Track_var::LAST_POS_REL:
			var.var[idx] = track.last_last_pos() - base_start;
			break;
		// Sequence-based functions are already handled above
		default:
			if (!TrackExpressionVars::is_sequence_based_function(var.val_func))
				verror("Internal error: unsupported function %d", var.val_func);
			break;
		}

		// Restore track state after processing single filtered part
		track.read_interval(base_interval);
	} else {
		// No filter or single unmasked part - use normal path
		switch (var.val_func) {
		case TrackExpressionVars::Track_var::REG:
		case TrackExpressionVars::Track_var::PV:
			var.var[idx] = track.last_avg();
			break;
		case TrackExpressionVars::Track_var::REG_MIN:
		case TrackExpressionVars::Track_var::PV_MIN:
			var.var[idx] = track.last_min();
			break;
		case TrackExpressionVars::Track_var::REG_MAX:
		case TrackExpressionVars::Track_var::PV_MAX:
			var.var[idx] = track.last_max();
			break;
		case TrackExpressionVars::Track_var::MAX_POS_ABS:
			var.var[idx] = track.last_max_pos();
			break;
		case TrackExpressionVars::Track_var::MAX_POS_REL:
			var.var[idx] = track.last_max_pos() - base_start;
			break;
		case TrackExpressionVars::Track_var::MIN_POS_ABS:
			var.var[idx] = track.last_min_pos();
			break;
		case TrackExpressionVars::Track_var::MIN_POS_REL:
			var.var[idx] = track.last_min_pos() - base_start;
			break;
		case TrackExpressionVars::Track_var::REG_NEAREST:
			var.var[idx] = track.last_nearest();
			break;
		case TrackExpressionVars::Track_var::STDDEV:
			var.var[idx] = track.last_stddev();
			break;
		case TrackExpressionVars::Track_var::SUM:
			var.var[idx] = track.last_sum();
			break;
		case TrackExpressionVars::Track_var::QUANTILE:
			var.var[idx] = track.last_quantile(var.percentile);
			break;
		case TrackExpressionVars::Track_var::EXISTS:
			var.var[idx] = track.last_exists();
			break;
		case TrackExpressionVars::Track_var::SIZE:
			var.var[idx] = track.last_size();
			break;
		case TrackExpressionVars::Track_var::SAMPLE:
			var.var[idx] = track.last_sample();
			break;
		case TrackExpressionVars::Track_var::SAMPLE_POS_ABS:
			var.var[idx] = track.last_sample_pos();
			break;
		case TrackExpressionVars::Track_var::SAMPLE_POS_REL:
			var.var[idx] = track.last_sample_pos() - base_start;
			break;
		case TrackExpressionVars::Track_var::FIRST:
			var.var[idx] = track.last_first();
			break;
		case TrackExpressionVars::Track_var::FIRST_POS_ABS:
			var.var[idx] = track.last_first_pos();
			break;
		case TrackExpressionVars::Track_var::FIRST_POS_REL:
			var.var[idx] = track.last_first_pos() - base_start;
			break;
		case TrackExpressionVars::Track_var::LAST:
			var.var[idx] = track.last_last();
			break;
		case TrackExpressionVars::Track_var::LAST_POS_ABS:
			var.var[idx] = track.last_last_pos();
			break;
		case TrackExpressionVars::Track_var::LAST_POS_REL:
			var.var[idx] = track.last_last_pos() - base_start;
			break;
		// Sequence-based functions are already handled above
		default:
			if (!TrackExpressionVars::is_sequence_based_function(var.val_func))
				verror("Internal error: unsupported function %d", var.val_func);
			break;
		}

		if (var.requires_pv) {
			double val = var.var[idx];
			if (!std::isnan(val)) {
				int bin = var.pv_binned.binfinder.val2bin(val);
				if (bin < 0) {
					if (val <= var.pv_binned.binfinder.get_breaks().front())
						var.var[idx] = var.pv_binned.bins[0];
					else
						var.var[idx] = 1.;
				} else
					var.var[idx] = var.pv_binned.bins[bin];
			}
		}
	}
}

void TrackVarProcessor::process_single_track_var_2d(
	TrackExpressionVars::Track_var &var,
	const GInterval2D &interval,
	const DiagonalBand &band,
	unsigned idx)
{
	GenomeTrack2D &track = *(GenomeTrack2D *)var.track_n_imdf->track;

	if (var.track_n_imdf->imdf2d && var.track_n_imdf->imdf2d->out_of_range) {
		var.var[idx] = numeric_limits<double>::quiet_NaN();
		return;
	}

	switch (var.val_func) {
	case TrackExpressionVars::Track_var::REG:
		var.var[idx] = track.last_avg();
		break;
	case TrackExpressionVars::Track_var::REG_MIN:
		var.var[idx] = track.last_min();
		break;
	case TrackExpressionVars::Track_var::REG_MAX:
		var.var[idx] = track.last_max();
		break;
	case TrackExpressionVars::Track_var::WEIGHTED_SUM:
		var.var[idx] = track.last_weighted_sum();
		break;
	case TrackExpressionVars::Track_var::OCCUPIED_AREA:
		var.var[idx] = track.last_occupied_area();
		break;
	default:
		verror("Internal error: unsupported function %d", var.val_func);
	}
}

// Filter aggregation helper methods
double TrackVarProcessor::aggregate_avg_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	long double total_weight = 0;
	long double total_weighted_sum = 0;

	for (const auto& part : parts) {
		track.read_interval(part);
		double part_avg = track.last_avg();
		if (!std::isnan(part_avg)) {
			double part_len = part.end - part.start;
			total_weighted_sum += part_avg * part_len;
			total_weight += part_len;
		}
	}

	return (total_weight > 0) ? (total_weighted_sum / total_weight) : numeric_limits<double>::quiet_NaN();
}

double TrackVarProcessor::aggregate_sum_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	double total_sum = 0;
	bool has_value = false;

	for (const auto& part : parts) {
		track.read_interval(part);
		double part_sum = track.last_sum();
		if (!std::isnan(part_sum)) {
			total_sum += part_sum;
			has_value = true;
		}
	}

	return has_value ? total_sum : numeric_limits<double>::quiet_NaN();
}

double TrackVarProcessor::aggregate_min_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	double min_val = numeric_limits<double>::infinity();

	for (const auto& part : parts) {
		track.read_interval(part);
		double part_min = track.last_min();
		if (!std::isnan(part_min)) {
			min_val = std::min(min_val, part_min);
		}
	}

	return (min_val != numeric_limits<double>::infinity()) ? min_val : numeric_limits<double>::quiet_NaN();
}

double TrackVarProcessor::aggregate_max_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	double max_val = -numeric_limits<double>::infinity();

	for (const auto& part : parts) {
		track.read_interval(part);
		double part_max = track.last_max();
		if (!std::isnan(part_max)) {
			max_val = std::max(max_val, part_max);
		}
	}

	return (max_val != -numeric_limits<double>::infinity()) ? max_val : numeric_limits<double>::quiet_NaN();
}

bool TrackVarProcessor::find_best_max_pos_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, double &best_pos)
{
	double best_val = -numeric_limits<double>::infinity();
	bool has_value = false;

	for (const auto &part : parts) {
		track.read_interval(part);
		double part_max = track.last_max();
		double part_pos = track.last_max_pos();

		if (std::isnan(part_max) || std::isnan(part_pos))
			continue;

		if (!has_value || part_max > best_val || (part_max == best_val && part_pos < best_pos)) {
			best_val = part_max;
			best_pos = part_pos;
			has_value = true;
		}
	}

	return has_value;
}

double TrackVarProcessor::aggregate_max_pos_abs_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	double best_pos = numeric_limits<double>::quiet_NaN();
	if (!find_best_max_pos_with_filter(track, parts, best_pos))
		return numeric_limits<double>::quiet_NaN();

	return best_pos;
}

double TrackVarProcessor::aggregate_max_pos_rel_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, int64_t base_start)
{
	double best_pos = numeric_limits<double>::quiet_NaN();
	if (!find_best_max_pos_with_filter(track, parts, best_pos))
		return numeric_limits<double>::quiet_NaN();

	return best_pos - base_start;
}

bool TrackVarProcessor::find_best_min_pos_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, double &best_pos)
{
	double best_val = numeric_limits<double>::infinity();
	bool has_value = false;

	for (const auto &part : parts) {
		track.read_interval(part);
		double part_min = track.last_min();
		double part_pos = track.last_min_pos();

		if (std::isnan(part_min) || std::isnan(part_pos))
			continue;

		if (!has_value || part_min < best_val || (part_min == best_val && part_pos < best_pos)) {
			best_val = part_min;
			best_pos = part_pos;
			has_value = true;
		}
	}

	return has_value;
}

double TrackVarProcessor::aggregate_min_pos_abs_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	double best_pos = numeric_limits<double>::quiet_NaN();
	if (!find_best_min_pos_with_filter(track, parts, best_pos))
		return numeric_limits<double>::quiet_NaN();

	return best_pos;
}

double TrackVarProcessor::aggregate_min_pos_rel_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, int64_t base_start)
{
	double best_pos = numeric_limits<double>::quiet_NaN();
	if (!find_best_min_pos_with_filter(track, parts, best_pos))
		return numeric_limits<double>::quiet_NaN();

	return best_pos - base_start;
}

double TrackVarProcessor::aggregate_stddev_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	long double total_weight = 0;
	long double total_weighted_sum = 0;
	long double M2 = 0;  // For variance calculation

	for (const auto& part : parts) {
		track.read_interval(part);
		double part_avg = track.last_avg();
		double part_stddev = track.last_stddev();
		double part_len = part.end - part.start;

		if (!std::isnan(part_avg) && !std::isnan(part_stddev)) {
			// Welford online algorithm for combining variances
			double delta = part_avg - (total_weight > 0 ? total_weighted_sum / total_weight : 0);
			total_weight += part_len;
			total_weighted_sum += part_avg * part_len;
			M2 += part_stddev * part_stddev * part_len + part_len * total_weight / (total_weight + part_len) * delta * delta;
		}
	}

	return (total_weight > 0) ? std::sqrt(M2 / total_weight) : numeric_limits<double>::quiet_NaN();
}

double TrackVarProcessor::aggregate_quantile_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, double percentile)
{
	// Create a combined percentiler and add all samples from all parts
	StreamPercentiler<float> combined_sp;
	combined_sp.init(m_iu.get_max_data_size(),
	                 m_iu.get_quantile_edge_data_size(),
	                 m_iu.get_quantile_edge_data_size());

	for (const auto& part : parts) {
		track.read_interval(part);

		const auto& sp = track.get_percentiler();
		const auto& samples = sp.samples();
		const auto& lowest = sp.lowest_vals();
		const auto& highest = sp.highest_vals();

		// Add all values to combined percentiler
		for (const auto& val : lowest) {
			if (!std::isnan(val)) {
				combined_sp.add(val, unif_rand);
			}
		}
		for (const auto& val : samples) {
			if (!std::isnan(val)) {
				combined_sp.add(val, unif_rand);
			}
		}
		for (const auto& val : highest) {
			if (!std::isnan(val)) {
				combined_sp.add(val, unif_rand);
			}
		}
	}

	if (combined_sp.stream_size() > 0) {
		bool is_estimated;
		return combined_sp.get_percentile(percentile, is_estimated);
	}

	return numeric_limits<double>::quiet_NaN();
}

double TrackVarProcessor::aggregate_exists_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	// Check if any value exists in any part
	for (const auto& part : parts) {
		track.read_interval(part);
		float exists = track.last_exists();
		if (exists == 1.0f) {
			return 1.0;
		}
	}
	return 0.0;
}

double TrackVarProcessor::aggregate_size_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	// Sum the sizes across all parts
	double total_size = 0.0;
	for (const auto& part : parts) {
		track.read_interval(part);
		total_size += track.last_size();
	}
	return total_size;
}

double TrackVarProcessor::aggregate_sample_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	// Collect all non-NaN values from all parts and sample one
	vector<float> all_values;
	for (const auto& part : parts) {
		track.read_interval(part);
		float val = track.last_sample();
		if (!std::isnan(val)) {
			all_values.push_back(val);
		}
	}

	if (all_values.empty()) {
		return numeric_limits<double>::quiet_NaN();
	}

	// Use R's random number generator for reproducibility
	GetRNGstate();
	int idx = (int)(unif_rand() * all_values.size());
	PutRNGstate();

	return all_values[idx];
}

double TrackVarProcessor::aggregate_sample_pos_abs_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	// Collect all non-NaN positions from all parts and sample one
	vector<double> all_positions;
	for (const auto& part : parts) {
		track.read_interval(part);
		double pos = track.last_sample_pos();
		if (!std::isnan(pos)) {
			all_positions.push_back(pos);
		}
	}

	if (all_positions.empty()) {
		return numeric_limits<double>::quiet_NaN();
	}

	// Use R's random number generator for reproducibility
	GetRNGstate();
	int idx = (int)(unif_rand() * all_positions.size());
	PutRNGstate();

	return all_positions[idx];
}

double TrackVarProcessor::aggregate_sample_pos_rel_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, int64_t base_start)
{
	double abs_pos = aggregate_sample_pos_abs_with_filter(track, parts);
	if (std::isnan(abs_pos)) {
		return abs_pos;
	}
	return abs_pos - base_start;
}

