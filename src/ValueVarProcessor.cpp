#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

#include "ValueVarProcessor.h"
#include "TrackExpressionVars.h"
#include "GenomeTrackInMemory.h"
#include "rdbutils.h"

using namespace rdb;
using namespace std;

void ValueVarProcessor::process_value_vars(
	TrackExpressionVars::Value_vars &value_vars,
	const GInterval &interval,
	unsigned idx)
{
	for (TrackExpressionVars::Value_vars::iterator ivar = value_vars.begin(); ivar != value_vars.end(); ++ivar) {
		process_single_value_var(*ivar, interval, idx);
	}
}

void ValueVarProcessor::process_single_value_var(
	TrackExpressionVars::Value_var &var,
	const GInterval &interval,
	unsigned idx)
{
	const GInterval &eval_interval = var.imdf1d ? var.imdf1d->interval : interval;

	if (var.imdf1d && var.imdf1d->out_of_range) {
		var.var[idx] = numeric_limits<double>::quiet_NaN();
		return;
	}

	// Apply filter if present
	if (var.filter) {
		std::vector<GInterval> eval_intervals;
		var.filter->subtract(eval_interval, eval_intervals);

		if (eval_intervals.empty()) {
			var.var[idx] = numeric_limits<double>::quiet_NaN();
			return;
		}

		// For filtered intervals, we need to combine results across multiple sub-intervals
		// This matches the pattern used for track variables with filters
		struct FilteredValueAggregate {
			double total_sum = 0;
			double total_sum_sq = 0;
			int64_t total_size = 0;
			double min_val = numeric_limits<double>::quiet_NaN();
			double max_val = numeric_limits<double>::quiet_NaN();
			double first_val = numeric_limits<double>::quiet_NaN();
			double last_val = numeric_limits<double>::quiet_NaN();
			double sample_val = numeric_limits<double>::quiet_NaN();
			double min_pos = numeric_limits<double>::quiet_NaN();
			double max_pos = numeric_limits<double>::quiet_NaN();
			double first_pos = numeric_limits<double>::quiet_NaN();
			double last_pos = numeric_limits<double>::quiet_NaN();
			double sample_pos = numeric_limits<double>::quiet_NaN();
			double nearest_val = numeric_limits<double>::quiet_NaN();
			double lse_val = -numeric_limits<double>::infinity();
			bool exists = false;
			bool has_value = false;

			void add_interval(GenomeTrackInMemory &track) {
				double part_avg = track.last_avg();
				double part_sum = track.last_sum();
				if (std::isnan(part_avg) || std::isnan(part_sum))
					return;

				double part_size = track.last_size();
				double part_stddev = track.last_stddev();
				double part_min = track.last_min();
				double part_max = track.last_max();
				double part_first = track.last_first();
				double part_last = track.last_last();
				double part_sample = track.last_sample();
				double part_min_pos = track.last_min_pos();
				double part_max_pos = track.last_max_pos();
				double part_first_pos = track.last_first_pos();
				double part_last_pos = track.last_last_pos();
				double part_sample_pos = track.last_sample_pos();
				double part_nearest = track.last_nearest();
				double part_exists = track.last_exists();

				total_sum += part_sum;
				total_size += static_cast<int64_t>(part_size);

				double part_lse = track.last_lse();
				if (!std::isnan(part_lse))
					lse_accumulate(lse_val, part_lse);

				if (part_size > 1 && !std::isnan(part_stddev)) {
					double var_term = part_stddev * part_stddev * (part_size - 1);
					total_sum_sq += var_term + part_avg * part_avg * part_size;
				} else
					total_sum_sq += part_avg * part_avg * part_size;

				if (!std::isnan(part_min)) {
					if (std::isnan(min_val) || part_min < min_val)
						min_val = part_min;
				}
				if (!std::isnan(part_max)) {
					if (std::isnan(max_val) || part_max > max_val)
						max_val = part_max;
				}

				if (!std::isnan(part_first) && std::isnan(first_val))
					first_val = part_first;
				if (!std::isnan(part_last))
					last_val = part_last;

				if (!std::isnan(part_sample) && std::isnan(sample_val))
					sample_val = part_sample;

				if (!std::isnan(part_min_pos)) {
					if (std::isnan(min_pos) || part_min_pos < min_pos)
						min_pos = part_min_pos;
				}
				if (!std::isnan(part_max_pos)) {
					if (std::isnan(max_pos) || part_max_pos > max_pos)
						max_pos = part_max_pos;
				}
				if (!std::isnan(part_first_pos) && std::isnan(first_pos))
					first_pos = part_first_pos;
				if (!std::isnan(part_last_pos))
					last_pos = part_last_pos;

				if (!std::isnan(part_sample_pos) && std::isnan(sample_pos))
					sample_pos = part_sample_pos;

				if (!std::isnan(part_nearest) && std::isnan(nearest_val))
					nearest_val = part_nearest;

				exists = exists || (!std::isnan(part_exists) && part_exists != 0);
				has_value = true;
			}
		} agg;

		for (const auto &eval_int : eval_intervals) {
			var.track->read_interval(eval_int);
			agg.add_interval(*var.track);
		}

		// Finalize result for filtered case
		if (!agg.has_value) {
			var.var[idx] = numeric_limits<double>::quiet_NaN();
		} else {
			switch (var.val_func) {
				case TrackExpressionVars::Value_var::AVG:
					var.var[idx] = agg.total_size > 0 ? agg.total_sum / agg.total_size : numeric_limits<double>::quiet_NaN();
					break;
				case TrackExpressionVars::Value_var::STDDEV:
					if (agg.total_size > 1) {
						double mean = agg.total_sum / agg.total_size;
						double variance = agg.total_sum_sq / (agg.total_size - 1) - mean * mean * ((double)agg.total_size / (agg.total_size - 1));
						var.var[idx] = sqrt(std::max(0.0, variance));
					} else
						var.var[idx] = numeric_limits<double>::quiet_NaN();
					break;
				case TrackExpressionVars::Value_var::SUM:
					var.var[idx] = agg.total_sum;
					break;
				case TrackExpressionVars::Value_var::LSE:
					var.var[idx] = std::isfinite(agg.lse_val) ? agg.lse_val : numeric_limits<double>::quiet_NaN();
					break;
				case TrackExpressionVars::Value_var::MIN:
					var.var[idx] = agg.min_val;
					break;
				case TrackExpressionVars::Value_var::MAX:
					var.var[idx] = agg.max_val;
					break;
				case TrackExpressionVars::Value_var::EXISTS:
					var.var[idx] = agg.exists ? 1.0 : 0.0;
					break;
				case TrackExpressionVars::Value_var::SIZE:
					var.var[idx] = agg.total_size;
					break;
				case TrackExpressionVars::Value_var::FIRST:
					var.var[idx] = agg.first_val;
					break;
				case TrackExpressionVars::Value_var::LAST:
					var.var[idx] = agg.last_val;
					break;
				case TrackExpressionVars::Value_var::SAMPLE:
					var.var[idx] = agg.sample_val;
					break;
				case TrackExpressionVars::Value_var::FIRST_POS_ABS:
					var.var[idx] = agg.first_pos;
					break;
				case TrackExpressionVars::Value_var::FIRST_POS_REL:
					var.var[idx] = std::isnan(agg.first_pos) ? numeric_limits<double>::quiet_NaN() : agg.first_pos - eval_interval.start;
					break;
				case TrackExpressionVars::Value_var::LAST_POS_ABS:
					var.var[idx] = agg.last_pos;
					break;
				case TrackExpressionVars::Value_var::LAST_POS_REL:
					var.var[idx] = std::isnan(agg.last_pos) ? numeric_limits<double>::quiet_NaN() : agg.last_pos - eval_interval.start;
					break;
				case TrackExpressionVars::Value_var::SAMPLE_POS_ABS:
					var.var[idx] = agg.sample_pos;
					break;
				case TrackExpressionVars::Value_var::SAMPLE_POS_REL:
					var.var[idx] = std::isnan(agg.sample_pos) ? numeric_limits<double>::quiet_NaN() : agg.sample_pos - eval_interval.start;
					break;
				case TrackExpressionVars::Value_var::MIN_POS_ABS:
					var.var[idx] = agg.min_pos;
					break;
				case TrackExpressionVars::Value_var::MIN_POS_REL:
					var.var[idx] = std::isnan(agg.min_pos) ? numeric_limits<double>::quiet_NaN() : agg.min_pos - eval_interval.start;
					break;
				case TrackExpressionVars::Value_var::MAX_POS_ABS:
					var.var[idx] = agg.max_pos;
					break;
				case TrackExpressionVars::Value_var::MAX_POS_REL:
					var.var[idx] = std::isnan(agg.max_pos) ? numeric_limits<double>::quiet_NaN() : agg.max_pos - eval_interval.start;
					break;
				case TrackExpressionVars::Value_var::NEAREST:
					var.var[idx] = agg.nearest_val;
					break;
				case TrackExpressionVars::Value_var::QUANTILE: {
					// For quantile with filter, just use first non-empty result
					float val = numeric_limits<float>::quiet_NaN();
					for (const auto &eval_int : eval_intervals) {
						var.track->read_interval(eval_int);
						val = var.track->last_quantile(var.percentile);
						if (!std::isnan(val))
							break;
					}
					var.var[idx] = val;
					break;
				}
				default:
					var.var[idx] = numeric_limits<double>::quiet_NaN();
					break;
			}
		}
	} else {
		// No filter - simple case, use track interface directly
		var.track->read_interval(eval_interval);

		// Extract result based on function
		switch (var.val_func) {
			case TrackExpressionVars::Value_var::AVG:
				var.var[idx] = var.track->last_avg();
				break;
			case TrackExpressionVars::Value_var::MIN:
				var.var[idx] = var.track->last_min();
				break;
			case TrackExpressionVars::Value_var::MAX:
				var.var[idx] = var.track->last_max();
				break;
			case TrackExpressionVars::Value_var::SUM:
				var.var[idx] = var.track->last_sum();
				break;
			case TrackExpressionVars::Value_var::LSE:
				var.var[idx] = var.track->last_lse();
				break;
			case TrackExpressionVars::Value_var::STDDEV:
				var.var[idx] = var.track->last_stddev();
				break;
			case TrackExpressionVars::Value_var::QUANTILE:
				var.var[idx] = var.track->last_quantile(var.percentile);
				break;
			case TrackExpressionVars::Value_var::NEAREST:
				var.var[idx] = var.track->last_nearest();
				break;
			case TrackExpressionVars::Value_var::EXISTS:
				var.var[idx] = var.track->last_exists();
				break;
			case TrackExpressionVars::Value_var::SIZE:
				var.var[idx] = var.track->last_size();
				break;
			case TrackExpressionVars::Value_var::FIRST:
				var.var[idx] = var.track->last_first();
				break;
			case TrackExpressionVars::Value_var::LAST:
				var.var[idx] = var.track->last_last();
				break;
			case TrackExpressionVars::Value_var::SAMPLE:
				var.var[idx] = var.track->last_sample();
				break;
			case TrackExpressionVars::Value_var::FIRST_POS_ABS:
				var.var[idx] = var.track->last_first_pos();
				break;
			case TrackExpressionVars::Value_var::FIRST_POS_REL:
				var.var[idx] = var.track->last_first_pos() - eval_interval.start;
				break;
			case TrackExpressionVars::Value_var::LAST_POS_ABS:
				var.var[idx] = var.track->last_last_pos();
				break;
			case TrackExpressionVars::Value_var::LAST_POS_REL:
				var.var[idx] = var.track->last_last_pos() - eval_interval.start;
				break;
			case TrackExpressionVars::Value_var::SAMPLE_POS_ABS:
				var.var[idx] = var.track->last_sample_pos();
				break;
			case TrackExpressionVars::Value_var::SAMPLE_POS_REL:
				var.var[idx] = var.track->last_sample_pos() - eval_interval.start;
				break;
			case TrackExpressionVars::Value_var::MIN_POS_ABS:
				var.var[idx] = var.track->last_min_pos();
				break;
			case TrackExpressionVars::Value_var::MIN_POS_REL:
				var.var[idx] = var.track->last_min_pos() - eval_interval.start;
				break;
			case TrackExpressionVars::Value_var::MAX_POS_ABS:
				var.var[idx] = var.track->last_max_pos();
				break;
			case TrackExpressionVars::Value_var::MAX_POS_REL:
				var.var[idx] = var.track->last_max_pos() - eval_interval.start;
				break;
			default:
				verror("Internal error: unsupported value variable function %d", var.val_func);
		}
	}
}

