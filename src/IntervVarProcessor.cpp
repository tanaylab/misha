#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>
#include <unordered_set>
#include <algorithm>

#include "IntervVarProcessor.h"
#include "TrackExpressionVars.h"
#include "GIntervals.h"
#include "rdbutils.h"

using namespace rdb;
using namespace std;

void IntervVarProcessor::process_interv_vars(
	TrackExpressionVars::Interv_vars &interv_vars,
	const GInterval &interval,
	unsigned idx)
{
	for (TrackExpressionVars::Interv_vars::iterator ivar = interv_vars.begin(); ivar != interv_vars.end(); ++ivar) {
		switch (ivar->val_func) {
		case TrackExpressionVars::Interv_var::DIST:
			process_distance(*ivar, interval, idx);
			break;
		case TrackExpressionVars::Interv_var::DIST_CENTER:
			process_distance_center(*ivar, interval, idx);
			break;
		case TrackExpressionVars::Interv_var::COVERAGE:
			process_coverage(*ivar, interval, idx);
			break;
		case TrackExpressionVars::Interv_var::NEIGHBOR_COUNT:
			process_neighbor_count(*ivar, interval, idx);
			break;
		default:
			verror("Internal error: unsupported interval function %d", ivar->val_func);
		}
	}
}

void IntervVarProcessor::process_distance(
	TrackExpressionVars::Interv_var &var,
	const GInterval &interval,
	unsigned idx)
{
	// if iterator modifier exists, iterator intervals might not come sorted => perform a binary search
	if (var.imdf1d) {
		const GInterval &eval_interval = var.imdf1d->interval;
		double min_dist = numeric_limits<double>::max();
		double dist;
		int64_t coord = (eval_interval.start + eval_interval.end) / 2;
		GIntervals::const_iterator iinterv = lower_bound(var.sintervs.begin(), var.sintervs.end(), eval_interval, GIntervals::compare_by_start_coord);

		if (iinterv != var.sintervs.end() && iinterv->chromid == eval_interval.chromid)
			min_dist = iinterv->dist2coord(coord, var.dist_margin);

		if (iinterv != var.sintervs.begin() && (iinterv - 1)->chromid == eval_interval.chromid) {
			dist = (iinterv - 1)->dist2coord(coord, var.dist_margin);
			if (fabs(min_dist) > fabs(dist))
				min_dist = dist;
		}

		// if min_dist == double_max then we haven't found an interval with the same chromosome as the iterator interval =>
		// we can skip the second binary search
		if (min_dist == numeric_limits<double>::max())
			var.var[idx] = numeric_limits<double>::quiet_NaN();
		else {
			iinterv = lower_bound(var.eintervs.begin(), var.eintervs.end(), eval_interval, GIntervals::compare_by_end_coord);

			if (iinterv != var.eintervs.end() && iinterv->chromid == eval_interval.chromid) {
				dist = iinterv->dist2coord(coord, var.dist_margin);
				if (fabs(min_dist) > fabs(dist))
					min_dist = dist;
			}

			if (iinterv != var.eintervs.begin() && (iinterv - 1)->chromid == eval_interval.chromid) {
				dist = (iinterv - 1)->dist2coord(coord, var.dist_margin);
				if (fabs(min_dist) > fabs(dist))
					min_dist = dist;
			}

			var.var[idx] = min_dist;
		}
	} else {
		const GIntervals *pintervs[2] = { &var.sintervs, &var.eintervs };
		GIntervals::const_iterator *piinterv[2] = { &var.siinterv, &var.eiinterv };
		double dist[2] = { 0, 0 };

		for (int i = 0; i < 2; ++i) {
			const GIntervals &intervs = *pintervs[i];
			GIntervals::const_iterator &iinterv = *piinterv[i];

			while (iinterv != intervs.end() && iinterv->chromid < interval.chromid)
				++iinterv;

			if (iinterv == intervs.end() || iinterv->chromid != interval.chromid)
				dist[i] = numeric_limits<double>::quiet_NaN();
			else {
				int64_t coord = (interval.start + interval.end) / 2;
				dist[i] = (double)iinterv->dist2coord(coord, var.dist_margin);
				GIntervals::const_iterator iinterv_next = iinterv + 1;

				while (iinterv_next != intervs.end() && iinterv_next->chromid == interval.chromid) {
					double dist_next = iinterv_next->dist2coord(coord, var.dist_margin);

					if (fabs(dist[i]) < fabs(dist_next))
						break;

					iinterv = iinterv_next;
					dist[i] = dist_next;
					++iinterv_next;
				}
			}
		}

		var.var[idx] = fabs(dist[0]) < fabs(dist[1]) ? dist[0] : dist[1];
	}
}

void IntervVarProcessor::process_distance_center(
	TrackExpressionVars::Interv_var &var,
	const GInterval &interval,
	unsigned idx)
{
	// if iterator modifier exists, iterator intervals might not come sorted => perform a binary search
	if (var.imdf1d) {
		int64_t coord = (var.imdf1d->interval.start + var.imdf1d->interval.end) / 2;
		GInterval eval_interval(var.imdf1d->interval.chromid, coord, coord + 1, 0);
		GIntervals::const_iterator iinterv = lower_bound(var.sintervs.begin(), var.sintervs.end(), eval_interval, GIntervals::compare_by_start_coord);
		double dist = numeric_limits<double>::quiet_NaN();

		var.var[idx] = numeric_limits<double>::quiet_NaN();

		if (iinterv != var.sintervs.end() && iinterv->chromid == eval_interval.chromid)
			dist = iinterv->dist2center(coord);

		if (dist != numeric_limits<double>::quiet_NaN() && iinterv != var.sintervs.begin() && (iinterv - 1)->chromid == eval_interval.chromid)
			dist = (iinterv - 1)->dist2center(coord);

		var.var[idx] = dist;
	} else {
		int64_t coord = (interval.start + interval.end) / 2;
		GIntervals::const_iterator &iinterv = var.siinterv;
		double dist = numeric_limits<double>::quiet_NaN();

		while (iinterv != var.sintervs.end() && var.siinterv->chromid < interval.chromid)
			++iinterv;

		while (iinterv != var.sintervs.end() && iinterv->chromid == interval.chromid && iinterv->start <= coord) {
			if (iinterv->end > coord)
				dist = iinterv->dist2center(coord);
			++iinterv;
		}

		var.var[idx] = dist;
	}
}

void IntervVarProcessor::process_neighbor_count(
	TrackExpressionVars::Interv_var &var,
	const GInterval &interval,
	unsigned idx)
{
	const GInterval &eval_interval = var.imdf1d ? var.imdf1d->interval : interval;

	if (var.imdf1d && var.imdf1d->out_of_range) {
		var.var[idx] = 0;
		return;
	}

	std::vector<GInterval> eval_intervals;
	if (var.filter) {
		var.filter->subtract(eval_interval, eval_intervals);
		if (eval_intervals.empty()) {
			var.var[idx] = numeric_limits<double>::quiet_NaN();
			return;
		}
	} else {
		eval_intervals.push_back(eval_interval);
	}

	size_t neighbor_count = 0;
	if (!var.imdf1d && !var.filter) {
		GIntervals::const_iterator &eiter = var.eiinterv;
		const GIntervals &expanded = var.eintervs;

		while (eiter != expanded.end() && (eiter->chromid < eval_interval.chromid ||
			   (eiter->chromid == eval_interval.chromid && eiter->end <= eval_interval.start)))
			++eiter;

		GIntervals::const_iterator scan = eiter;
		while (scan != expanded.end() && scan->chromid == eval_interval.chromid && scan->start < eval_interval.end) {
			if (scan->end > eval_interval.start)
				++neighbor_count;
			++scan;
		}
	} else {
		const GIntervals &expanded = var.eintervs;
		std::unordered_set<size_t> counted;
		counted.reserve(eval_intervals.size() * 2);

		for (const auto &eval_int : eval_intervals) {
			auto it = lower_bound(expanded.begin(), expanded.end(), eval_int, GIntervals::compare_by_start_coord);

			// Walk backward to the first interval on this chromosome
			// Since intervals are only sorted by start (not end), we can't make assumptions
			// about whether earlier intervals might have large spans that overlap the query
			while (it != expanded.begin()) {
				auto prev = it - 1;
				if (prev->chromid != eval_int.chromid)
					break;
				--it;
			}

			// Scan forward checking all intervals on this chromosome for overlaps
			for (; it != expanded.end() && it->chromid == eval_int.chromid; ++it) {
				if (it->end <= eval_int.start)
					continue;
				if (it->start >= eval_int.end)
					break;

				size_t expanded_idx = it - expanded.begin();
				if (counted.insert(expanded_idx).second)
					++neighbor_count;
			}
		}
	}

	var.var[idx] = static_cast<double>(neighbor_count);
}

void IntervVarProcessor::process_coverage(
	TrackExpressionVars::Interv_var &var,
	const GInterval &interval,
	unsigned idx)
{
	const GInterval &eval_interval = var.imdf1d ? var.imdf1d->interval : interval;

	if (var.imdf1d && var.imdf1d->out_of_range) {
		var.var[idx] = 0;
		return;
	}

	// Apply filter if present
	std::vector<GInterval> eval_intervals;
	if (var.filter) {
		var.filter->subtract(eval_interval, eval_intervals);
		if (eval_intervals.empty()) {
			var.var[idx] = numeric_limits<double>::quiet_NaN();
			return;
		}
	} else {
		eval_intervals.push_back(eval_interval);
	}

	int64_t total_overlap = 0;
	int64_t total_unmasked_length = 0;

	// Calculate coverage for each unmasked sub-interval
	for (const auto &eval_int : eval_intervals) {
		total_unmasked_length += eval_int.range();
		int64_t sub_overlap = 0;
		GIntervals::const_iterator iinterv;

		// For non-sequential access or first access
		if (var.imdf1d || var.siinterv == var.sintervs.end()) {
			iinterv = lower_bound(var.sintervs.begin(), var.sintervs.end(), eval_int,
								  GIntervals::compare_by_start_coord);

			// Check previous interval too
			if (iinterv != var.sintervs.begin()) {
				auto prev = iinterv - 1;
				if (prev->chromid == eval_int.chromid && prev->end > eval_int.start) {
					int64_t overlap_start = max(eval_int.start, prev->start);
					int64_t overlap_end = min(eval_int.end, prev->end);
					sub_overlap += overlap_end - overlap_start;
				}
			}
		} else {
			// For sequential access, start from last position
			iinterv = var.siinterv;

			// Skip past intervals from previous chromosomes
			while (iinterv != var.sintervs.end() && iinterv->chromid < eval_int.chromid)
				++iinterv;

			// But check if we need to back up
			while (iinterv != var.sintervs.begin() &&
				   (iinterv - 1)->chromid == eval_int.chromid &&
				   (iinterv - 1)->end > eval_int.start) {
				--iinterv;
			}
		}

		// Check forward intervals
		while (iinterv != var.sintervs.end() &&
			   iinterv->chromid == eval_int.chromid &&
			   iinterv->start < eval_int.end) {
			if (iinterv->end > eval_int.start) {
				int64_t overlap_start = max(eval_int.start, iinterv->start);
				int64_t overlap_end = min(eval_int.end, iinterv->end);
				sub_overlap += overlap_end - overlap_start;
			}
			++iinterv;
		}

		if (!var.imdf1d) {
			var.siinterv = iinterv;
		}

		total_overlap += sub_overlap;
	}

	var.var[idx] = (double)total_overlap / total_unmasked_length;
}

