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
		case TrackExpressionVars::Interv_var::DIST_EDGE:
			process_distance_edge(*ivar, interval, idx);
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

namespace {
// Build (or rebuild on chromosome change) the per-chromosome nearest-neighbor
// index for `var`, populated from its start-sorted source intervals. Returns the
// finder ready to query the given chromosome. An empty chromosome yields a finder
// with no objects, so NNIterator::begin returns false and callers emit NaN.
void ensure_finder(TrackExpressionVars::Interv_var &var, int chromid, const GenomeChromKey &chromkey)
{
	if (var.finder && var.finder_chromid == chromid)
		return;

	if (!var.finder)
		var.finder = std::make_unique<SegmentFinder<GInterval>>();

	var.finder->reset(0, (int64_t)chromkey.get_chrom_size(chromid));
	for (var.sintervs.begin_chrom_iter(chromid); !var.sintervs.isend_chrom(); var.sintervs.next_in_chrom())
		var.finder->insert(var.sintervs.cur_interval());

	var.finder_chromid = chromid;
}
}

void IntervVarProcessor::process_distance(
	TrackExpressionVars::Interv_var &var,
	const GInterval &interval,
	unsigned idx)
{
	const GInterval &q = var.imdf1d ? var.imdf1d->interval : interval;

	if (var.imdf1d && var.imdf1d->out_of_range) {
		var.var[idx] = numeric_limits<double>::quiet_NaN();
		return;
	}

	ensure_finder(var, q.chromid, m_iu.get_chromkey());

	int64_t coord = (q.start + q.end) / 2;
	SegmentFinder<GInterval>::NNIterator nn(var.finder.get());

	// No source intervals on this chromosome => NaN.
	if (!nn.begin(Segment(coord, coord + 1))) {
		var.var[idx] = numeric_limits<double>::quiet_NaN();
		return;
	}

	// The nearest interval to the query-bin center, measured at its edges
	// (dist2coord returns the signed distance, 0 / normalized fraction when inside).
	var.var[idx] = nn->dist2coord(coord, var.dist_margin);
}

void IntervVarProcessor::process_distance_center(
	TrackExpressionVars::Interv_var &var,
	const GInterval &interval,
	unsigned idx)
{
	const GInterval &q = var.imdf1d ? var.imdf1d->interval : interval;

	if (var.imdf1d && var.imdf1d->out_of_range) {
		var.var[idx] = numeric_limits<double>::quiet_NaN();
		return;
	}

	ensure_finder(var, q.chromid, m_iu.get_chromkey());

	int64_t coord = (q.start + q.end) / 2;
	double best = numeric_limits<double>::quiet_NaN();

	// distance.center is defined only when the query-bin center sits inside a
	// source interval. NNIterator yields candidates by increasing proximity, so
	// every interval containing the point precedes the first one that doesn't;
	// among the containing intervals we keep the nearest center. With
	// non-overlapping sources there is at most one such interval, preserving the
	// historical behavior; overlapping sources now resolve to the nearest center
	// instead of erroring.
	SegmentFinder<GInterval>::NNIterator nn(var.finder.get());
	for (bool ok = nn.begin(Segment(coord, coord + 1)); ok; ok = nn.next()) {
		const GInterval &iv = *nn;
		if (coord < iv.start || coord >= iv.end)
			break;
		double d = iv.dist2center(coord);
		if (std::isnan(best) || fabs(d) < fabs(best))
			best = d;
	}

	var.var[idx] = best;
}

void IntervVarProcessor::process_distance_edge(
	TrackExpressionVars::Interv_var &var,
	const GInterval &interval,
	unsigned idx)
{
	const GInterval &q = var.imdf1d ? var.imdf1d->interval : interval;

	if (var.imdf1d && var.imdf1d->out_of_range) {
		var.var[idx] = numeric_limits<double>::quiet_NaN();
		return;
	}

	ensure_finder(var, q.chromid, m_iu.get_chromkey());

	SegmentFinder<GInterval>::NNIterator nn(var.finder.get());

	// No source intervals on this chromosome => NaN.
	if (!nn.begin(Segment(q.start, q.end))) {
		var.var[idx] = numeric_limits<double>::quiet_NaN();
		return;
	}

	// Edge-to-edge distance to the nearest source interval, signed by the source
	// strand and 0 on overlap - identical to gintervals.neighbors, since both
	// consume the same NNIterator over start-sorted intervals.
	var.var[idx] = (double)q.dist2interv(*nn);
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

		// Scan backward to handle non-monotone access (overlapping regions)
		while (eiter != expanded.begin()) {
			GIntervals::const_iterator prev = eiter - 1;
			if (prev->chromid != eval_interval.chromid)
				break;
			if (prev->end > eval_interval.start)
				--eiter;
			else
				break;
		}

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

