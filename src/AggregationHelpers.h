#ifndef AGGREGATION_HELPERS_H
#define AGGREGATION_HELPERS_H

#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <cstring>
#include "rdbutils.h"

using namespace std;

enum class AggregationType {
	MEAN,
	MEDIAN,
	SUM,
	MIN,
	MAX,
	COUNT,
	FIRST,
	LAST,
	NTH,
	MAX_COV_LEN,
	MIN_COV_LEN,
	MAX_COV_FRAC,
	MIN_COV_FRAC,
	WEIGHTED_MEAN,
	WEIGHTED_SUM,
	WEIGHTED_MEDIAN,
	COVERAGE
};

struct Contribution {
	double value;
	double overlap_len;
	double coverage_frac;
	int64_t start;
	int64_t end;
	bool is_na;
	int64_t chain_id;
};

struct AggregationConfig {
	AggregationType type;
	bool na_rm;
	int min_n;      // -1 means disabled
	int nth_index;  // only used for NTH; 1-based; -1 means unset
};

struct AggregationState {
	vector<Contribution> contributions;
	bool has_na;

	AggregationState() : has_na(false) {}

	void reset() {
		contributions.clear();
		has_na = false;
	}
};

inline AggregationType parse_aggregation_type(const char *agg) {
	if (!strcmp(agg, "mean"))
		return AggregationType::MEAN;
	if (!strcmp(agg, "median"))
		return AggregationType::MEDIAN;
	if (!strcmp(agg, "sum"))
		return AggregationType::SUM;
	if (!strcmp(agg, "min"))
		return AggregationType::MIN;
	if (!strcmp(agg, "max"))
		return AggregationType::MAX;
	if (!strcmp(agg, "count"))
		return AggregationType::COUNT;
	if (!strcmp(agg, "first"))
		return AggregationType::FIRST;
	if (!strcmp(agg, "last"))
		return AggregationType::LAST;
	if (!strcmp(agg, "nth"))
		return AggregationType::NTH;
	if (!strcmp(agg, "max.coverage_len"))
		return AggregationType::MAX_COV_LEN;
	if (!strcmp(agg, "min.coverage_len"))
		return AggregationType::MIN_COV_LEN;
	if (!strcmp(agg, "max.coverage_frac"))
		return AggregationType::MAX_COV_FRAC;
	if (!strcmp(agg, "min.coverage_frac"))
		return AggregationType::MIN_COV_FRAC;
	if (!strcmp(agg, "weighted.mean"))
		return AggregationType::WEIGHTED_MEAN;
	if (!strcmp(agg, "weighted.sum"))
		return AggregationType::WEIGHTED_SUM;
	if (!strcmp(agg, "weighted.median"))
		return AggregationType::WEIGHTED_MEDIAN;
	if (!strcmp(agg, "coverage"))
		return AggregationType::COVERAGE;
	rdb::verror("Unknown multi_target_agg value: %s", agg);
	return AggregationType::MEAN;  // unreachable
}

inline void aggregation_state_add(AggregationState &state,
                                         double value,
                                         double overlap_len,
                                         double locus_len,
                                         int64_t overlap_start,
                                         int64_t overlap_end,
                                         int64_t chain_id)
{
	Contribution contrib;
	contrib.value = value;
	contrib.overlap_len = std::max(0.0, overlap_len);

	double denom = locus_len > 0.0 ? locus_len : 0.0;
	contrib.coverage_frac = (denom > 0.0) ? contrib.overlap_len / denom : 0.0;
	contrib.start = overlap_start;
	contrib.end = overlap_end;
	contrib.is_na = std::isnan(value);
	contrib.chain_id = chain_id;

	if (contrib.is_na)
		state.has_na = true;

	state.contributions.push_back(contrib);
}

inline double aggregate_values(const AggregationConfig &cfg, const AggregationState &state) {
	// Under na.rm = FALSE, any NaN contribution makes the whole locus NaN.
	if (!cfg.na_rm) {
		for (const auto &contrib : state.contributions) {
			if (contrib.is_na)
				return numeric_limits<double>::quiet_NaN();
		}
	}

	// Merge contributions by chain_id so a single chain (possibly split into several
	// pieces over this locus - parallel slices under "agg", or several source bins
	// mapping in) counts once. NaN pieces are dropped FIRST (na.rm) so that a chain
	// mapping both finite and NaN source bins into this locus keeps its finite value
	// instead of being discarded wholesale by an is_na carried over the whole chain.
	vector<Contribution> merged;
	merged.reserve(state.contributions.size());

	for (const auto &contrib : state.contributions) {
		if (contrib.is_na)
			continue; // na.rm == true here (the na.rm == false case already returned)
		bool found = false;
		for (auto &existing : merged) {
			if (existing.chain_id == contrib.chain_id) {
				existing.overlap_len += contrib.overlap_len;
				existing.coverage_frac += contrib.coverage_frac;
				existing.start = std::min(existing.start, contrib.start);
				existing.end = std::max(existing.end, contrib.end);
				found = true;
				break;
			}
		}
		if (!found)
			merged.push_back(contrib);
	}

	vector<const Contribution *> valid;
	valid.reserve(merged.size());
	for (const auto &contrib : merged)
		valid.push_back(&contrib);

	if (cfg.min_n >= 0 && (int)valid.size() < cfg.min_n)
		return numeric_limits<double>::quiet_NaN();

	if (cfg.type == AggregationType::COUNT) {
		if (valid.empty())
			return 0.0;
		return static_cast<double>(valid.size());
	}

	if (valid.empty())
		return numeric_limits<double>::quiet_NaN();

	switch (cfg.type) {
	case AggregationType::MEAN: {
		double sum = 0.0;
		for (const Contribution *c : valid)
			sum += c->value;
		return sum / static_cast<double>(valid.size());
	}
	case AggregationType::MEDIAN: {
		vector<double> values;
		values.reserve(valid.size());
		for (const Contribution *c : valid)
			values.push_back(c->value);
		std::sort(values.begin(), values.end());
		size_t mid = values.size() / 2;
		if (values.size() % 2 == 0)
			return (values[mid - 1] + values[mid]) / 2.0;
		return values[mid];
	}
	case AggregationType::SUM: {
		double sum = 0.0;
		for (const Contribution *c : valid)
			sum += c->value;
		return sum;
	}
	case AggregationType::MIN: {
		double result = valid.front()->value;
		for (size_t i = 1; i < valid.size(); ++i)
			result = std::min(result, valid[i]->value);
		return result;
	}
	case AggregationType::MAX: {
		double result = valid.front()->value;
		for (size_t i = 1; i < valid.size(); ++i)
			result = std::max(result, valid[i]->value);
		return result;
	}
	case AggregationType::FIRST: {
		const Contribution *first = *std::min_element(valid.begin(), valid.end(),
			                                               [](const Contribution *a, const Contribution *b) {
				if (a->start != b->start)
					return a->start < b->start;
				if (a->end != b->end)
					return a->end < b->end;
				return a->value > b->value;
		});
	return first->value;
}
case AggregationType::LAST: {
	const Contribution *last = *std::max_element(valid.begin(), valid.end(),
		                                              [](const Contribution *a, const Contribution *b) {
				if (a->start != b->start)
					return a->start < b->start;
				if (a->end != b->end)
					return a->end < b->end;
				return a->value > b->value;
	});
	return last->value;
}
case AggregationType::NTH: {
	std::sort(valid.begin(), valid.end(),
	          [](const Contribution *a, const Contribution *b) {
				  if (a->start != b->start)
					  return a->start < b->start;
				  if (a->end != b->end)
					  return a->end < b->end;
				  return a->value > b->value;
	          });
	// nth
	if (cfg.nth_index <= 0)
		return numeric_limits<double>::quiet_NaN();
	size_t idx = static_cast<size_t>(cfg.nth_index - 1);
	if (idx >= valid.size())
		return numeric_limits<double>::quiet_NaN();
	return valid[idx]->value;
}
	case AggregationType::MAX_COV_LEN:
	case AggregationType::MIN_COV_LEN:
	case AggregationType::MAX_COV_FRAC:
	case AggregationType::MIN_COV_FRAC: {
		const Contribution *best = valid.front();
		for (size_t i = 1; i < valid.size(); ++i) {
			const Contribution *curr = valid[i];
			bool better = false;
			switch (cfg.type) {
			case AggregationType::MAX_COV_LEN:
				if (curr->overlap_len > best->overlap_len)
					better = true;
				else if (curr->overlap_len == best->overlap_len && curr->value > best->value)
					better = true;
				break;
			case AggregationType::MIN_COV_LEN:
				if (curr->overlap_len < best->overlap_len)
					better = true;
				else if (curr->overlap_len == best->overlap_len && curr->value > best->value)
					better = true;
				break;
			case AggregationType::MAX_COV_FRAC:
				if (curr->coverage_frac > best->coverage_frac)
					better = true;
				else if (curr->coverage_frac == best->coverage_frac && curr->value > best->value)
					better = true;
				break;
			case AggregationType::MIN_COV_FRAC:
				if (curr->coverage_frac < best->coverage_frac)
					better = true;
				else if (curr->coverage_frac == best->coverage_frac && curr->value > best->value)
					better = true;
				break;
			default:
				break;
			}
			if (better)
				best = curr;
		}
		return best->value;
	}
	case AggregationType::WEIGHTED_MEAN: {
		double sum = 0.0;
		double w = 0.0;
		for (const Contribution *c : valid) {
			sum += c->value * c->overlap_len;
			w += c->overlap_len;
		}
		if (w <= 0.0)
			return numeric_limits<double>::quiet_NaN();
		return sum / w;
	}
	case AggregationType::WEIGHTED_SUM: {
		double sum = 0.0;
		for (const Contribution *c : valid)
			sum += c->value * c->overlap_len;
		return sum;
	}
	case AggregationType::COVERAGE: {
		// sum(v_i * overlap_i / locus_len). For ChIP-seq pileup with v=1:
		// bin value = average reads-per-base over the bin.
		double sum = 0.0;
		for (const Contribution *c : valid)
			sum += c->value * c->coverage_frac;
		return sum;
	}
	case AggregationType::WEIGHTED_MEDIAN: {
		// Lower weighted median: sort by value asc, accumulate overlap_len;
		// return the value of the first contribution whose running cumulative
		// overlap reaches or exceeds total_overlap / 2.
		vector<const Contribution *> sorted(valid.begin(), valid.end());
		std::sort(sorted.begin(), sorted.end(),
			[](const Contribution *a, const Contribution *b) {
				return a->value < b->value;
			});
		double total = 0.0;
		for (const Contribution *c : sorted)
			total += c->overlap_len;
		if (total <= 0.0)
			return numeric_limits<double>::quiet_NaN();
		double half = total / 2.0;
		double acc = 0.0;
		for (const Contribution *c : sorted) {
			acc += c->overlap_len;
			if (acc >= half)
				return c->value;
		}
		return sorted.back()->value;  // unreachable; defensive
	}
	case AggregationType::COUNT:  // handled earlier
		break;
	}

	return numeric_limits<double>::quiet_NaN();
}

#endif // AGGREGATION_HELPERS_H
