/*
 * GenomeTrackInMemory.cpp
 *
 * In-memory track implementation with count-based statistics (matching GenomeTrackSparse)
 */

#include "GenomeTrackInMemory.h"
#include "TGLException.h"
#include <algorithm>

using namespace std;

GenomeTrackInMemory::GenomeTrackInMemory() :
	GenomeTrack1D(COMPUTED),
	m_last_min_pos(numeric_limits<double>::quiet_NaN())
{
	// COMPUTED type since it's not a file-based track
}

void GenomeTrackInMemory::init_from_data(
	const GIntervals &intervals,
	const vector<float> &vals,
	int chromid)
{
	if (intervals.size() != vals.size()) {
		TGLError<GenomeTrackInMemory>(
			"Intervals and values must have the same size: %zu vs %zu",
			intervals.size(), vals.size());
	}

	m_intervals = intervals;
	m_vals = vals;
	m_chromid = chromid;
	m_icur_interval = m_intervals.begin();
}

void GenomeTrackInMemory::read_interval(const GInterval &interval)
{
	// Initialize all results to NaN
	m_last_avg = m_last_nearest = m_last_min = m_last_max =
		m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();

	if (has_function(LSE))
		m_last_lse = numeric_limits<float>::quiet_NaN();
	if (has_function(MAX_POS))
		m_last_max_pos = numeric_limits<double>::quiet_NaN();
	if (has_function(MIN_POS))
		m_last_min_pos = numeric_limits<double>::quiet_NaN();
	if (has_function(EXISTS))
		m_last_exists = 0;
	if (has_function(SIZE))
		m_last_size = 0;
	if (has_function(SAMPLE))
		m_last_sample = numeric_limits<float>::quiet_NaN();
	if (has_function(SAMPLE_POS))
		m_last_sample_pos = numeric_limits<double>::quiet_NaN();
	if (has_function(FIRST))
		m_last_first = numeric_limits<float>::quiet_NaN();
	if (has_function(FIRST_POS))
		m_last_first_pos = numeric_limits<double>::quiet_NaN();
	if (has_function(LAST))
		m_last_last = numeric_limits<float>::quiet_NaN();
	if (has_function(LAST_POS))
		m_last_last_pos = numeric_limits<double>::quiet_NaN();

	if (m_use_quantile)
		m_sp.reset();

	if (m_intervals.empty())
		return;

	// Early exit checks (from GenomeTrackSparse)
	if (m_intervals.front().start >= interval.end) {
		m_last_nearest = m_vals.front();
		return;
	}

	if (m_intervals.back().end <= interval.start) {
		m_last_nearest = m_vals.back();
		return;
	}

	// Try sequential access first (optimization from GenomeTrackSparse)
	if (check_first_overlap(m_icur_interval, interval)) {
		calc_vals_coverage_weighted(interval);
	} else if (m_icur_interval + 1 < m_intervals.end() &&
	           check_first_overlap(m_icur_interval + 1, interval)) {
		++m_icur_interval;
		calc_vals_coverage_weighted(interval);
	} else {
		// Binary search (pattern from GenomeTrackSparse)
		GIntervals::const_iterator istart_interval = m_intervals.begin();
		GIntervals::const_iterator iend_interval = m_intervals.end();

		while (iend_interval - istart_interval > 1) {
			GIntervals::const_iterator imid_interval =
				istart_interval + (iend_interval - istart_interval) / 2;

			if (check_first_overlap(imid_interval, interval)) {
				m_icur_interval = imid_interval;
				calc_vals_coverage_weighted(interval);
				break;
			}

			if (GIntervals::compare_by_start_coord(*imid_interval, interval))
				istart_interval = imid_interval;
			else
				iend_interval = imid_interval;
		}

		if (iend_interval - istart_interval == 1 &&
		    check_first_overlap(istart_interval, interval)) {
			m_icur_interval = istart_interval;
			calc_vals_coverage_weighted(interval);
		}

		// Handle nearest value
		if (iend_interval - istart_interval == 1)
			m_last_nearest = iend_interval == m_intervals.end() ||
			                 interval.dist2interv(*istart_interval) <= interval.dist2interv(*iend_interval) ?
				m_vals[istart_interval - m_intervals.begin()] : m_vals[iend_interval - m_intervals.begin()];
	}
}

double GenomeTrackInMemory::last_max_pos() const
{
	return m_last_max_pos;
}

double GenomeTrackInMemory::last_min_pos() const
{
	return m_last_min_pos;
}

void GenomeTrackInMemory::calc_vals_coverage_weighted(const GInterval &interval)
{
	// Use count-based statistics (exactly like GenomeTrackSparse), NOT coverage-weighted
	uint64_t num_vs = 0;
	double stddev_mean = 0;
	double stddev_m2 = 0;
	double sum_accum = 0;
	float v;

	// For sampling
	vector<float> all_values;
	vector<double> all_positions;
	if (has_function(SAMPLE) || has_function(SAMPLE_POS))
		all_values.reserve(100);
	if (has_function(SAMPLE_POS))
		all_positions.reserve(100);

	m_last_min = numeric_limits<float>::max();
	m_last_max = -numeric_limits<float>::max();

	if (has_function(MAX_POS))
		m_last_max_pos = numeric_limits<double>::quiet_NaN();
	if (has_function(MIN_POS))
		m_last_min_pos = numeric_limits<double>::quiet_NaN();
	if (has_function(LSE))
		m_last_lse = -numeric_limits<float>::infinity();

	// Iterate through all overlapping intervals
	for (GIntervals::const_iterator iinterv = m_icur_interval;
	     iinterv != m_intervals.end();
	     ++iinterv) {

		if (!iinterv->do_overlap(interval))
			break;

		v = m_vals[iinterv - m_intervals.begin()];

		if (!std::isnan(v)) {
			sum_accum += v;

			// Min/Max tracking
			if (v < m_last_min) {
				m_last_min = v;
				if (has_function(MIN_POS))
					m_last_min_pos = iinterv->start;
			} else if (has_function(MIN_POS) && v == m_last_min) {
				if (std::isnan(m_last_min_pos) || iinterv->start < m_last_min_pos)
					m_last_min_pos = iinterv->start;
			}

			if (v > m_last_max) {
				m_last_max = v;
				if (has_function(MAX_POS))
					m_last_max_pos = iinterv->start;
			}

			++num_vs;
			if (has_function(STDDEV)) {
				const double delta = v - stddev_mean;
				stddev_mean += delta / static_cast<double>(num_vs);
				const double delta2 = v - stddev_mean;
				stddev_m2 += delta * delta2;
			}

			if (has_function(LSE))
				lse_accumulate(m_last_lse, v);

			if (m_use_quantile)
				m_sp.add(v, s_rnd_func);

			// Exists
			if (has_function(EXISTS))
				m_last_exists = 1;

			// First
			if (has_function(FIRST) && std::isnan(m_last_first))
				m_last_first = v;

			if (has_function(FIRST_POS) && std::isnan(m_last_first_pos))
				m_last_first_pos = iinterv->start;

			// Last
			if (has_function(LAST))
				m_last_last = v;

			if (has_function(LAST_POS))
				m_last_last_pos = iinterv->start;

			// Sampling
			if (has_function(SAMPLE))
				all_values.push_back(v);
			if (has_function(SAMPLE_POS))
				all_positions.push_back(iinterv->start);
		}
	}

	// Finalize size
	if (has_function(SIZE))
		m_last_size = num_vs;

	// Sample from collected values
	if (has_function(SAMPLE) && !all_values.empty()) {
		int idx = (int)(s_rnd_func() * all_values.size());
		if (idx >= (int)all_values.size()) idx = (int)all_values.size() - 1;
		if (idx < 0) idx = 0;
		m_last_sample = all_values[idx];
	}

	if (has_function(SAMPLE_POS) && !all_positions.empty()) {
		int idx = (int)(s_rnd_func() * all_positions.size());
		if (idx >= (int)all_positions.size()) idx = (int)all_positions.size() - 1;
		if (idx < 0) idx = 0;
		m_last_sample_pos = all_positions[idx];
	}

	m_last_sum = (float)sum_accum;
	if (num_vs > 0)
		m_last_avg = m_last_nearest = (float)(sum_accum / num_vs);
	else {
		m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_sum = numeric_limits<float>::quiet_NaN();
		if (has_function(LSE))
			m_last_lse = numeric_limits<float>::quiet_NaN();
		if (has_function(MIN_POS))
			m_last_min_pos = numeric_limits<double>::quiet_NaN();
	}

	// Handle min/max edge cases
	if (m_last_min == numeric_limits<float>::max())
		m_last_min = numeric_limits<float>::quiet_NaN();
	if (m_last_max == -numeric_limits<float>::max())
		m_last_max = numeric_limits<float>::quiet_NaN();

	// Unbiased sample standard deviation via Welford's stable algorithm.
	if (has_function(STDDEV))
		m_last_stddev = num_vs > 1 ? sqrt(stddev_m2 / static_cast<double>(num_vs - 1)) : numeric_limits<float>::quiet_NaN();
}
