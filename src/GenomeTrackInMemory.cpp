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

	if (m_functions[LSE])
		m_last_lse = numeric_limits<float>::quiet_NaN();
	if (m_functions[MAX_POS])
		m_last_max_pos = numeric_limits<double>::quiet_NaN();
	if (m_functions[MIN_POS])
		m_last_min_pos = numeric_limits<double>::quiet_NaN();
	if (m_functions[EXISTS])
		m_last_exists = 0;
	if (m_functions[SIZE])
		m_last_size = 0;
	if (m_functions[SAMPLE])
		m_last_sample = numeric_limits<float>::quiet_NaN();
	if (m_functions[SAMPLE_POS])
		m_last_sample_pos = numeric_limits<double>::quiet_NaN();
	if (m_functions[FIRST])
		m_last_first = numeric_limits<float>::quiet_NaN();
	if (m_functions[FIRST_POS])
		m_last_first_pos = numeric_limits<double>::quiet_NaN();
	if (m_functions[LAST])
		m_last_last = numeric_limits<float>::quiet_NaN();
	if (m_functions[LAST_POS])
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
	float num_vs = 0;
	double mean_square_sum = 0;
	float v;

	// For sampling
	vector<float> all_values;
	vector<double> all_positions;
	if (m_functions[SAMPLE] || m_functions[SAMPLE_POS])
		all_values.reserve(100);
	if (m_functions[SAMPLE_POS])
		all_positions.reserve(100);

	m_last_sum = 0;
	m_last_min = numeric_limits<float>::max();
	m_last_max = -numeric_limits<float>::max();

	if (m_functions[MAX_POS])
		m_last_max_pos = numeric_limits<double>::quiet_NaN();
	if (m_functions[MIN_POS])
		m_last_min_pos = numeric_limits<double>::quiet_NaN();
	if (m_functions[LSE])
		m_last_lse = -numeric_limits<float>::infinity();

	// Iterate through all overlapping intervals
	for (GIntervals::const_iterator iinterv = m_icur_interval;
	     iinterv != m_intervals.end();
	     ++iinterv) {

		if (!iinterv->do_overlap(interval))
			break;

		v = m_vals[iinterv - m_intervals.begin()];

		if (!std::isnan(v)) {
			m_last_sum += v;

			// Min/Max tracking
			if (v < m_last_min) {
				m_last_min = v;
				if (m_functions[MIN_POS])
					m_last_min_pos = iinterv->start;
			} else if (m_functions[MIN_POS] && v == m_last_min) {
				if (std::isnan(m_last_min_pos) || iinterv->start < m_last_min_pos)
					m_last_min_pos = iinterv->start;
			}

			if (v > m_last_max) {
				m_last_max = v;
				if (m_functions[MAX_POS])
					m_last_max_pos = iinterv->start;
			}

			if (m_functions[STDDEV])
				mean_square_sum += v * v;

			if (m_functions[LSE])
				lse_accumulate(m_last_lse, v);

			if (m_use_quantile)
				m_sp.add(v, s_rnd_func);

			// Exists
			if (m_functions[EXISTS])
				m_last_exists = 1;

			// First
			if (m_functions[FIRST] && std::isnan(m_last_first))
				m_last_first = v;

			if (m_functions[FIRST_POS] && std::isnan(m_last_first_pos))
				m_last_first_pos = iinterv->start;

			// Last
			if (m_functions[LAST])
				m_last_last = v;

			if (m_functions[LAST_POS])
				m_last_last_pos = iinterv->start;

			// Sampling
			if (m_functions[SAMPLE])
				all_values.push_back(v);
			if (m_functions[SAMPLE_POS])
				all_positions.push_back(iinterv->start);

			++num_vs;
		}
	}

	// Finalize size
	if (m_functions[SIZE])
		m_last_size = num_vs;

	// Sample from collected values
	if (m_functions[SAMPLE] && !all_values.empty()) {
		int idx = (int)(s_rnd_func() * all_values.size());
		m_last_sample = all_values[idx];
	}

	if (m_functions[SAMPLE_POS] && !all_positions.empty()) {
		int idx = (int)(s_rnd_func() * all_positions.size());
		m_last_sample_pos = all_positions[idx];
	}

	if (num_vs > 0)
		m_last_avg = m_last_nearest = m_last_sum / num_vs;
	else {
		m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_sum = numeric_limits<float>::quiet_NaN();
		if (m_functions[LSE])
			m_last_lse = numeric_limits<float>::quiet_NaN();
		if (m_functions[MIN_POS])
			m_last_min_pos = numeric_limits<double>::quiet_NaN();
	}

	// Handle min/max edge cases
	if (m_last_min == numeric_limits<float>::max())
		m_last_min = numeric_limits<float>::quiet_NaN();
	if (m_last_max == -numeric_limits<float>::max())
		m_last_max = numeric_limits<float>::quiet_NaN();

	// Calculate unbiased standard deviation (exactly like GenomeTrackSparse)
	// sqrt(sum((x-mean)^2) / (N-1)) = sqrt(sum(x^2)/(N-1) - N*(mean^2)/(N-1))
	if (m_functions[STDDEV])
		m_last_stddev = num_vs > 1 ? sqrt(mean_square_sum / (num_vs - 1) - (m_last_avg * (double)m_last_avg) * (num_vs / (num_vs - 1))) : numeric_limits<float>::quiet_NaN();
}
