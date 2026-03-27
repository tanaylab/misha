#include <errno.h>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "TGLException.h"
#include "GenomeTrackFixedBin.h"
#include "TrackIndex.h"

namespace {
enum ReducerBits : uint32_t {
	RBIT_SUM = 1u << 0,
	RBIT_LSE = 1u << 1,
	RBIT_EXISTS = 1u << 2,
	RBIT_SIZE = 1u << 3
};
}

void GenomeTrackFixedBin::sync_master_state_from_dependent()
{
	if (!m_master_obj)
		return;

	bool changed = false;
	uint32_t new_bits = m_func_mask & ~m_master_obj->m_func_mask;
	if (new_bits) {
		m_master_obj->m_func_mask |= new_bits;
		changed = true;
	}

	if (m_use_quantile && !m_master_obj->m_use_quantile) {
		m_master_obj->m_use_quantile = true;
		m_master_obj->m_sp = m_sp;
		changed = true;
	}

	if (changed) {
		// Recompute path classification on the next call with the merged function set.
		m_master_obj->m_fast_path_mode = 0;
		m_master_obj->m_fast_reducer_bits = 0;
	}
}

void GenomeTrackFixedBin::copy_state_from_master()
{
	if (!m_master_obj)
		return;

	if (has_function(AVG))
		m_last_avg = m_master_obj->m_last_avg;
	if (has_function(MIN))
		m_last_min = m_master_obj->m_last_min;
	if (has_function(MAX))
		m_last_max = m_master_obj->m_last_max;
	if (has_function(MAX_POS))
		m_last_max_pos = m_master_obj->m_last_max_pos;
	if (has_function(MIN_POS))
		m_last_min_pos = m_master_obj->m_last_min_pos;
	if (has_function(NEAREST))
		m_last_nearest = m_master_obj->m_last_nearest;
	if (has_function(STDDEV))
		m_last_stddev = m_master_obj->m_last_stddev;
	if (has_function(SUM))
		m_last_sum = m_master_obj->m_last_sum;
	if (has_function(LSE))
		m_last_lse = m_master_obj->m_last_lse;
	if (has_function(EXISTS))
		m_last_exists = m_master_obj->m_last_exists;
	if (has_function(SIZE))
		m_last_size = m_master_obj->m_last_size;
	if (has_function(SAMPLE))
		m_last_sample = m_master_obj->m_last_sample;
	if (has_function(SAMPLE_POS))
		m_last_sample_pos = m_master_obj->m_last_sample_pos;
	if (has_function(FIRST))
		m_last_first = m_master_obj->m_last_first;
	if (has_function(FIRST_POS))
		m_last_first_pos = m_master_obj->m_last_first_pos;
	if (has_function(LAST))
		m_last_last = m_master_obj->m_last_last;
	if (has_function(LAST_POS))
		m_last_last_pos = m_master_obj->m_last_last_pos;
	if (m_use_quantile)
		m_sp = m_master_obj->m_sp;
}

void GenomeTrackFixedBin::classify_fast_path_mode()
{
	// Single-function fast path: if exactly one function is registered and no quantile,
	// we can run a tight loop for just that function.
	if (!m_use_quantile && __builtin_popcount(m_func_mask) == 1) {
		m_fast_path_mode = 3;
		m_single_func = static_cast<Functions>(__builtin_ctz(m_func_mask));
		return;
	}

	bool reducer_only = !m_use_quantile;
	bool avg_nearest_only = !m_use_quantile;
	bool has_avg_nearest = false;
	uint32_t bits = 0;

	for (int i = 0; i < NUM_FUNCS; ++i) {
		if (!(m_func_mask & (1u << i)))
			continue;

		const int func = i;
		if (func == AVG || func == NEAREST)
			has_avg_nearest = true;
		else
			avg_nearest_only = false;

		if (!reducer_only)
			continue;

		switch (func) {
			case SUM: bits |= RBIT_SUM; break;
			case LSE: bits |= RBIT_LSE; break;
			case EXISTS: bits |= RBIT_EXISTS; break;
			case SIZE: bits |= RBIT_SIZE; break;
			default:
				reducer_only = false;
				break;
		}
	}

	if (reducer_only && bits != 0)
		m_fast_path_mode = 1;
	else if (avg_nearest_only && has_avg_nearest)
		m_fast_path_mode = 2;
	else
		m_fast_path_mode = -1;
	m_fast_reducer_bits = bits;
}

void GenomeTrackFixedBin::assign_single_bin_value(float value, double overlap_start)
{
	m_last_avg = value;
	m_last_min = m_last_max = m_last_nearest = m_last_sum = m_last_avg;
	m_last_stddev = numeric_limits<float>::quiet_NaN();

	if (has_function(LSE))
		m_last_lse = m_last_avg;
	if (has_function(MAX_POS))
		m_last_max_pos = overlap_start;
	if (has_function(MIN_POS))
		m_last_min_pos = overlap_start;
	if (has_function(EXISTS))
		m_last_exists = 1;
	if (has_function(SIZE))
		m_last_size = 1;
	if (has_function(SAMPLE))
		m_last_sample = m_last_avg;
	if (has_function(SAMPLE_POS))
		m_last_sample_pos = overlap_start;
	if (has_function(FIRST))
		m_last_first = m_last_avg;
	if (has_function(FIRST_POS))
		m_last_first_pos = overlap_start;
	if (has_function(LAST))
		m_last_last = m_last_avg;
	if (has_function(LAST_POS))
		m_last_last_pos = overlap_start;
	if (m_use_quantile && !std::isnan(m_last_avg))
		m_sp.add(m_last_avg, s_rnd_func);
}

void GenomeTrackFixedBin::assign_single_bin_missing()
{
	m_last_min = m_last_max = m_last_nearest = m_last_avg = m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();
	if (has_function(LSE))
		m_last_lse = numeric_limits<float>::quiet_NaN();
	if (has_function(MAX_POS))
		m_last_max_pos = numeric_limits<double>::quiet_NaN();
	if (has_function(MIN_POS))
		m_last_min_pos = numeric_limits<double>::quiet_NaN();
}

void GenomeTrackFixedBin::reset_sliding_window_state()
{
	m_lse_sliding_valid = false;
	m_sliding_sum = 0;
	m_sliding_sum_comp = 0;
	m_sliding_num_vs = 0;
	m_running_lse_initialized = false;
}

void GenomeTrackFixedBin::read_interval_reducers_only(const GInterval &interval)
{
	const bool need_sum = (m_fast_reducer_bits & RBIT_SUM) != 0;
	const bool need_lse = (m_fast_reducer_bits & RBIT_LSE) != 0;
	const bool need_exists = (m_fast_reducer_bits & RBIT_EXISTS) != 0;
	const bool need_size = (m_fast_reducer_bits & RBIT_SIZE) != 0;
	if (!need_lse)
		m_running_lse_initialized = false;

	auto assign_from_state = [&]() {
		if (m_sliding_num_vs > 0) {
			if (need_sum)
				m_last_sum = (float)m_sliding_sum;
			if (need_lse)
				m_last_lse = m_running_lse.window.empty() ? numeric_limits<float>::quiet_NaN() : (float)m_running_lse.value();
			if (need_exists)
				m_last_exists = 1;
		} else {
			if (need_sum)
				m_last_sum = numeric_limits<float>::quiet_NaN();
			if (need_lse)
				m_last_lse = numeric_limits<float>::quiet_NaN();
			if (need_exists)
				m_last_exists = 0;
		}
		if (need_size)
			m_last_size = m_sliding_num_vs;
	};

	auto assign_single_value = [&](bool have_value, float v) {
		bool has_num = have_value && !std::isnan(v);
		if (need_sum)
			m_last_sum = has_num ? v : numeric_limits<float>::quiet_NaN();
		if (need_lse)
			m_last_lse = has_num ? v : numeric_limits<float>::quiet_NaN();
		if (need_exists)
			m_last_exists = has_num ? 1 : 0;
		if (need_size)
			m_last_size = has_num ? 1 : 0;
	};

	// Common case: iterator advances exactly one dense bin.
	if (interval.start == m_cur_coord && interval.end == m_cur_coord + m_bin_size) {
		float v = numeric_limits<float>::quiet_NaN();
		bool have_value = read_next_bin(v);
		if (have_value) {
			m_cached_bin_idx = (int64_t)(interval.start / m_bin_size);
			m_cached_bin_val = v;
			m_cache_valid = true;
		}
		assign_single_value(have_value, v);
		m_lse_sliding_valid = false;
		m_sliding_sum = 0;
		m_sliding_sum_comp = 0;
		m_sliding_num_vs = 0;
		m_running_lse_initialized = false;
		return;
	}

	int64_t sbin = (int64_t)(interval.start / m_bin_size);
	int64_t ebin = (int64_t)ceil(interval.end / (double)m_bin_size);

	if (ebin == sbin + 1) {
		float v = numeric_limits<float>::quiet_NaN();
		bool have_value = false;

		if (m_cache_valid && m_cached_bin_idx == sbin) {
			v = m_cached_bin_val;
			m_cur_coord = (sbin + 1) * m_bin_size;
			have_value = true;
		} else {
			if (m_cur_coord != sbin * m_bin_size)
				goto_bin(sbin);
			if (read_next_bin(v)) {
				have_value = true;
				m_cached_bin_idx = sbin;
				m_cached_bin_val = v;
				m_cache_valid = true;
			}
		}

		assign_single_value(have_value, v);
		m_lse_sliding_valid = false;
		m_sliding_sum = 0;
		m_sliding_sum_comp = 0;
		m_sliding_num_vs = 0;
		m_running_lse_initialized = false;
		return;
	}

	const int64_t window_size = ebin - sbin;

	// Sliding reducers update when window shifts forward by <= previous window.
	if (m_lse_sliding_valid && window_size > 0 && (!need_lse || m_running_lse_initialized)) {
		int64_t step = sbin - m_lse_prev_sbin;
		int64_t prev_window = m_lse_prev_ebin - m_lse_prev_sbin;
		if (step > 0 && step <= prev_window && window_size == prev_window &&
			(int64_t)m_lse_window_bins.size() == prev_window) {
			int64_t appended = 0;

			if (step == 1) {
				float new_val = numeric_limits<float>::quiet_NaN();
				if (m_cur_coord != m_lse_prev_ebin * m_bin_size)
					goto_bin(m_lse_prev_ebin);
				if (read_next_bin(new_val) && !m_lse_window_bins.empty()) {
					float old_val = m_lse_window_bins.front();
					m_lse_window_bins.pop_front();
					if (!std::isnan(old_val)) {
						kahan_sub_from_sliding_sum(old_val);
						--m_sliding_num_vs;
						if (need_lse)
							m_running_lse.pop_front();
					}

					m_lse_window_bins.push_back(new_val);
					if (!std::isnan(new_val)) {
						kahan_add_to_sliding_sum(new_val);
						++m_sliding_num_vs;
						if (need_lse)
							m_running_lse.push(new_val);
					}
					appended = 1;
				}
			} else {
				vector<float> new_vals;
				appended = read_bins_bulk(m_lse_prev_ebin, step, new_vals);
				if (appended == step) {
					for (int64_t i = 0; i < step && !m_lse_window_bins.empty(); ++i) {
						float old_val = m_lse_window_bins.front();
						m_lse_window_bins.pop_front();
						if (!std::isnan(old_val)) {
							kahan_sub_from_sliding_sum(old_val);
							--m_sliding_num_vs;
							if (need_lse)
								m_running_lse.pop_front();
						}
					}

					for (int64_t i = 0; i < step; ++i) {
						float new_val = new_vals[i];
						m_lse_window_bins.push_back(new_val);
						if (!std::isnan(new_val)) {
							kahan_add_to_sliding_sum(new_val);
							++m_sliding_num_vs;
							if (need_lse)
								m_running_lse.push(new_val);
						}
					}
				}
			}

			if (appended == step) {
				assign_from_state();
				m_lse_prev_sbin = sbin;
				m_lse_prev_ebin = ebin;
				m_lse_sliding_valid = true;
				m_cached_bin_idx = ebin - 1;
				m_cached_bin_val = m_lse_window_bins.back();
				m_cache_valid = true;
				return;
			}
		}
	}

	// Fallback: full window read.
	vector<float> bin_vals;
	int64_t bins_read = read_bins_bulk(sbin, window_size, bin_vals);

	if (need_lse)
		m_running_lse.clear();
	m_sliding_sum = 0;
	m_sliding_sum_comp = 0;
	m_sliding_num_vs = 0;
	m_lse_window_bins.clear();
	for (int64_t i = 0; i < bins_read; ++i) {
		float v = bin_vals[i];
		m_lse_window_bins.push_back(v);
		if (!std::isnan(v)) {
			kahan_add_to_sliding_sum(v);
			++m_sliding_num_vs;
			if (need_lse)
				m_running_lse.push(v);
		}
	}
	m_running_lse_initialized = need_lse;

	if (bins_read > 0) {
		m_cached_bin_idx = sbin + bins_read - 1;
		m_cached_bin_val = bin_vals[bins_read - 1];
		m_cache_valid = true;
	}

	assign_from_state();
	m_lse_prev_sbin = sbin;
	m_lse_prev_ebin = ebin;
	m_lse_sliding_valid = bins_read == window_size && window_size > 0;
}

void GenomeTrackFixedBin::read_interval_avg_nearest_only(const GInterval &interval)
{
	if (interval.start == m_cur_coord && interval.end == m_cur_coord + m_bin_size) {
		if (read_next_bin(m_last_avg)) {
			m_last_min = m_last_max = m_last_nearest = m_last_sum = m_last_avg;
			m_last_stddev = numeric_limits<float>::quiet_NaN();
		} else {
			m_last_min = m_last_max = m_last_nearest = m_last_avg = m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();
		}
		return;
	}

	int64_t sbin = (int64_t)(interval.start / m_bin_size);
	int64_t ebin = (int64_t)ceil(interval.end / (double)m_bin_size);

	if (ebin == sbin + 1) {
		goto_bin(sbin);
		if (read_next_bin(m_last_avg)) {
			m_last_min = m_last_max = m_last_nearest = m_last_sum = m_last_avg;
			m_last_stddev = numeric_limits<float>::quiet_NaN();
		} else {
			m_last_min = m_last_max = m_last_nearest = m_last_avg = m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();
		}
	} else {
		uint64_t num_vs = 0;
		double stddev_mean = 0;
		double stddev_m2 = 0;
		float v;

		double sum_accum = 0;
		m_last_min = numeric_limits<float>::max();
		m_last_max = -numeric_limits<float>::max();

		goto_bin(sbin);
		for (int64_t bin = sbin; bin < ebin; ++bin) {
			if (read_next_bin(v) && !std::isnan(v)) {
				sum_accum += v;
				m_last_min = min(m_last_min, v);
				m_last_max = max(m_last_max, v);
				++num_vs;
				const double delta = v - stddev_mean;
				stddev_mean += delta / static_cast<double>(num_vs);
				const double delta2 = v - stddev_mean;
				stddev_m2 += delta * delta2;
			}
		}

		m_last_sum = (float)sum_accum;
		if (num_vs > 0) {
			m_last_avg = m_last_nearest = (float)(sum_accum / num_vs);
			m_last_stddev = num_vs > 1 ? sqrt(stddev_m2 / static_cast<double>(num_vs - 1)) : numeric_limits<float>::quiet_NaN();
		} else {
			m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_sum = numeric_limits<float>::quiet_NaN();
			m_last_stddev = numeric_limits<float>::quiet_NaN();
		}
	}
}

void GenomeTrackFixedBin::read_interval_single_function(const GInterval &interval)
{
	int64_t sbin = (int64_t)(interval.start / m_bin_size);
	int64_t ebin = (int64_t)ceil(interval.end / (double)m_bin_size);

	// Read all bins in the interval
	vector<float> bin_vals;
	int64_t bins_read = read_bins_bulk(sbin, ebin - sbin, bin_vals);

	if (bins_read > 0) {
		m_cached_bin_idx = sbin + bins_read - 1;
		m_cached_bin_val = bin_vals[bins_read - 1];
		m_cache_valid = true;
	}

	switch (m_single_func) {
	case AVG:
	case NEAREST: {
		uint64_t num_vs = 0;
		double sum_accum = 0;
		for (int64_t i = 0; i < bins_read; ++i) {
			if (!std::isnan(bin_vals[i])) {
				sum_accum += bin_vals[i];
				++num_vs;
			}
		}
		if (num_vs > 0) {
			m_last_avg = m_last_nearest = (float)(sum_accum / num_vs);
			m_last_sum = (float)sum_accum;
		} else {
			m_last_avg = m_last_nearest = m_last_sum = numeric_limits<float>::quiet_NaN();
		}
		m_last_min = m_last_max = m_last_stddev = numeric_limits<float>::quiet_NaN();
		break;
	}
	case SUM: {
		double sum_accum = 0;
		uint64_t num_vs = 0;
		for (int64_t i = 0; i < bins_read; ++i) {
			if (!std::isnan(bin_vals[i])) {
				sum_accum += bin_vals[i];
				++num_vs;
			}
		}
		m_last_sum = num_vs > 0 ? (float)sum_accum : numeric_limits<float>::quiet_NaN();
		m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_stddev = numeric_limits<float>::quiet_NaN();
		break;
	}
	case LSE: {
		m_last_lse = -numeric_limits<float>::infinity();
		uint64_t num_vs = 0;
		for (int64_t i = 0; i < bins_read; ++i) {
			if (!std::isnan(bin_vals[i])) {
				lse_accumulate(m_last_lse, bin_vals[i]);
				++num_vs;
			}
		}
		if (num_vs == 0)
			m_last_lse = numeric_limits<float>::quiet_NaN();
		m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_sum = m_last_stddev = numeric_limits<float>::quiet_NaN();
		break;
	}
	case MIN: {
		m_last_min = numeric_limits<float>::max();
		uint64_t num_vs = 0;
		for (int64_t i = 0; i < bins_read; ++i) {
			if (!std::isnan(bin_vals[i])) {
				if (bin_vals[i] < m_last_min)
					m_last_min = bin_vals[i];
				++num_vs;
			}
		}
		if (num_vs == 0)
			m_last_min = numeric_limits<float>::quiet_NaN();
		// Also compute avg/sum since the generic path does
		m_last_avg = m_last_nearest = m_last_max = m_last_sum = m_last_stddev = numeric_limits<float>::quiet_NaN();
		break;
	}
	case MAX: {
		m_last_max = -numeric_limits<float>::max();
		uint64_t num_vs = 0;
		for (int64_t i = 0; i < bins_read; ++i) {
			if (!std::isnan(bin_vals[i])) {
				if (bin_vals[i] > m_last_max)
					m_last_max = bin_vals[i];
				++num_vs;
			}
		}
		if (num_vs == 0)
			m_last_max = numeric_limits<float>::quiet_NaN();
		m_last_avg = m_last_nearest = m_last_min = m_last_sum = m_last_stddev = numeric_limits<float>::quiet_NaN();
		break;
	}
	case EXISTS: {
		m_last_exists = 0;
		for (int64_t i = 0; i < bins_read; ++i) {
			if (!std::isnan(bin_vals[i])) {
				m_last_exists = 1;
				break;
			}
		}
		m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_sum = m_last_stddev = numeric_limits<float>::quiet_NaN();
		break;
	}
	case SIZE: {
		uint64_t num_vs = 0;
		for (int64_t i = 0; i < bins_read; ++i) {
			if (!std::isnan(bin_vals[i]))
				++num_vs;
		}
		m_last_size = num_vs;
		m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_sum = m_last_stddev = numeric_limits<float>::quiet_NaN();
		break;
	}
	case STDDEV: {
		uint64_t num_vs = 0;
		double stddev_mean = 0;
		double stddev_m2 = 0;
		for (int64_t i = 0; i < bins_read; ++i) {
			if (!std::isnan(bin_vals[i])) {
				++num_vs;
				const double delta = bin_vals[i] - stddev_mean;
				stddev_mean += delta / static_cast<double>(num_vs);
				const double delta2 = bin_vals[i] - stddev_mean;
				stddev_m2 += delta * delta2;
			}
		}
		m_last_stddev = num_vs > 1 ? sqrt(stddev_m2 / static_cast<double>(num_vs - 1)) : numeric_limits<float>::quiet_NaN();
		m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_sum = numeric_limits<float>::quiet_NaN();
		break;
	}
	default:
		// For less common functions (MAX_POS, MIN_POS, SAMPLE, SAMPLE_POS, FIRST, FIRST_POS,
		// LAST, LAST_POS), fall back to the generic path. This should not happen since
		// classify_fast_path_mode only sets mode 3 for popcount==1 and these functions
		// are handled, but as a safety net:
		m_fast_path_mode = -1;
		read_interval(interval);
		return;
	}

	// Reset sliding window state since single-function path doesn't maintain it
	m_lse_sliding_valid = false;
	m_running_lse_initialized = false;
}

void GenomeTrackFixedBin::read_interval(const GInterval &interval)
{
	if (m_master_obj) {
		if (!m_master_synced) {
			sync_master_state_from_dependent();
			m_master_synced = true;
		}
		m_master_obj->read_interval(interval);
		copy_state_from_master();
		return;
	}

	if (m_fast_path_mode == 0)
		classify_fast_path_mode();

	if (m_fast_path_mode == 1) {
		read_interval_reducers_only(interval);
		return;
	}

	if (m_fast_path_mode == 2) {
		read_interval_avg_nearest_only(interval);
		return;
	}

	if (m_fast_path_mode == 3) {
		read_interval_single_function(interval);
		return;
	}

	if (m_use_quantile)
		m_sp.reset();
	if (!has_function(LSE))
		m_running_lse_initialized = false;

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

	// optimization of the most common case when the expression iterator starts at 0 and steps by bin_size
	if (interval.start == m_cur_coord && interval.end == m_cur_coord + m_bin_size) {
		if (read_next_bin(m_last_avg)) {
			m_cached_bin_idx = (int64_t)(interval.start / m_bin_size);
			m_cached_bin_val = m_last_avg;
			m_cache_valid = true;
			assign_single_bin_value(m_last_avg, interval.start);
		} else {
			assign_single_bin_missing();
		}
		reset_sliding_window_state();
		return;
	}

	int64_t sbin = (int64_t)(interval.start / m_bin_size);
	int64_t ebin = (int64_t)ceil(interval.end / (double)m_bin_size);

	const bool single_bin = ebin == sbin + 1;
	float cached_val = numeric_limits<float>::quiet_NaN();
	bool use_cache = false;
	bool have_value = false;

	if (single_bin && m_cache_valid && m_cached_bin_idx == sbin) {
		cached_val = m_cached_bin_val;
		use_cache = true;
		have_value = true;
	}

	if (single_bin) {
		if (!use_cache) {
			if (m_cur_coord != sbin * m_bin_size)
				goto_bin(sbin);
			if (read_next_bin(m_last_avg)) {
				m_cached_bin_idx = sbin;
				m_cached_bin_val = m_last_avg;
				m_cache_valid = true;
				have_value = true;
			}
		} else {
			m_last_avg = cached_val;
			// Keep virtual cursor at the end of this bin to match read_next_bin behaviour
			m_cur_coord = (sbin + 1) * m_bin_size;
		}

		if (have_value) {
			double overlap_start = std::max(static_cast<double>(sbin * m_bin_size), static_cast<double>(interval.start));
			assign_single_bin_value(m_last_avg, overlap_start);
		} else {
			assign_single_bin_missing();
		}
		reset_sliding_window_state();
	} else {
		uint64_t num_vs = 0;
		double stddev_mean = 0;
		double stddev_m2 = 0;
		double sum_accum = 0;
		const int64_t window_size = ebin - sbin;
		const bool simple_sliding_compatible =
			!m_use_quantile &&
			!has_function(MIN) &&
			!has_function(MAX) &&
			!has_function(MAX_POS) &&
			!has_function(MIN_POS) &&
			!has_function(STDDEV) &&
			!has_function(SAMPLE) &&
			!has_function(SAMPLE_POS) &&
			!has_function(FIRST) &&
			!has_function(FIRST_POS) &&
			!has_function(LAST) &&
			!has_function(LAST_POS);

		// Reuse scratch buffers for sampling (avoids per-call allocation)
		m_scratch_all_values.clear();
		m_scratch_all_positions.clear();

		const bool need_min = has_function(MIN) || has_function(MIN_POS);
		const bool need_max = has_function(MAX) || has_function(MAX_POS);
		m_last_min = need_min ? numeric_limits<float>::max() : numeric_limits<float>::quiet_NaN();
		m_last_max = need_max ? -numeric_limits<float>::max() : numeric_limits<float>::quiet_NaN();
		if (has_function(MAX_POS))
			m_last_max_pos = numeric_limits<double>::quiet_NaN();
		if (has_function(MIN_POS))
			m_last_min_pos = numeric_limits<double>::quiet_NaN();
		if (has_function(LSE))
			m_last_lse = -numeric_limits<float>::infinity();

		bool simple_sliding_used = false;

		if (simple_sliding_compatible && m_lse_sliding_valid && window_size > 0 &&
			(!has_function(LSE) || m_running_lse_initialized)) {
			int64_t step = sbin - m_lse_prev_sbin;
			int64_t prev_window = m_lse_prev_ebin - m_lse_prev_sbin;

				if (step > 0 && step <= prev_window && window_size == prev_window &&
					(int64_t)m_lse_window_bins.size() == prev_window) {
					int64_t appended = 0;
					if (step == 1) {
						float new_val = numeric_limits<float>::quiet_NaN();
					if (m_cur_coord != m_lse_prev_ebin * m_bin_size)
						goto_bin(m_lse_prev_ebin);
					appended = read_next_bin(new_val) ? 1 : 0;
					if (appended == 1 && !m_lse_window_bins.empty()) {
						float old_val = m_lse_window_bins.front();
						m_lse_window_bins.pop_front();
						if (!std::isnan(old_val)) {
							kahan_sub_from_sliding_sum(old_val);
							--m_sliding_num_vs;
							if (has_function(LSE))
								m_running_lse.pop_front();
						}

						m_lse_window_bins.push_back(new_val);
							if (!std::isnan(new_val)) {
								kahan_add_to_sliding_sum(new_val);
								++m_sliding_num_vs;
								if (has_function(LSE))
									m_running_lse.push(new_val);
							}
						}
					} else {
						vector<float> new_vals;
						appended = read_bins_bulk(m_lse_prev_ebin, step, new_vals);
						if (appended == step) {
							for (int64_t i = 0; i < step && !m_lse_window_bins.empty(); ++i) {
								float old_val = m_lse_window_bins.front();
								m_lse_window_bins.pop_front();
								if (!std::isnan(old_val)) {
									kahan_sub_from_sliding_sum(old_val);
									--m_sliding_num_vs;
									if (has_function(LSE))
										m_running_lse.pop_front();
								}
							}
							for (int64_t i = 0; i < appended; ++i) {
								float new_val = new_vals[i];
								m_lse_window_bins.push_back(new_val);
							if (!std::isnan(new_val)) {
								kahan_add_to_sliding_sum(new_val);
								++m_sliding_num_vs;
								if (has_function(LSE))
									m_running_lse.push(new_val);
								}
							}
						}
					}

				if (appended == step) {
					if (!m_lse_window_bins.empty()) {
						m_cached_bin_idx = ebin - 1;
						m_cached_bin_val = m_lse_window_bins.back();
						m_cache_valid = true;
					}

					if (m_sliding_num_vs > 0) {
						m_last_sum = (float)m_sliding_sum;
						if (has_function(AVG) || has_function(NEAREST))
							m_last_avg = m_last_nearest = (float)(m_sliding_sum / m_sliding_num_vs);
						if (has_function(LSE))
							m_last_lse = m_running_lse.window.empty()
								? numeric_limits<float>::quiet_NaN()
								: (float)m_running_lse.value();
						if (has_function(EXISTS))
							m_last_exists = 1;
					} else {
						m_last_sum = numeric_limits<float>::quiet_NaN();
						if (has_function(AVG) || has_function(NEAREST))
							m_last_avg = m_last_nearest = numeric_limits<float>::quiet_NaN();
						if (has_function(LSE))
							m_last_lse = numeric_limits<float>::quiet_NaN();
						if (has_function(EXISTS))
							m_last_exists = 0;
					}
					if (has_function(SIZE))
						m_last_size = m_sliding_num_vs;

					m_lse_prev_sbin = sbin;
					m_lse_prev_ebin = ebin;
					m_lse_sliding_valid = true;
					simple_sliding_used = true;
				}
			}
		}

		if (!simple_sliding_used) {
			// Bulk read all bins at once instead of one-by-one
			vector<float> bin_vals;
			int64_t bins_read = read_bins_bulk(sbin, window_size, bin_vals);
			bool lse_sliding_used = false;

			if (has_function(LSE) && m_lse_sliding_valid && bins_read > 0 &&
				m_running_lse_initialized && !simple_sliding_compatible) {
				int64_t step = sbin - m_lse_prev_sbin;
				int64_t prev_window = m_lse_prev_ebin - m_lse_prev_sbin;

				if (step > 0 && step <= prev_window && window_size == prev_window &&
					bins_read == window_size) {
					// Pop 'step' old bins from front
					for (int64_t i = 0; i < step && !m_lse_window_bins.empty(); i++) {
						float old_val = m_lse_window_bins.front();
						m_lse_window_bins.pop_front();
						if (!std::isnan(old_val))
							m_running_lse.pop_front();
					}
					// Push 'step' new bins at back
					for (int64_t i = bins_read - step; i < bins_read; i++) {
						float new_val = bin_vals[i];
						m_lse_window_bins.push_back(new_val);
						if (!std::isnan(new_val))
							m_running_lse.push(new_val);
					}
					// Extract result
					m_last_lse = m_running_lse.window.empty()
						? numeric_limits<float>::quiet_NaN()
						: (float)m_running_lse.value();

					lse_sliding_used = true;
				}
			}

			for (int64_t i = 0; i < bins_read; ++i) {
				int64_t bin = sbin + i;
				float v = bin_vals[i];

				m_cached_bin_idx = bin;
				m_cached_bin_val = v;
				m_cache_valid = true;

				if (!std::isnan(v)) {
					sum_accum += v;
					double bin_start = static_cast<double>(bin * m_bin_size);
					double overlap_start = std::max(bin_start, static_cast<double>(interval.start));
					if (need_min) {
						if (v < m_last_min) {
							m_last_min = v;
							if (has_function(MIN_POS))
								m_last_min_pos = overlap_start;
						} else if (has_function(MIN_POS) && v == m_last_min) {
							double candidate_pos = overlap_start;
							if (std::isnan(m_last_min_pos) || candidate_pos < m_last_min_pos)
								m_last_min_pos = candidate_pos;
						}
					}
					if (need_max) {
						if (v > m_last_max) {
							m_last_max = v;
							if (has_function(MAX_POS))
								m_last_max_pos = overlap_start;
						}
					}

					if (has_function(STDDEV)) {
						const double delta = v - stddev_mean;
						stddev_mean += delta / static_cast<double>(num_vs + 1);
						const double delta2 = v - stddev_mean;
						stddev_m2 += delta * delta2;
					}

					if (has_function(LSE) && !lse_sliding_used)
						lse_accumulate(m_last_lse, v);

					if (m_use_quantile && !std::isnan(v))
						m_sp.add(v, s_rnd_func);

					// New virtual track computations
					if (has_function(EXISTS))
						m_last_exists = 1;

					if (has_function(FIRST) && std::isnan(m_last_first))
						m_last_first = v;

					if (has_function(FIRST_POS) && std::isnan(m_last_first_pos))
						m_last_first_pos = overlap_start;

					if (has_function(LAST))
						m_last_last = v;

					if (has_function(LAST_POS))
						m_last_last_pos = overlap_start;

					if (has_function(SAMPLE))
						m_scratch_all_values.push_back(v);
					if (has_function(SAMPLE_POS))
						m_scratch_all_positions.push_back(overlap_start);

					++num_vs;
				}
			}

			// Finalize size
			if (has_function(SIZE))
				m_last_size = num_vs;

			// Sample from collected values
			if (has_function(SAMPLE) && !m_scratch_all_values.empty()) {
				int idx = (int)(s_rnd_func() * m_scratch_all_values.size());
				if (idx >= (int)m_scratch_all_values.size()) idx = (int)m_scratch_all_values.size() - 1;
				if (idx < 0) idx = 0;
				m_last_sample = m_scratch_all_values[idx];
			}

			if (has_function(SAMPLE_POS) && !m_scratch_all_positions.empty()) {
				int idx = (int)(s_rnd_func() * m_scratch_all_positions.size());
				if (idx >= (int)m_scratch_all_positions.size()) idx = (int)m_scratch_all_positions.size() - 1;
				if (idx < 0) idx = 0;
				m_last_sample_pos = m_scratch_all_positions[idx];
			}

			m_last_sum = (float)sum_accum;
			if (num_vs > 0)
				m_last_avg = m_last_nearest = (float)(sum_accum / num_vs);
			else {
				m_last_avg = m_last_nearest = m_last_sum = numeric_limits<float>::quiet_NaN();
				if (need_min)
					m_last_min = numeric_limits<float>::quiet_NaN();
				if (need_max)
					m_last_max = numeric_limits<float>::quiet_NaN();
				if (has_function(LSE))
					m_last_lse = numeric_limits<float>::quiet_NaN();
				if (has_function(MIN_POS))
					m_last_min_pos = numeric_limits<double>::quiet_NaN();
			}

			// Unbiased sample standard deviation via Welford's stable algorithm.
			if (has_function(STDDEV))
				m_last_stddev = num_vs > 1 ? sqrt(stddev_m2 / static_cast<double>(num_vs - 1)) : numeric_limits<float>::quiet_NaN();

			if (has_function(LSE) || simple_sliding_compatible) {
				if (!lse_sliding_used) {
					// Full computation was done; initialize sliding state from bin_vals
					m_running_lse.clear();
					m_lse_window_bins.clear();
					m_sliding_sum = 0;
					m_sliding_sum_comp = 0;
					m_sliding_num_vs = 0;
					for (int64_t i = 0; i < bins_read; i++) {
						m_lse_window_bins.push_back(bin_vals[i]);
						if (!std::isnan(bin_vals[i])) {
							kahan_add_to_sliding_sum(bin_vals[i]);
							++m_sliding_num_vs;
							if (has_function(LSE))
								m_running_lse.push(bin_vals[i]);
						}
					}
					m_running_lse_initialized = has_function(LSE);
				}
				m_lse_prev_sbin = sbin;
				m_lse_prev_ebin = ebin;
				m_lse_sliding_valid = bins_read == window_size && window_size > 0;
			} else {
				m_lse_sliding_valid = false;
				m_running_lse_initialized = false;
			}
			}
		}
	}

double GenomeTrackFixedBin::last_max_pos() const
{
	return m_last_max_pos;
}

double GenomeTrackFixedBin::last_min_pos() const
{
	return m_last_min_pos;
}

int64_t GenomeTrackFixedBin::read_bins_bulk(int64_t start_bin, int64_t num_bins, std::vector<float> &vals)
{
	if (num_bins <= 0) {
		vals.clear();
		return 0;
	}

	// Clamp to available samples
	int64_t available = m_num_samples - start_bin;
	if (available <= 0) {
		vals.clear();
		return 0;
	}
	int64_t to_read = std::min(num_bins, available);

	vals.resize(to_read);

	if (m_mmap_data) {
		// mmap path: direct memcpy from mapped region
		memcpy(vals.data(), m_mmap_data + start_bin, to_read * sizeof(float));
		m_cur_bin = start_bin + to_read;
	} else {
		if (m_cur_coord != start_bin * m_bin_size)
			goto_bin(start_bin);

		// Bulk read all bins in one syscall
		size_t bytes_to_read = to_read * sizeof(float);
		uint64_t bytes_read = m_bfile.read(vals.data(), bytes_to_read);

		if (bytes_read != bytes_to_read) {
			if (m_bfile.error())
				TGLError<GenomeTrackFixedBin>("Failed to read a dense track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
			to_read = bytes_read / sizeof(float);
			vals.resize(to_read);
		}
	}

	// Convert infinity to NaN (matching read_next_bin behavior)
	for (int64_t i = 0; i < to_read; ++i) {
		if (std::isinf(vals[i]))
			vals[i] = numeric_limits<float>::quiet_NaN();
	}

	// Update cursor position
	m_cur_coord = (start_bin + to_read) * m_bin_size;

	return to_read;
}

void GenomeTrackFixedBin::read_header_at_current_pos_(BufferedFile &bf)
{
	int32_t signature = 0;
	if (bf.read(&signature, sizeof(signature)) != sizeof(signature) || signature <= 0)
		TGLError<GenomeTrackFixedBin>("Invalid fixed-bin header in %s", bf.file_name().c_str());
	if (bf.read(&m_bin_size, sizeof(m_bin_size)) != sizeof(m_bin_size))
		TGLError<GenomeTrackFixedBin>("Invalid fixed-bin header in %s", bf.file_name().c_str());
}

void GenomeTrackFixedBin::init_read(const char *filename, const char *mode, int chromid)
{
	m_base_offset = 0; // Reset for per-chromosome
	m_cur_coord = 0;
	uint64_t header_start = 0;
	uint64_t total_bytes = 0;
	m_cached_bin_idx = -1;
	m_cached_bin_val = numeric_limits<float>::quiet_NaN();
	m_cache_valid = false;
	m_sliding_sum = 0;
	m_sliding_sum_comp = 0;
	m_sliding_num_vs = 0;
	m_lse_sliding_valid = false;
	m_running_lse_initialized = false;

	// Check for indexed format FIRST
	const std::string track_dir = GenomeTrack::get_track_dir(filename);
	const std::string idx_path = track_dir + "/track.idx";

	struct stat idx_st;
	if (stat(idx_path.c_str(), &idx_st) == 0) {
		// --- INDEXED PATH ---
		const std::string dat_path  = track_dir + "/track.dat";

		// Reopen file if: not open, path changed, or mode changed
		if (!m_dat_open || m_dat_path != dat_path || m_dat_mode != mode) {
			m_bfile.close();
			if (m_bfile.open(dat_path.c_str(), mode))
				TGLError<GenomeTrackFixedBin>("Cannot open %s: %s", dat_path.c_str(), strerror(errno));
			m_dat_open = true;
			m_dat_path = dat_path;
			m_dat_mode = mode;
		}

		auto idx   = get_track_index(track_dir);
		if (!idx)
			TGLError<GenomeTrackFixedBin>("Failed to load track index for %s", track_dir.c_str());

		auto entry = idx->get_entry(chromid);
		if (!entry || entry->length == 0) {
			// Chromosome not in index or empty contig - treat as empty
			m_num_samples = 0;
			m_chromid = chromid;
			return;
		}

		if (m_bfile.seek(entry->offset, SEEK_SET))
			TGLError<GenomeTrackFixedBin>("Failed to seek to offset %llu in %s",
				(unsigned long long)entry->offset, dat_path.c_str());

		header_start = entry->offset;
		// For indexed format, read just bin_size (no signature)
		// The data was copied as-is from per-chromosome files which have: bin_size + values
		if (m_bfile.read(&m_bin_size, sizeof(m_bin_size)) != sizeof(m_bin_size))
			TGLError<GenomeTrackFixedBin>("Invalid fixed-bin header in %s", dat_path.c_str());

		m_base_offset = entry->offset; 
		total_bytes = entry->length;
	} else {
		// --- PER-CHROMOSOME PATH ---
		m_bfile.close();
		m_dat_open = false;

		if (m_bfile.open(filename, mode))
			TGLError<GenomeTrackFixedBin>("%s", strerror(errno));

		if (m_bfile.read(&m_bin_size, sizeof(m_bin_size)) != sizeof(m_bin_size)) {
			if (m_bfile.error())
				TGLError<GenomeTrackFixedBin>("Failed to read a dense track file %s: %s", filename, strerror(errno));
			TGLError<GenomeTrackFixedBin>("Invalid format of a dense track file %s", filename);
		}

		header_start = 0;
		total_bytes = m_bfile.file_size();
	}

	// --- COMMON LOGIC ---
	const uint64_t header_size = m_bfile.tell() - header_start;
	if (total_bytes < header_size)
		TGLError<GenomeTrackFixedBin>("Invalid format of a dense track file %s", filename);
	const uint64_t data_bytes = total_bytes - header_size;

	if (m_bin_size <= 0 || (data_bytes % sizeof(float)) != 0)
		TGLError<GenomeTrackFixedBin>("Invalid format of a dense track file %s", filename);

	m_num_samples = (int64_t)(data_bytes / sizeof(float));
	m_chromid = chromid;

	// Set up mmap for read-only mode (naryn pattern)
	// For indexed format: reuse existing mmap if same file (avoid re-mmap per chromosome)
	m_mmap_data = nullptr;
	m_mmap_num_bins = 0;
	m_cur_bin = 0;

	if (strcmp(mode, "rb") == 0 && m_num_samples > 0) {
		const std::string file_path = m_dat_open ? m_dat_path : std::string(filename);

		// Only re-mmap if file changed (indexed format reuses same track.dat)
		if (!m_mmap.is_open() || m_mmap_path != file_path) {
			m_mmap.close();
			m_mmap.open(file_path, true /* sequential */);
			m_mmap_path = file_path;
		}

		if (m_mmap.is_open()) {
			const uint64_t data_offset = m_base_offset + sizeof(m_bin_size);
			if (data_offset + m_num_samples * sizeof(float) <= m_mmap.size()) {
				m_mmap_data = reinterpret_cast<const float *>(m_mmap.data() + data_offset);
				m_mmap_num_bins = m_num_samples;
			} else {
				m_mmap.close();  // file too small, fall back to BufferedFile
				m_mmap_path.clear();
			}
		}
	} else {
		// Write/update mode or empty: close any existing mmap
		if (m_mmap.is_open()) {
			m_mmap.close();
			m_mmap_path.clear();
		}
	}
}

void GenomeTrackFixedBin::init_write(const char *filename, unsigned bin_size, int chromid)
{
	m_num_samples = 0;
	m_cur_coord = 0;

	mode_t old_umask = umask(07);

	if (m_bfile.open(filename, "wb")) {
		umask(old_umask);
		TGLError<GenomeTrackFixedBin>("Opening a dense track file %s: %s", filename, strerror(errno));
	}

	umask(old_umask);

	m_bin_size = bin_size;
	if (m_bfile.write(&m_bin_size, sizeof(m_bin_size)) != sizeof(m_bin_size)) {
		if (m_bfile.error())
			TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s: %s", filename, strerror(errno));
		TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s", filename);
	}

	m_chromid = chromid;
}
