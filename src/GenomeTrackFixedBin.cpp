#include <errno.h>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "TGLException.h"
#include "GenomeTrackFixedBin.h"
#include "TrackIndex.h"

void GenomeTrackFixedBin::read_interval_sum_only(const GInterval &interval)
{
	// Common case: iterator advances exactly one dense bin.
	if (interval.start == m_cur_coord && interval.end == m_cur_coord + m_bin_size) {
		float v = numeric_limits<float>::quiet_NaN();
		if (read_next_bin(v) && !std::isnan(v)) {
			m_last_sum = v;
			m_cached_bin_idx = (int64_t)(interval.start / m_bin_size);
			m_cached_bin_val = v;
			m_cache_valid = true;
		} else
			m_last_sum = numeric_limits<float>::quiet_NaN();
		m_lse_sliding_valid = false;
		m_sliding_sum = 0;
		m_sliding_num_vs = 0;
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

		m_last_sum = (have_value && !std::isnan(v)) ? v : numeric_limits<float>::quiet_NaN();
		m_lse_sliding_valid = false;
		m_sliding_sum = 0;
		m_sliding_num_vs = 0;
		return;
	}

	const int64_t window_size = ebin - sbin;

	// Sliding sum update for one-bin steps.
	if (m_lse_sliding_valid && window_size > 0) {
		int64_t step = sbin - m_lse_prev_sbin;
		int64_t prev_window = m_lse_prev_ebin - m_lse_prev_sbin;
		if (step > 0 && step <= prev_window && window_size == prev_window &&
			(int64_t)m_lse_window_bins.size() == prev_window) {
			if (step == 1) {
				float new_val = numeric_limits<float>::quiet_NaN();
				if (m_cur_coord != m_lse_prev_ebin * m_bin_size)
					goto_bin(m_lse_prev_ebin);
				if (read_next_bin(new_val) && !m_lse_window_bins.empty()) {
					float old_val = m_lse_window_bins.front();
					m_lse_window_bins.pop_front();
					if (!std::isnan(old_val)) {
						m_sliding_sum -= old_val;
						--m_sliding_num_vs;
					}

					m_lse_window_bins.push_back(new_val);
					if (!std::isnan(new_val)) {
						m_sliding_sum += new_val;
						++m_sliding_num_vs;
					}

					m_last_sum = m_sliding_num_vs > 0 ? (float)m_sliding_sum : numeric_limits<float>::quiet_NaN();
					m_lse_prev_sbin = sbin;
					m_lse_prev_ebin = ebin;
					m_lse_sliding_valid = true;
					m_cached_bin_idx = ebin - 1;
					m_cached_bin_val = new_val;
					m_cache_valid = true;
					return;
				}
			}
		}
	}

	// Fallback: full window read.
	vector<float> bin_vals;
	int64_t bins_read = read_bins_bulk(sbin, window_size, bin_vals);

	m_sliding_sum = 0;
	m_sliding_num_vs = 0;
	m_lse_window_bins.clear();
	for (int64_t i = 0; i < bins_read; ++i) {
		float v = bin_vals[i];
		m_lse_window_bins.push_back(v);
		if (!std::isnan(v)) {
			m_sliding_sum += v;
			++m_sliding_num_vs;
		}
	}

	if (bins_read > 0) {
		m_cached_bin_idx = sbin + bins_read - 1;
		m_cached_bin_val = bin_vals[bins_read - 1];
		m_cache_valid = true;
	}

	m_last_sum = m_sliding_num_vs > 0 ? (float)m_sliding_sum : numeric_limits<float>::quiet_NaN();
	m_lse_prev_sbin = sbin;
	m_lse_prev_ebin = ebin;
	m_lse_sliding_valid = bins_read == window_size && window_size > 0;
}

void GenomeTrackFixedBin::read_interval(const GInterval &interval)
{
	if (m_fast_path_mode == 0) {
		bool sum_only = !m_use_quantile && m_functions[SUM];
		for (size_t i = 0; i < m_functions.size() && sum_only; ++i) {
			if ((int)i != SUM && m_functions[i])
				sum_only = false;
		}
		m_fast_path_mode = sum_only ? 1 : -1;
	}

	if (m_fast_path_mode == 1) {
		read_interval_sum_only(interval);
		return;
	}

	if (m_use_quantile)
		m_sp.reset();

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

	// optimization of the most common case when the expression iterator starts at 0 and steps by bin_size
	if (interval.start == m_cur_coord && interval.end == m_cur_coord + m_bin_size) {
		if (read_next_bin(m_last_avg)) {
			m_cached_bin_idx = (int64_t)(interval.start / m_bin_size);
			m_cached_bin_val = m_last_avg;
			m_cache_valid = true;

			m_last_min = m_last_max = m_last_nearest = m_last_sum = m_last_avg;
			m_last_stddev = numeric_limits<float>::quiet_NaN();
			if (m_functions[LSE])
				m_last_lse = m_last_avg;
			if (m_functions[MAX_POS])
				m_last_max_pos = interval.start;
			if (m_functions[MIN_POS])
				m_last_min_pos = interval.start;
			if (m_functions[EXISTS])
				m_last_exists = 1;
			if (m_functions[SIZE])
				m_last_size = 1;
			if (m_functions[SAMPLE])
				m_last_sample = m_last_avg;
			if (m_functions[SAMPLE_POS])
				m_last_sample_pos = interval.start;
			if (m_functions[FIRST])
				m_last_first = m_last_avg;
			if (m_functions[FIRST_POS])
				m_last_first_pos = interval.start;
			if (m_functions[LAST])
				m_last_last = m_last_avg;
			if (m_functions[LAST_POS])
				m_last_last_pos = interval.start;
			if (m_use_quantile && !std::isnan(m_last_avg))
				m_sp.add(m_last_avg, s_rnd_func);
		} else {
			m_last_min = m_last_max = m_last_nearest = m_last_avg = m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();
			if (m_functions[LSE])
				m_last_lse = numeric_limits<float>::quiet_NaN();
			if (m_functions[MAX_POS])
				m_last_max_pos = numeric_limits<double>::quiet_NaN();
			if (m_functions[MIN_POS])
				m_last_min_pos = numeric_limits<double>::quiet_NaN();
		}
		m_lse_sliding_valid = false;
		m_sliding_sum = 0;
		m_sliding_num_vs = 0;
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
			m_last_min = m_last_max = m_last_nearest = m_last_sum = m_last_avg;
			m_last_stddev = numeric_limits<float>::quiet_NaN();
			if (m_functions[LSE])
				m_last_lse = m_last_avg;
			double overlap_start = std::max(static_cast<double>(sbin * m_bin_size), static_cast<double>(interval.start));
			if (m_functions[MAX_POS])
				m_last_max_pos = overlap_start;
			if (m_functions[MIN_POS])
				m_last_min_pos = overlap_start;
			if (m_functions[EXISTS])
				m_last_exists = 1;
			if (m_functions[SIZE])
				m_last_size = 1;
			if (m_functions[SAMPLE])
				m_last_sample = m_last_avg;
			if (m_functions[SAMPLE_POS])
				m_last_sample_pos = overlap_start;
			if (m_functions[FIRST])
				m_last_first = m_last_avg;
			if (m_functions[FIRST_POS])
				m_last_first_pos = overlap_start;
			if (m_functions[LAST])
				m_last_last = m_last_avg;
			if (m_functions[LAST_POS])
				m_last_last_pos = overlap_start;
			if (m_use_quantile && !std::isnan(m_last_avg))
				m_sp.add(m_last_avg, s_rnd_func);
		} else {
			m_last_min = m_last_max = m_last_nearest = m_last_avg = m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();
			if (m_functions[LSE])
				m_last_lse = numeric_limits<float>::quiet_NaN();
			if (m_functions[MAX_POS])
				m_last_max_pos = numeric_limits<double>::quiet_NaN();
			if (m_functions[MIN_POS])
				m_last_min_pos = numeric_limits<double>::quiet_NaN();
		}
		m_lse_sliding_valid = false;
		m_sliding_sum = 0;
		m_sliding_num_vs = 0;
	} else {
		float num_vs = 0;
		double mean_square_sum = 0;
		const int64_t window_size = ebin - sbin;
		const bool simple_sliding_compatible =
			!m_use_quantile &&
			!m_functions[MIN] &&
			!m_functions[MAX] &&
			!m_functions[MAX_POS] &&
			!m_functions[MIN_POS] &&
			!m_functions[STDDEV] &&
			!m_functions[SAMPLE] &&
			!m_functions[SAMPLE_POS] &&
			!m_functions[FIRST] &&
			!m_functions[FIRST_POS] &&
			!m_functions[LAST] &&
			!m_functions[LAST_POS];

		// For sampling, collect all values/positions
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

		bool simple_sliding_used = false;

		if (simple_sliding_compatible && m_lse_sliding_valid && window_size > 0) {
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
							m_sliding_sum -= old_val;
							--m_sliding_num_vs;
							if (m_functions[LSE])
								m_running_lse.pop_front();
						}

						m_lse_window_bins.push_back(new_val);
						if (!std::isnan(new_val)) {
							m_sliding_sum += new_val;
							++m_sliding_num_vs;
							if (m_functions[LSE])
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
								m_sliding_sum -= old_val;
								--m_sliding_num_vs;
								if (m_functions[LSE])
									m_running_lse.pop_front();
							}
						}
						for (int64_t i = 0; i < appended; ++i) {
							float new_val = new_vals[i];
							m_lse_window_bins.push_back(new_val);
							if (!std::isnan(new_val)) {
								m_sliding_sum += new_val;
								++m_sliding_num_vs;
								if (m_functions[LSE])
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
						if (m_functions[AVG] || m_functions[NEAREST])
							m_last_avg = m_last_nearest = (float)(m_sliding_sum / m_sliding_num_vs);
						if (m_functions[LSE])
							m_last_lse = m_running_lse.window.empty()
								? numeric_limits<float>::quiet_NaN()
								: (float)m_running_lse.value();
						if (m_functions[EXISTS])
							m_last_exists = 1;
					} else {
						m_last_sum = numeric_limits<float>::quiet_NaN();
						if (m_functions[AVG] || m_functions[NEAREST])
							m_last_avg = m_last_nearest = numeric_limits<float>::quiet_NaN();
						if (m_functions[LSE])
							m_last_lse = numeric_limits<float>::quiet_NaN();
						if (m_functions[EXISTS])
							m_last_exists = 0;
					}
					if (m_functions[SIZE])
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

			if (m_functions[LSE] && m_lse_sliding_valid && bins_read > 0 && !simple_sliding_compatible) {
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
					m_last_sum += v;
					double bin_start = static_cast<double>(bin * m_bin_size);
					double overlap_start = std::max(bin_start, static_cast<double>(interval.start));
					if (v < m_last_min) {
						m_last_min = v;
						if (m_functions[MIN_POS])
							m_last_min_pos = overlap_start;
					} else if (m_functions[MIN_POS] && v == m_last_min) {
						double candidate_pos = overlap_start;
						if (std::isnan(m_last_min_pos) || candidate_pos < m_last_min_pos)
							m_last_min_pos = candidate_pos;
					}
					if (v > m_last_max) {
						m_last_max = v;
						if (m_functions[MAX_POS])
							m_last_max_pos = overlap_start;
					}

					if (m_functions[STDDEV])
						mean_square_sum += v * v;

					if (m_functions[LSE] && !lse_sliding_used)
						lse_accumulate(m_last_lse, v);

					if (m_use_quantile && !std::isnan(v))
						m_sp.add(v, s_rnd_func);

					// New virtual track computations
					if (m_functions[EXISTS])
						m_last_exists = 1;

					if (m_functions[FIRST] && std::isnan(m_last_first))
						m_last_first = v;

					if (m_functions[FIRST_POS] && std::isnan(m_last_first_pos))
						m_last_first_pos = overlap_start;

					if (m_functions[LAST])
						m_last_last = v;

					if (m_functions[LAST_POS])
						m_last_last_pos = overlap_start;

					if (m_functions[SAMPLE])
						all_values.push_back(v);
					if (m_functions[SAMPLE_POS])
						all_positions.push_back(overlap_start);

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

			// we are calaculating unbiased standard deviation:
			// sqrt(sum((x-mean)^2) / (N-1)) = sqrt(sum(x^2)/(N-1) - N*(mean^2)/(N-1))
			if (m_functions[STDDEV])
				m_last_stddev = num_vs > 1 ? sqrt(mean_square_sum / (num_vs - 1) - (m_last_avg * (double)m_last_avg) * (num_vs / (num_vs - 1))) : numeric_limits<float>::quiet_NaN();

			if (m_functions[LSE] || simple_sliding_compatible) {
				if (!lse_sliding_used) {
					// Full computation was done; initialize sliding state from bin_vals
					m_running_lse.clear();
					m_lse_window_bins.clear();
					m_sliding_sum = 0;
					m_sliding_num_vs = 0;
					for (int64_t i = 0; i < bins_read; i++) {
						m_lse_window_bins.push_back(bin_vals[i]);
						if (!std::isnan(bin_vals[i])) {
							m_sliding_sum += bin_vals[i];
							++m_sliding_num_vs;
							if (m_functions[LSE])
								m_running_lse.push(bin_vals[i]);
						}
					}
				}
				m_lse_prev_sbin = sbin;
				m_lse_prev_ebin = ebin;
				m_lse_sliding_valid = bins_read == window_size && window_size > 0;
			} else {
				m_lse_sliding_valid = false;
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
	if (m_cur_coord != start_bin * m_bin_size)
		goto_bin(start_bin);

	// Bulk read all bins in one syscall
	size_t bytes_to_read = to_read * sizeof(float);
	uint64_t bytes_read = m_bfile.read(vals.data(), bytes_to_read);

	if (bytes_read != bytes_to_read) {
		if (m_bfile.error())
			TGLError<GenomeTrackFixedBin>("Failed to read a dense track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		// Partial read - adjust size
		to_read = bytes_read / sizeof(float);
		vals.resize(to_read);
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
	m_sliding_num_vs = 0;
	m_lse_sliding_valid = false;

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
	double num_samples = (total_bytes - header_size) / (double)sizeof(float);

	if (m_bin_size <= 0 || num_samples != (int64_t)num_samples)
		TGLError<GenomeTrackFixedBin>("Invalid format of a dense track file %s", filename);

	m_num_samples = (int64_t)num_samples;
	m_chromid = chromid;
}

void GenomeTrackFixedBin::init_write(const char *filename, unsigned bin_size, int chromid)
{
	m_num_samples = 0;
	m_cur_coord = 0;

	umask(07);

	if (m_bfile.open(filename, "wb"))
		TGLError<GenomeTrackFixedBin>("Opening a dense track file %s: %s", filename, strerror(errno));

	m_bin_size = bin_size;
	if (m_bfile.write(&m_bin_size, sizeof(m_bin_size)) != sizeof(m_bin_size)) {
		if (m_bfile.error())
			TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s: %s", filename, strerror(errno));
		TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s", filename);
	}

	m_chromid = chromid;
}
