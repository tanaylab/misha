/*
 * GenomeTrackSparse.h
 *
 *  Created on: May 24, 2011
 *      Author: hoichman
 */

#ifndef GENOMETRACKSPARSE_H_
#define GENOMETRACKSPARSE_H_

#include <cmath>
#include <limits>
#include <vector>
#include <string>

#include "GenomeTrack1D.h"
#include "GIntervals.h"

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrackSparse : public GenomeTrack1D {
public:
	GenomeTrackSparse();

	virtual void read_interval(const GInterval &interval);
	virtual double last_max_pos() const;
	virtual double last_min_pos() const;
	void set_master_obj(GenomeTrackSparse *master_obj) { m_master_obj = master_obj; m_master_synced = false; }

	void init_read(const char *filename, int chromid);
	void init_write(const char *filename, int chromid);

	void write_next_interval(const GInterval &interval, float val);

	const GIntervals &get_intervals();
	const vector<float> &get_vals();

protected:
	static const int RECORD_SIZE;

	GIntervals    m_intervals;
	vector<float> m_vals;
	bool          m_loaded;
	int64_t       m_num_records;
	GIntervals::const_iterator m_icur_interval;
	double        m_last_min_pos;
	GenomeTrackSparse *m_master_obj{NULL};
	bool          m_master_synced{false};

	// State for indexed "smart handle"
	std::string m_dat_path;
	std::string m_dat_mode;
	bool        m_dat_open{false};

	void read_file_into_mem();
	void calc_vals(const GInterval &interval);
	bool check_first_overlap(const GIntervals::const_iterator &iinterval1, const GInterval &interval2);
	void sync_master_state_from_dependent();
	void copy_state_from_master();

	// Helper to parse header at current file position
	void read_header_at_current_pos_(BufferedFile &bf);

	// On-disk record size constant (matches RECORD_SIZE)
	static constexpr size_t kSparseRecBytes = sizeof(int64_t) * 2 + sizeof(float);
};


//------------------------------------ IMPLEMENTATION --------------------------------

inline bool GenomeTrackSparse::check_first_overlap(const GIntervals::const_iterator &iinterval1, const GInterval &interval2)
{
	return iinterval1->do_overlap(interval2) && (iinterval1 == m_intervals.begin() || !(iinterval1 - 1)->do_overlap(interval2));
}

inline void GenomeTrackSparse::calc_vals(const GInterval &interval)
{
	float num_vs = 0;
	double mean_square_sum = 0;
	float v;

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

	for (GIntervals::const_iterator iinterv = m_icur_interval; iinterv != m_intervals.end(); ++iinterv) {
		if (!iinterv->do_overlap(interval))
			break;

		v = m_vals[iinterv - m_intervals.begin()];
		if (!std::isnan(v)) {
			m_last_sum += v;
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

			// New virtual track computations
			if (m_functions[EXISTS])
				m_last_exists = 1;

			if (m_functions[FIRST] && std::isnan(m_last_first))
				m_last_first = v;

			if (m_functions[FIRST_POS] && std::isnan(m_last_first_pos))
				m_last_first_pos = iinterv->start;

			if (m_functions[LAST])
				m_last_last = v;

			if (m_functions[LAST_POS])
				m_last_last_pos = iinterv->start;

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

	// we are calaculating unbiased standard deviation:
	// sqrt(sum((x-mean)^2) / (N-1)) = sqrt(sum(x^2)/(N-1) - N*(mean^2)/(N-1))
	if (m_functions[STDDEV])
		m_last_stddev = num_vs > 1 ? sqrt(mean_square_sum / (num_vs - 1) - (m_last_avg * (double)m_last_avg) * (num_vs / (num_vs - 1))) : numeric_limits<float>::quiet_NaN();
}

#endif /* GENOMETRACKSPARSE_H_ */
