/*
 * GenomeTrack1D.h
 *
 *  Created on: Jan 25, 2012
 *      Author: hoichman
 */

#ifndef GENOMETRACK1D_H_
#define GENOMETRACK1D_H_

#include <cstdint>
#include <limits>
#include "GenomeTrack.h"
#include "GInterval.h"
#include "StreamPercentiler.h"

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrack1D : public GenomeTrack {
public:
	enum Functions { AVG, MIN, MAX, NEAREST, STDDEV, SUM, MAX_POS, MIN_POS, EXISTS, SIZE, SAMPLE, SAMPLE_POS, FIRST, FIRST_POS, LAST, LAST_POS, NUM_FUNCS };

	virtual ~GenomeTrack1D() {}

	int get_chrom_id() const { return m_chromid; }

	virtual void read_interval(const GInterval &interval) = 0;

	void register_function(Functions func) { m_functions[func] = true; }
	void register_quantile(uint64_t rnd_sampling_buf_size, uint64_t lowest_vals_buf_size, uint64_t highest_vals_buf_size);

	float last_avg() const { return m_last_avg; }
	float last_min() const { return m_last_min; }
	float last_max() const { return m_last_max; }
	virtual double last_max_pos() const = 0;
	virtual double last_min_pos() const = 0;
	float last_nearest() const { return m_last_nearest; }
	float last_stddev() const { return m_last_stddev; }
	float last_sum() const { return m_last_sum; }
	float last_quantile(double percentile);
	float last_exists() const { return m_last_exists; }
	float last_size() const { return m_last_size; }
	float last_sample() const { return m_last_sample; }
	double last_sample_pos() const { return m_last_sample_pos; }
	float last_first() const { return m_last_first; }
	double last_first_pos() const { return m_last_first_pos; }
	float last_last() const { return m_last_last; }
	double last_last_pos() const { return m_last_last_pos; }

	// Access to stream percentiler for combining quantiles across intervals
	const StreamPercentiler<float>& get_percentiler() const { return m_sp; }

	const string &file_name() const { return m_bfile.file_name(); }

protected:
	vector<bool> m_functions;
	bool         m_use_quantile;
	int          m_chromid;

	float        m_last_avg;
	float        m_last_min;
	float        m_last_max;
	double       m_last_max_pos;
	double       m_last_min_pos;
	float        m_last_nearest;
	float        m_last_stddev;
	float        m_last_sum;
	float        m_last_exists;
	float        m_last_size;
	float        m_last_sample;
	double       m_last_sample_pos;
	float        m_last_first;
	double       m_last_first_pos;
	float        m_last_last;
	double       m_last_last_pos;
	StreamPercentiler<float> m_sp;

	GenomeTrack1D(Type type) : GenomeTrack(type), m_use_quantile(false) {
		m_functions.resize(NUM_FUNCS, false);
		m_last_max_pos = numeric_limits<double>::quiet_NaN();
		m_last_min_pos = numeric_limits<double>::quiet_NaN();
		m_last_sample_pos = numeric_limits<double>::quiet_NaN();
		m_last_first_pos = numeric_limits<double>::quiet_NaN();
		m_last_last_pos = numeric_limits<double>::quiet_NaN();
	}
};


//--------------------------------------------- IMPLEMENTATION -----------------------------------------------------------

inline void GenomeTrack1D::register_quantile(uint64_t rnd_sampling_buf_size, uint64_t lowest_vals_buf_size, uint64_t highest_vals_buf_size)
{
	m_sp.init(rnd_sampling_buf_size, lowest_vals_buf_size, highest_vals_buf_size);
	m_use_quantile = true;
}

inline float GenomeTrack1D::last_quantile(double percentile)
{
	if (m_sp.stream_size()) {
		bool is_estimated;
		return m_sp.get_percentile(percentile, is_estimated);
	}
	return numeric_limits<float>::quiet_NaN();
}

#endif /* GENOMETRACK1D_H_ */
