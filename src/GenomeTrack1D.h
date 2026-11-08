/*
 * GenomeTrack1D.h
 *
 *  Created on: Jan 25, 2012
 *      Author: hoichman
 */

#ifndef GENOMETRACK1D_H_
#define GENOMETRACK1D_H_

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>
#include "GenomeTrack.h"
#include "GInterval.h"
#include "StreamPercentiler.h"

// Numerically stable log-sum-exp helper: l1 = log(exp(l1) + exp(l2))
inline void lse_accumulate(float &l1, float l2) {
	if (l1 > l2) {
		if (!std::isinf(l2))
			l1 += std::log(1.0f + std::exp(l2 - l1));
	} else {
		if (std::isinf(l1))
			l1 = l2;
		else
			l1 = l2 + std::log(1.0f + std::exp(l1 - l2));
	}
}

inline void lse_accumulate(double &l1, double l2) {
	if (l1 > l2) {
		if (!std::isinf(l2))
			l1 += std::log(1.0 + std::exp(l2 - l1));
	} else {
		if (std::isinf(l1))
			l1 = l2;
		else
			l1 = l2 + std::log(1.0 + std::exp(l1 - l2));
	}
}

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrack1D : public GenomeTrack {
public:
	// WEIGHTED_SUM: sum of value x the bp of the item inside the interval (each other function counts an item in full)
	enum Functions { AVG, MIN, MAX, NEAREST, STDDEV, SUM, LSE, MAX_POS, MIN_POS, EXISTS, SIZE, SAMPLE, SAMPLE_POS, FIRST, FIRST_POS, LAST, LAST_POS, WEIGHTED_SUM, NUM_FUNCS };

	virtual ~GenomeTrack1D() {}

	int get_chrom_id() const { return m_chromid; }

	virtual void read_interval(const GInterval &interval) = 0;

	void register_function(Functions func) { m_func_mask |= (1u << func); }
	bool has_function(Functions func) const { return (m_func_mask & (1u << func)) != 0; }
	void register_quantile(uint64_t rnd_sampling_buf_size, uint64_t lowest_vals_buf_size, uint64_t highest_vals_buf_size);

	float last_avg() const { return m_last_avg; }
	float last_min() const { return m_last_min; }
	float last_max() const { return m_last_max; }
	virtual double last_max_pos() const = 0;
	virtual double last_min_pos() const = 0;
	float last_nearest() const { return m_last_nearest; }
	float last_stddev() const { return m_last_stddev; }
	float last_sum() const { return m_last_sum; }
	float last_lse() const { return m_last_lse; }
	float last_quantile(double percentile);
	float last_exists() const { return m_last_exists; }
	float last_size() const { return m_last_size; }
	float last_sample() const { return m_last_sample; }
	double last_sample_pos() const { return m_last_sample_pos; }
	float last_first() const { return m_last_first; }
	double last_first_pos() const { return m_last_first_pos; }
	float last_last() const { return m_last_last; }
	double last_last_pos() const { return m_last_last_pos; }
	double last_weighted_sum() const { return m_last_weighted_sum; }

	// Access to stream percentiler for combining quantiles across intervals
	const StreamPercentiler<float>& get_percentiler() const { return m_sp; }

	// The track item (dense bin, or sparse, array or in-memory interval) that holds coord on chromid, if any
	virtual bool item_at(int chromid, int64_t coord, int64_t &start, int64_t &end) { return false; }

	// The filter's unmasked parts (sorted, disjoint), each with the item it shares with the part before it cut off,
	// so that reading the parts one by one reads an item that the filter splits only once.
	std::vector<GInterval> parts_reading_items_once(const std::vector<GInterval> &parts);

	const string &file_name() const { return m_bfile.file_name(); }

protected:
	uint32_t     m_func_mask{0};
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
	float        m_last_lse;
	float        m_last_exists;
	float        m_last_size;
	float        m_last_sample;
	double       m_last_sample_pos;
	float        m_last_first;
	double       m_last_first_pos;
	float        m_last_last;
	double       m_last_last_pos;
	double       m_last_weighted_sum;
	StreamPercentiler<float> m_sp;

	GenomeTrack1D(Type type) : GenomeTrack(type), m_use_quantile(false) {
		m_last_max_pos = numeric_limits<double>::quiet_NaN();
		m_last_min_pos = numeric_limits<double>::quiet_NaN();
		m_last_sample_pos = numeric_limits<double>::quiet_NaN();
		m_last_first_pos = numeric_limits<double>::quiet_NaN();
		m_last_last_pos = numeric_limits<double>::quiet_NaN();
		m_last_weighted_sum = numeric_limits<double>::quiet_NaN();
	}

	// item_at over intervals sorted by start that do not overlap each other
	static bool interval_at(std::vector<GInterval>::const_iterator begin, std::vector<GInterval>::const_iterator end,
	                        int64_t coord, int64_t &start, int64_t &finish);
};


//--------------------------------------------- IMPLEMENTATION -----------------------------------------------------------

inline std::vector<GInterval> GenomeTrack1D::parts_reading_items_once(const std::vector<GInterval> &parts)
{
	std::vector<GInterval> res;
	res.reserve(parts.size());
	for (size_t i = 0; i < parts.size(); ++i) {
		GInterval part = parts[i];
		int64_t item_start, item_end;
		// An item holding this part's start that begins before the previous part ends overlaps that part as well
		// and was read with it (or, if that part was cut away entirely, with the part before it).
		if (i && item_at(part.chromid, part.start, item_start, item_end) && item_start < parts[i - 1].end)
			part.start = item_end;
		if (part.start < part.end)
			res.push_back(part);
	}
	return res;
}

inline bool GenomeTrack1D::interval_at(std::vector<GInterval>::const_iterator begin, std::vector<GInterval>::const_iterator end,
                                       int64_t coord, int64_t &start, int64_t &finish)
{
	auto it = std::upper_bound(begin, end, coord, [](int64_t c, const GInterval &interv) { return c < interv.start; });
	if (it == begin || (--it)->end <= coord)
		return false;
	start = it->start;
	finish = it->end;
	return true;
}

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
