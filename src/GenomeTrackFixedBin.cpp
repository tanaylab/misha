#include <errno.h>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "TGLException.h"
#include "GenomeTrackFixedBin.h"

void GenomeTrackFixedBin::read_interval(const GInterval &interval)
{
	if (m_use_quantile)
		m_sp.reset();

	if (m_functions[MIN_POS])
		m_last_min_pos = numeric_limits<double>::quiet_NaN();

	// optimization of the most common case when the expression iterator starts at 0 and steps by bin_size
	if (interval.start == m_cur_coord && interval.end == m_cur_coord + m_bin_size) {
		if (read_next_bin(m_last_avg)) {
			m_last_min = m_last_max = m_last_nearest = m_last_sum = m_last_avg;
			m_last_stddev = numeric_limits<float>::quiet_NaN();
			if (m_functions[MAX_POS])
				m_last_max_pos = interval.start;
			if (m_functions[MIN_POS])
				m_last_min_pos = interval.start;
			if (m_use_quantile && !std::isnan(m_last_avg))
				m_sp.add(m_last_avg, s_rnd_func);
		} else {
			m_last_min = m_last_max = m_last_nearest = m_last_avg = m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();
			if (m_functions[MAX_POS])
				m_last_max_pos = numeric_limits<double>::quiet_NaN();
			if (m_functions[MIN_POS])
				m_last_min_pos = numeric_limits<double>::quiet_NaN();
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
			double overlap_start = std::max(static_cast<double>(sbin * m_bin_size), static_cast<double>(interval.start));
			if (m_functions[MAX_POS])
				m_last_max_pos = overlap_start;
			if (m_functions[MIN_POS])
				m_last_min_pos = overlap_start;
			if (m_use_quantile && !std::isnan(m_last_avg))
				m_sp.add(m_last_avg, s_rnd_func);
		} else {
			m_last_min = m_last_max = m_last_nearest = m_last_avg = m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();
			if (m_functions[MAX_POS])
				m_last_max_pos = numeric_limits<double>::quiet_NaN();
			if (m_functions[MIN_POS])
				m_last_min_pos = numeric_limits<double>::quiet_NaN();
		}
	} else {
		float num_vs = 0;
		double mean_square_sum = 0;
		float v;

		m_last_sum = 0;
		m_last_min = numeric_limits<float>::max();
		m_last_max = -numeric_limits<float>::max();
		if (m_functions[MAX_POS])
			m_last_max_pos = numeric_limits<double>::quiet_NaN();
		if (m_functions[MIN_POS])
			m_last_min_pos = numeric_limits<double>::quiet_NaN();

		goto_bin(sbin);
		for (int64_t bin = sbin; bin < ebin; ++bin) {
			if (read_next_bin(v) && !std::isnan(v)) {
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

				if (m_use_quantile && !std::isnan(v))
					m_sp.add(v, s_rnd_func);

				++num_vs;
			}
		}

		if (num_vs > 0)
			m_last_avg = m_last_nearest = m_last_sum / num_vs;
		else {
			m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_sum = numeric_limits<float>::quiet_NaN();
			if (m_functions[MIN_POS])
				m_last_min_pos = numeric_limits<double>::quiet_NaN();
		}

		// we are calaculating unbiased standard deviation:
		// sqrt(sum((x-mean)^2) / (N-1)) = sqrt(sum(x^2)/(N-1) - N*(mean^2)/(N-1))
		if (m_functions[STDDEV])
			m_last_stddev = num_vs > 1 ? sqrt(mean_square_sum / (num_vs - 1) - (m_last_avg * (double)m_last_avg) * (num_vs / (num_vs - 1))) : numeric_limits<float>::quiet_NaN();
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

void GenomeTrackFixedBin::init_read(const char *filename, const char *mode, int chromid)
{
	m_cur_coord = 0;

    if (m_bfile.open(filename, mode))
        TGLError<GenomeTrackFixedBin>("%s", strerror(errno));

    if (m_bfile.read(&m_bin_size, sizeof(m_bin_size)) != sizeof(m_bin_size)) {
        if (m_bfile.error())
            TGLError<GenomeTrackFixedBin>("Failed to read a dense track file %s: %s", filename, strerror(errno));
        TGLError<GenomeTrackFixedBin>("Invalid format of a dense track file %s", filename);
    }

    // determine the number of samples in the file
    double num_samples = (m_bfile.file_size() - m_bfile.tell()) / (double)sizeof(float);

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
