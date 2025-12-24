/*
 * GenomeTrackFixedBin.h
 *
 *  Created on: May 15, 2011
 *      Author: hoichman
 */

#ifndef GENOMETRACKFIXEDBIN_H_
#define GENOMETRACKFIXEDBIN_H_

#include <cstdint>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "GenomeTrack1D.h"

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrackFixedBin : public GenomeTrack1D {
public:
	GenomeTrackFixedBin() : GenomeTrack1D(FIXED_BIN), m_bin_size(0), m_num_samples(0), m_cur_coord(0), m_last_min_pos(numeric_limits<double>::quiet_NaN()) {}

	void read_interval(const GInterval &interval) override;
	double last_max_pos() const override;
	double last_min_pos() const override;

	void init_read(const char *filename, int chromid) { init_read(filename, "rb", chromid); }
	void init_write(const char *filename, unsigned bin_size, int chromid);

	void init_update(const char *filename, int chromid) { init_read(filename, "rb+", chromid); }

	unsigned get_bin_size() const { return m_bin_size; }
	int64_t  get_num_samples() const { return m_num_samples; }

	void goto_bin(uint64_t bin);

	bool read_next_bin(float &val);
	// Bulk read multiple bins into buffer, returns number of bins actually read
	int64_t read_bins_bulk(int64_t start_bin, int64_t num_bins, std::vector<float> &vals);
	void write_next_bin(float val);
	void write_next_bins(float *vals, uint64_t num_vals);

protected:
	unsigned  m_bin_size;
	int64_t   m_num_samples;
	int64_t   m_cur_coord;
	double    m_last_min_pos;
	int64_t   m_base_offset{0};
	int64_t   m_cached_bin_idx{-1};
	float     m_cached_bin_val{numeric_limits<float>::quiet_NaN()};
	bool      m_cache_valid{false};

	// State for indexed "smart handle"
	std::string m_dat_path;
	std::string m_dat_mode;
	bool        m_dat_open{false};

	void init_read(const char *filename, const char *mode, int chromid);

	// Helper to parse header at current file position
	void read_header_at_current_pos_(BufferedFile &bf);
};


//------------------------------ IMPLEMENTATION ------------------------------------

inline void GenomeTrackFixedBin::goto_bin(uint64_t bin)
{
	// Add m_base_offset to the absolute seek for indexed format support
	if (m_bfile.seek((long)(m_base_offset + sizeof(m_bin_size) + (uint64_t)bin * sizeof(float)), SEEK_SET))
		TGLError<GenomeTrackFixedBin>("Failed to seek a dense track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
	m_cur_coord = bin * m_bin_size;
}


inline bool GenomeTrackFixedBin::read_next_bin(float &val)
{
	if (m_bfile.read(&val, sizeof(val)) != sizeof(val)) {
		if (m_bfile.error())
			TGLError<GenomeTrackFixedBin>("Failed to read a dense track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		return false;
	}

	if (isinf(val))
		val = numeric_limits<float>::quiet_NaN();

	m_cur_coord += m_bin_size;
	return true;
}

inline void GenomeTrackFixedBin::write_next_bin(float val)
{
	if (m_bfile.write(&val, sizeof(val)) != sizeof(val)) {
		if (m_bfile.error())
			TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s", m_bfile.file_name().c_str());
	}
	m_num_samples++;
	m_cur_coord += m_bin_size;
}

inline void GenomeTrackFixedBin::write_next_bins(float *vals, uint64_t num_vals)
{
	uint64_t size = sizeof(vals[0]) * num_vals;
	if (m_bfile.write(vals, size) != size) {
		if (m_bfile.error())
			TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s", m_bfile.file_name().c_str());
	}
	m_num_samples += num_vals;
	m_cur_coord += m_bin_size * num_vals;
}

#endif /* GENOMETRACKFIXEDBIN_H_ */
