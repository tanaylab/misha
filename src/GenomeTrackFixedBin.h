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
#include <deque>
#include <limits>
#include <string>
#include <vector>

#include "GenomeTrack1D.h"
#include "MmapFile.h"
#include "utils/RunningLogSumExp.h"

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrackFixedBin : public GenomeTrack1D {
public:
	GenomeTrackFixedBin() : GenomeTrack1D(FIXED_BIN), m_bin_size(0), m_num_samples(0), m_cur_coord(0), m_last_min_pos(numeric_limits<double>::quiet_NaN()) {}

	void read_interval(const GInterval &interval) override;
	double last_max_pos() const override;
	double last_min_pos() const override;
	void set_master_obj(GenomeTrackFixedBin *master_obj) { m_master_obj = master_obj; m_master_synced = false; }

	void init_read(const char *filename, int chromid) { init_read(filename, "rb", chromid, true); }
	void init_write(const char *filename, unsigned bin_size, int chromid);

	void init_update(const char *filename, int chromid) { init_read(filename, "rb+", chromid, true); }

	// Metadata-only init for callers that just need bin_size / num_samples
	// (e.g. the create_expr_iterator validation loop). Skips the mmap setup
	// and avoids re-reading the header on every chromosome of the same track —
	// bin_size is identical across chromosomes, so we cache it on the object
	// and only stat() subsequent files for size. With many tracks × many
	// chromosomes this dominated gextract setup time even after MAP_POPULATE
	// was removed.
	void init_read_metadata(const char *filename, int chromid);

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

	void init_read(const char *filename, const char *mode, int chromid, bool setup_mmap = true);

	// Helper to parse header at current file position
	void read_header_at_current_pos_(BufferedFile &bf);

	// Sliding window LSE state
	RunningLogSumExp m_running_lse;
	std::deque<float> m_lse_window_bins;
	int64_t m_lse_prev_sbin{-1};
	int64_t m_lse_prev_ebin{-1};
	double m_sliding_sum{0.0};
	double m_sliding_sum_comp{0.0};
	int64_t m_sliding_num_vs{0};
	bool m_lse_sliding_valid{false};
	bool m_running_lse_initialized{false};
	int m_fast_path_mode{0}; // 0=unknown, 1=reducer-only fast path, 2=avg/nearest-only fast path, 3=single-function fast path, -1=generic path
	Functions m_single_func{AVG}; // used when m_fast_path_mode == 3
	uint32_t m_fast_reducer_bits{0};
	GenomeTrackFixedBin *m_master_obj{NULL};
	bool m_master_synced{false};

	// mmap-backed read path (naryn pattern): pointer dereference instead of fread
	MmapFile m_mmap;
	std::string m_mmap_path;  // track which file is mmap'd (avoid re-mmap on chrom switch)
	const float *m_mmap_data{nullptr};  // points to first bin value in mmap'd region
	int64_t m_mmap_num_bins{0};
	int64_t m_cur_bin{0};  // current bin index for mmap path

	// Scratch buffers reused across read_interval calls to avoid per-call heap allocation
	std::vector<float> m_scratch_all_values;
	std::vector<double> m_scratch_all_positions;
	std::vector<float> m_scratch_bin_vals;

	void read_interval_reducers_only(const GInterval &interval);
	void read_interval_avg_nearest_only(const GInterval &interval);
	void read_interval_single_function(const GInterval &interval);
	void sync_master_state_from_dependent();
	void copy_state_from_master();
	void classify_fast_path_mode();
	void assign_single_bin_value(float value, double overlap_start);
	void assign_single_bin_missing();
	void reset_sliding_window_state();

	// Kahan compensated summation helpers for m_sliding_sum
	inline void kahan_add_to_sliding_sum(double value) {
		double y = value - m_sliding_sum_comp;
		double t = m_sliding_sum + y;
		m_sliding_sum_comp = (t - m_sliding_sum) - y;
		m_sliding_sum = t;
	}
	inline void kahan_sub_from_sliding_sum(double value) {
		kahan_add_to_sliding_sum(-value);
	}
};


//------------------------------ IMPLEMENTATION ------------------------------------

inline void GenomeTrackFixedBin::goto_bin(uint64_t bin)
{
	if (m_mmap_data) {
		m_cur_bin = bin;
	} else {
		if (m_bfile.seek((long)(m_base_offset + sizeof(m_bin_size) + (uint64_t)bin * sizeof(float)), SEEK_SET))
			TGLError<GenomeTrackFixedBin>("Failed to seek a dense track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
	}
	m_cur_coord = bin * m_bin_size;
}


inline bool GenomeTrackFixedBin::read_next_bin(float &val)
{
	if (m_mmap_data) {
		if (m_cur_bin >= m_mmap_num_bins)
			return false;
		val = m_mmap_data[m_cur_bin++];
		if (isinf(val))
			val = numeric_limits<float>::quiet_NaN();
		m_cur_coord += m_bin_size;
		return true;
	}

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
