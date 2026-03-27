/*
 * GenomeTrack2D.h
 *
 *  Created on: Jan 12, 2012
 *      Author: hoichman
 */

#ifndef GENOMETRACK2D_H_
#define GENOMETRACK2D_H_

#include <cstdint>
#include <string>
#include <utility>

#include "BufferedFile.h"
#include "DiagonalBand.h"
#include "Rectangle.h"

#include "GenomeTrack.h"
#include "GInterval2D.h"

using namespace std;

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrack2D : public GenomeTrack {
public:
	enum Functions { WEIGHTED_SUM, OCCUPIED_AREA, AVG, MIN, MAX, NUM_FUNCS };

	GenomeTrack2D();

	virtual ~GenomeTrack2D() {}

	virtual void init_read(const char *filename, int chromid1, int chromid2);
	virtual void init_write(const char *filename, int chromid1, int chromid2);

	pair<int, int> get_chrom_ids() const { return pair<int, int>(m_chromid1, m_chromid2); }

	virtual void read_interval(const Rectangle &rectangle, const DiagonalBand &band) = 0;

	virtual bool begin_interval() = 0;
	virtual bool next_interval() = 0;
	virtual bool is_end_interval() = 0;
	const GInterval2D &cur_interval() const { return m_interval; }

	void register_function(Functions func) { m_func_mask |= (1u << func); }
	bool has_function(Functions func) const { return (m_func_mask & (1u << func)) != 0; }

	float last_weighted_sum() const { return m_last_weighted_sum; }
	float last_occupied_area() const { return m_last_occupied_area; }
	float last_avg() const { return m_last_weighted_sum / m_last_occupied_area; }
	float last_min() const { return m_last_min; }
	float last_max() const { return m_last_max; }

	const std::string &file_name() const { return m_bfile.file_name(); }

	// Returns true if the current chromosome pair has data available.
	// In indexed mode, m_bfile may stay open even for empty pairs,
	// so opened() alone is not sufficient.
	bool has_data_for_pair() const { return m_pair_has_data; }

protected:
	uint32_t     m_func_mask{0};
	int          m_chromid1;
	int          m_chromid2;
	bool         m_loaded;

	int64_t      m_last_occupied_area;
	double       m_last_weighted_sum;
	float        m_last_min;
	float        m_last_max;

	GInterval2D  m_interval;

	// State for indexed "smart handle"
	std::string m_dat_path;
	bool        m_dat_open{false};
	bool        m_pair_has_data{false};  // true if current pair has data (replaces opened() check for 2D)
	int64_t     m_base_offset_2d{0};     // base offset for quad tree in indexed format

	GenomeTrack2D(Type type) : GenomeTrack(type), m_loaded(false), m_dat_open(false), m_pair_has_data(false), m_base_offset_2d(0) {}
};

#endif /* GENOMETRACK2D_H_ */
