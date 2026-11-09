/*
 * GenomeTrackInMemory.h
 *
 * In-memory track with count-based statistics (matching GenomeTrackSparse)
 * Used for value-based virtual tracks
 */

#ifndef GENOMETRACKINMEMORY_H_
#define GENOMETRACKINMEMORY_H_

#include <cmath>
#include <limits>
#include <vector>

#include "GenomeTrack1D.h"
#include "GIntervals.h"

using namespace std;

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrackInMemory : public GenomeTrack1D {
public:
	GenomeTrackInMemory();

	// Initialize from in-memory data (no file I/O)
	void init_from_data(const GIntervals &intervals, const vector<float> &vals, int chromid);

	virtual void read_interval(const GInterval &interval);
	virtual double last_max_pos() const;
	virtual double last_min_pos() const;

	virtual bool item_at(int chromid, int64_t coord, int64_t &start, int64_t &end) {
		select_chrom(chromid);
		return interval_at(m_chrom_begin, m_chrom_end, coord, start, end);
	}

	// The intervals and values of every chromosome; [begin, end) are those of chromid
	const GIntervals &get_intervals() const { return m_intervals; }
	const vector<float> &get_vals() const { return m_vals; }
	void chrom_range(int chromid, size_t &begin, size_t &end) {
		select_chrom(chromid);
		begin = m_chrom_begin - m_intervals.cbegin();
		end = m_chrom_end - m_intervals.cbegin();
	}

protected:
	GIntervals    m_intervals;
	vector<float> m_vals;
	GIntervals::const_iterator m_icur_interval;
	// m_intervals holds every chromosome (sorted by chrom, then start); queries see only [m_chrom_begin, m_chrom_end)
	GIntervals::const_iterator m_chrom_begin;
	GIntervals::const_iterator m_chrom_end;
	int           m_cur_chromid;
	double        m_last_min_pos;

	void select_chrom(int chromid);

	// Count-based calculation (matching GenomeTrackSparse exactly)
	void calc_vals_coverage_weighted(const GInterval &interval);

	// Helper for finding first overlapping interval
	bool check_first_overlap(const GIntervals::const_iterator &iinterval1, const GInterval &interval2);
};


//------------------------------------ IMPLEMENTATION --------------------------------

inline bool GenomeTrackInMemory::check_first_overlap(
	const GIntervals::const_iterator &iinterval1,
	const GInterval &interval2)
{
	return iinterval1->do_overlap(interval2) &&
	       (iinterval1 == m_intervals.begin() || !(iinterval1 - 1)->do_overlap(interval2));
}

#endif /* GENOMETRACKINMEMORY_H_ */
