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

protected:
	GIntervals    m_intervals;
	vector<float> m_vals;
	GIntervals::const_iterator m_icur_interval;
	double        m_last_min_pos;

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
