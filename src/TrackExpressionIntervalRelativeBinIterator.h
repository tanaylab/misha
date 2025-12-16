#ifndef TRACKEXPRESSIONINTERVALRELATIVEBINITERATOR_H_
#define TRACKEXPRESSIONINTERVALRELATIVEBINITERATOR_H_

#include "TrackExpressionIterator.h"
#include "GIntervals.h"

//------------------------------ TrackExpressionIntervalRelativeBinIterator ---------------------------------
// Creates bins aligned to each source interval's start position rather than chromosome position 0.
// Each source interval generates bins independently, and the iterator tracks which source interval
// spawned each output bin via get_cur_src_interval_idx().

class TrackExpressionIntervalRelativeBinIterator : public TrackExpression1DIterator {
public:
	enum PartialBinsMode { CLIP = 0, EXACT = 1 };

	TrackExpressionIntervalRelativeBinIterator() :
		TrackExpression1DIterator(INTERVALS1D),
		m_src_intervals(NULL),
		m_binsize(0),
		m_partial_bins_mode(CLIP),
		m_cur_bin(0),
		m_num_bins(0),
		m_cur_src_interval_idx(-1) {}

	bool begin(const GIntervals &src_intervals, int64_t binsize,
	           PartialBinsMode partial_bins_mode, GIntervalsFetcher1D &scope);
	virtual bool next();

	int64_t get_bin_size() const { return m_binsize; }
	int64_t get_cur_src_interval_idx() const { return m_cur_src_interval_idx; }
	const GInterval &get_cur_src_interval() const { return *m_icur_src_interval; }

private:
	GIntervals                 *m_src_intervals;
	GIntervals::const_iterator  m_icur_src_interval;
	int64_t                     m_binsize;
	PartialBinsMode             m_partial_bins_mode;
	int64_t                     m_cur_bin;     // Current bin within source interval (0-based)
	int64_t                     m_num_bins;    // Total bins for current source interval
	int64_t                     m_cur_src_interval_idx;  // Index of current source interval (0-based)

	bool advance_to_next_src_interval();
	void compute_bins_for_current_interval();
};

#endif /* TRACKEXPRESSIONINTERVALRELATIVEBINITERATOR_H_ */
